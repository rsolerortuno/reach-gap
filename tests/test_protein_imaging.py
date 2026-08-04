from __future__ import annotations

import io
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from PIL import Image

from reach_gap import protein_imaging
from reach_gap.protein_imaging import (
    calibrate_image_threshold,
    distance_band_summary,
    distance_to_mask,
    parse_external_ome_channel,
    read_pyramidal_plane,
    sample_local_mean,
    structural_mask,
)
from reach_gap.real_rcc_imaging import discover_rcc_protein_image_inputs


def _ome_xml(file_name: str = "ch0028_cd31.ome.tif") -> str:
    channels = "".join(
        f'<Channel ID="Channel:0:{index}" Name="channel_{index}" SamplesPerPixel="1"/>'
        for index in range(28)
    )
    channels += '<Channel ID="Channel:0:28" Name="CD31" SamplesPerPixel="1"/>'
    return (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
        '<Image ID="Image:0"><Pixels ID="Pixels:0" DimensionOrder="XYCZT" '
        'Type="uint16" SizeX="6" SizeY="6" SizeC="29" SizeZ="1" SizeT="1" '
        'PhysicalSizeX="0.2125" PhysicalSizeY="0.2125">'
        f"{channels}"
        f'<TiffData FirstC="28"><UUID FileName="{file_name}">urn:uuid:test</UUID></TiffData>'
        "</Pixels></Image></OME>"
    )


def test_parse_external_ome_channel() -> None:
    observed = parse_external_ome_channel(_ome_xml(), "ch0028_cd31.ome.tif", level=4)
    assert observed.channel_index == 28
    assert observed.channel_name == "CD31"
    assert observed.size_x == 6
    assert observed.level_pixel_size_x_um == pytest.approx(3.4)


def test_decode_tiled_jpeg2000_without_imagecodecs(tmp_path: Path, monkeypatch) -> None:
    payloads: list[bytes] = []
    expected = np.zeros((6, 6), dtype=np.uint16)
    for tile_index in range(4):
        tile = np.full((4, 4), tile_index + 1, dtype=np.uint16)
        buffer = io.BytesIO()
        Image.fromarray(tile).save(buffer, format="JPEG2000")
        payloads.append(buffer.getvalue())
        y_start = (tile_index // 2) * 4
        x_start = (tile_index % 2) * 4
        expected[y_start : min(y_start + 4, 6), x_start : min(x_start + 4, 6)] = tile[
            : min(4, 6 - y_start), : min(4, 6 - x_start)
        ]
    source_path = tmp_path / "ch0028_cd31.ome.tif"
    offsets: list[int] = []
    with source_path.open("wb") as handle:
        for payload in payloads:
            offsets.append(handle.tell())
            handle.write(payload)

    page = SimpleNamespace(
        tilewidth=4,
        tilelength=4,
        shape=(6, 6),
        dtype=np.dtype(np.uint16),
        dataoffsets=offsets,
        databytecounts=[len(payload) for payload in payloads],
        pages=[],
    )

    class FakeTiffFile:
        def __init__(self, _path: Path) -> None:
            self.pages = [page]
            self.ome_metadata = _ome_xml()

        def __enter__(self) -> FakeTiffFile:
            return self

        def __exit__(self, *_args: object) -> None:
            return None

    fake_tifffile = SimpleNamespace(TiffFile=FakeTiffFile)
    monkeypatch.setattr(
        protein_imaging,
        "_require_imaging_dependencies",
        lambda: (fake_tifffile, Image),
    )
    observed, metadata = read_pyramidal_plane(source_path, level=0)
    assert np.array_equal(observed, expected)
    assert metadata.channel_name == "CD31"


def test_calibration_structure_and_distance_are_deterministic() -> None:
    samples = np.array([0.0, 0.1, 0.2, 1.0, 1.2, 2.0])
    reference = np.array([False, False, False, True, True, True])
    calibration = calibrate_image_threshold(samples, reference, negative_quantile=0.9)
    assert 0.1 < calibration.image_threshold <= 0.2
    assert calibration.false_positive_rate <= 1 / 3
    assert calibration.true_positive_rate == 1.0

    image = np.zeros((9, 9), dtype=np.float64)
    image[4, 3:6] = 5.0
    mask, diagnostics = structural_mask(
        image,
        threshold=1.0,
        pixel_size_um=2.0,
        minimum_component_area_um2=8.0,
        closing_radius_um=0.0,
    )
    assert diagnostics["retained_components"] == 1
    distance = distance_to_mask(mask, pixel_size_um=2.0)
    assert distance[4, 4] == 0.0
    assert distance[0, 4] == pytest.approx(8.0)

    rows = distance_band_summary(distance, {"signal": image}, bands_um=((0, 5), (5, None)))
    assert sum(int(row["pixels"]) for row in rows) == image.size


def test_local_mean_and_input_discovery(tmp_path: Path) -> None:
    image = np.arange(25, dtype=np.float64).reshape(5, 5)
    observed = sample_local_mean(
        image,
        [2.0],
        [2.0],
        pixel_size_x_um=1.0,
        pixel_size_y_um=1.0,
        radius_pixels=1,
    )
    assert observed[0] == pytest.approx(12.0)

    morphology = tmp_path / "morphology"
    qc = tmp_path / "qc"
    background = tmp_path / "background"
    morphology.mkdir()
    qc.mkdir()
    background.mkdir()
    names = {
        "CD31": "ch0028_cd31.ome.tif",
        "alphaSMA": "ch0032_alphasma.ome.tif",
        "Vimentin": "ch0031_vimentin.ome.tif",
        "PanCK": "ch0030_panck.ome.tif",
        "PD-L1": "ch0021_pd-l1.ome.tif",
        "VISTA": "ch0020_vista.ome.tif",
    }
    for file_name in names.values():
        (morphology / file_name).write_bytes(b"image")
        (qc / file_name).write_bytes(b"mask")
    for colour in ("blu", "grn", "yel", "red", "nuv"):
        (background / f"background_02_{colour}.tiff").write_bytes(b"background")

    channels, masks, backgrounds = discover_rcc_protein_image_inputs(
        morphology,
        qc_mask_dir=qc,
        background_dir=background,
    )
    assert set(channels) == set(names)
    assert set(masks) == set(names)
    assert set(backgrounds) == {"blu", "grn", "yel", "red", "nuv"}
