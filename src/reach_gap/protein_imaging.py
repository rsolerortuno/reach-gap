"""Low-memory utilities for Xenium morphology-focus protein images.

The public RCC images are tiled pyramidal OME-TIFF files whose tiles use JPEG 2000.
This module decodes one pyramid level at a time without loading the native plane and
without relying on the optional ``imagecodecs`` package. It uses Pillow only for the
individual JPEG 2000 tile payloads.

All thresholds in this module are measurement-calibration thresholds. They may use
matched Xenium cell-level protein measurements, but never clinical outcomes or drug
response labels.
"""

from __future__ import annotations

import io
import math
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeAlias, cast

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import ndimage as ndi

FloatArray: TypeAlias = NDArray[np.float64]
BoolArray: TypeAlias = NDArray[np.bool_]


@dataclass(frozen=True)
class OmeChannelMetadata:
    """Metadata for one externally stored OME channel plane."""

    channel_index: int
    channel_name: str
    size_x: int
    size_y: int
    physical_size_x_um: float
    physical_size_y_um: float
    pyramid_level: int

    @property
    def level_pixel_size_x_um(self) -> float:
        return float(self.physical_size_x_um * (2**self.pyramid_level))

    @property
    def level_pixel_size_y_um(self) -> float:
        return float(self.physical_size_y_um * (2**self.pyramid_level))


@dataclass(frozen=True)
class ImageCalibration:
    """Transparent calibration of image signal against matched cell measurements."""

    negative_quantile: float
    image_threshold: float
    reference_positive_fraction: float
    image_positive_fraction: float
    false_positive_rate: float
    true_positive_rate: float


def _require_imaging_dependencies() -> tuple[Any, Any]:
    try:
        import tifffile
        from PIL import Image
    except ImportError as exc:  # pragma: no cover - depends on optional installation
        raise RuntimeError(
            "Protein-image analysis requires the 'xenium' optional dependencies "
            "(tifffile and Pillow)."
        ) from exc
    return tifffile, Image


def parse_external_ome_channel(ome_xml: str, file_name: str, level: int) -> OmeChannelMetadata:
    """Resolve the channel represented by one externally stored OME-TIFF file."""

    root = ET.fromstring(ome_xml)
    namespace = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    pixels = root.find(".//ome:Pixels", namespace)
    if pixels is None:
        raise ValueError("OME Pixels element is missing")
    channels = pixels.findall("ome:Channel", namespace)
    first_channel: int | None = None
    for data in pixels.findall("ome:TiffData", namespace):
        uuid_element = data.find("ome:UUID", namespace)
        if uuid_element is None:
            continue
        candidate = uuid_element.attrib.get("FileName")
        if candidate == file_name:
            first_channel = int(data.attrib["FirstC"])
            break
    if first_channel is None:
        raise ValueError(f"OME metadata does not reference external file {file_name!r}")
    if first_channel >= len(channels):
        raise ValueError("OME FirstC exceeds the declared channel count")
    name = channels[first_channel].attrib.get("Name", f"channel_{first_channel}")
    return OmeChannelMetadata(
        channel_index=first_channel,
        channel_name=name,
        size_x=int(pixels.attrib["SizeX"]),
        size_y=int(pixels.attrib["SizeY"]),
        physical_size_x_um=float(pixels.attrib["PhysicalSizeX"]),
        physical_size_y_um=float(pixels.attrib["PhysicalSizeY"]),
        pyramid_level=level,
    )


def read_pyramidal_plane(
    path: Path, *, level: int = 3, ome_file_name: str | None = None
) -> tuple[NDArray[Any], OmeChannelMetadata]:
    """Decode one pyramid level from a tiled single-channel external OME-TIFF.

    The function reads compressed tile payloads directly and decodes each with Pillow.
    Peak memory is therefore approximately the requested level plus one 1024x1024 tile.
    """

    if level < 0:
        raise ValueError("Pyramid level must be non-negative")
    tifffile, image_module = _require_imaging_dependencies()
    with tifffile.TiffFile(path) as handle:
        base = handle.pages[0]
        page = base if level == 0 else base.pages[level - 1]
        if page.tilewidth is None or page.tilelength is None:
            raise ValueError("Only tiled OME-TIFF planes are supported")
        height, width = (int(page.shape[-2]), int(page.shape[-1]))
        tile_width = int(page.tilewidth)
        tile_height = int(page.tilelength)
        tiles_x = math.ceil(width / tile_width)
        tiles_y = math.ceil(height / tile_height)
        if len(page.dataoffsets) != tiles_x * tiles_y:
            raise ValueError("Unexpected tile count for image dimensions")
        output = np.zeros((height, width), dtype=page.dtype)
        with path.open("rb") as source:
            for tile_index, (offset, byte_count) in enumerate(
                zip(page.dataoffsets, page.databytecounts, strict=True)
            ):
                if byte_count == 0:
                    continue
                source.seek(int(offset))
                payload = source.read(int(byte_count))
                with image_module.open(io.BytesIO(payload)) as tile_image:
                    tile = np.asarray(tile_image)
                if tile.ndim != 2:
                    raise ValueError("Expected a single-channel JPEG 2000 tile")
                y_start = (tile_index // tiles_x) * tile_height
                x_start = (tile_index % tiles_x) * tile_width
                y_size = min(tile_height, height - y_start)
                x_size = min(tile_width, width - x_start)
                output[y_start : y_start + y_size, x_start : x_start + x_size] = tile[
                    :y_size, :x_size
                ]
        if not handle.ome_metadata:
            raise ValueError("OME metadata is missing")
        metadata = parse_external_ome_channel(
            handle.ome_metadata, ome_file_name or path.name, level
        )
    return output, metadata


def image_mode(values: NDArray[Any], valid_mask: BoolArray | None = None) -> float:
    """Return the most common integer intensity among valid pixels."""

    array = np.asarray(values)
    selected = array.ravel() if valid_mask is None else array[np.asarray(valid_mask, dtype=bool)]
    selected = selected[np.isfinite(selected)]
    if selected.size == 0:
        raise ValueError("No valid image values")
    selected_int = np.asarray(selected, dtype=np.int64)
    if selected_int.min() < 0:
        raise ValueError("Image mode expects non-negative integer intensities")
    return float(np.argmax(np.bincount(selected_int)))


def sample_local_mean(
    image: NDArray[Any],
    x_um: ArrayLike,
    y_um: ArrayLike,
    *,
    pixel_size_x_um: float,
    pixel_size_y_um: float,
    radius_pixels: int = 1,
) -> FloatArray:
    """Sample a square local mean around cell centroids."""

    if radius_pixels < 0:
        raise ValueError("radius_pixels must be non-negative")
    array = np.asarray(image)
    x = np.rint(np.asarray(x_um, dtype=np.float64) / pixel_size_x_um).astype(np.int64)
    y = np.rint(np.asarray(y_um, dtype=np.float64) / pixel_size_y_um).astype(np.int64)
    if np.any(x < 0) or np.any(x >= array.shape[1]) or np.any(y < 0) or np.any(y >= array.shape[0]):
        raise ValueError("At least one cell coordinate lies outside the image")
    samples: list[NDArray[Any]] = []
    for y_offset in range(-radius_pixels, radius_pixels + 1):
        yy = np.clip(y + y_offset, 0, array.shape[0] - 1)
        for x_offset in range(-radius_pixels, radius_pixels + 1):
            xx = np.clip(x + x_offset, 0, array.shape[1] - 1)
            samples.append(array[yy, xx])
    return np.asarray(np.mean(np.stack(samples).astype(np.float32), axis=0), dtype=np.float64)


def calibrate_image_threshold(
    image_samples: ArrayLike,
    reference_positive: ArrayLike,
    *,
    negative_quantile: float = 0.995,
) -> ImageCalibration:
    """Set an image threshold from the declared tail of matched negative cells.

    The threshold is not optimized for accuracy. It fixes a transparent empirical
    false-positive target among cells negative by the matched Xenium aggregate signal.
    """

    if not 0.5 < negative_quantile < 1.0:
        raise ValueError("negative_quantile must be between 0.5 and 1")
    samples = np.asarray(image_samples, dtype=np.float64)
    positive = np.asarray(reference_positive, dtype=bool)
    if samples.shape != positive.shape:
        raise ValueError("Image samples and reference labels must have equal shape")
    if not np.any(positive) or not np.any(~positive):
        raise ValueError("Calibration requires both positive and negative reference cells")
    threshold = float(np.quantile(samples[~positive], negative_quantile))
    image_positive = samples > threshold
    return ImageCalibration(
        negative_quantile=negative_quantile,
        image_threshold=threshold,
        reference_positive_fraction=float(np.mean(positive)),
        image_positive_fraction=float(np.mean(image_positive)),
        false_positive_rate=float(np.mean(image_positive[~positive])),
        true_positive_rate=float(np.mean(image_positive[positive])),
    )


def structural_mask(
    image: NDArray[Any],
    *,
    threshold: float,
    pixel_size_um: float,
    minimum_component_area_um2: float,
    closing_radius_um: float = 1.7,
    qc_mask: BoolArray | None = None,
) -> tuple[BoolArray, dict[str, float | int]]:
    """Create a conservative connected structural mask from image signal."""

    if pixel_size_um <= 0 or minimum_component_area_um2 <= 0:
        raise ValueError("Physical sizes must be positive")
    array = np.asarray(image, dtype=np.float64)
    candidate = array > threshold
    if qc_mask is not None:
        candidate &= ~np.asarray(qc_mask, dtype=bool)
    labels_raw, raw_components = cast(tuple[NDArray[Any], int], ndi.label(candidate))
    labels = np.asarray(labels_raw, dtype=np.int64)
    sizes = np.bincount(labels.ravel())
    minimum_pixels = max(1, math.ceil(minimum_component_area_um2 / pixel_size_um**2))
    keep = sizes >= minimum_pixels
    keep[0] = False
    retained = keep[labels]
    radius_pixels = max(0, round(closing_radius_um / pixel_size_um))
    if radius_pixels:
        yy, xx = np.ogrid[-radius_pixels : radius_pixels + 1, -radius_pixels : radius_pixels + 1]
        footprint = xx * xx + yy * yy <= radius_pixels * radius_pixels
        retained = ndi.binary_closing(retained, structure=footprint)
    retained_labels_raw, retained_components = cast(tuple[NDArray[Any], int], ndi.label(retained))
    retained_labels = np.asarray(retained_labels_raw, dtype=np.int64)
    retained_sizes = np.bincount(retained_labels.ravel())[1:]
    diagnostics: dict[str, float | int] = {
        "threshold": float(threshold),
        "minimum_component_area_um2": float(minimum_component_area_um2),
        "minimum_component_pixels": minimum_pixels,
        "raw_components": int(raw_components),
        "retained_components": int(retained_components),
        "pixel_fraction": float(np.mean(retained)),
        "median_component_area_um2": (
            float(np.median(retained_sizes) * pixel_size_um**2) if retained_sizes.size else 0.0
        ),
    }
    return np.asarray(retained, dtype=bool), diagnostics


def distance_to_mask(mask: BoolArray, *, pixel_size_um: float) -> FloatArray:
    """Compute Euclidean distance to the nearest positive mask pixel in microns."""

    array = np.asarray(mask, dtype=bool)
    if not np.any(array):
        return np.full(array.shape, np.nan, dtype=np.float64)
    distance = np.asarray(ndi.distance_transform_edt(~array), dtype=np.float64)
    return distance * float(pixel_size_um)


def excess_signal(image: NDArray[Any], *, baseline: float) -> FloatArray:
    """Return non-negative signal above a declared image baseline."""

    return np.maximum(np.asarray(image, dtype=np.float64) - baseline, 0.0)


def local_gaussian_signal(
    image: NDArray[Any], *, pixel_size_um: float, sigma_um: float, baseline: float
) -> NDArray[np.float32]:
    """Smooth baseline-subtracted image signal over a physical neighbourhood."""

    if sigma_um <= 0:
        raise ValueError("sigma_um must be positive")
    signal = excess_signal(image, baseline=baseline).astype(np.float32)
    return np.asarray(
        ndi.gaussian_filter(signal, sigma=sigma_um / pixel_size_um, output=np.float32),
        dtype=np.float32,
    )


def sample_pixels(
    image: NDArray[Any],
    x_um: ArrayLike,
    y_um: ArrayLike,
    *,
    pixel_size_x_um: float,
    pixel_size_y_um: float,
) -> FloatArray:
    """Sample nearest image pixels at physical coordinates."""

    array = np.asarray(image)
    x = np.rint(np.asarray(x_um, dtype=np.float64) / pixel_size_x_um).astype(np.int64)
    y = np.rint(np.asarray(y_um, dtype=np.float64) / pixel_size_y_um).astype(np.int64)
    if np.any(x < 0) or np.any(x >= array.shape[1]) or np.any(y < 0) or np.any(y >= array.shape[0]):
        raise ValueError("At least one coordinate lies outside the image")
    return np.asarray(array[y, x], dtype=np.float64)


def distance_band_summary(
    distance_um: NDArray[Any],
    signals: dict[str, NDArray[Any]],
    *,
    valid_mask: BoolArray | None = None,
    bands_um: Sequence[tuple[float, float | None]] = (
        (0.0, 10.0),
        (10.0, 25.0),
        (25.0, 50.0),
        (50.0, 100.0),
        (100.0, None),
    ),
) -> list[dict[str, float | int | str]]:
    """Summarise image signals in preregistered distance bands."""

    distance = np.asarray(distance_um, dtype=np.float64)
    valid = np.isfinite(distance)
    if valid_mask is not None:
        valid &= np.asarray(valid_mask, dtype=bool)
    rows: list[dict[str, float | int | str]] = []
    for low, high in bands_um:
        active = valid & (distance >= low)
        if high is not None:
            active &= distance < high
        label = f"{low:g}-{high:g}" if high is not None else f">={low:g}"
        row: dict[str, float | int | str] = {
            "distance_band_um": label,
            "pixels": int(np.sum(active)),
        }
        for name, values in signals.items():
            selected = np.asarray(values, dtype=np.float64)[active]
            row[f"{name}_median"] = float(np.median(selected)) if selected.size else float("nan")
            row[f"{name}_q90"] = (
                float(np.quantile(selected, 0.90)) if selected.size else float("nan")
            )
        rows.append(row)
    return rows
