from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from tifffile import imwrite

from reach_gap.he_pathology import (
    detect_lumen_candidates,
    prepare_he_pathology_rcc,
    reassemble_split_file,
    tissue_mask,
)


def test_reassemble_split_file(tmp_path: Path) -> None:
    payload = (b"reach-gap-split-data-" * 1000) + b"end"
    source_md5 = hashlib.md5(payload).hexdigest()
    parts: list[dict[str, object]] = []
    chunk_size = 7000
    for index, start in enumerate(range(0, len(payload), chunk_size), start=1):
        chunk = payload[start : start + chunk_size]
        name = f"source.bin.part{index:03d}-of-004"
        (tmp_path / name).write_bytes(chunk)
        parts.append(
            {
                "index": index,
                "name": name,
                "size": len(chunk),
                "sha256": hashlib.sha256(chunk).hexdigest(),
            }
        )
    manifest = {
        "source": {"size": len(payload), "md5Checksum": source_md5},
        "parts": parts,
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    output = tmp_path / "reassembled.bin"
    report = reassemble_split_file(manifest_path, tmp_path, output)
    assert output.read_bytes() == payload
    assert report["status"] == "VERIFIED"
    assert report["md5"] == source_md5


def test_lumen_candidate_is_explicitly_unvalidated() -> None:
    rgb = np.full((64, 64, 3), 180, dtype=np.uint8)
    rgb[24:40, 24:40] = 250
    tissue = np.ones((64, 64), dtype=np.bool_)
    tissue[24:40, 24:40] = False
    hematoxylin = np.full((64, 64), 0.35, dtype=np.float64)
    lumen_mask, candidates = detect_lumen_candidates(
        rgb,
        tissue,
        hematoxylin,
        pixel_size_x_um=2.0,
        pixel_size_y_um=2.0,
    )
    assert lumen_mask.any()
    assert candidates
    assert all(candidate.confidence <= 1.0 for candidate in candidates)


def test_prepare_he_pathology_abstains(tmp_path: Path) -> None:
    height, width = 128, 96
    rgb = np.full((height, width, 3), 250, dtype=np.uint8)
    rgb[10:118, 8:88] = np.array([205, 135, 175], dtype=np.uint8)
    rgb[50:65, 42:56] = 250
    he_path = tmp_path / "synthetic.ome.tif"
    imwrite(
        he_path,
        rgb,
        photometric="rgb",
        metadata={
            "axes": "YXS",
            "PhysicalSizeX": 0.5,
            "PhysicalSizeXUnit": "µm",
            "PhysicalSizeY": 0.5,
            "PhysicalSizeYUnit": "µm",
        },
        ome=True,
    )
    annotations = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "id": "tumour-1",
                "properties": {"name": "Tumor"},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[[8, 10], [88, 10], [88, 118], [8, 118], [8, 10]]],
                },
            }
        ],
    }
    annotation_path = tmp_path / "annotations.geojson"
    annotation_path.write_text(json.dumps(annotations), encoding="utf-8")
    alignment_path = tmp_path / "alignment.csv"
    np.savetxt(alignment_path, np.eye(3), delimiter=",")

    result = prepare_he_pathology_rcc(
        he_path=he_path,
        annotation_path=annotation_path,
        alignment_path=alignment_path,
        output_dir=tmp_path / "output",
        analysis_level=0,
    )
    assert result["status"].endswith("TARGET_INDEX_NOT_COMPUTED")
    assert result["absolute_index"]["status"] == "NOT_COMPUTED"
    reasons = set(result["absolute_index"]["reasons"])
    assert "TARGET_EXPRESSION_NOT_AVAILABLE" in reasons
    claims = json.loads((tmp_path / "output" / "claims.json").read_text())
    assert "Reachable fraction." in claims["unsupported"]
    assert (tmp_path / "output" / "qc" / "pathology_overlay.png").exists()


def test_tissue_mask_rejects_white_background() -> None:
    rgb = np.full((32, 32, 3), 255, dtype=np.uint8)
    rgb[8:24, 8:24] = np.array([190, 110, 160], dtype=np.uint8)
    mask = tissue_mask(rgb)
    assert not mask[0, 0]
    assert mask[16, 16]
