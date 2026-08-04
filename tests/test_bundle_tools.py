from __future__ import annotations

import json
import zipfile
from pathlib import Path

from reach_gap.bundle_tools import extract_xenium_essentials_low_memory
from reach_gap.he_pathology import reassemble_split_file


def test_streaming_essential_extraction_and_reassembly(tmp_path: Path) -> None:
    payload = (b"0123456789abcdef" * 5000) + b"tail"
    small = b"cell_id,x_centroid,y_centroid,cell_area,nucleus_area\n"
    zip_path = tmp_path / "bundle.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("outs/cell_feature_matrix.h5", payload)
        archive.writestr("outs/cells.csv.gz", small)

    output_dir = tmp_path / "package"
    result = extract_xenium_essentials_low_memory(
        zip_path,
        output_dir,
        part_size_bytes=20_000,
        chunk_size=4096,
        selected_members=["outs/cell_feature_matrix.h5", "outs/cells.csv.gz"],
    )
    assert result["status"] == "ESSENTIAL_XENIUM_PACKAGE_EXTRACTED"
    assert result["selected_member_count"] == 2
    split = next(item for item in result["members"] if item["mode"] == "split_parts")
    manifest_path = Path(split["manifest"])
    manifest = json.loads(manifest_path.read_text())
    assert manifest["source"]["size"] == len(payload)
    reconstructed = tmp_path / "reconstructed.h5"
    reassemble_split_file(manifest_path, Path(split["parts_dir"]), reconstructed)
    assert reconstructed.read_bytes() == payload
    single = next(item for item in result["members"] if item["mode"] == "single_file")
    assert Path(single["path"]).read_bytes() == small
