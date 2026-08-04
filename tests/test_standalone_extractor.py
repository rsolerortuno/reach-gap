from __future__ import annotations

import importlib.util
import json
import sys
import zipfile
from pathlib import Path


def _load_script():
    path = Path(__file__).parents[1] / "scripts" / "extract_xenium_essentials_standalone.py"
    spec = importlib.util.spec_from_file_location("standalone_xenium_extractor", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_standalone_streams_selected_members_to_parts(tmp_path: Path) -> None:
    module = _load_script()
    archive_path = tmp_path / "bundle.zip"
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("outs/cells.csv.gz", b"cell-table")
        archive.writestr("outs/cell_feature_matrix.h5", b"0123456789" * 5)
        archive.writestr("outs/gene_panel.json", b"{}")
        archive.writestr("outs/morphology_focus/morphology.ome.tif", b"not-selected")

    output = tmp_path / "essential"
    result = module.extract(
        archive_path,
        output,
        part_size=19,
        chunk_size=7,
        list_only=False,
    )

    assert result["status"] == "ESSENTIAL_XENIUM_PACKAGE_EXTRACTED"
    assert result["selected_member_count"] == 3
    manifest = json.loads((output / "essential_package_manifest.json").read_text())
    matrix = next(row for row in manifest["members"] if row["member"].endswith(".h5"))
    assert matrix["mode"] == "split_parts"
    assert matrix["part_count"] == 3
    matrix_dir = output / "split" / "cell_feature_matrix.h5"
    reconstructed = b"".join(path.read_bytes() for path in sorted(matrix_dir.glob("*.part*")))
    assert reconstructed == b"0123456789" * 5
    assert not any("morphology" in row["member"] for row in manifest["members"])


def test_standalone_list_only_does_not_extract_bytes(tmp_path: Path) -> None:
    module = _load_script()
    archive_path = tmp_path / "bundle.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("outs/cells.parquet", b"cells")
        archive.writestr("outs/cell_feature_matrix.h5", b"matrix")

    output = tmp_path / "listed"
    result = module.extract(
        archive_path,
        output,
        part_size=10,
        chunk_size=4,
        list_only=True,
    )

    assert result["status"] == "ESSENTIAL_MEMBERS_LISTED_NOT_EXTRACTED"
    assert (output / "selected_essential_members.json").exists()
    assert not (output / "files").exists()
    assert not (output / "split").exists()
