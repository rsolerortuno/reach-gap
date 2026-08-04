from __future__ import annotations

import json
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix

from reach_gap.xenium import (
    DEFAULT_MARKERS,
    DEFAULT_TARGET_ALIASES,
    canonical_feature_name,
    extract_members,
    extract_selected_h5_features,
    geojson_vertices,
    infer_annotation_transform,
    inspect_zip,
    read_10x_h5_catalog,
    resolve_feature_indices,
    score_cells,
    select_essential_members,
    summarise_h5_features,
    validate_cell_barcode_identity,
)


def _write_h5(path: Path) -> None:
    names = np.array([b"PECAM1", b"CD31", b"COL1A1", b"PD-L1", b"PanCK"])
    ids = np.array([b"ENSG1", b"TXP1", b"ENSG2", b"TXP2", b"TXP3"])
    types = np.array(
        [
            b"Gene Expression",
            b"Protein Expression",
            b"Gene Expression",
            b"Protein Expression",
            b"Protein Expression",
        ]
    )
    barcodes = np.array([b"cell-1", b"cell-2", b"cell-3", b"cell-4"])
    dense = np.array(
        [
            [8, 0, 1, 0],
            [12, 0, 1, 0],
            [0, 5, 8, 1],
            [0, 1, 9, 0],
            [0, 8, 7, 1],
        ],
        dtype=np.float32,
    )
    matrix = csc_matrix(dense)
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset("barcodes", data=barcodes)
        group.create_dataset("data", data=matrix.data)
        group.create_dataset("indices", data=matrix.indices)
        group.create_dataset("indptr", data=matrix.indptr)
        group.create_dataset("shape", data=np.array(matrix.shape, dtype=np.int64))
        features = group.create_group("features")
        features.create_dataset("id", data=ids)
        features.create_dataset("name", data=names)
        features.create_dataset("feature_type", data=types)


def test_h5_catalog_summary_and_selected_extraction(tmp_path: Path) -> None:
    h5_path = tmp_path / "cell_feature_matrix.h5"
    _write_h5(h5_path)
    catalog = read_10x_h5_catalog(h5_path)
    assert catalog.matrix_shape == (5, 4)
    selected, resolution = resolve_feature_indices(catalog)
    assert len(selected) == 5
    assert resolution["selected_count"] == 5
    table = extract_selected_h5_features(h5_path, catalog, selected, chunk_cells=2)
    assert list(table["cell_id"]) == ["cell-1", "cell-2", "cell-3", "cell-4"]
    assert float(table.loc[0, "protein__CD31"]) == 12.0
    assert float(table.loc[2, "protein__PD_L1"]) == 9.0
    summary = summarise_h5_features(h5_path, catalog, chunk_nnz=3)
    pd_l1 = summary.loc[summary["feature_name"] == "PD-L1"].iloc[0]
    assert pd_l1["total_signal"] == 10.0
    assert pd_l1["positive_cells"] == 2


def test_score_cells_generates_targets_and_vascular_geometry(tmp_path: Path) -> None:
    h5_path = tmp_path / "cell_feature_matrix.h5"
    _write_h5(h5_path)
    catalog = read_10x_h5_catalog(h5_path)
    selected, resolution = resolve_feature_indices(
        catalog, markers=DEFAULT_MARKERS, targets=DEFAULT_TARGET_ALIASES
    )
    expression = extract_selected_h5_features(h5_path, catalog, selected)
    cells = pd.DataFrame(
        {
            "cell_id": ["cell-1", "cell-2", "cell-3", "cell-4"],
            "x_centroid": [0.0, 20.0, 40.0, 60.0],
            "y_centroid": [0.0, 0.0, 0.0, 0.0],
            "x_um": [0.0, 20.0, 40.0, 60.0],
            "y_um": [0.0, 0.0, 0.0, 0.0],
            "cell_area": [50.0] * 4,
            "nucleus_area": [20.0] * 4,
        }
    )
    scored, diagnostics = score_cells(cells, expression, resolution, local_neighbours=2)
    assert "target__PD_L1__signal" in scored
    assert "distance_to_vessel_um" in scored
    assert diagnostics["vessel_positive_cells"] >= 1
    assert diagnostics["targets"]["PD-L1"]["measurement"] == "protein_intensity"
    assert np.isfinite(scored["distance_to_vessel_um"]).all()


def test_zip_inventory_and_selective_extraction(tmp_path: Path) -> None:
    archive_path = tmp_path / "bundle.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("outs/cells.csv.gz", b"fake")
        archive.writestr("outs/cell_feature_matrix.h5", b"fake-h5")
        archive.writestr("outs/gene_panel.json", b"{}")
        archive.writestr("outs/analysis/clustering/graphclust/clusters.csv", b"cell,cluster\n")
        archive.writestr("outs/transcripts.parquet", b"huge-not-selected")
    inventory = inspect_zip(archive_path)
    members = select_essential_members(inventory)
    assert "outs/cells.csv.gz" in members
    assert "outs/cell_feature_matrix.h5" in members
    assert "outs/transcripts.parquet" not in members
    output_dir = tmp_path / "extracted"
    extracted = extract_members(archive_path, members, output_dir)
    assert Path(extracted["outs/gene_panel.json"]).exists()


def test_annotation_transform_inference(tmp_path: Path) -> None:
    geojson_path = tmp_path / "annotation.geojson"
    geojson_path.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [[[0, 0], [100, 0], [100, 100], [0, 100], [0, 0]]],
                        },
                        "properties": {"name": "Tumor"},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    vertices = geojson_vertices(geojson_path)
    cells = pd.DataFrame({"x_um": [0.0, 100.0], "y_um": [0.0, 100.0]})
    best, candidates = infer_annotation_transform(
        vertices,
        cells,
        np.eye(3),
        candidate_pixel_sizes_um=(1.0, 0.5),
    )
    assert candidates
    assert best is not None
    assert best.scale_x == 1.0


def test_antibody_capture_is_treated_as_protein() -> None:
    assert canonical_feature_name("PD-L1", "Antibody Capture") == "protein__PD_L1"
    assert canonical_feature_name("PD-L1", "Protein Expression") == "protein__PD_L1"


def test_cell_barcode_identity_rejects_same_length_mismatch(tmp_path: Path) -> None:
    h5_path = tmp_path / "cell_feature_matrix.h5"
    _write_h5(h5_path)
    catalog = read_10x_h5_catalog(h5_path)
    cells = pd.DataFrame({"cell_id": ["cell-1", "cell-2", "cell-3", "wrong-cell"]})
    try:
        validate_cell_barcode_identity(cells, catalog)
    except ValueError as error:
        assert "barcode identities differ" in str(error)
    else:
        raise AssertionError("Expected barcode mismatch to be rejected")


def test_prepare_rcc_xenium_end_to_end(tmp_path: Path, monkeypatch) -> None:
    from reach_gap import xenium

    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    source_h5 = tmp_path / "source.h5"
    _write_h5(source_h5)

    cells = pd.DataFrame(
        {
            "cell_id": ["cell-1", "cell-2", "cell-3", "cell-4"],
            "x_centroid": [0.0, 20.0, 40.0, 60.0],
            "y_centroid": [0.0, 0.0, 0.0, 0.0],
            "cell_area": [50.0] * 4,
            "nucleus_area": [20.0] * 4,
        }
    )
    cells_path = tmp_path / "cells.csv.gz"
    cells.to_csv(cells_path, index=False, compression="gzip")

    zip_name = "Xenium_V1_Human_Kidney_FFPE_Protein_updated_outs.zip"
    zip_path = raw_dir / zip_name
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.write(cells_path, "outs/cells.csv.gz")
        archive.write(source_h5, "outs/cell_feature_matrix.h5")
        archive.writestr("outs/gene_panel.json", "{}")
        archive.writestr("outs/protein_panel.json", "{}")
        archive.writestr("outs/metrics_summary.csv", "metric,value\n")

    he_name = "Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_image.ome.tif"
    alignment_name = "Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_imagealignment.csv"
    annotation_name = "Xenium_V1_Human_Kidney_FFPE_Protein_updated_annotation.geojson"
    (raw_dir / he_name).write_bytes(b"synthetic-he-placeholder")
    np.savetxt(raw_dir / alignment_name, np.eye(3), delimiter=",")
    (raw_dir / annotation_name).write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [[[-5, -5], [65, -5], [65, 5], [-5, 5], [-5, -5]]],
                        },
                        "properties": {"name": "Tumor"},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    for name in (zip_name, he_name, alignment_name, annotation_name):
        path = raw_dir / name
        monkeypatch.setitem(
            xenium.EXPECTED_RCC_FILES,
            name,
            {"size": path.stat().st_size, "md5": xenium.md5_file(path)},
        )

    result = xenium.prepare_rcc_xenium(
        raw_dir=raw_dir,
        output_dir=raw_dir / "reach-gap-analysis",
        verify_large_md5=True,
    )
    assert result["cells"] == 4
    assert result["absolute_index"]["status"] == "NOT_COMPUTED"
    assert (raw_dir / "reach-gap-analysis" / "processing_manifest.json").exists()
    assert result["targets"]["PD-L1"]
