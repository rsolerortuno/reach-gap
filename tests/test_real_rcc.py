from __future__ import annotations

import numpy as np
import pandas as pd

from reach_gap.real_rcc import rna_protein_concordance, target_spatial_summary


def test_rna_protein_concordance_reports_low_overlap() -> None:
    table = pd.DataFrame(
        {
            "protein__PD_L1": [1.0, 0.0, 1.0, 0.0],
            "rna__CD274": [0.0, 1.0, 1.0, 0.0],
        }
    )
    observed = rna_protein_concordance(table)
    assert observed.loc[0, "target"] == "PD-L1"
    assert np.isclose(observed.loc[0, "nonzero_jaccard"], 1.0 / 3.0)


def test_target_summary_uses_tumour_region() -> None:
    scored = pd.DataFrame(
        {
            "pathology_region": ["Tumor", "UNANNOTATED"],
            "in_molecular_tumour_neighbourhood": [True, False],
            "cell_is_malignant_proxy": [True, False],
            "immune_score": [0.0, 1.0],
            "distance_to_vessel_um": [10.0, 60.0],
            "local_caf_score": [0.2, 0.3],
            "local_ecm_score": [0.0, 0.0],
            "target__PD_L1__signal": [0.7, 0.8],
            "target__PD_L1__positive": [True, True],
        }
    )
    diagnostics = {"targets": {"PD-L1": {"measurement": "protein_intensity"}}}
    observed = target_spatial_summary(scored, diagnostics)
    tumour = observed.loc[observed["subset"] == "tumour_region"].iloc[0]
    assert tumour["cells"] == 1
    assert tumour["median_distance_to_vessel_um_positive"] == 10.0


def test_prepare_essential_workflow_small_fixture(tmp_path, monkeypatch) -> None:
    import h5py
    from scipy.sparse import csc_matrix

    from reach_gap import real_rcc

    input_dir = tmp_path / "input"
    input_dir.mkdir()
    cells = pd.DataFrame(
        {
            "cell_id": ["cell-1", "cell-2", "cell-3", "cell-4"],
            "x_centroid": [0.0, 20.0, 40.0, 60.0],
            "y_centroid": [0.0, 0.0, 0.0, 0.0],
            "cell_area": [50.0] * 4,
            "nucleus_area": [20.0] * 4,
        }
    )
    cells.to_csv(input_dir / "cells.csv.gz", index=False, compression="gzip")
    names = np.array([b"PECAM1", b"CD31", b"PD-L1", b"PanCK"])
    feature_types = np.array(
        [
            b"Gene Expression",
            b"Protein Expression",
            b"Protein Expression",
            b"Protein Expression",
        ]
    )
    matrix = csc_matrix(
        np.array(
            [[8, 0, 1, 0], [12, 0, 1, 0], [0, 1, 9, 0], [0, 8, 7, 1]],
            dtype=np.float32,
        )
    )
    with h5py.File(input_dir / "cell_feature_matrix.h5", "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset(
            "barcodes", data=np.array([b"cell-1", b"cell-2", b"cell-3", b"cell-4"])
        )
        group.create_dataset("data", data=matrix.data)
        group.create_dataset("indices", data=matrix.indices)
        group.create_dataset("indptr", data=matrix.indptr)
        group.create_dataset("shape", data=np.array(matrix.shape, dtype=np.int64))
        features = group.create_group("features")
        features.create_dataset("id", data=np.array([b"1", b"2", b"3", b"4"]))
        features.create_dataset("name", data=names)
        features.create_dataset("feature_type", data=feature_types)
    (input_dir / "gene_panel.json").write_text("{}", encoding="utf-8")
    (input_dir / "protein_panel.json").write_text("{}", encoding="utf-8")
    (input_dir / "metrics_summary.csv").write_text("metric,value\n", encoding="utf-8")
    (input_dir / "experiment.xenium").write_text("{}", encoding="utf-8")

    fake_vessel = pd.DataFrame(
        {
            "definition": ["fixture"],
            "vessel_fraction": [0.25],
            "target": ["PD-L1"],
            "target_median_um": [20.0],
            "within_50um": [1.0],
            "within_100um": [1.0],
        }
    )
    monkeypatch.setattr(
        real_rcc,
        "vessel_definition_sensitivity",
        lambda scored, expression: (fake_vessel, {"status": "fixture"}),
    )
    output = tmp_path / "output"
    result = real_rcc.prepare_rcc_xenium_essentials(input_dir, output)
    assert result["cells"] == 4
    assert result["absolute_index"]["status"] == "NOT_COMPUTED"
    assert (output / "target_spatial_summary.csv").exists()
