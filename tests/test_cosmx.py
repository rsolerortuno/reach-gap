from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from reach_gap.cosmx import (
    polygon_metrics,
    prepare_cosmx_external_validation,
    relative_rna_geometry,
    segmentation_summary,
)


def test_polygon_metrics_and_segmentation(tmp_path: Path) -> None:
    polygon = pd.DataFrame(
        {
            "fov": [1, 1, 1, 1, 1, 1, 1, 1],
            "cellID": [1, 1, 1, 1, 2, 2, 2, 2],
            "x_local_px": [0, 2, 2, 0, 10, 12, 12, 10],
            "y_local_px": [0, 0, 2, 2, 0, 0, 2, 2],
        }
    )
    path = tmp_path / "polygons.csv.gz"
    polygon.to_csv(path, index=False, compression="gzip")
    metrics = polygon_metrics(path)
    assert np.allclose(metrics["polygon_area_px2"], [4.0, 4.0])
    assert np.allclose(metrics["polygon_centroid_x_px"], [1.0, 11.0])
    metadata = pd.DataFrame(
        {
            "fov": [1, 1],
            "cell_ID": [1, 2],
            "Area": [4.0, 4.0],
            "CenterX_local_px": [1.0, 11.0],
            "CenterY_local_px": [1.0, 1.0],
            "CenterX_global_px": [1.0, 11.0],
            "CenterY_global_px": [1.0, 1.0],
        }
    )
    summary = segmentation_summary(metadata, metrics)
    assert summary["matched_fraction"] == 1.0
    assert summary["within_10pct_area"] == 1.0
    assert summary["within_2px_centroid"] == 1.0


def test_relative_geometry_is_in_pixels_and_not_perfusion() -> None:
    metadata = pd.DataFrame(
        {
            "fov": [1, 1, 1],
            "cell_ID": [1, 2, 3],
            "Area": [4.0] * 3,
            "CenterX_local_px": [0.0, 10.0, 20.0],
            "CenterY_local_px": [0.0, 0.0, 0.0],
            "CenterX_global_px": [0.0, 10.0, 20.0],
            "CenterY_global_px": [0.0, 0.0, 0.0],
        }
    )
    expression = pd.DataFrame(
        {
            "fov": [1, 1, 1],
            "cell_ID": [1, 2, 3],
            "PECAM1": [1, 0, 0],
            "VWF": [1, 0, 0],
            "KDR": [1, 0, 0],
            "CDH5": [1, 0, 0],
            "ESAM": [1, 0, 0],
            "ENG": [1, 0, 0],
            "RAMP2": [1, 0, 0],
            "ACKR1": [1, 0, 0],
            "ERBB2": [0, 1, 1],
        }
    )
    targets, definitions = relative_rna_geometry(metadata, expression)
    balanced = targets[
        (targets["definition"] == "balanced_4_of_8") & (targets["target"] == "ERBB2")
    ]
    assert balanced.iloc[0]["median_distance_px"] == 15.0
    assert balanced.iloc[0]["distance_unit"] == "px"
    assert definitions.iloc[0]["source_semantics"] == "RNA_ENDOTHELIAL_PROXY_NOT_PERFUSION"


def test_prepare_cosmx_workflow_small_fixture(tmp_path: Path) -> None:
    sample = tmp_path / "sample"
    sample.mkdir()
    metadata = pd.DataFrame(
        {
            "fov": [1, 1],
            "cell_ID": [1, 2],
            "Area": [4.0, 4.0],
            "CenterX_local_px": [1.0, 11.0],
            "CenterY_local_px": [1.0, 1.0],
            "CenterX_global_px": [1.0, 11.0],
            "CenterY_global_px": [1.0, 1.0],
        }
    )
    expression = pd.DataFrame(
        {
            "fov": [1, 1],
            "cell_ID": [1, 2],
            "PECAM1": [1, 0],
            "VWF": [1, 0],
            "KDR": [1, 0],
            "CDH5": [1, 0],
            "ESAM": [1, 0],
            "ENG": [1, 0],
            "RAMP2": [1, 0],
            "ACKR1": [1, 0],
            "ERBB2": [0, 1],
        }
    )
    polygon = pd.DataFrame(
        {
            "fov": [1] * 8,
            "cellID": [1] * 4 + [2] * 4,
            "x_local_px": [0, 2, 2, 0, 10, 12, 12, 10],
            "y_local_px": [0, 0, 2, 2, 0, 0, 2, 2],
        }
    )
    metadata_path = sample / "metadata.csv.gz"
    expression_path = sample / "expression.csv.gz"
    polygon_path = sample / "polygons.csv.gz"
    metadata.to_csv(metadata_path, index=False, compression="gzip")
    expression.to_csv(expression_path, index=False, compression="gzip")
    polygon.to_csv(polygon_path, index=False, compression="gzip")
    output = tmp_path / "out"
    result = prepare_cosmx_external_validation(
        {
            "fixture": {
                "metadata": metadata_path,
                "expression": expression_path,
                "polygons": polygon_path,
            }
        },
        output,
    )
    assert result["absolute_index"]["status"] == "NOT_COMPUTED"
    assert result["vascular_source"] == "RNA_ENDOTHELIAL_PROXY_NOT_PERFUSION"
    assert (output / "cosmx_target_proxy_geometry.csv").exists()


def test_geometry_sensitivity_audit_flags_unstable_definitions() -> None:
    from reach_gap.cosmx import geometry_sensitivity_audit

    table = pd.DataFrame(
        {
            "sample": ["a"] * 6 + ["b"] * 6,
            "definition": ["inclusive_2_of_8", "balanced_4_of_8", "strict_6_of_8"] * 4,
            "target": ["ERBB2"] * 3 + ["EGFR"] * 3 + ["ERBB2"] * 3 + ["EGFR"] * 3,
            "median_distance_px": [
                10.0,
                20.0,
                100.0,
                20.0,
                30.0,
                80.0,
                15.0,
                25.0,
                90.0,
                10.0,
                35.0,
                120.0,
            ],
            "within_100px": [0.9, 0.7, 0.4, 0.8, 0.6, 0.3, 0.85, 0.65, 0.35, 0.9, 0.55, 0.2],
        }
    )
    sensitivity, ranks, summary = geometry_sensitivity_audit(table)
    assert sensitivity["definition_distance_ratio"].max() >= 10.0
    assert len(ranks) == 1
    assert summary["interpretation"] == "DEFINITION_SENSITIVE_AND_TARGET_RANKING_NOT_STABLE"
