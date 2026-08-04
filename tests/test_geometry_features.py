from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from reach_gap.features import extract_features, ingest_cell_table, load_features
from reach_gap.geometry import GeometryConfig, simulate_geometry


def test_geometry_is_deterministic() -> None:
    config = GeometryConfig(size=24, cell_count=100, seed=4)
    first = simulate_geometry(config)
    second = simulate_geometry(config)
    assert np.array_equal(first.vessel_mask, second.vessel_mask)
    assert np.allclose(first.antigen, second.antigen)
    assert np.array_equal(first.cell_rows, second.cell_rows)


def test_geometry_rejects_invalid_settings() -> None:
    with pytest.raises(ValueError):
        simulate_geometry(GeometryConfig(size=8))
    with pytest.raises(ValueError):
        simulate_geometry(GeometryConfig(vessel_count=0))


def test_extract_features_calibration_and_distance() -> None:
    geometry = simulate_geometry(GeometryConfig(size=20, cell_count=80, seed=2))
    features = extract_features(geometry, antigen_calibration_nM_per_signal=100.0)
    assert features.antigen_calibrated
    assert np.all(features.antigen_nM >= 0.0)
    assert np.all(features.vessel_distance_um[features.vessel_mask] == 0.0)
    uncalibrated = extract_features(geometry, antigen_calibration_nM_per_signal=None)
    assert not uncalibrated.antigen_calibrated


def test_ingest_cell_table(tmp_path: Path) -> None:
    table = pd.DataFrame(
        {
            "cell_id": ["a", "b", "c"],
            "x_um": [0.0, 10.0, 20.0],
            "y_um": [0.0, 10.0, 20.0],
            "is_tumour": [True, True, False],
            "target_signal": [0.8, 0.2, 0.6],
            "vessel_signal": [1.0, 0.0, 0.0],
            "ecm_score": [0.1, 0.5, 0.2],
            "caf_score": [0.2, 0.6, 0.1],
        }
    )
    cells = tmp_path / "cells.csv"
    table.to_csv(cells, index=False)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "source_id": "unit-test",
                "platform": "Xenium",
                "segmentation_version": "v1",
                "vessel_definition": "marker",
                "antigen_calibration_nM_per_signal": 50.0,
            }
        ),
        encoding="utf-8",
    )
    feature_path, manifest_out = ingest_cell_table(cells, manifest, tmp_path / "out", grid_size=8)
    loaded = load_features(feature_path)
    assert loaded.antigen_calibrated
    assert np.any(loaded.vessel_mask)
    assert manifest_out.exists()


def test_ingest_missing_column(tmp_path: Path) -> None:
    cells = tmp_path / "cells.csv"
    pd.DataFrame({"cell_id": ["a"]}).to_csv(cells, index=False)
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError):
        ingest_cell_table(cells, manifest, tmp_path / "out")
