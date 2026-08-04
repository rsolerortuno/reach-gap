from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from reach_gap.relative_accessibility import (
    _minmax,
    _rank_from_scores,
    _simulate,
    compute_relative_accessibility,
    run_relative_accessibility,
)


def repository_root() -> Path:
    return Path(__file__).resolve().parents[1]


def test_relative_accessibility_has_stable_vista_top() -> None:
    result = compute_relative_accessibility(repository_root(), draws=20_000, seed=17)
    assert result.status.endswith("ABSOLUTE_REACHABILITY_NOT_COMPUTED")
    assert result.target_count == 4
    assert result.vessel_definition_count == 6
    assert result.stable_top_target == "VISTA"
    assert result.stable_top_probability >= 0.97
    summaries = {row.target: row for row in result.target_summaries}
    assert summaries["VISTA"].median_rank == 1.0
    assert summaries["PD-L1"].tumour_positive_fraction == pytest.approx(0.0464339213017698)
    assert result.pairwise_win_probability["VISTA"]["PD-1"] >= 0.98
    assert set(result.leave_one_component_out_top_probability) == {
        "tumour_positive_fraction",
        "median_structural_proximity",
        "within_50um",
    }


def test_relative_accessibility_is_deterministic() -> None:
    one = compute_relative_accessibility(repository_root(), draws=2_000, seed=4)
    two = compute_relative_accessibility(repository_root(), draws=2_000, seed=4)
    assert one.model_dump() == two.model_dump()


def test_write_relative_accessibility_package(tmp_path: Path) -> None:
    payload = run_relative_accessibility(repository_root(), tmp_path, draws=2_000, seed=7)
    assert payload["summary"]["stable_top_target"] == "VISTA"
    result = json.loads((tmp_path / "relative_accessibility_v0.8.json").read_text())
    assert result["target_count"] == 4
    with (tmp_path / "target_rank_summary_v0.8.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["target"] == "VISTA"
    assert (tmp_path / "target_top_rank_probability_v0.8.svg").is_file()


def test_minmax_and_rank_helpers() -> None:
    assert np.allclose(_minmax(np.array([2.0, 4.0])), [0.0, 1.0])
    assert np.allclose(_minmax(np.array([2.0, 4.0]), higher_is_better=False), [1.0, 0.0])
    assert np.allclose(_minmax(np.array([3.0, 3.0])), [0.5, 0.5])
    assert _rank_from_scores(np.array([0.2, 0.8, 0.5])).tolist() == [3, 1, 2]


def test_simulate_rejects_invalid_draws() -> None:
    matrix = np.eye(3, dtype=float)
    with pytest.raises(ValueError):
        _simulate([matrix], draws=99, seed=1)
    with pytest.raises(ValueError):
        _simulate([matrix], draws=100, seed=1, active_components=())


def test_missing_relative_input_columns_fail(tmp_path: Path) -> None:
    result_dir = tmp_path / "results/real_rcc_xenium"
    result_dir.mkdir(parents=True)
    pd.DataFrame({"target": ["A"]}).to_csv(result_dir / "target_spatial_summary.csv", index=False)
    pd.DataFrame({"definition": ["d"]}).to_csv(
        result_dir / "vessel_calling_sensitivity.csv", index=False
    )
    with pytest.raises(ValueError, match="Missing target columns"):
        compute_relative_accessibility(tmp_path, draws=100)
