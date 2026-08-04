from __future__ import annotations

import pandas as pd

from reach_gap.features import extract_features
from reach_gap.geometry import GeometryConfig, simulate_geometry
from reach_gap.retrospective import run_retrospective
from reach_gap.schemas import ModelParameters, ParameterRange
from reach_gap.sensitivity import sobol_sensitivity


def test_small_sobol_analysis() -> None:
    geometry = simulate_geometry(GeometryConfig(size=16, cell_count=100, seed=6))
    features = extract_features(geometry)
    ranges = {
        "diffusion_um2_s": ParameterRange(8.0, 20.0, "loguniform"),
        "internalisation_s": ParameterRange(2.0e-5, 2.0e-4, "loguniform"),
    }
    result = sobol_sensitivity(features, ModelParameters(), ranges=ranges, sample_power=2)
    assert result["status"] in {"OK", "DEGENERATE_OUTPUT"}
    assert set(result["indices"]) == set(ranges)


def _retrospective_table() -> pd.DataFrame:
    rows = []
    for index in range(8):
        rows.append(
            {
                "programme_id": f"P{index}",
                "target": "T",
                "molecule": f"M{index}",
                "format": "IgG" if index < 4 else "ADC",
                "payload": "none",
                "dose": "verified",
                "indication": "A",
                "line_of_therapy": "1",
                "target_expression_evidence": "verified",
                "outcome": "success" if index % 2 == 0 else "failure",
                "reason_for_discontinuation": "verified",
                "source_identifiers": "verified",
                "expression_reach_gap": 0.1 + 0.05 * index,
            }
        )
    return pd.DataFrame(rows)


def test_retrospective_evaluation_only() -> None:
    result = run_retrospective(_retrospective_table(), permutations=50)
    assert result.status == "EXPLORATORY_ONLY"
    assert result.n_programmes == 8
    assert result.permutation_p is not None


def test_synthetic_rows_excluded_by_default() -> None:
    table = _retrospective_table()
    table.loc[0, "programme_id"] = "SYNTHETIC_0"
    result = run_retrospective(table)
    assert result.status == "NOT_COMPUTED"


def test_clinical_labels_cannot_enter_parameter_setting() -> None:
    import inspect

    import reach_gap.config as parameter_config

    signature = inspect.signature(parameter_config.sample_parameters)
    assert set(signature.parameters) == {"rng", "base", "ranges"}
    source = inspect.getsource(parameter_config).lower()
    assert "outcome" not in source
    assert "retrospective" not in source
    assert "success" not in source
    assert "failure" not in source
