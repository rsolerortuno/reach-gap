from __future__ import annotations

from dataclasses import replace

from reach_gap.features import extract_features
from reach_gap.geometry import GeometryConfig, simulate_geometry
from reach_gap.indexing import compute_index
from reach_gap.schemas import ModelParameters


def test_index_has_interval_and_claims() -> None:
    geometry = simulate_geometry(GeometryConfig(size=20, cell_count=150, seed=3))
    features = extract_features(geometry, antigen_calibration_nM_per_signal=200.0)
    output, claims = compute_index(
        features,
        ModelParameters(vessel_concentration_nM=300.0),
        ranges={},
        draws=4,
        reachable_decision_threshold=0.4,
    )
    assert output.reachable_fraction.lower <= output.reachable_fraction.median
    assert output.reachable_fraction.median <= output.reachable_fraction.upper
    assert output.expression_reach_gap.median >= 0.0
    assert claims.unsupported


def test_uncalibrated_features_abstain() -> None:
    geometry = simulate_geometry(GeometryConfig(size=18, cell_count=100, seed=7))
    features = extract_features(geometry, antigen_calibration_nM_per_signal=None)
    output, _ = compute_index(features, ModelParameters(), ranges={}, draws=4)
    assert output.status == "INSUFFICIENT_EVIDENCE"
    assert "ANTIGEN_DENSITY_UNCALIBRATED" in output.abstention_reasons


def test_no_vessels_abstain() -> None:
    geometry = simulate_geometry(GeometryConfig(size=18, cell_count=100, seed=9))
    features = extract_features(geometry)
    no_vessels = replace(features, vessel_mask=features.vessel_mask & False)
    output, _ = compute_index(no_vessels, ModelParameters(), ranges={}, draws=4)
    assert output.status == "INSUFFICIENT_EVIDENCE"
    assert "NO_VASCULAR_STRUCTURES" in output.abstention_reasons
