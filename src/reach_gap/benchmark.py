"""First-class deterministic simulation benchmark and required ablations."""

from __future__ import annotations

import json
import time
import tracemalloc
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.ndimage import binary_dilation, binary_erosion, distance_transform_edt, gaussian_filter

from reach_gap.features import extract_features
from reach_gap.geometry import GeometryConfig, simulate_geometry
from reach_gap.indexing import compute_index
from reach_gap.retrospective import run_retrospective
from reach_gap.schemas import FloatArray, ModelParameters, SpatialFeatures
from reach_gap.sensitivity import sobol_sensitivity
from reach_gap.solver import solve_transport_robust


def _rmse(observed: FloatArray, predicted: FloatArray) -> float:
    return float(np.sqrt(np.mean((observed - predicted) ** 2)))


def _classification_accuracy(
    observed: FloatArray, predicted: FloatArray, threshold: float
) -> float:
    return float(np.mean((observed >= threshold) == (predicted >= threshold)))


def _penetration_depth(features: SpatialFeatures, bound: FloatArray, threshold: float) -> float:
    rows = features.cell_rows[features.cell_target_positive]
    cols = features.cell_cols[features.cell_target_positive]
    if len(rows) == 0:
        return 0.0
    values = bound[rows, cols]
    distances = features.vessel_distance_um[rows, cols]
    reached = values >= threshold
    return float(np.max(distances[reached])) if np.any(reached) else 0.0


def _reachable_fraction(features: SpatialFeatures, bound: FloatArray, threshold: float) -> float:
    target = features.cell_is_tumour & features.cell_target_positive
    if not np.any(target):
        return 0.0
    values = bound[features.cell_rows[target], features.cell_cols[target]]
    return float(np.mean(values >= threshold))


def _replace_vessels(features: SpatialFeatures, vessel_mask: NDArray[np.bool_]) -> SpatialFeatures:
    return replace(
        features,
        vessel_mask=vessel_mask.astype(np.bool_),
        vessel_distance_um=np.asarray(distance_transform_edt(~vessel_mask), dtype=np.float64)
        * features.dx_um,
    )


def _synthetic_retrospective() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for index in range(12):
        rows.append(
            {
                "programme_id": f"SYNTHETIC_{index:02d}",
                "target": f"T{index % 3}",
                "molecule": f"M{index}",
                "format": "IgG" if index % 2 == 0 else "ADC",
                "payload": "none" if index % 2 == 0 else "synthetic_payload",
                "dose": f"synthetic_{index % 3}",
                "indication": "A" if index % 4 < 2 else "B",
                "line_of_therapy": "synthetic",
                "target_expression_evidence": "SYNTHETIC_SCHEMA_EXAMPLE",
                "outcome": "success" if index in {0, 3, 4, 7, 9, 10} else "failure",
                "reason_for_discontinuation": "SYNTHETIC_SCHEMA_EXAMPLE",
                "source_identifiers": "SYNTHETIC_SCHEMA_EXAMPLE",
                "expression_reach_gap": float(0.2 + 0.04 * ((index * 7) % 5)),
            }
        )
    return pd.DataFrame(rows)


def run_benchmark(*, quick: bool = False, seed: int = 17) -> dict[str, Any]:
    """Run all preregistered simulation ablations and theory checks."""

    tracemalloc.start()
    started = time.perf_counter()
    size = 30 if quick else 42
    cell_count = 500 if quick else 900
    geometry = simulate_geometry(
        GeometryConfig(
            size=size,
            cell_count=cell_count,
            vessel_count=5,
            stroma_level=0.65,
            antigen_level=0.80,
            seed=seed,
        )
    )
    features = extract_features(geometry, antigen_calibration_nM_per_signal=300.0)
    base = ModelParameters()
    truth = solve_transport_robust(features, base)
    if not truth.converged:
        raise RuntimeError("Reference simulation failed to converge")

    rng_observation = np.random.default_rng(seed + 100)
    observed_vessels = binary_dilation(features.vessel_mask, iterations=1)
    observed = replace(
        features,
        vessel_mask=observed_vessels,
        vessel_distance_um=np.asarray(distance_transform_edt(~observed_vessels), dtype=np.float64)
        * features.dx_um,
        ecm=np.clip(
            gaussian_filter(features.ecm, 0.8)
            + rng_observation.normal(0.0, 0.03, size=features.ecm.shape),
            0.0,
            1.0,
        ),
        caf=np.clip(
            gaussian_filter(features.caf, 0.8)
            + rng_observation.normal(0.0, 0.03, size=features.caf.shape),
            0.0,
            1.0,
        ),
        antigen_nM=np.clip(
            gaussian_filter(features.antigen_nM, 0.8)
            * (1.0 + rng_observation.normal(0.0, 0.08, size=features.antigen_nM.shape)),
            0.0,
            None,
        ),
    )
    mechanistic_solution = solve_transport_robust(observed, base)
    if not mechanistic_solution.converged:
        raise RuntimeError("Observed-feature mechanistic solve failed to converge")
    true_bound = truth.bound_fraction[features.tumour_mask]
    distance = observed.vessel_distance_um[features.tumour_mask]
    ecm = observed.ecm[features.tumour_mask]
    caf = observed.caf[features.tumour_mask]
    antigen = observed.antigen_nM[features.tumour_mask]
    distance_only = np.exp(-distance / 80.0)
    naive_linear = 1.8 - 0.015 * distance - 1.1 * ecm - 0.7 * caf - 0.002 * antigen
    naive = 1.0 / (1.0 + np.exp(-naive_linear))
    mechanistic = mechanistic_solution.bound_fraction[features.tumour_mask]
    comparison = {
        "mechanistic": {
            "rmse": _rmse(true_bound, mechanistic),
            "classification_accuracy": _classification_accuracy(
                true_bound, mechanistic, base.engagement_threshold
            ),
        },
        "naive_weighted_sum": {
            "rmse": _rmse(true_bound, naive),
            "classification_accuracy": _classification_accuracy(
                true_bound, naive, base.engagement_threshold
            ),
            "label": "EXPLICITLY_NAIVE_UNFITTED_BASELINE",
        },
        "distance_only": {
            "rmse": _rmse(true_bound, distance_only),
            "classification_accuracy": _classification_accuracy(
                true_bound, distance_only, base.engagement_threshold
            ),
        },
        "primary_metric": "rmse",
        "physics_buys_lower_rmse": bool(
            _rmse(true_bound, mechanistic) < _rmse(true_bound, naive)
            and _rmse(true_bound, mechanistic) < _rmse(true_bound, distance_only)
        ),
        "classification_winner": max(
            {
                "mechanistic": _classification_accuracy(
                    true_bound, mechanistic, base.engagement_threshold
                ),
                "naive_weighted_sum": _classification_accuracy(
                    true_bound, naive, base.engagement_threshold
                ),
                "distance_only": _classification_accuracy(
                    true_bound, distance_only, base.engagement_threshold
                ),
            },
            key=lambda name: {
                "mechanistic": _classification_accuracy(
                    true_bound, mechanistic, base.engagement_threshold
                ),
                "naive_weighted_sum": _classification_accuracy(
                    true_bound, naive, base.engagement_threshold
                ),
                "distance_only": _classification_accuracy(
                    true_bound, distance_only, base.engagement_threshold
                ),
            }[name],
        ),
        "qualification": (
            "Ground truth and prediction share the same equation family, but prediction uses "
            "blurred, noisy and morphologically perturbed observed features. This validates "
            "implementation robustness, not biology."
        ),
    }

    format_diffusion = {"IgG": 6.0, "Fab": 20.0, "scFv": 60.0}
    size_depths: dict[str, float] = {}
    size_series_parameters = replace(
        base,
        vessel_concentration_nM=8.0,
        internalisation_s=2.0e-4,
    )
    for name, diffusion in format_diffusion.items():
        result = solve_transport_robust(
            features,
            replace(size_series_parameters, diffusion_um2_s=diffusion),
        )
        size_depths[name] = _penetration_depth(
            features, result.bound_fraction, base.engagement_threshold
        )

    antigen_low = replace(features, antigen_nM=features.antigen_nM * 0.25)
    antigen_high = replace(features, antigen_nM=features.antigen_nM * 2.0)
    low_antigen_solution = solve_transport_robust(antigen_low, base)
    high_antigen_solution = solve_transport_robust(antigen_high, base)
    antigen_reach = {
        "low": _reachable_fraction(
            antigen_low, low_antigen_solution.bound_fraction, base.engagement_threshold
        ),
        "high": _reachable_fraction(
            antigen_high, high_antigen_solution.bound_fraction, base.engagement_threshold
        ),
    }

    dose_reach: dict[str, float] = {}
    for dose in (10.0, 20.0, 100.0):
        solution = solve_transport_robust(features, replace(base, vessel_concentration_nM=dose))
        dose_reach[str(int(dose))] = _reachable_fraction(
            features, solution.bound_fraction, base.engagement_threshold
        )

    internalisation_reach: dict[str, float] = {}
    for rate in (1.0e-5, 3.0e-4):
        solution = solve_transport_robust(features, replace(base, internalisation_s=rate))
        internalisation_reach[f"{rate:.1e}"] = _reachable_fraction(
            features, solution.bound_fraction, base.engagement_threshold
        )

    low_barrier_features = replace(
        features,
        ecm=np.clip(features.ecm * 0.25, 0.0, 1.0),
        caf=np.clip(features.caf * 0.25, 0.0, 1.0),
        antigen_nM=features.antigen_nM * 0.25,
    )
    high_barrier_features = replace(
        features,
        ecm=np.clip(features.ecm * 1.35, 0.0, 1.0),
        caf=np.clip(features.caf * 1.35, 0.0, 1.0),
        antigen_nM=features.antigen_nM * 1.75,
    )
    low_barrier_solution = solve_transport_robust(
        low_barrier_features, replace(base, internalisation_s=1.0e-5)
    )
    high_barrier_solution = solve_transport_robust(
        high_barrier_features, replace(base, internalisation_s=3.0e-4)
    )
    target_fraction = float(np.mean(features.cell_target_positive[features.cell_is_tumour]))
    low_barrier_gap = target_fraction * (
        1.0
        - _reachable_fraction(
            low_barrier_features,
            low_barrier_solution.bound_fraction,
            base.engagement_threshold,
        )
    )
    high_barrier_gap = target_fraction * (
        1.0
        - _reachable_fraction(
            high_barrier_features,
            high_barrier_solution.bound_fraction,
            base.engagement_threshold,
        )
    )

    uncertainty_scenarios = {
        "low_barrier": low_barrier_features,
        "reference": features,
        "high_barrier": high_barrier_features,
    }
    uncertainty_results: dict[str, dict[str, object]] = {}
    for scenario_index, (scenario_name, scenario_features) in enumerate(
        uncertainty_scenarios.items()
    ):
        index_output, _ = compute_index(
            scenario_features,
            base,
            draws=8 if quick else 20,
            seed=seed + 20 + scenario_index,
        )
        uncertainty_results[scenario_name] = {
            "status": index_output.status,
            "decision_stability": index_output.decision_stability,
            "reachable_fraction_interval": index_output.reachable_fraction.model_dump(),
            "abstention_reasons": index_output.abstention_reasons,
        }
    abstention_share = float(
        np.mean(
            [result["status"] == "INSUFFICIENT_EVIDENCE" for result in uncertainty_results.values()]
        )
    )

    sensitivity_geometry = simulate_geometry(
        GeometryConfig(
            size=24 if quick else 28,
            cell_count=350,
            vessel_count=4,
            stroma_level=0.65,
            antigen_level=0.80,
            seed=seed + 2,
        )
    )
    sensitivity_features = extract_features(
        sensitivity_geometry, antigen_calibration_nM_per_signal=300.0
    )
    sensitivity = sobol_sensitivity(
        sensitivity_features, base, sample_power=3 if quick else 4, seed=seed
    )

    rng = np.random.default_rng(seed)
    target = features.cell_is_tumour & features.cell_target_positive
    jitter_reach: list[float] = []
    for _ in range(25):
        jitter_rows = np.clip(
            features.cell_rows + rng.integers(-1, 2, size=len(features.cell_rows)),
            0,
            features.vessel_mask.shape[0] - 1,
        )
        jitter_cols = np.clip(
            features.cell_cols + rng.integers(-1, 2, size=len(features.cell_cols)),
            0,
            features.vessel_mask.shape[1] - 1,
        )
        values = truth.bound_fraction[jitter_rows[target], jitter_cols[target]]
        jitter_reach.append(float(np.mean(values >= base.engagement_threshold)))
    vessel_sensitivity: dict[str, float] = {}
    for label, mask in {
        "eroded": binary_erosion(features.vessel_mask),
        "original": features.vessel_mask,
        "dilated": binary_dilation(features.vessel_mask),
    }.items():
        if not np.any(mask):
            vessel_sensitivity[label] = float("nan")
            continue
        modified = _replace_vessels(features, mask)
        solution = solve_transport_robust(modified, base)
        vessel_sensitivity[label] = _reachable_fraction(
            modified, solution.bound_fraction, base.engagement_threshold
        )

    synthetic_table = _synthetic_retrospective()
    retrospective = run_retrospective(
        synthetic_table, permutations=200 if quick else 1000, seed=seed, allow_synthetic=True
    )

    elapsed = time.perf_counter() - started
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    theory_checks = {
        "size_series_increasing_penetration": (
            size_depths["IgG"] < size_depths["Fab"] < size_depths["scFv"]
        ),
        "higher_antigen_reduces_reach": antigen_reach["high"] < antigen_reach["low"],
        "higher_dose_rescues_reach": dose_reach["10"] < dose_reach["20"] < dose_reach["100"],
        "higher_internalisation_reduces_reach": internalisation_reach["3.0e-04"]
        <= internalisation_reach["1.0e-05"],
        "high_barrier_scenario_increases_gap": high_barrier_gap > low_barrier_gap,
    }
    status = "PASS" if all(theory_checks.values()) else "FAIL"
    indices = sensitivity.get("indices", {})
    top_total = None
    if isinstance(indices, dict) and indices:
        top_total = max(
            indices,
            key=lambda key: float(indices[key]["total_order"]),
        )

    return {
        "status": status,
        "seed": seed,
        "simulation_only": True,
        "model_comparison": comparison,
        "global_sensitivity": sensitivity,
        "top_total_sensitivity_parameter": top_total,
        "size_series_penetration_depth_um": size_depths,
        "dose_sweep_reachable_fraction": dose_reach,
        "antigen_density_reachable_fraction": antigen_reach,
        "internalisation_reachable_fraction": internalisation_reach,
        "barrier_scenario_expression_reach_gap": {
            "low_barrier": low_barrier_gap,
            "high_barrier": high_barrier_gap,
        },
        "uncertainty_abstention": {
            "scenarios": uncertainty_results,
            "abstention_share": abstention_share,
            "substantial_threshold": 1.0 / 3.0,
            "prediction_supported": abstention_share >= 1.0 / 3.0,
        },
        "theory_checks": theory_checks,
        "segmentation_robustness": {
            "cell_jitter_reachable_fraction_5_95": [
                float(np.quantile(jitter_reach, 0.05)),
                float(np.quantile(jitter_reach, 0.95)),
            ],
            "vessel_mask_sensitivity": vessel_sensitivity,
        },
        "retrospective": {
            "status": "NOT_COMPUTED_REAL_DATA",
            "synthetic_schema_test": {
                "status": retrospective.status,
                "n_programmes": retrospective.n_programmes,
                "median_difference": retrospective.median_difference,
                "permutation_p": retrospective.permutation_p,
                "stratified_difference": retrospective.stratified_difference,
                "label": "SYNTHETIC_SCHEMA_EXAMPLE_NOT_SCIENTIFIC_EVIDENCE",
            },
            "distance_matched_null": "NOT_COMPUTED_REAL_DATA",
        },
        "preregistered_predictions": {
            "high_density_fast_internalisation_stroma_gap": (
                "SUPPORTED_IN_SIMULATION"
                if high_barrier_gap > low_barrier_gap
                else "NOT_SUPPORTED_IN_SIMULATION"
            ),
            "uncertainty_forces_substantial_abstention": (
                "SUPPORTED_IN_SIMULATION"
                if abstention_share >= 1.0 / 3.0
                else "NOT_SUPPORTED_IN_SIMULATION"
            ),
            "clinical_retrospective_null_or_weak": "NOT_COMPUTED_REAL_DATA",
        },
        "runtime_seconds": elapsed,
        "peak_memory_mb": peak / (1024.0 * 1024.0),
    }


def write_benchmark(result: dict[str, Any], output: Path) -> Path:
    """Write a benchmark result document."""

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return output
