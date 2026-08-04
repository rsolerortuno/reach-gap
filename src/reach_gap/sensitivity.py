"""Global variance-based parameter sensitivity analysis."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.stats import qmc

from reach_gap.config import DEFAULT_RANGES
from reach_gap.schemas import ModelParameters, ParameterRange, SpatialFeatures
from reach_gap.solver import solve_transport_robust


def _scale_unit(value: float, bounds: ParameterRange) -> float:
    if bounds.distribution == "loguniform":
        return float(
            np.exp(np.log(bounds.low) + value * (np.log(bounds.high) - np.log(bounds.low)))
        )
    return float(bounds.low + value * (bounds.high - bounds.low))


def _evaluate_gap(
    vector: NDArray[np.float64],
    names: list[str],
    ranges: dict[str, ParameterRange],
    base: ModelParameters,
    features: SpatialFeatures,
) -> float:
    values = base.as_dict()
    for name, unit_value in zip(names, vector, strict=True):
        values[name] = _scale_unit(float(unit_value), ranges[name])
    parameters = ModelParameters(**values)
    solution = solve_transport_robust(features, parameters)
    if not solution.converged:
        return float("nan")
    tumour = features.cell_is_tumour
    target = tumour & features.cell_target_positive
    target_count = int(np.sum(target))
    tumour_count = int(np.sum(tumour))
    if target_count == 0 or tumour_count == 0:
        return 0.0
    bound = solution.bound_fraction[features.cell_rows[target], features.cell_cols[target]]
    reachable = float(np.mean(bound >= parameters.engagement_threshold))
    target_fraction = target_count / tumour_count
    return target_fraction * (1.0 - reachable)


def sobol_sensitivity(
    features: SpatialFeatures,
    base: ModelParameters,
    *,
    ranges: dict[str, ParameterRange] | None = None,
    sample_power: int = 4,
    seed: int = 17,
) -> dict[str, object]:
    """Estimate first-order and total Sobol indices using hybrid quasi-random matrices."""

    active = DEFAULT_RANGES if ranges is None else ranges
    names = list(active)
    dimension = len(names)
    if sample_power < 2:
        raise ValueError("sample_power must be at least 2")
    sampler_a = qmc.Sobol(d=dimension, scramble=True, seed=seed)
    sampler_b = qmc.Sobol(d=dimension, scramble=True, seed=seed + 1)
    matrix_a = sampler_a.random_base2(m=sample_power)
    matrix_b = sampler_b.random_base2(m=sample_power)
    output_a = np.array(
        [_evaluate_gap(row, names, active, base, features) for row in matrix_a], dtype=np.float64
    )
    output_b = np.array(
        [_evaluate_gap(row, names, active, base, features) for row in matrix_b], dtype=np.float64
    )
    if np.any(~np.isfinite(output_a)) or np.any(~np.isfinite(output_b)):
        raise RuntimeError("Sensitivity solve did not converge")
    variance = float(np.var(np.concatenate([output_a, output_b]), ddof=1))
    if variance <= 1.0e-15:
        return {
            "status": "DEGENERATE_OUTPUT",
            "sample_count": len(output_a),
            "variance": variance,
            "indices": {name: {"first_order": 0.0, "total_order": 0.0} for name in names},
        }
    indices: dict[str, dict[str, float]] = {}
    for column, name in enumerate(names):
        hybrid = matrix_a.copy()
        hybrid[:, column] = matrix_b[:, column]
        output_hybrid = np.array(
            [_evaluate_gap(row, names, active, base, features) for row in hybrid],
            dtype=np.float64,
        )
        first = float(np.mean(output_b * (output_hybrid - output_a)) / variance)
        total = float(0.5 * np.mean((output_a - output_hybrid) ** 2) / variance)
        indices[name] = {
            "first_order": float(np.clip(first, -1.0, 1.0)),
            "total_order": float(np.clip(total, 0.0, 1.5)),
        }
    return {
        "status": "OK",
        "sample_count": len(output_a),
        "variance": variance,
        "indices": indices,
        "seed": seed,
        "method": "Saltelli-style hybrid matrices from scrambled Sobol sequences",
    }
