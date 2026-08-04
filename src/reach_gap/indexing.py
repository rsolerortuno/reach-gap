"""Index computation, uncertainty propagation, barrier attribution and claims."""

from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any, Literal, cast

import numpy as np
import scipy

from reach_gap.config import DEFAULT_RANGES, sample_parameters, serialise_ranges
from reach_gap.schemas import (
    ClaimsDocument,
    FloatArray,
    IndexOutput,
    Interval,
    ModelParameters,
    ParameterRange,
    SpatialFeatures,
)
from reach_gap.solver import solve_transport_robust


def _interval(values: list[float]) -> Interval:
    array = np.asarray(values, dtype=np.float64)
    return Interval(
        median=float(np.median(array)),
        lower=float(np.quantile(array, 0.05)),
        upper=float(np.quantile(array, 0.95)),
    )


def _cell_values(field: FloatArray, features: SpatialFeatures) -> FloatArray:
    return np.asarray(field[features.cell_rows, features.cell_cols], dtype=np.float64)


def _barrier_summary(
    feature_draws: list[tuple[float, float, float]],
) -> tuple[str, dict[str, float]]:
    names = ("vascular_distance", "matrix_fibroblast", "binding_site")
    if not feature_draws:
        return "BARRIER_INDISTINGUISHABLE", {name: 0.0 for name in names}
    weights = np.asarray(feature_draws, dtype=np.float64)
    row_sums = weights.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    normalised = weights / row_sums
    means = normalised.mean(axis=0)
    winners = np.argmax(normalised, axis=1)
    winner_counts = np.bincount(winners, minlength=3) / float(len(winners))
    order = np.argsort(means)[::-1]
    top = int(order[0])
    runner_up = int(order[1])
    identifiable = winner_counts[top] >= 0.60 and means[top] - means[runner_up] >= 0.15
    dominant = names[top] if identifiable else "BARRIER_INDISTINGUISHABLE"
    return dominant, {name: float(value) for name, value in zip(names, means, strict=True)}


def compute_index(
    features: SpatialFeatures,
    base_parameters: ModelParameters,
    *,
    ranges: dict[str, ParameterRange] | None = None,
    draws: int = 32,
    seed: int = 17,
    reachable_decision_threshold: float = 0.5,
    maximum_valid_vessel_distance_um: float = 600.0,
) -> tuple[IndexOutput, ClaimsDocument]:
    """Propagate parameter uncertainty through the full solver and compute the reach gap."""

    if draws < 4:
        raise ValueError("At least four uncertainty draws are required")
    active_ranges = DEFAULT_RANGES if ranges is None else ranges
    reasons: list[str] = []
    if not np.any(features.vessel_mask):
        reasons.append("NO_VASCULAR_STRUCTURES")
    if not features.antigen_calibrated:
        reasons.append("ANTIGEN_DENSITY_UNCALIBRATED")
    max_distance = float(np.max(features.vessel_distance_um[features.tumour_mask]))
    if max_distance > maximum_valid_vessel_distance_um:
        reasons.append("GEOMETRY_OUTSIDE_VALIDATED_REGIME")

    tumour_cells = features.cell_is_tumour
    target_cells = tumour_cells & features.cell_target_positive
    target_count = int(np.sum(target_cells))
    tumour_count = int(np.sum(tumour_cells))
    target_fraction = 0.0 if tumour_count == 0 else target_count / tumour_count
    if target_count == 0:
        reasons.append("NO_TARGET_POSITIVE_CELLS")

    reachable_values: list[float] = []
    population_values: list[float] = []
    gap_values: list[float] = []
    penetration_values: list[float] = []
    barrier_draws: list[tuple[float, float, float]] = []
    rng = np.random.default_rng(seed)

    if not reasons:
        for _ in range(draws):
            parameters = sample_parameters(rng, base_parameters, active_ranges)
            solution = solve_transport_robust(features, parameters)
            if not solution.converged:
                reasons.append("SOLVER_DID_NOT_CONVERGE")
                break
            bound_cells = _cell_values(solution.bound_fraction, features)
            target_bound = bound_cells[target_cells]
            reachable_mask = target_bound >= parameters.engagement_threshold
            reachable = float(np.mean(reachable_mask))
            population_reachable = target_fraction * reachable
            gap = target_fraction * (1.0 - reachable)
            reachable_values.append(reachable)
            population_values.append(population_reachable)
            gap_values.append(gap)
            if np.any(reachable_mask):
                target_distances = _cell_values(features.vessel_distance_um, features)[target_cells]
                penetration_values.append(float(np.max(target_distances[reachable_mask])))
            else:
                penetration_values.append(0.0)

            target_rows = features.cell_rows[target_cells]
            target_cols = features.cell_cols[target_cells]
            median_diffusion = float(
                np.median(solution.effective_diffusion_um2_s[target_rows, target_cols])
            )
            median_antigen = (
                float(np.median(features.antigen_nM[target_rows, target_cols]))
                * parameters.antigen_calibration_factor
            )
            local_sink = parameters.clearance_s + (
                parameters.internalisation_s
                * median_antigen
                / (parameters.kd_nM + parameters.vessel_concentration_nM)
            )
            penetration_length = np.sqrt(median_diffusion / max(local_sink, 1.0e-12))
            median_distance = float(
                np.median(features.vessel_distance_um[target_rows, target_cols])
            )
            vascular_penalty = (median_distance / max(penetration_length, features.dx_um)) ** 2
            median_effective_diffusion = float(
                np.median(solution.effective_diffusion_um2_s[target_rows, target_cols])
            )
            matrix_penalty = np.log(
                max(parameters.diffusion_um2_s / max(median_effective_diffusion, 1.0e-12), 1.0)
            )
            binding_penalty = (
                parameters.internalisation_s
                * median_antigen
                * median_distance**2
                / (
                    max(median_effective_diffusion, 1.0e-12)
                    * (parameters.kd_nM + parameters.vessel_concentration_nM)
                )
            )
            barrier_draws.append((vascular_penalty, matrix_penalty, binding_penalty))

    status: Literal["OK", "INSUFFICIENT_EVIDENCE"] = (
        "OK" if not reasons else "INSUFFICIENT_EVIDENCE"
    )
    if not reachable_values:
        reachable_values = [0.0]
        population_values = [0.0]
        gap_values = [target_fraction]
        penetration_values = [0.0]
    reachable_interval = _interval(reachable_values)
    stable = (
        reachable_interval.upper < reachable_decision_threshold
        or reachable_interval.lower > reachable_decision_threshold
    )
    if not stable and "PARAMETER_INTERVAL_SPANS_DECISION_THRESHOLD" not in reasons:
        reasons.append("PARAMETER_INTERVAL_SPANS_DECISION_THRESHOLD")
        status = "INSUFFICIENT_EVIDENCE"
    dominant, barrier_weights = _barrier_summary(barrier_draws)

    output = IndexOutput(
        status=status,
        abstention_reasons=reasons,
        target_positive_fraction=target_fraction,
        reachable_fraction=reachable_interval,
        population_reachable_fraction=_interval(population_values),
        expression_reach_gap=_interval(gap_values),
        penetration_depth_um=_interval(penetration_values) if penetration_values else None,
        dominant_barrier=dominant,
        barrier_weights=barrier_weights,
        decision_stability=stable,
        engagement_threshold=base_parameters.engagement_threshold,
        reachable_decision_threshold=reachable_decision_threshold,
        parameter_draws=draws,
        seed=seed,
        parameters=base_parameters.as_dict(),
        parameter_ranges=serialise_ranges(active_ranges),
        library_versions={
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
    )
    claims = build_claims(output)
    return output, claims


def build_claims(output: IndexOutput) -> ClaimsDocument:
    """Derive a bounded claim set from exactly what one run computed."""

    permitted = [
        (
            "The reported values are outputs of a two-dimensional mechanistic simulation "
            "or calibrated spatial run."
        ),
        "The expression-reach gap is reported with propagated parameter uncertainty.",
    ]
    conditional = [
        (
            "Reachability is conditional on the stated vessel concentration, engagement "
            "threshold, geometry and parameter ranges."
        ),
        (
            "Barrier attribution is qualitative and is emitted only when uncertainty draws "
            "separate one mechanism."
        ),
    ]
    unsupported = [
        "Clinical efficacy or failure prediction.",
        "Whole-body pharmacokinetics.",
        "Three-dimensional tumour penetration.",
        "A claim that changing the dominant barrier will improve patient outcome.",
    ]
    if output.status != "OK":
        permitted.append("The run abstained rather than issuing a stable reachability conclusion.")
    return ClaimsDocument(
        permitted=permitted,
        conditional=conditional,
        unsupported=unsupported,
        abstention_reasons=output.abstention_reasons,
        interval_basis="Empirical 5th-95th percentiles across full-pipeline parameter draws.",
    )


def write_index_outputs(
    output: IndexOutput, claims: ClaimsDocument, output_dir: Path
) -> tuple[Path, Path]:
    """Write index and claims JSON documents."""

    output_dir.mkdir(parents=True, exist_ok=True)
    index_path = output_dir / "index.json"
    claims_path = output_dir / "claims.json"
    index_path.write_text(output.model_dump_json(indent=2), encoding="utf-8")
    claims_path.write_text(claims.model_dump_json(indent=2), encoding="utf-8")
    return index_path, claims_path


def output_as_dict(output: IndexOutput) -> dict[str, Any]:
    """Return an ordinary mapping for composite result documents."""

    return cast(dict[str, Any], json.loads(output.model_dump_json()))
