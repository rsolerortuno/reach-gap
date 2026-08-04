"""Finite-volume steady reaction-diffusion solver."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix, diags, lil_matrix
from scipy.sparse.linalg import spsolve

from reach_gap.schemas import FloatArray, ModelParameters, SolverResult, SpatialFeatures


def effective_diffusion(features: SpatialFeatures, parameters: ModelParameters) -> FloatArray:
    """Compute positive scalar diffusivity from matrix and fibroblast scores."""

    attenuation = np.exp(-parameters.beta_ecm * features.ecm - parameters.beta_caf * features.caf)
    return np.asarray(parameters.diffusion_um2_s * attenuation, dtype=np.float64)


def _harmonic(left: float, right: float) -> float:
    if left <= 0.0 or right <= 0.0:
        raise ValueError("Diffusivity must be positive")
    return 2.0 * left * right / (left + right)


def _assemble_diffusion(
    diffusion: FloatArray,
    vessel_mask: NDArray[np.bool_],
    dx_um: float,
    vessel_concentration_nM: float,
) -> tuple[csr_matrix, FloatArray]:
    rows, cols = diffusion.shape
    size = rows * cols
    matrix = lil_matrix((size, size), dtype=np.float64)
    rhs = np.zeros(size, dtype=np.float64)
    inv_dx2 = 1.0 / (dx_um * dx_um)
    for row in range(rows):
        for col in range(cols):
            index = row * cols + col
            if vessel_mask[row, col]:
                matrix[index, index] = 1.0
                rhs[index] = vessel_concentration_nM
                continue
            diagonal = 0.0
            for d_row, d_col in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                n_row = row + d_row
                n_col = col + d_col
                if n_row < 0 or n_row >= rows or n_col < 0 or n_col >= cols:
                    continue
                weight = (
                    _harmonic(float(diffusion[row, col]), float(diffusion[n_row, n_col])) * inv_dx2
                )
                diagonal += weight
                matrix[index, n_row * cols + n_col] = -weight
            matrix[index, index] = diagonal
    return matrix.tocsr(), rhs


def solve_transport(
    features: SpatialFeatures,
    parameters: ModelParameters,
    *,
    max_iterations: int = 80,
    tolerance: float = 1.0e-7,
    relaxation: float = 0.75,
) -> SolverResult:
    """Solve the nonlinear steady transport equation by Picard iteration."""

    if not np.any(features.vessel_mask):
        raise ValueError("No vascular structures were identified")
    if parameters.kd_nM <= 0 or parameters.vessel_concentration_nM < 0:
        raise ValueError("Invalid concentration parameter")
    if not 0.0 < relaxation <= 1.0:
        raise ValueError("Relaxation must lie in (0, 1]")

    diffusion = effective_diffusion(features, parameters)
    characteristic_sink = parameters.clearance_s + (
        parameters.internalisation_s
        * parameters.antigen_calibration_factor
        * np.maximum(features.antigen_nM, 0.0)
        / (parameters.kd_nM + max(parameters.vessel_concentration_nM, 1.0e-12))
    )
    median_sink = max(float(np.median(characteristic_sink)), 1.0e-12)
    length = np.sqrt(float(np.median(diffusion)) / median_sink)
    initial = parameters.vessel_concentration_nM * np.exp(
        -features.vessel_distance_um / max(length, features.dx_um)
    )
    initial[features.vessel_mask] = parameters.vessel_concentration_nM
    concentration = initial
    residual = float("inf")
    diffusion_matrix, rhs = _assemble_diffusion(
        diffusion,
        features.vessel_mask,
        features.dx_um,
        parameters.vessel_concentration_nM,
    )
    vessel_flat = features.vessel_mask.ravel()

    for iteration in range(1, max_iterations + 1):
        sink = parameters.clearance_s + (
            parameters.internalisation_s
            * parameters.antigen_calibration_factor
            * np.maximum(features.antigen_nM, 0.0)
            / (parameters.kd_nM + np.maximum(concentration, 0.0))
        )
        sink_flat = np.asarray(sink.ravel(), dtype=np.float64)
        sink_flat[vessel_flat] = 0.0
        matrix = diffusion_matrix + diags(sink_flat, offsets=0, format="csr")
        solved = np.asarray(spsolve(matrix, rhs), dtype=np.float64).reshape(
            features.vessel_mask.shape
        )
        solved = np.clip(solved, 0.0, parameters.vessel_concentration_nM)
        updated = relaxation * solved + (1.0 - relaxation) * concentration
        denominator = max(float(np.max(np.abs(updated))), 1.0e-12)
        residual = float(np.max(np.abs(updated - concentration)) / denominator)
        concentration = updated
        if residual < tolerance:
            bound = concentration / (parameters.kd_nM + concentration)
            return SolverResult(
                concentration_nM=concentration,
                bound_fraction=bound,
                effective_diffusion_um2_s=diffusion,
                iterations=iteration,
                converged=True,
                residual=residual,
            )

    bound = concentration / (parameters.kd_nM + concentration)
    return SolverResult(
        concentration_nM=concentration,
        bound_fraction=bound,
        effective_diffusion_um2_s=diffusion,
        iterations=max_iterations,
        converged=False,
        residual=residual,
    )


def solve_transport_robust(
    features: SpatialFeatures,
    parameters: ModelParameters,
) -> SolverResult:
    """Retry the same equation with three damping schedules before declaring failure."""

    attempts = ((100, 0.75), (200, 0.5), (360, 0.25))
    last: SolverResult | None = None
    for max_iterations, relaxation in attempts:
        last = solve_transport(
            features,
            parameters,
            max_iterations=max_iterations,
            relaxation=relaxation,
        )
        if last.converged:
            return last
    if last is None:
        raise RuntimeError("No solver attempt was executed")
    return last


def analytical_slab_profile(
    x_um: FloatArray,
    length_um: float,
    diffusion_um2_s: float,
    sink_s: float,
    boundary_nM: float,
) -> FloatArray:
    """Analytical 1D profile with a Dirichlet source and distal no-flux boundary."""

    if diffusion_um2_s <= 0 or sink_s < 0 or length_um <= 0:
        raise ValueError("Invalid analytical slab parameter")
    if sink_s == 0:
        return np.full_like(x_um, boundary_nM, dtype=np.float64)
    lam = np.sqrt(sink_s / diffusion_um2_s)
    return np.asarray(
        boundary_nM * np.cosh(lam * (length_um - x_um)) / np.cosh(lam * length_um),
        dtype=np.float64,
    )


def save_solution(result: SolverResult, path: Path) -> None:
    """Save solver fields and convergence metadata."""

    np.savez_compressed(
        path,
        concentration_nM=result.concentration_nM,
        bound_fraction=result.bound_fraction,
        effective_diffusion_um2_s=result.effective_diffusion_um2_s,
        iterations=np.array(result.iterations),
        converged=np.array(result.converged),
        residual=np.array(result.residual),
    )


def load_solution(path: Path) -> SolverResult:
    """Load solver fields from disk."""

    with np.load(path, allow_pickle=False) as data:
        return SolverResult(
            concentration_nM=data["concentration_nM"].astype(np.float64),
            bound_fraction=data["bound_fraction"].astype(np.float64),
            effective_diffusion_um2_s=data["effective_diffusion_um2_s"].astype(np.float64),
            iterations=int(data["iterations"]),
            converged=bool(data["converged"]),
            residual=float(data["residual"]),
        )
