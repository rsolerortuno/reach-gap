from __future__ import annotations

from dataclasses import replace

import numpy as np

from reach_gap.features import extract_features
from reach_gap.geometry import GeometryConfig, simulate_geometry
from reach_gap.schemas import ModelParameters, SpatialFeatures
from reach_gap.solver import analytical_slab_profile, effective_diffusion, solve_transport


def _uniform_slab() -> SpatialFeatures:
    rows, cols = 8, 31
    vessel = np.zeros((rows, cols), dtype=np.bool_)
    vessel[:, 0] = True
    zeros = np.zeros((rows, cols), dtype=np.float64)
    tumour = np.ones((rows, cols), dtype=np.bool_)
    cell_rows, cell_cols = np.indices((rows, cols))
    return SpatialFeatures(
        vessel_mask=vessel,
        vessel_distance_um=np.tile(np.arange(cols, dtype=np.float64) * 10.0, (rows, 1)),
        ecm=zeros,
        caf=zeros,
        antigen_nM=zeros,
        tumour_mask=tumour,
        cell_rows=cell_rows.ravel().astype(np.int64),
        cell_cols=cell_cols.ravel().astype(np.int64),
        cell_is_tumour=np.ones(rows * cols, dtype=np.bool_),
        cell_target_positive=np.ones(rows * cols, dtype=np.bool_),
        dx_um=10.0,
        antigen_calibrated=True,
        seed=1,
    )


def test_solver_matches_linear_analytical_slab() -> None:
    features = _uniform_slab()
    parameters = ModelParameters(
        diffusion_um2_s=12.0,
        clearance_s=1.0e-4,
        internalisation_s=0.0,
        vessel_concentration_nM=100.0,
    )
    result = solve_transport(features, parameters)
    assert result.converged
    x = np.arange(features.shape[1], dtype=np.float64) * features.dx_um
    expected = analytical_slab_profile(
        x, x[-1], parameters.diffusion_um2_s, parameters.clearance_s, 100.0
    )
    assert np.max(np.abs(result.concentration_nM[3] - expected)) / 100.0 < 0.015


def test_zero_sink_is_constant() -> None:
    features = _uniform_slab()
    parameters = ModelParameters(clearance_s=0.0, internalisation_s=0.0)
    result = solve_transport(features, parameters)
    assert result.converged
    assert np.allclose(result.concentration_nM, parameters.vessel_concentration_nM, atol=1.0e-5)


def test_effective_diffusion_decreases_with_stroma() -> None:
    geometry = simulate_geometry(GeometryConfig(size=20, cell_count=100, seed=5))
    features = extract_features(geometry)
    parameters = ModelParameters()
    diffusion = effective_diffusion(features, parameters)
    assert np.all(diffusion > 0.0)
    assert float(diffusion.max()) <= parameters.diffusion_um2_s


def test_theory_monotonicity() -> None:
    geometry = simulate_geometry(GeometryConfig(size=24, cell_count=200, seed=8))
    features = extract_features(geometry, antigen_calibration_nM_per_signal=250.0)
    base = ModelParameters()
    igg = solve_transport(features, replace(base, diffusion_um2_s=8.0))
    scfv = solve_transport(features, replace(base, diffusion_um2_s=60.0))
    target = features.cell_target_positive
    rows, cols = features.cell_rows[target], features.cell_cols[target]
    assert np.mean(scfv.bound_fraction[rows, cols]) >= np.mean(igg.bound_fraction[rows, cols])

    low_antigen = replace(features, antigen_nM=features.antigen_nM * 0.2)
    high_antigen = replace(features, antigen_nM=features.antigen_nM * 2.0)
    low = solve_transport(low_antigen, base)
    high = solve_transport(high_antigen, base)
    assert np.mean(low.bound_fraction[rows, cols]) >= np.mean(high.bound_fraction[rows, cols])

    low_dose = solve_transport(features, replace(base, vessel_concentration_nM=20.0))
    high_dose = solve_transport(features, replace(base, vessel_concentration_nM=500.0))
    assert np.mean(high_dose.bound_fraction[rows, cols]) >= np.mean(
        low_dose.bound_fraction[rows, cols]
    )
