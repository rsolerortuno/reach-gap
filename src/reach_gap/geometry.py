"""Deterministic synthetic tumour geometry generation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter

from reach_gap.schemas import FloatArray, TissueGeometry


@dataclass(frozen=True)
class GeometryConfig:
    """Configuration for a synthetic two-dimensional tumour section."""

    size: int = 48
    dx_um: float = 10.0
    vessel_count: int = 7
    cell_count: int = 900
    stroma_level: float = 0.45
    antigen_level: float = 0.65
    antigen_perivascular_bias: float = 0.25
    seed: int = 17


def _normalise(field: FloatArray) -> FloatArray:
    minimum = float(field.min())
    maximum = float(field.max())
    if maximum <= minimum:
        return np.zeros_like(field)
    return (field - minimum) / (maximum - minimum)


def simulate_geometry(config: GeometryConfig) -> TissueGeometry:
    """Generate vessels, stroma, matrix, antigen and segmented cell centroids."""

    if config.size < 12:
        raise ValueError("Geometry size must be at least 12")
    if config.vessel_count < 1:
        raise ValueError("At least one vessel is required for simulation")
    rng = np.random.default_rng(config.seed)
    shape: tuple[int, int] = (config.size, config.size)
    yy, xx = np.indices(shape)

    vessel_mask = np.zeros(shape, dtype=np.bool_)
    margin = max(2, config.size // 12)
    centres = rng.integers(margin, config.size - margin, size=(config.vessel_count, 2))
    radii = rng.integers(1, max(2, config.size // 20 + 1), size=config.vessel_count)
    for (row, col), radius in zip(centres, radii, strict=True):
        vessel_mask |= (yy - row) ** 2 + (xx - col) ** 2 <= int(radius) ** 2

    raw_caf: FloatArray = np.asarray(
        gaussian_filter(rng.random(shape), sigma=max(1.0, config.size / 12.0)), dtype=np.float64
    )
    caf = np.clip(config.stroma_level * (0.35 + 0.9 * _normalise(raw_caf)), 0.0, 1.0)
    raw_ecm: FloatArray = np.asarray(
        gaussian_filter(0.65 * caf + 0.35 * rng.random(shape), sigma=2.0), dtype=np.float64
    )
    ecm = np.clip(config.stroma_level * (0.3 + _normalise(raw_ecm)), 0.0, 1.0)

    vessel_proximity: FloatArray = np.asarray(
        gaussian_filter(vessel_mask.astype(np.float64), sigma=max(1.5, config.size / 18.0)),
        dtype=np.float64,
    )
    vessel_proximity = _normalise(vessel_proximity)
    antigen_noise = _normalise(gaussian_filter(rng.random(shape), sigma=2.2))
    antigen = config.antigen_level * (
        0.55 + 0.45 * antigen_noise + config.antigen_perivascular_bias * vessel_proximity
    )
    antigen = np.clip(antigen, 0.0, 1.0)

    tumour_radius = config.size * 0.47
    tumour_mask = (yy - (config.size - 1) / 2.0) ** 2 + (
        xx - (config.size - 1) / 2.0
    ) ** 2 <= tumour_radius**2
    tumour_mask &= ~vessel_mask
    valid = np.argwhere(tumour_mask)
    replace = config.cell_count > len(valid)
    chosen = valid[rng.choice(len(valid), size=config.cell_count, replace=replace)]
    rows = chosen[:, 0].astype(np.int64)
    cols = chosen[:, 1].astype(np.int64)
    signal = np.clip(antigen[rows, cols] + rng.normal(0.0, 0.04, size=config.cell_count), 0.0, 1.0)

    return TissueGeometry(
        vessel_mask=vessel_mask,
        ecm=ecm.astype(np.float64),
        caf=caf.astype(np.float64),
        antigen=antigen.astype(np.float64),
        tumour_mask=tumour_mask,
        cell_rows=rows,
        cell_cols=cols,
        cell_is_tumour=np.ones(config.cell_count, dtype=np.bool_),
        cell_target_signal=signal.astype(np.float64),
        dx_um=config.dx_um,
        seed=config.seed,
    )
