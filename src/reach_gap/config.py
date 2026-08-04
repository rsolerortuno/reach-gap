"""Default scenario ranges and parameter sampling."""

from __future__ import annotations

from dataclasses import fields

import numpy as np

from reach_gap.schemas import ModelParameters, ParameterRange

DEFAULT_RANGES: dict[str, ParameterRange] = {
    "diffusion_um2_s": ParameterRange(5.0, 30.0, "loguniform"),
    "kd_nM": ParameterRange(0.1, 10.0, "loguniform"),
    "antigen_calibration_factor": ParameterRange(0.25, 4.0, "loguniform"),
    "internalisation_s": ParameterRange(1.0e-5, 5.0e-4, "loguniform"),
    "clearance_s": ParameterRange(1.0e-8, 2.0e-5, "loguniform"),
    "beta_ecm": ParameterRange(0.2, 2.0, "uniform"),
    "beta_caf": ParameterRange(0.0, 1.5, "uniform"),
}


def sample_parameters(
    rng: np.random.Generator,
    base: ModelParameters,
    ranges: dict[str, ParameterRange] | None = None,
) -> ModelParameters:
    """Sample one parameter set while preserving operational settings."""

    active = DEFAULT_RANGES if ranges is None else ranges
    values = base.as_dict()
    valid_names = {field.name for field in fields(ModelParameters)}
    for name, bounds in active.items():
        if name not in valid_names:
            raise KeyError(f"Unknown model parameter: {name}")
        if bounds.low > bounds.high:
            raise ValueError(f"Invalid range for {name}")
        if bounds.distribution == "loguniform":
            if bounds.low <= 0:
                raise ValueError(f"Log-uniform range for {name} must be positive")
            draw = float(np.exp(rng.uniform(np.log(bounds.low), np.log(bounds.high))))
        else:
            draw = float(rng.uniform(bounds.low, bounds.high))
        values[name] = draw
    return ModelParameters(**values)


def serialise_ranges(
    ranges: dict[str, ParameterRange] | None = None,
) -> dict[str, dict[str, object]]:
    """Return ranges as JSON-safe nested mappings."""

    active = DEFAULT_RANGES if ranges is None else ranges
    return {
        name: {"low": value.low, "high": value.high, "distribution": value.distribution}
        for name, value in active.items()
    }
