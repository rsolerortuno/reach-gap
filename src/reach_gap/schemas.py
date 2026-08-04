"""Typed public schemas used across the package."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Literal, TypeAlias

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, ConfigDict, Field

FloatArray: TypeAlias = NDArray[np.float64]
BoolArray: TypeAlias = NDArray[np.bool_]
IntArray: TypeAlias = NDArray[np.int64]


@dataclass(frozen=True)
class ModelParameters:
    """Physical and operational parameters for one transport solve."""

    diffusion_um2_s: float = 12.0
    kd_nM: float = 5.0
    antigen_calibration_factor: float = 1.0
    internalisation_s: float = 1.0e-4
    clearance_s: float = 5.0e-6
    beta_ecm: float = 1.2
    beta_caf: float = 0.8
    vessel_concentration_nM: float = 20.0
    engagement_threshold: float = 0.5

    def as_dict(self) -> dict[str, float]:
        """Return a JSON-serialisable representation."""

        return asdict(self)


@dataclass(frozen=True)
class ParameterRange:
    """Closed uncertainty range for one scalar parameter."""

    low: float
    high: float
    distribution: Literal["uniform", "loguniform"] = "uniform"


@dataclass(frozen=True)
class TissueGeometry:
    """Raster fields and cell coordinates for a synthetic or ingested section."""

    vessel_mask: BoolArray
    ecm: FloatArray
    caf: FloatArray
    antigen: FloatArray
    tumour_mask: BoolArray
    cell_rows: IntArray
    cell_cols: IntArray
    cell_is_tumour: BoolArray
    cell_target_signal: FloatArray
    dx_um: float
    seed: int

    @property
    def shape(self) -> tuple[int, int]:
        """Return raster shape."""

        return self.antigen.shape


@dataclass(frozen=True)
class SpatialFeatures:
    """Feature arrays consumed by the mechanistic solver and index."""

    vessel_mask: BoolArray
    vessel_distance_um: FloatArray
    ecm: FloatArray
    caf: FloatArray
    antigen_nM: FloatArray
    tumour_mask: BoolArray
    cell_rows: IntArray
    cell_cols: IntArray
    cell_is_tumour: BoolArray
    cell_target_positive: BoolArray
    dx_um: float
    antigen_calibrated: bool
    seed: int

    @property
    def shape(self) -> tuple[int, int]:
        """Return raster shape."""

        return self.antigen_nM.shape


@dataclass(frozen=True)
class SolverResult:
    """Numerical fields and convergence metadata."""

    concentration_nM: FloatArray
    bound_fraction: FloatArray
    effective_diffusion_um2_s: FloatArray
    iterations: int
    converged: bool
    residual: float


class ProvenanceManifest(BaseModel):
    """Input and execution provenance manifest."""

    model_config = ConfigDict(extra="forbid")

    schema_version: str = "1.0"
    source_id: str
    source_url: str | None = None
    licence: str | None = None
    checksum_sha256: str | None = None
    platform: str
    coordinate_unit: Literal["um"] = "um"
    segmentation_version: str
    vessel_definition: str
    perfusion_measured: bool = False
    antigen_calibration_nM_per_signal: float | None = Field(default=None, gt=0)
    seed: int = 17


class Interval(BaseModel):
    """Median and central uncertainty interval."""

    median: float
    lower: float
    upper: float


class IndexOutput(BaseModel):
    """Public index output with first-class uncertainty and abstention."""

    model_config = ConfigDict(extra="forbid")

    status: Literal["OK", "INSUFFICIENT_EVIDENCE"]
    abstention_reasons: list[str]
    target_positive_fraction: float
    reachable_fraction: Interval
    population_reachable_fraction: Interval
    expression_reach_gap: Interval
    penetration_depth_um: Interval | None
    dominant_barrier: str
    barrier_weights: dict[str, float]
    decision_stability: bool
    engagement_threshold: float
    reachable_decision_threshold: float
    parameter_draws: int
    seed: int
    parameters: dict[str, float]
    parameter_ranges: dict[str, dict[str, Any]]
    library_versions: dict[str, str]


class ClaimsDocument(BaseModel):
    """Statements permitted by one run and statements that remain unsupported."""

    model_config = ConfigDict(extra="forbid")

    permitted: list[str]
    conditional: list[str]
    unsupported: list[str]
    abstention_reasons: list[str]
    interval_basis: str
