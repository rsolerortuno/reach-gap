"""Spatial feature extraction and conservative flat-table ingestion."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt

from reach_gap.schemas import ProvenanceManifest, SpatialFeatures, TissueGeometry

REQUIRED_COLUMNS = {
    "cell_id",
    "x_um",
    "y_um",
    "is_tumour",
    "target_signal",
    "vessel_signal",
    "ecm_score",
    "caf_score",
}


def extract_features(
    geometry: TissueGeometry,
    *,
    antigen_calibration_nM_per_signal: float | None = 200.0,
    target_positive_threshold: float = 0.5,
) -> SpatialFeatures:
    """Convert geometry fields into model-ready spatial features."""

    if antigen_calibration_nM_per_signal is not None and antigen_calibration_nM_per_signal <= 0:
        raise ValueError("Antigen calibration must be positive")
    distances = (
        np.asarray(distance_transform_edt(~geometry.vessel_mask), dtype=np.float64) * geometry.dx_um
    )
    calibrated = antigen_calibration_nM_per_signal is not None
    scale = 1.0 if antigen_calibration_nM_per_signal is None else antigen_calibration_nM_per_signal
    antigen_nM = geometry.antigen * scale
    target_positive = geometry.cell_target_signal >= target_positive_threshold
    return SpatialFeatures(
        vessel_mask=geometry.vessel_mask,
        vessel_distance_um=distances.astype(np.float64),
        ecm=np.clip(geometry.ecm, 0.0, 1.0),
        caf=np.clip(geometry.caf, 0.0, 1.0),
        antigen_nM=antigen_nM.astype(np.float64),
        tumour_mask=geometry.tumour_mask,
        cell_rows=geometry.cell_rows,
        cell_cols=geometry.cell_cols,
        cell_is_tumour=geometry.cell_is_tumour,
        cell_target_positive=target_positive.astype(np.bool_),
        dx_um=geometry.dx_um,
        antigen_calibrated=calibrated,
        seed=geometry.seed,
    )


def ingest_cell_table(
    cell_csv: Path,
    manifest_path: Path,
    output_dir: Path,
    *,
    grid_size: int = 64,
) -> tuple[Path, Path]:
    """Rasterise a generic spatial cell table without inferring biological markers."""

    table = pd.read_csv(cell_csv)
    missing = REQUIRED_COLUMNS.difference(table.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    manifest = ProvenanceManifest.model_validate_json(manifest_path.read_text(encoding="utf-8"))
    if grid_size < 8:
        raise ValueError("grid_size must be at least 8")

    x = table["x_um"].to_numpy(dtype=np.float64)
    y = table["y_um"].to_numpy(dtype=np.float64)
    x_min, x_max = float(x.min()), float(x.max())
    y_min, y_max = float(y.min()), float(y.max())
    span = max(x_max - x_min, y_max - y_min, 1.0)
    dx_um = span / float(grid_size - 1)
    cols = np.clip(np.rint((x - x_min) / dx_um), 0, grid_size - 1).astype(np.int64)
    rows = np.clip(np.rint((y - y_min) / dx_um), 0, grid_size - 1).astype(np.int64)

    count = np.zeros((grid_size, grid_size), dtype=np.float64)
    target = np.zeros_like(count)
    vessel = np.zeros_like(count)
    ecm = np.zeros_like(count)
    caf = np.zeros_like(count)
    tumour = np.zeros((grid_size, grid_size), dtype=np.bool_)
    for row, col, target_value, vessel_value, ecm_value, caf_value, is_tumour in zip(
        rows,
        cols,
        table["target_signal"].to_numpy(dtype=np.float64),
        table["vessel_signal"].to_numpy(dtype=np.float64),
        table["ecm_score"].to_numpy(dtype=np.float64),
        table["caf_score"].to_numpy(dtype=np.float64),
        table["is_tumour"].to_numpy(dtype=np.bool_),
        strict=True,
    ):
        count[row, col] += 1.0
        target[row, col] += target_value
        vessel[row, col] += vessel_value
        ecm[row, col] += ecm_value
        caf[row, col] += caf_value
        tumour[row, col] |= bool(is_tumour)
    occupied = count > 0
    for field in (target, vessel, ecm, caf):
        field[occupied] /= count[occupied]
    vessel_mask = vessel >= 0.5
    calibration = manifest.antigen_calibration_nM_per_signal
    antigen_scale = 1.0 if calibration is None else calibration
    features = SpatialFeatures(
        vessel_mask=vessel_mask,
        vessel_distance_um=(
            np.asarray(distance_transform_edt(~vessel_mask), dtype=np.float64) * dx_um
        ).astype(np.float64),
        ecm=np.clip(ecm, 0.0, 1.0),
        caf=np.clip(caf, 0.0, 1.0),
        antigen_nM=(np.clip(target, 0.0, None) * antigen_scale).astype(np.float64),
        tumour_mask=tumour,
        cell_rows=rows,
        cell_cols=cols,
        cell_is_tumour=table["is_tumour"].to_numpy(dtype=np.bool_),
        cell_target_positive=(table["target_signal"].to_numpy(dtype=np.float64) >= 0.5),
        dx_um=dx_um,
        antigen_calibrated=calibration is not None,
        seed=manifest.seed,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    features_path = output_dir / "features.npz"
    save_features(features, features_path)
    manifest_out = output_dir / "manifest.json"
    manifest_out.write_text(
        json.dumps(manifest.model_dump(mode="json"), indent=2), encoding="utf-8"
    )
    return features_path, manifest_out


def save_features(features: SpatialFeatures, path: Path) -> None:
    """Save features in a compact, deterministic NumPy archive."""

    np.savez_compressed(
        path,
        vessel_mask=features.vessel_mask,
        vessel_distance_um=features.vessel_distance_um,
        ecm=features.ecm,
        caf=features.caf,
        antigen_nM=features.antigen_nM,
        tumour_mask=features.tumour_mask,
        cell_rows=features.cell_rows,
        cell_cols=features.cell_cols,
        cell_is_tumour=features.cell_is_tumour,
        cell_target_positive=features.cell_target_positive,
        dx_um=np.array(features.dx_um),
        antigen_calibrated=np.array(features.antigen_calibrated),
        seed=np.array(features.seed),
    )


def load_features(path: Path) -> SpatialFeatures:
    """Load a saved feature archive."""

    with np.load(path, allow_pickle=False) as data:
        return SpatialFeatures(
            vessel_mask=data["vessel_mask"].astype(np.bool_),
            vessel_distance_um=data["vessel_distance_um"].astype(np.float64),
            ecm=data["ecm"].astype(np.float64),
            caf=data["caf"].astype(np.float64),
            antigen_nM=data["antigen_nM"].astype(np.float64),
            tumour_mask=data["tumour_mask"].astype(np.bool_),
            cell_rows=data["cell_rows"].astype(np.int64),
            cell_cols=data["cell_cols"].astype(np.int64),
            cell_is_tumour=data["cell_is_tumour"].astype(np.bool_),
            cell_target_positive=data["cell_target_positive"].astype(np.bool_),
            dx_um=float(data["dx_um"]),
            antigen_calibrated=bool(data["antigen_calibrated"]),
            seed=int(data["seed"]),
        )
