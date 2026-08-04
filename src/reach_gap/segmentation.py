"""Streaming segmentation-boundary robustness checks for Xenium outputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from reach_gap.parquet_lite import iter_flat_parquet_row_groups


def _polygon_metrics(cell_id: str, x: np.ndarray, y: np.ndarray) -> dict[str, object]:
    if len(x) < 3:
        return {
            "cell_id": cell_id,
            "boundary_vertices": len(x),
            "boundary_area": np.nan,
            "boundary_perimeter": np.nan,
            "boundary_compactness": np.nan,
        }
    x_next = np.roll(x, -1)
    y_next = np.roll(y, -1)
    area = 0.5 * abs(float(np.sum(x * y_next - y * x_next, dtype=np.float64)))
    perimeter = float(np.sum(np.hypot(x_next - x, y_next - y), dtype=np.float64))
    compactness = float(4.0 * np.pi * area / (perimeter * perimeter)) if perimeter > 0 else np.nan
    return {
        "cell_id": cell_id,
        "boundary_vertices": len(x),
        "boundary_area": area,
        "boundary_perimeter": perimeter,
        "boundary_compactness": compactness,
    }


def _complete_groups(
    ids: np.ndarray, x: np.ndarray, y: np.ndarray, *, retain_last: bool
) -> tuple[pd.DataFrame, tuple[np.ndarray, np.ndarray, np.ndarray] | None]:
    if len(ids) == 0:
        return pd.DataFrame(), None
    starts = np.r_[0, np.flatnonzero(ids[1:] != ids[:-1]) + 1]
    stops = np.r_[starts[1:], len(ids)]
    process_count = len(starts) - 1 if retain_last else len(starts)
    if process_count == 0:
        return pd.DataFrame(), (ids.copy(), x.copy(), y.copy())

    process_stop = int(stops[process_count - 1])
    p_ids = ids[:process_stop]
    p_x = x[:process_stop]
    p_y = y[:process_stop]
    p_starts = starts[:process_count]
    p_stops = stops[:process_count]

    x_next = np.roll(p_x, -1)
    y_next = np.roll(p_y, -1)
    last_indices = p_stops - 1
    x_next[last_indices] = p_x[p_starts]
    y_next[last_indices] = p_y[p_starts]
    cross = p_x * y_next - p_y * x_next
    edge = np.hypot(x_next - p_x, y_next - p_y)
    area = 0.5 * np.abs(np.add.reduceat(cross, p_starts))
    perimeter = np.add.reduceat(edge, p_starts)
    vertices = p_stops - p_starts
    compactness = np.divide(
        4.0 * np.pi * area,
        perimeter * perimeter,
        out=np.full_like(area, np.nan),
        where=perimeter > 0,
    )
    frame = pd.DataFrame(
        {
            "cell_id": p_ids[p_starts].astype(str),
            "boundary_vertices": vertices.astype(np.int32),
            "boundary_area": area,
            "boundary_perimeter": perimeter,
            "boundary_compactness": compactness,
        }
    )
    carry = None
    if retain_last:
        carry_start = int(starts[-1])
        carry = (ids[carry_start:].copy(), x[carry_start:].copy(), y[carry_start:].copy())
    return frame, carry


def boundary_polygon_metrics(path: Path) -> pd.DataFrame:
    """Compute polygon metrics by streaming ordered Xenium boundary vertices."""

    parts: list[pd.DataFrame] = []
    carry: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
    groups = iter_flat_parquet_row_groups(path, ["cell_id", "vertex_x", "vertex_y"])
    for frame in groups:
        ids = frame["cell_id"].astype(str).to_numpy()
        x = frame["vertex_x"].to_numpy(dtype=np.float64)
        y = frame["vertex_y"].to_numpy(dtype=np.float64)
        if carry is not None:
            ids = np.concatenate([carry[0], ids])
            x = np.concatenate([carry[1], x])
            y = np.concatenate([carry[2], y])
        complete, carry = _complete_groups(ids, x, y, retain_last=True)
        if not complete.empty:
            parts.append(complete)
    if carry is not None:
        complete, _ = _complete_groups(*carry, retain_last=False)
        if not complete.empty:
            parts.append(complete)
    result = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if result["cell_id"].duplicated().any():
        raise ValueError("Boundary table produced duplicate polygon cell IDs")
    return result


def segmentation_robustness_summary(
    cells: pd.DataFrame,
    cell_metrics: pd.DataFrame,
    nucleus_metrics: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Compare boundary-derived polygon areas with Xenium's cell summary areas."""

    merged = cells[["cell_id", "cell_area", "nucleus_area"]].merge(
        cell_metrics, on="cell_id", how="left", validate="one_to_one"
    )
    nucleus = nucleus_metrics.rename(
        columns={
            "boundary_vertices": "nucleus_boundary_vertices",
            "boundary_area": "nucleus_boundary_area",
            "boundary_perimeter": "nucleus_boundary_perimeter",
            "boundary_compactness": "nucleus_boundary_compactness",
        }
    )
    merged = merged.merge(nucleus, on="cell_id", how="left", validate="one_to_one")
    merged["cell_area_relative_error"] = np.abs(
        merged["boundary_area"] - merged["cell_area"]
    ) / np.maximum(merged["cell_area"], 1.0e-9)
    merged["nucleus_area_relative_error"] = np.abs(
        merged["nucleus_boundary_area"] - merged["nucleus_area"]
    ) / np.maximum(merged["nucleus_area"], 1.0e-9)

    def metrics(prefix: str, area_error: str, vertices: str, compactness: str) -> dict[str, object]:
        values = merged[area_error].to_numpy(dtype=np.float64)
        finite = np.isfinite(values)
        return {
            f"{prefix}_polygons": int(finite.sum()),
            f"{prefix}_coverage_fraction": float(finite.mean()),
            f"{prefix}_median_relative_area_error": float(np.nanmedian(values)),
            f"{prefix}_q95_relative_area_error": float(np.nanquantile(values, 0.95)),
            f"{prefix}_within_1pct": float(np.nanmean(values <= 0.01)),
            f"{prefix}_within_5pct": float(np.nanmean(values <= 0.05)),
            f"{prefix}_within_10pct": float(np.nanmean(values <= 0.10)),
            f"{prefix}_median_vertices": float(np.nanmedian(merged[vertices])),
            f"{prefix}_median_compactness": float(np.nanmedian(merged[compactness])),
        }

    summary: dict[str, object] = {
        "cells_total": len(merged),
        **metrics("cell", "cell_area_relative_error", "boundary_vertices", "boundary_compactness"),
        **metrics(
            "nucleus",
            "nucleus_area_relative_error",
            "nucleus_boundary_vertices",
            "nucleus_boundary_compactness",
        ),
    }
    return merged, summary
