"""Cross-platform CosMx validation with explicit pixel units and abstention.

This module validates cell-table, polygon, and relative RNA geometry adapters.  It
never interprets RNA-positive endothelial cells as perfused vessels and never
converts pixel distances to physical length without an externally supplied scale.
"""

from __future__ import annotations

import json
import time
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

_ENDOTHELIAL_MARKERS = ("PECAM1", "VWF", "KDR", "CDH5", "ESAM", "ENG", "RAMP2", "ACKR1")
_TARGETS = ("ERBB2", "EGFR", "CD274", "VSIR", "PDCD1", "LAG3")
_STROMA_MARKERS = ("COL1A1", "COL3A1", "ACTA2", "VIM")
_EPITHELIAL_MARKERS = ("EPCAM", "KRT8", "KRT18")
_PROXY_DEFINITIONS: Mapping[str, int] = {
    "inclusive_2_of_8": 2,
    "balanced_4_of_8": 4,
    "strict_6_of_8": 6,
}


def _normalise_columns(table: pd.DataFrame) -> pd.DataFrame:
    unnamed = [column for column in table.columns if str(column).startswith("Unnamed:")]
    return table.drop(columns=unnamed, errors="ignore")


def read_cosmx_metadata(path: Path) -> pd.DataFrame:
    """Read the CosMx cell metadata and enforce a unique `(fov, cell_ID)` key."""

    table = _normalise_columns(pd.read_csv(path))
    required = {
        "fov",
        "cell_ID",
        "Area",
        "CenterX_local_px",
        "CenterY_local_px",
        "CenterX_global_px",
        "CenterY_global_px",
    }
    missing = sorted(required.difference(table.columns))
    if missing:
        raise KeyError(f"CosMx metadata is missing required columns: {missing}")
    table = table.copy()
    table["fov"] = table["fov"].astype(np.int32)
    table["cell_ID"] = table["cell_ID"].astype(np.int64)
    if table.duplicated(["fov", "cell_ID"]).any():
        raise ValueError("CosMx metadata contains duplicate (fov, cell_ID) keys")
    return table


def read_cosmx_expression(path: Path, genes: Iterable[str] | None = None) -> pd.DataFrame:
    """Read selected expression columns without loading the complete panel into memory."""

    header = _normalise_columns(pd.read_csv(path, nrows=0))
    available = set(header.columns)
    selected = list(_ENDOTHELIAL_MARKERS + _TARGETS + _STROMA_MARKERS + _EPITHELIAL_MARKERS)
    if genes is not None:
        selected = list(dict.fromkeys(genes))
    usecols = [column for column in ["fov", "cell_ID", *selected] if column in available]
    if "fov" not in usecols or "cell_ID" not in usecols:
        raise KeyError("CosMx expression matrix must contain fov and cell_ID")
    table = _normalise_columns(pd.read_csv(path, usecols=usecols))
    table["fov"] = table["fov"].astype(np.int32)
    table["cell_ID"] = table["cell_ID"].astype(np.int64)
    if table.duplicated(["fov", "cell_ID"]).any():
        raise ValueError("CosMx expression contains duplicate (fov, cell_ID) keys")
    return table


def polygon_metrics(path: Path) -> pd.DataFrame:
    """Compute cell polygon area and centroid using a vectorised shoelace reduction."""

    polygon = _normalise_columns(pd.read_csv(path))
    required = {"fov", "cellID", "x_local_px", "y_local_px"}
    missing = sorted(required.difference(polygon.columns))
    if missing:
        raise KeyError(f"CosMx polygon table is missing required columns: {missing}")
    polygon = polygon[["fov", "cellID", "x_local_px", "y_local_px"]].copy()
    polygon["fov"] = polygon["fov"].astype(np.int32)
    polygon["cellID"] = polygon["cellID"].astype(np.int64)
    polygon = polygon.sort_values(["fov", "cellID"], kind="stable").reset_index(drop=True)

    fov = polygon["fov"].to_numpy(dtype=np.int64)
    cell = polygon["cellID"].to_numpy(dtype=np.int64)
    x = polygon["x_local_px"].to_numpy(dtype=np.float64)
    y = polygon["y_local_px"].to_numpy(dtype=np.float64)
    start = np.r_[True, (fov[1:] != fov[:-1]) | (cell[1:] != cell[:-1])]
    starts = np.flatnonzero(start)
    ends = np.r_[starts[1:] - 1, len(polygon) - 1]

    next_x = np.roll(x, -1)
    next_y = np.roll(y, -1)
    next_x[ends] = x[starts]
    next_y[ends] = y[starts]
    cross = x * next_y - next_x * y
    edge_length = np.hypot(next_x - x, next_y - y)

    twice_area = np.add.reduceat(cross, starts)
    signed_area = 0.5 * twice_area
    area = np.abs(signed_area)
    perimeter = np.add.reduceat(edge_length, starts)
    vertices = np.diff(np.r_[starts, len(polygon)])
    sum_cx = np.add.reduceat((x + next_x) * cross, starts)
    sum_cy = np.add.reduceat((y + next_y) * cross, starts)
    mean_x = np.add.reduceat(x, starts) / vertices
    mean_y = np.add.reduceat(y, starts) / vertices
    nonzero = np.abs(signed_area) > np.finfo(np.float64).eps
    centroid_x = np.where(nonzero, sum_cx / (6.0 * signed_area), mean_x)
    centroid_y = np.where(nonzero, sum_cy / (6.0 * signed_area), mean_y)

    return pd.DataFrame(
        {
            "fov": fov[starts].astype(np.int32),
            "cell_ID": cell[starts].astype(np.int64),
            "polygon_area_px2": area,
            "polygon_perimeter_px": perimeter,
            "polygon_centroid_x_px": centroid_x,
            "polygon_centroid_y_px": centroid_y,
            "polygon_vertices": vertices.astype(np.int32),
        }
    )


def segmentation_summary(metadata: pd.DataFrame, metrics: pd.DataFrame) -> dict[str, Any]:
    """Compare supplied metadata geometry to independently reconstructed polygons."""

    merged = metadata.merge(metrics, on=["fov", "cell_ID"], how="left", validate="one_to_one")
    matched = merged["polygon_area_px2"].notna().to_numpy()
    metadata_area = merged.loc[matched, "Area"].to_numpy(dtype=np.float64)
    polygon_area = merged.loc[matched, "polygon_area_px2"].to_numpy(dtype=np.float64)
    relative_error = np.abs(polygon_area - metadata_area) / np.maximum(metadata_area, 1.0)
    centroid_error = np.hypot(
        merged.loc[matched, "polygon_centroid_x_px"].to_numpy(dtype=np.float64)
        - merged.loc[matched, "CenterX_local_px"].to_numpy(dtype=np.float64),
        merged.loc[matched, "polygon_centroid_y_px"].to_numpy(dtype=np.float64)
        - merged.loc[matched, "CenterY_local_px"].to_numpy(dtype=np.float64),
    )
    vertices = merged.loc[matched, "polygon_vertices"].to_numpy(dtype=np.float64)
    return {
        "cells": len(metadata),
        "matched_polygon_cells": int(matched.sum()),
        "matched_fraction": float(matched.mean()),
        "median_relative_area_error": float(np.median(relative_error)),
        "q95_relative_area_error": float(np.quantile(relative_error, 0.95)),
        "within_10pct_area": float(np.mean(relative_error <= 0.10)),
        "median_centroid_error_px": float(np.median(centroid_error)),
        "q95_centroid_error_px": float(np.quantile(centroid_error, 0.95)),
        "within_2px_centroid": float(np.mean(centroid_error <= 2.0)),
        "median_polygon_vertices": float(np.median(vertices)),
    }


def _distances_within_fov(table: pd.DataFrame, proxy: np.ndarray) -> np.ndarray:
    coordinates = table[["CenterX_local_px", "CenterY_local_px"]].to_numpy(dtype=np.float64)
    fovs = table["fov"].to_numpy(dtype=np.int32)
    output: np.ndarray = np.full(len(table), np.nan, dtype=np.float64)
    for fov in np.unique(fovs):
        in_fov = fovs == fov
        source = in_fov & proxy
        if not np.any(source):
            continue
        tree = cKDTree(coordinates[source])
        output[in_fov] = tree.query(coordinates[in_fov], k=1, workers=-1)[0]
    return output


def relative_rna_geometry(
    metadata: pd.DataFrame, expression: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarise target-to-endothelial-proxy geometry in uncalibrated pixels."""

    data = metadata.merge(expression, on=["fov", "cell_ID"], how="inner", validate="one_to_one")
    marker_columns = [marker for marker in _ENDOTHELIAL_MARKERS if marker in data.columns]
    if len(marker_columns) < 4:
        raise ValueError("At least four endothelial RNA markers are required")
    marker_count = (data[marker_columns].to_numpy(dtype=np.float64) > 0).sum(axis=1)
    fov_values = data["fov"].to_numpy(dtype=np.int32)
    target_rows: list[dict[str, Any]] = []
    definition_rows: list[dict[str, Any]] = []
    for definition, required_count in _PROXY_DEFINITIONS.items():
        effective_count = min(required_count, len(marker_columns))
        proxy = marker_count >= effective_count
        distance = _distances_within_fov(data, proxy)
        definition_rows.append(
            {
                "definition": definition,
                "required_markers": effective_count,
                "available_markers": len(marker_columns),
                "proxy_cells": int(proxy.sum()),
                "proxy_fraction": float(proxy.mean()),
                "fovs_with_proxy": int(np.unique(fov_values[proxy]).size),
                "distance_unit": "px",
                "source_semantics": "RNA_ENDOTHELIAL_PROXY_NOT_PERFUSION",
            }
        )
        for target in _TARGETS:
            if target not in data.columns:
                continue
            signal = data[target].to_numpy(dtype=np.float64)
            positive = signal > 0
            valid = positive & np.isfinite(distance)
            target_rows.append(
                {
                    "definition": definition,
                    "target": target,
                    "cells": len(data),
                    "positive_cells": int(positive.sum()),
                    "positive_fraction": float(positive.mean()),
                    "positive_with_distance": int(valid.sum()),
                    "median_distance_px": (
                        float(np.median(distance[valid])) if np.any(valid) else np.nan
                    ),
                    "q90_distance_px": (
                        float(np.quantile(distance[valid], 0.90)) if np.any(valid) else np.nan
                    ),
                    "within_25px": (
                        float(np.mean(distance[valid] <= 25.0)) if np.any(valid) else np.nan
                    ),
                    "within_100px": (
                        float(np.mean(distance[valid] <= 100.0)) if np.any(valid) else np.nan
                    ),
                    "distance_unit": "px",
                    "source_semantics": "RNA_ENDOTHELIAL_PROXY_NOT_PERFUSION",
                }
            )
    return pd.DataFrame(target_rows), pd.DataFrame(definition_rows)


def geometry_sensitivity_audit(
    targets: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Quantify threshold sensitivity and cross-sample target-rank stability.

    The audit is deliberately descriptive.  Pixel distances from different
    samples are not pooled as physical measurements.
    """

    required = {"sample", "definition", "target", "median_distance_px", "within_100px"}
    missing = sorted(required.difference(targets.columns))
    if missing:
        raise KeyError(f"CosMx geometry audit is missing columns: {missing}")

    sensitivity = targets.groupby(["sample", "target"], as_index=False).agg(
        minimum_median_distance_px=("median_distance_px", "min"),
        maximum_median_distance_px=("median_distance_px", "max"),
        minimum_within_100px=("within_100px", "min"),
        maximum_within_100px=("within_100px", "max"),
    )
    sensitivity["definition_distance_ratio"] = sensitivity[
        "maximum_median_distance_px"
    ] / sensitivity["minimum_median_distance_px"].replace(0.0, np.nan)

    balanced = targets[targets["definition"] == "balanced_4_of_8"].pivot(
        index="target", columns="sample", values="median_distance_px"
    )
    rank_rows: list[dict[str, Any]] = []
    samples = [str(sample) for sample in balanced.columns]
    for left_index, left_sample in enumerate(samples):
        for right_sample in samples[left_index + 1 :]:
            paired = balanced[[left_sample, right_sample]].dropna()
            statistic = (
                float(spearmanr(paired[left_sample], paired[right_sample]).statistic)
                if len(paired) >= 3
                else np.nan
            )
            rank_rows.append(
                {
                    "sample_left": left_sample,
                    "sample_right": right_sample,
                    "targets_compared": len(paired),
                    "spearman_target_rank": statistic,
                }
            )
    rank_table = pd.DataFrame(rank_rows)
    ratios = sensitivity["definition_distance_ratio"].dropna().to_numpy(dtype=np.float64)
    rank_values = (
        rank_table["spearman_target_rank"].dropna().to_numpy(dtype=np.float64)
        if not rank_table.empty
        else np.asarray([], dtype=np.float64)
    )
    summary = {
        "median_definition_distance_ratio": float(np.median(ratios)) if ratios.size else np.nan,
        "maximum_definition_distance_ratio": float(np.max(ratios)) if ratios.size else np.nan,
        "median_pairwise_balanced_target_rank_spearman": (
            float(np.median(rank_values)) if rank_values.size else np.nan
        ),
        "target_ranking_generalises_across_samples": bool(
            rank_values.size > 0 and np.all(rank_values >= 0.7)
        ),
        "interpretation": (
            "DEFINITION_SENSITIVE_AND_TARGET_RANKING_NOT_STABLE"
            if ratios.size and (float(np.median(ratios)) > 2.0 or not np.all(rank_values >= 0.7))
            else "RELATIVELY_STABLE_WITHIN_TESTED_DEFINITIONS"
        ),
    }
    return sensitivity, rank_table, summary


def _plot_cosmx(segmentation: pd.DataFrame, targets: pd.DataFrame, output_dir: Path) -> list[str]:
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover
        return []
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []

    figure = plt.figure(figsize=(8, 5))
    plt.bar(segmentation["sample"], segmentation["median_relative_area_error"] * 100.0)
    plt.ylabel("Median polygon area error (%)")
    plt.title("CosMx segmentation reconstruction")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    path = output_dir / "cosmx_segmentation_area_error.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    paths.append(str(path))

    selected = targets[targets["target"].isin(["ERBB2", "EGFR", "CD274", "VSIR"])]
    selected = selected[selected["definition"] == "balanced_4_of_8"]
    figure = plt.figure(figsize=(9, 5))
    for target, group in selected.groupby("target"):
        plt.plot(group["sample"], group["median_distance_px"], marker="o", label=target)
    plt.ylabel("Median distance to RNA endothelial proxy (px)")
    plt.title("Cross-sample relative RNA geometry; not perfusion")
    plt.xticks(rotation=30, ha="right")
    plt.legend()
    plt.tight_layout()
    path = output_dir / "cosmx_target_proxy_distance.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    paths.append(str(path))
    return paths


def prepare_cosmx_external_validation(
    samples: Mapping[str, Mapping[str, Path]], output_dir: Path
) -> dict[str, Any]:
    """Run four-sample CosMx adapter, segmentation, and relative geometry validation."""

    started = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)
    segmentation_rows: list[dict[str, Any]] = []
    target_tables: list[pd.DataFrame] = []
    definition_tables: list[pd.DataFrame] = []
    for sample, paths in samples.items():
        metadata = read_cosmx_metadata(Path(paths["metadata"]))
        expression = read_cosmx_expression(Path(paths["expression"]))
        metrics = polygon_metrics(Path(paths["polygons"]))
        segmentation = segmentation_summary(metadata, metrics)
        expression_keys = expression[["fov", "cell_ID"]]
        matched_expression = metadata.merge(
            expression_keys,
            on=["fov", "cell_ID"],
            how="inner",
            validate="one_to_one",
        )
        segmentation_rows.append(
            {
                "sample": sample,
                **segmentation,
                "expression_cells": len(expression),
                "metadata_expression_matched_cells": len(matched_expression),
                "metadata_expression_matched_fraction": len(matched_expression) / len(metadata),
            }
        )
        targets, definitions = relative_rna_geometry(metadata, expression)
        targets.insert(0, "sample", sample)
        definitions.insert(0, "sample", sample)
        target_tables.append(targets)
        definition_tables.append(definitions)

    segmentation_table = pd.DataFrame(segmentation_rows)
    target_table = pd.concat(target_tables, ignore_index=True)
    definition_table = pd.concat(definition_tables, ignore_index=True)
    segmentation_table.to_csv(output_dir / "cosmx_segmentation_summary.csv", index=False)
    target_table.to_csv(output_dir / "cosmx_target_proxy_geometry.csv", index=False)
    definition_table.to_csv(output_dir / "cosmx_proxy_definitions.csv", index=False)
    sensitivity_table, rank_table, sensitivity_summary = geometry_sensitivity_audit(target_table)
    sensitivity_table.to_csv(output_dir / "cosmx_definition_sensitivity.csv", index=False)
    rank_table.to_csv(output_dir / "cosmx_target_rank_stability.csv", index=False)
    figures = _plot_cosmx(segmentation_table, target_table, output_dir / "figures")
    status = {
        "status": "CROSS_PLATFORM_RELATIVE_GEOMETRY_VALIDATED",
        "samples": len(samples),
        "cells": int(segmentation_table["cells"].sum()),
        "polygon_match_min": float(segmentation_table["matched_fraction"].min()),
        "expression_match_min": float(
            segmentation_table["metadata_expression_matched_fraction"].min()
        ),
        "distance_unit": "px",
        "vascular_source": "RNA_ENDOTHELIAL_PROXY_NOT_PERFUSION",
        "geometry_sensitivity": sensitivity_summary,
        "absolute_index": {
            "status": "NOT_COMPUTED",
            "reasons": [
                "CosMx pixel size was not supplied in the flat files used here",
                "RNA endothelial proxies do not identify functionally perfused vessels",
                "RNA counts are not calibrated surface-antigen density",
                "No administered antibody distribution is measured in these samples",
            ],
        },
        "runtime_seconds": time.time() - started,
        "figures": figures,
    }
    (output_dir / "cosmx_external_validation.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    return status
