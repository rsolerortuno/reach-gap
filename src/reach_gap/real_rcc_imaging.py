"""Real RCC Xenium protein-image validation and relative geometry analysis.

This module adds image-derived vascular and stromal geometry to the compact RCC
Xenium analysis. It intentionally stops before the mechanistic reachability index:
CD31-positive structures are not equivalent to perfused vessels, and image intensity
is not an absolute surface-antigen density.
"""

from __future__ import annotations

import json
import math
import platform
import resource
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from reach_gap.protein_imaging import (
    calibrate_image_threshold,
    distance_band_summary,
    distance_to_mask,
    excess_signal,
    image_mode,
    local_gaussian_signal,
    read_pyramidal_plane,
    sample_local_mean,
    sample_pixels,
    structural_mask,
)
from reach_gap.xenium import (
    extract_selected_h5_features,
    otsu_threshold,
    read_10x_h5_catalog,
    robust_scale,
    sha256_file,
)

_REQUIRED_CHANNELS = ("CD31", "alphaSMA", "Vimentin", "PanCK", "PD-L1", "VISTA")
_H5_COLUMNS = {
    "CD31": "protein__CD31",
    "alphaSMA": "protein__alphaSMA",
    "Vimentin": "protein__Vimentin",
    "PanCK": "protein__PanCK",
    "PD-L1": "protein__PD_L1",
    "VISTA": "protein__VISTA",
}
_TARGET_POSITIVE_COLUMNS = {
    "PD-L1": "target__PD_L1__positive",
    "VISTA": "target__VISTA__positive",
}
_BACKGROUND_COLOUR = {
    "CD31": "grn",
    "alphaSMA": "grn",
    "Vimentin": "blu",
    "PanCK": "red",
    "PD-L1": "yel",
    "VISTA": "grn",
}
_VESSEL_DEFINITIONS = {
    "inclusive_fpr_1pct": 0.99,
    "balanced_fpr_0_5pct": 0.995,
    "strict_fpr_0_1pct": 0.999,
}

_CHANNEL_FILE_NAMES = {
    "CD31": "ch0028_cd31.ome.tif",
    "alphaSMA": "ch0032_alphasma.ome.tif",
    "Vimentin": "ch0031_vimentin.ome.tif",
    "PanCK": "ch0030_panck.ome.tif",
    "PD-L1": "ch0021_pd-l1.ome.tif",
    "VISTA": "ch0020_vista.ome.tif",
}
_BACKGROUND_FILE_NAMES = {
    "blu": "background_02_blu.tiff",
    "grn": "background_02_grn.tiff",
    "yel": "background_02_yel.tiff",
    "red": "background_02_red.tiff",
    "nuv": "background_02_nuv.tiff",
}


def discover_rcc_protein_image_inputs(
    morphology_dir: Path,
    *,
    qc_mask_dir: Path | None = None,
    background_dir: Path | None = None,
) -> tuple[dict[str, Path], dict[str, Path], dict[str, Path]]:
    """Resolve the declared RCC image channels and optional QC inputs by exact name."""

    channel_paths: dict[str, Path] = {}
    missing: list[str] = []
    for channel, file_name in _CHANNEL_FILE_NAMES.items():
        path = morphology_dir / file_name
        if path.is_file():
            channel_paths[channel] = path
        else:
            missing.append(file_name)
    if missing:
        raise FileNotFoundError(
            "Required morphology-focus files are missing: " + ", ".join(sorted(missing))
        )

    qc_paths: dict[str, Path] = {}
    if qc_mask_dir is not None:
        for channel, file_name in _CHANNEL_FILE_NAMES.items():
            path = qc_mask_dir / file_name
            if path.is_file():
                qc_paths[channel] = path

    background_paths: dict[str, Path] = {}
    if background_dir is not None:
        for colour, file_name in _BACKGROUND_FILE_NAMES.items():
            path = background_dir / file_name
            if path.is_file():
                background_paths[colour] = path
    return channel_paths, qc_paths, background_paths


def _load_scored_cells(directory: Path) -> pd.DataFrame:
    paths = sorted(directory.glob("rcc_cells_scored.part*.csv.gz"))
    if not paths:
        paths = sorted(directory.rglob("rcc_cells_scored.part*.csv.gz"))
    if not paths:
        raise FileNotFoundError(f"No scored RCC cell-table parts found under {directory}")
    table = pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
    if table["cell_id"].duplicated().any():
        raise ValueError("Scored cell table contains duplicate cell IDs")
    return table


def _extract_raw_proteins(h5_path: Path, cell_ids: np.ndarray) -> pd.DataFrame:
    catalog = read_10x_h5_catalog(h5_path)
    wanted = set(_H5_COLUMNS)
    indices = [
        index
        for index, (name, feature_type) in enumerate(
            zip(catalog.feature_names, catalog.feature_types, strict=True)
        )
        if name in wanted and feature_type == "Protein Expression"
    ]
    raw = extract_selected_h5_features(h5_path, catalog, indices)
    if not np.array_equal(raw["cell_id"].to_numpy(), cell_ids):
        raise ValueError("HDF5 barcodes are not identical to scored cell-table IDs")
    missing = sorted(set(_H5_COLUMNS.values()).difference(raw.columns))
    if missing:
        raise KeyError(f"Required protein features are missing from HDF5: {missing}")
    return raw


def _load_qc_mask(path: Path | None, level: int, expected_shape: tuple[int, int]) -> np.ndarray:
    if path is None:
        return np.zeros(expected_shape, dtype=bool)
    ome_name = path.name.removesuffix(".tiff") if path.name.endswith(".ome.tif.tiff") else None
    image, _ = read_pyramidal_plane(path, level=level, ome_file_name=ome_name)
    if image.shape != expected_shape:
        raise ValueError(f"QC mask shape {image.shape} does not match image {expected_shape}")
    return np.asarray(image > 0, dtype=bool)


def _reference_positive(raw_values: np.ndarray) -> tuple[np.ndarray, float]:
    scaled = robust_scale(raw_values)
    threshold = float(otsu_threshold(scaled))
    return np.asarray(scaled >= threshold, dtype=bool), threshold


def _safe_spearman(left: np.ndarray, right: np.ndarray) -> float | None:
    """Return Spearman correlation or None for empty/constant inputs."""

    x = np.asarray(left, dtype=np.float64)
    y = np.asarray(right, dtype=np.float64)
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    if x.size < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    value = float(spearmanr(x, y).statistic)
    return value if np.isfinite(value) else None


def _background_correlation(
    channel_path: Path,
    background_path: Path | None,
    *,
    level: int = 4,
) -> dict[str, float | str | None]:
    if background_path is None:
        return {
            "status": "NOT_COMPUTED",
            "reason": "MATCHED_BACKGROUND_QC_IMAGE_NOT_PROVIDED",
            "spearman_all": None,
            "spearman_low_signal": None,
        }
    channel, _ = read_pyramidal_plane(channel_path, level=level)
    try:
        import tifffile
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("tifffile is required for background QC") from exc
    background = np.asarray(tifffile.imread(background_path), dtype=np.float64)
    height = min(channel.shape[0], background.shape[0])
    width = min(channel.shape[1], background.shape[1])
    channel = np.asarray(channel[:height, :width], dtype=np.float64)
    background = background[:height, :width]
    baseline = image_mode(channel, channel > 0)
    residual = np.maximum(channel - baseline, 0.0).ravel()[::50]
    background_sample = background.ravel()[::50]
    all_corr = _safe_spearman(residual, background_sample)
    low = residual <= np.quantile(residual, 0.90)
    low_corr = _safe_spearman(residual[low], background_sample[low])
    return {
        "status": "COMPUTED_RELATIVE_QC",
        "reason": None,
        "spearman_all": all_corr,
        "spearman_low_signal": low_corr,
    }


def _write_cell_h5(table: pd.DataFrame, output_dir: Path) -> list[dict[str, Any]]:
    """Write all per-cell image derivatives to one portable HDF5 file."""

    try:
        import h5py
    except ImportError as exc:  # pragma: no cover - optional real-data dependency
        raise RuntimeError("Writing full image cell tables requires h5py") from exc
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "rcc_cells_image_geometry.h5"
    with h5py.File(str(path), "w") as handle:
        group = handle.create_group("cells")
        group.attrs["rows"] = len(table)
        group.attrs["columns"] = json.dumps(list(table.columns))
        for column in table.columns:
            values = table[column].to_numpy()
            if values.dtype.kind in {"O", "U"}:
                strings = np.asarray([str(value) for value in values])
                width = max(1, max(len(value.encode("utf-8")) for value in strings))
                encoded = np.asarray(
                    [value.encode("utf-8") for value in strings], dtype=f"S{width}"
                )
                group.create_dataset(column, data=encoded, chunks=True)
            else:
                group.create_dataset(column, data=values, chunks=True)
    return [{"path": str(path), "rows": len(table), "sha256": sha256_file(path)}]


def _plot_outputs(
    channel_summary: pd.DataFrame,
    vessel_summary: pd.DataFrame,
    stroma_summary: pd.DataFrame,
    output_dir: Path,
) -> list[Path]:
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover
        return []
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    figure = plt.figure(figsize=(8, 5))
    plt.bar(channel_summary["channel"], channel_summary["cell_image_spearman"])
    plt.axhline(0.0, linewidth=1)
    plt.ylabel("Spearman correlation")
    plt.title("Cell aggregate versus local morphology-focus signal")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    path = output_dir / "cell_image_concordance.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    paths.append(path)

    target_rows = vessel_summary[vessel_summary["target"].isin(["PD-L1", "VISTA"])]
    figure = plt.figure(figsize=(8, 5))
    for target, group in target_rows.groupby("target"):
        plt.plot(group["definition"], group["target_median_distance_um"], marker="o", label=target)
    plt.ylabel("Median distance to image CD31 structure (µm)")
    plt.title("Image-vessel definition sensitivity")
    plt.xticks(rotation=25, ha="right")
    plt.legend()
    plt.tight_layout()
    path = output_dir / "image_vessel_definition_sensitivity.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    paths.append(path)

    reference = stroma_summary[stroma_summary["definition"] == "balanced_fpr_0_5pct"]
    figure = plt.figure(figsize=(8, 5))
    plt.plot(
        reference["distance_band_um"],
        reference["alphaSMA_median"],
        marker="o",
        label="alphaSMA",
    )
    plt.plot(
        reference["distance_band_um"],
        reference["Vimentin_median"],
        marker="o",
        label="Vimentin",
    )
    plt.ylabel("Median excess image signal")
    plt.xlabel("Distance from image CD31 structure (µm)")
    plt.title("Perivascular stromal profile")
    plt.legend()
    plt.tight_layout()
    path = output_dir / "perivascular_stroma_profile.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    paths.append(path)
    return paths


def prepare_rcc_protein_imaging(
    *,
    scored_cells_dir: Path,
    h5_path: Path,
    channel_paths: Mapping[str, Path],
    output_dir: Path,
    qc_mask_paths: Mapping[str, Path] | None = None,
    background_paths: Mapping[str, Path] | None = None,
    level: int = 3,
    write_cell_tables: bool = False,
    source_checksums: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Run real protein-image validation and relative vascular/stromal geometry."""

    started = time.time()
    missing_channels = sorted(set(_REQUIRED_CHANNELS).difference(channel_paths))
    if missing_channels:
        raise KeyError(f"Required morphology channels are missing: {missing_channels}")
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    scored = _load_scored_cells(scored_cells_dir)
    raw = _extract_raw_proteins(h5_path, scored["cell_id"].to_numpy())
    x_um = scored["x_um"].to_numpy(dtype=np.float64)
    y_um = scored["y_um"].to_numpy(dtype=np.float64)

    images: dict[str, np.ndarray] = {}
    masks: dict[str, np.ndarray] = {}
    metadata: dict[str, Any] = {}
    samples: dict[str, np.ndarray] = {}
    positives: dict[str, np.ndarray] = {}
    baselines: dict[str, float] = {}
    channel_rows: list[dict[str, Any]] = []

    for channel in _REQUIRED_CHANNELS:
        path = Path(channel_paths[channel])
        image, meta = read_pyramidal_plane(path, level=level)
        if meta.channel_name.casefold().replace("-", "") != channel.casefold().replace("-", ""):
            raise ValueError(
                f"Channel filename/metadata mismatch for {channel}: OME reports {meta.channel_name}"
            )
        mask_path = None if qc_mask_paths is None else qc_mask_paths.get(channel)
        qc_mask = _load_qc_mask(
            Path(mask_path) if mask_path is not None else None, level, image.shape
        )
        valid = (image > 0) & ~qc_mask
        baseline = image_mode(image, valid)
        local = sample_local_mean(
            image,
            x_um,
            y_um,
            pixel_size_x_um=meta.level_pixel_size_x_um,
            pixel_size_y_um=meta.level_pixel_size_y_um,
            radius_pixels=1,
        )
        raw_values = raw[_H5_COLUMNS[channel]].to_numpy(dtype=np.float64)
        reference_positive, h5_threshold = _reference_positive(raw_values)
        calibration = calibrate_image_threshold(local, reference_positive, negative_quantile=0.995)
        correlation = spearmanr(local, raw_values).statistic
        background_path = None
        if background_paths is not None:
            background_path = background_paths.get(_BACKGROUND_COLOUR[channel])
        background_qc = _background_correlation(path, background_path)
        channel_rows.append(
            {
                "channel": channel,
                "ome_channel_index": meta.channel_index,
                "pyramid_level": level,
                "level_pixel_size_um": meta.level_pixel_size_x_um,
                "baseline_mode": baseline,
                "qc_mask_fraction": float(np.mean(qc_mask)),
                "h5_robust_otsu_threshold": h5_threshold,
                "h5_positive_fraction": float(np.mean(reference_positive)),
                "image_threshold_at_0_5pct_negative_tail": calibration.image_threshold,
                "image_positive_fraction": calibration.image_positive_fraction,
                "image_false_positive_rate": calibration.false_positive_rate,
                "image_true_positive_rate": calibration.true_positive_rate,
                "cell_image_spearman": float(correlation),
                "background_spearman_all": background_qc["spearman_all"],
                "background_spearman_low_signal": background_qc["spearman_low_signal"],
                "source_sha256": sha256_file(path),
            }
        )
        images[channel] = image
        masks[channel] = qc_mask
        metadata[channel] = meta
        samples[channel] = local
        positives[channel] = reference_positive
        baselines[channel] = baseline

    channel_summary = pd.DataFrame(channel_rows)
    channel_summary.to_csv(output_dir / "channel_cell_image_concordance.csv", index=False)

    pixel_size = float(metadata["CD31"].level_pixel_size_x_um)
    if not math.isclose(pixel_size, metadata["CD31"].level_pixel_size_y_um):
        raise ValueError("Anisotropic pixels are not supported for vessel distance")
    cd31_samples = samples["CD31"]
    cd31_negative = ~positives["CD31"]
    cell_x = np.rint(x_um / pixel_size).astype(np.int64)
    cell_y = np.rint(y_um / pixel_size).astype(np.int64)
    tumour = (
        scored["pathology_region"]
        .astype(str)
        .str.contains("Tumor", case=False, regex=False)
        .to_numpy(dtype=bool)
    )
    tissue_mask = np.zeros_like(images["CD31"], dtype=bool)
    for channel in ("CD31", "alphaSMA", "Vimentin", "PanCK"):
        tissue_mask |= images[channel] > 0

    vessel_rows: list[dict[str, Any]] = []
    stroma_rows: list[dict[str, Any]] = []
    cell_image = pd.DataFrame(
        {
            "cell_id": scored["cell_id"],
            "x_um": x_um,
            "y_um": y_um,
            "pathology_region": scored["pathology_region"],
        }
    )
    for channel in _REQUIRED_CHANNELS:
        safe = channel.replace("-", "_")
        cell_image[f"image__{safe}__local_mean"] = samples[channel]
        cell_image[f"image__{safe}__excess"] = np.maximum(
            samples[channel] - baselines[channel], 0.0
        )
        cell_image[f"image__{safe}__qc_masked"] = sample_pixels(
            masks[channel],
            x_um,
            y_um,
            pixel_size_x_um=metadata[channel].level_pixel_size_x_um,
            pixel_size_y_um=metadata[channel].level_pixel_size_y_um,
        ).astype(bool)

    smoothed_stroma: dict[str, np.ndarray] = {}
    for channel in ("alphaSMA", "Vimentin"):
        smoothed_stroma[channel] = local_gaussian_signal(
            images[channel],
            pixel_size_um=pixel_size,
            sigma_um=15.0,
            baseline=baselines[channel],
        )
        cell_image[f"image__{channel}__local_15um"] = sample_pixels(
            smoothed_stroma[channel],
            x_um,
            y_um,
            pixel_size_x_um=pixel_size,
            pixel_size_y_um=pixel_size,
        )

    reference_distance: np.ndarray | None = None
    reference_mask: np.ndarray | None = None
    for definition, negative_quantile in _VESSEL_DEFINITIONS.items():
        threshold = float(np.quantile(cd31_samples[cd31_negative], negative_quantile))
        vessel_mask, diagnostics = structural_mask(
            images["CD31"],
            threshold=threshold,
            pixel_size_um=pixel_size,
            minimum_component_area_um2=25.0,
            closing_radius_um=pixel_size,
            qc_mask=masks["CD31"],
        )
        distance_image = distance_to_mask(vessel_mask, pixel_size_um=pixel_size)
        cell_distance = distance_image[cell_y, cell_x]
        cell_image[f"distance_to_image_vessel__{definition}_um"] = cell_distance
        if definition == "balanced_fpr_0_5pct":
            reference_distance = cell_distance
            reference_mask = vessel_mask
        base = {
            "definition": definition,
            "negative_quantile": negative_quantile,
            **diagnostics,
            "all_cell_median_distance_um": float(np.median(cell_distance)),
            "tumour_cell_median_distance_um": float(np.median(cell_distance[tumour])),
            "h5_cd31_positive_within_3_4um": float(
                np.mean(cell_distance[positives["CD31"]] <= 3.4)
            ),
            "h5_cd31_negative_within_3_4um": float(
                np.mean(cell_distance[~positives["CD31"]] <= 3.4)
            ),
        }
        for target, positive_column in _TARGET_POSITIVE_COLUMNS.items():
            positive = scored[positive_column].to_numpy(dtype=bool) & tumour
            row = dict(base)
            row.update(
                {
                    "target": target,
                    "target_positive_tumour_cells": int(np.sum(positive)),
                    "target_median_distance_um": float(np.median(cell_distance[positive])),
                    "target_q90_distance_um": float(np.quantile(cell_distance[positive], 0.90)),
                    "target_within_25um": float(np.mean(cell_distance[positive] <= 25.0)),
                    "target_within_50um": float(np.mean(cell_distance[positive] <= 50.0)),
                    "target_within_100um": float(np.mean(cell_distance[positive] <= 100.0)),
                }
            )
            vessel_rows.append(row)
        profiles = distance_band_summary(
            distance_image,
            {
                "alphaSMA": excess_signal(images["alphaSMA"], baseline=baselines["alphaSMA"]),
                "Vimentin": excess_signal(images["Vimentin"], baseline=baselines["Vimentin"]),
            },
            valid_mask=tissue_mask,
        )
        for profile in profiles:
            profile["definition"] = definition
            stroma_rows.append(profile)

    if reference_distance is None or reference_mask is None:
        raise AssertionError("Balanced vessel definition was not generated")
    cell_image["distance_to_image_vessel_um"] = reference_distance
    cell_image["image_vessel_reference"] = reference_distance <= 3.4

    target_rows: list[dict[str, Any]] = []
    for target, positive_column in _TARGET_POSITIVE_COLUMNS.items():
        reference_positive = scored[positive_column].to_numpy(dtype=bool)
        image_threshold = float(
            channel_summary.loc[
                channel_summary["channel"] == target,
                "image_threshold_at_0_5pct_negative_tail",
            ].iloc[0]
        )
        image_positive = samples[target] > image_threshold
        cell_image[f"image__{target.replace('-', '_')}__positive"] = image_positive
        for subset_name, subset in {
            "all_cells": np.ones(len(scored), dtype=bool),
            "tumour_region": tumour,
            "immune_high": scored["immune_score"].to_numpy(dtype=float) >= 0.5,
            "malignant_proxy": scored["cell_is_malignant_proxy"].to_numpy(dtype=bool),
        }.items():
            active = subset
            both = active & reference_positive & image_positive
            h5_only = active & reference_positive & ~image_positive
            image_only = active & ~reference_positive & image_positive
            target_rows.append(
                {
                    "target": target,
                    "subset": subset_name,
                    "cells": int(np.sum(active)),
                    "h5_positive_fraction": float(np.mean(reference_positive[active])),
                    "image_positive_fraction": float(np.mean(image_positive[active])),
                    "both_positive_fraction": float(np.sum(both) / max(np.sum(active), 1)),
                    "h5_only_fraction": float(np.sum(h5_only) / max(np.sum(active), 1)),
                    "image_only_fraction": float(np.sum(image_only) / max(np.sum(active), 1)),
                    "both_positive_median_distance_um": (
                        float(np.median(reference_distance[both])) if np.any(both) else float("nan")
                    ),
                    "median_alphaSMA_15um_both_positive": (
                        float(
                            np.median(
                                cast(
                                    "pd.Series[Any]",
                                    cell_image.loc[both, "image__alphaSMA__local_15um"],
                                ).to_numpy(dtype=np.float64)
                            )
                        )
                        if np.any(both)
                        else float("nan")
                    ),
                    "median_vimentin_15um_both_positive": (
                        float(
                            np.median(
                                cast(
                                    "pd.Series[Any]",
                                    cell_image.loc[both, "image__Vimentin__local_15um"],
                                ).to_numpy(dtype=np.float64)
                            )
                        )
                        if np.any(both)
                        else float("nan")
                    ),
                }
            )

    vessel_summary = pd.DataFrame(vessel_rows)
    stroma_summary = pd.DataFrame(stroma_rows)
    target_summary = pd.DataFrame(target_rows)
    vessel_summary.to_csv(output_dir / "image_vessel_definition_sensitivity.csv", index=False)
    stroma_summary.to_csv(output_dir / "perivascular_stroma_profile.csv", index=False)
    target_summary.to_csv(output_dir / "target_image_assignment_summary.csv", index=False)
    cell_manifest = (
        _write_cell_h5(cell_image, output_dir / "cell_tables") if write_cell_tables else []
    )
    figure_paths = _plot_outputs(channel_summary, vessel_summary, stroma_summary, figures_dir)

    claims = {
        "permitted": [
            "The full RCC section was analysed at a declared OME pyramid level.",
            "CD31 image structures provide a relative vascular-geometry proxy.",
            "alphaSMA and Vimentin provide relative perivascular stromal signal.",
            "Matched morphology-focus and cell-aggregate protein measurements can be compared.",
            "Vessel-definition sensitivity is propagated into target-distance summaries.",
        ],
        "conditional": [
            "CD31 structures may represent vessel walls, but functional perfusion is unmeasured.",
            "Local image intensity is relative and depends on acquisition and correction settings.",
        ],
        "unsupported": [
            "The CD31 mask identifies functionally perfused vessels.",
            "Xenium image intensity equals surface-antigen molecules per cell.",
            "A real reachable_fraction or expression_reach_gap has been computed.",
            "The image geometry predicts clinical response or programme outcome.",
        ],
        "abstention_reasons": [
            "FUNCTIONALLY_PERFUSED_VESSELS_NOT_IDENTIFIED",
            "SURFACE_ANTIGEN_CALIBRATION_NOT_AVAILABLE",
            "MATRIX_TRANSPORT_COEFFICIENT_NOT_CALIBRATED",
            "DRUG_DISTRIBUTION_OR_ENGAGEMENT_NOT_MEASURED",
        ],
    }
    (output_dir / "claims.json").write_text(json.dumps(claims, indent=2), encoding="utf-8")

    input_hashes: dict[str, str] = {"h5": f"sha256:{sha256_file(h5_path)}"}
    for name, path in channel_paths.items():
        key = f"channel_{name}"
        if source_checksums is not None and key in source_checksums:
            input_hashes[key] = str(source_checksums[key])
        else:
            input_hashes[key] = f"sha256:{sha256_file(Path(path))}"
    if qc_mask_paths:
        input_hashes.update(
            {
                f"qc_mask_{name}": f"sha256:{sha256_file(Path(path))}"
                for name, path in qc_mask_paths.items()
            }
        )
    if background_paths:
        input_hashes.update(
            {
                f"background_{name}": f"sha256:{sha256_file(Path(path))}"
                for name, path in background_paths.items()
            }
        )
    summary = {
        "status": "REAL_PROTEIN_IMAGE_GEOMETRY_PREPARED_ABSOLUTE_INDEX_NOT_COMPUTED",
        "cells": len(scored),
        "channels": list(_REQUIRED_CHANNELS),
        "pyramid_level": level,
        "level_pixel_size_um": pixel_size,
        "runtime_seconds": time.time() - started,
        "peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
        "python": platform.python_version(),
        "primary_vessel_definition": "balanced_fpr_0_5pct",
        "input_sha256": input_hashes,
        "cell_table_parts": cell_manifest,
        "figures": [str(path) for path in figure_paths],
        "abstention_reasons": claims["abstention_reasons"],
    }
    (output_dir / "real_imaging_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary
