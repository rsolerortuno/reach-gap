"""Relative in-vivo perfusion validation from CD31/Hoechst confocal composites.

The BioStudies S-BIAD3159 TIFF exports are RGB composites rather than raw
single-channel stacks.  This module therefore validates a *relative* relation
between the Hoechst perfusion proxy and distance to CD31-positive structures.
It does not identify perfused vessels in the RCC Xenium section and it does not
turn CD31 into a pharmacological source boundary.
"""

from __future__ import annotations

import itertools
import json
import math
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.ndimage import distance_transform_edt
from scipy.stats import spearmanr

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]

DEFAULT_DISTANCE_BINS_UM: tuple[float, ...] = (0.0, 10.0, 25.0, 50.0, 100.0, math.inf)


def infer_tiff_pixel_size_um(path: Path) -> float | None:
    """Infer micrometres per pixel from an ImageJ TIFF when metadata permit it."""

    try:
        import tifffile
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Perfusion TIFF validation requires tifffile") from exc

    with tifffile.TiffFile(path) as tiff:
        page = cast(Any, tiff.pages[0])
        description = str(page.description or "").casefold()
        if "unit=micron" not in description and "unit=µm" not in description:
            return None
        tag = page.tags.get("XResolution")
        if tag is None:
            return None
        numerator, denominator = tag.value
        pixels_per_micron = float(numerator) / float(denominator)
        if pixels_per_micron <= 0:
            return None
        return 1.0 / pixels_per_micron


def read_rgb_tiff(path: Path, downsample: int = 1) -> NDArray[np.uint8]:
    """Read an RGB TIFF and optionally stride-downsample it."""

    if downsample < 1:
        raise ValueError("downsample must be >= 1")
    try:
        import tifffile
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Perfusion TIFF validation requires tifffile") from exc
    image = np.asarray(tifffile.imread(path))
    if image.ndim != 3 or image.shape[-1] < 3:
        raise ValueError(f"Expected an RGB TIFF, observed shape {image.shape}")
    return np.asarray(image[::downsample, ::downsample, :3], dtype=np.uint8)


def _otsu_threshold(values: FloatArray) -> float:
    if values.size == 0 or np.all(values == values[0]):
        return float(values[0]) if values.size else 1.0
    try:
        from skimage.filters import threshold_otsu
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Perfusion validation requires scikit-image") from exc
    return float(threshold_otsu(values))  # type: ignore[no-untyped-call]


def _clean_small_objects(mask: BoolArray, minimum_pixels: int) -> BoolArray:
    if minimum_pixels <= 1:
        return mask
    try:
        from scipy import ndimage
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Perfusion validation requires scipy") from exc
    labels_raw, count = cast(tuple[NDArray[Any], int], ndimage.label(mask))
    if count == 0:
        return mask
    labels = np.asarray(labels_raw, dtype=np.int64)
    sizes = np.bincount(labels.ravel())
    keep = sizes >= minimum_pixels
    keep[0] = False
    return np.asarray(keep[labels], dtype=np.bool_)


def segment_cd31(
    green: FloatArray,
    *,
    quantile_scale: float = 0.995,
    minimum_component_pixels: int = 20,
) -> tuple[BoolArray, dict[str, float]]:
    """Segment bright CD31 structures using an auditable within-image threshold."""

    if green.ndim != 2:
        raise ValueError("green channel must be two-dimensional")
    positive_values = green[green > 0]
    if positive_values.size == 0:
        raise ValueError("green channel contains no positive signal")
    scale = float(np.quantile(positive_values, quantile_scale))
    if scale <= 0:
        raise ValueError("green-channel robust scale must be positive")
    scaled = np.clip(green / scale, 0.0, 1.0)
    threshold = max(_otsu_threshold(scaled[scaled > 0]), 0.05)
    mask = np.asarray(scaled >= threshold, dtype=np.bool_)
    mask = _clean_small_objects(mask, minimum_component_pixels)
    return mask, {
        "green_scale_quantile": quantile_scale,
        "green_scale_value": scale,
        "otsu_threshold_relative": threshold,
    }


def perfusion_profile_from_rgb(
    rgb: NDArray[np.uint8],
    *,
    pixel_size_um: float,
    distance_bins_um: Sequence[float] = DEFAULT_DISTANCE_BINS_UM,
    red_correction_alpha: float = 0.0,
    minimum_component_pixels: int = 20,
) -> tuple[pd.DataFrame, dict[str, Any], BoolArray, FloatArray, FloatArray]:
    """Compute Hoechst intensity as a function of distance to CD31 structures."""

    if pixel_size_um <= 0:
        raise ValueError("pixel_size_um must be positive")
    bins = tuple(float(value) for value in distance_bins_um)
    if len(bins) < 2 or any(right <= left for left, right in itertools.pairwise(bins)):
        raise ValueError("distance bins must be strictly increasing")

    image = np.asarray(rgb, dtype=np.float64)
    red = image[..., 0]
    green = image[..., 1]
    blue = image[..., 2]
    corrected_blue = np.clip(blue - red_correction_alpha * red, 0.0, None)

    vessel_mask, threshold_metadata = segment_cd31(
        green, minimum_component_pixels=minimum_component_pixels
    )
    distance_raw = cast(FloatArray, distance_transform_edt(~vessel_mask))
    distance_um = np.asarray(distance_raw, dtype=np.float64) * pixel_size_um
    tissue_mask = np.asarray(np.max(image, axis=2) > 0, dtype=np.bool_)
    if not np.any(tissue_mask):
        raise ValueError("RGB image contains no positive signal")

    scale = float(np.quantile(corrected_blue[tissue_mask], 0.99))
    if scale <= 0:
        raise ValueError("Hoechst channel contains no positive signal")
    normalized_hoechst = np.clip(corrected_blue / scale, 0.0, 1.0)

    rows: list[dict[str, Any]] = []
    for lower, upper in itertools.pairwise(bins):
        selected = tissue_mask & (distance_um >= lower) & (distance_um < upper)
        values = normalized_hoechst[selected]
        rows.append(
            {
                "distance_lower_um": lower,
                "distance_upper_um": None if math.isinf(upper) else upper,
                "pixels": int(values.size),
                "hoechst_mean_relative": float(np.mean(values)) if values.size else math.nan,
                "hoechst_median_relative": float(np.median(values)) if values.size else math.nan,
                "hoechst_q90_relative": float(np.quantile(values, 0.90))
                if values.size
                else math.nan,
            }
        )

    distances = distance_um[tissue_mask]
    intensities = normalized_hoechst[tissue_mask]
    stride = max(1, distances.size // 200_000)
    correlation = cast(tuple[float, float], spearmanr(distances[::stride], intensities[::stride]))
    correlation_statistic = correlation[0]
    correlation_pvalue = correlation[1]
    near = tissue_mask & (distance_um < 10.0)
    far = tissue_mask & (distance_um >= 50.0) & (distance_um < 100.0)
    near_mean = float(np.mean(normalized_hoechst[near])) if np.any(near) else math.nan
    far_mean = float(np.mean(normalized_hoechst[far])) if np.any(far) else math.nan
    ratio = near_mean / far_mean if far_mean > 0 else math.inf

    summary: dict[str, Any] = {
        "pixel_size_um": pixel_size_um,
        "red_correction_alpha": red_correction_alpha,
        "tissue_signal_fraction": float(np.mean(tissue_mask)),
        "cd31_fraction_of_field": float(np.mean(vessel_mask)),
        "hoechst_q99_scale": scale,
        "distance_hoechst_spearman": correlation_statistic,
        "distance_hoechst_spearman_pvalue": correlation_pvalue,
        "near_lt10um_mean_relative": near_mean,
        "far_50_100um_mean_relative": far_mean,
        "near_to_far_mean_ratio": ratio,
        **threshold_metadata,
    }
    return pd.DataFrame(rows), summary, vessel_mask, distance_um, normalized_hoechst


def _plot_profiles(table: pd.DataFrame, output_dir: Path) -> list[str]:
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover
        return []
    output_dir.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(8, 5))
    for label, group in table.groupby("image_label", sort=True):
        midpoint = []
        for lower, upper in zip(
            group["distance_lower_um"], group["distance_upper_um"], strict=True
        ):
            midpoint.append(
                float(lower) + 25.0 if pd.isna(upper) else (float(lower) + float(upper)) / 2
            )
        plt.plot(midpoint, group["hoechst_mean_relative"], marker="o", label=str(label))
    plt.xlabel("Distance to CD31 structure (µm)")
    plt.ylabel("Mean relative Hoechst intensity")
    plt.title("Independent in-vivo perfusion-proxy profile")
    plt.legend()
    plt.tight_layout()
    path = output_dir / "hoechst_vs_cd31_distance.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return [str(path)]


def benchmark_perfusion_tiffs(
    images: Mapping[str, Path],
    output_dir: Path,
    *,
    downsample: int = 4,
    default_pixel_size_um: float = 0.3,
    red_correction_by_label: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    """Benchmark relative perfusion profiles across RGB composite TIFFs."""

    started = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles: list[pd.DataFrame] = []
    summaries: list[dict[str, Any]] = []
    correction = dict(red_correction_by_label or {})

    for label, path in images.items():
        native_pixel_size = infer_tiff_pixel_size_um(path)
        pixel_size = (native_pixel_size or default_pixel_size_um) * downsample
        profile, summary, _, _, _ = perfusion_profile_from_rgb(
            read_rgb_tiff(path, downsample=downsample),
            pixel_size_um=pixel_size,
            red_correction_alpha=float(correction.get(label, 0.0)),
        )
        profile.insert(0, "image_label", label)
        profile.insert(1, "image", path.name)
        profiles.append(profile)
        summaries.append(
            {
                "image_label": label,
                "image": path.name,
                "native_pixel_size_um": native_pixel_size,
                "pixel_size_source": "TIFF_XRESOLUTION_IMAGEJ_MICRON"
                if native_pixel_size is not None
                else "CONFIGURED_DEFAULT",
                **summary,
            }
        )

    if not profiles:
        raise ValueError("No perfusion TIFFs were supplied")
    profile_table = pd.concat(profiles, ignore_index=True)
    summary_table = pd.DataFrame(summaries)
    profile_table.to_csv(output_dir / "perfusion_distance_profiles.csv", index=False)
    summary_table.to_csv(output_dir / "perfusion_image_summary.csv", index=False)
    figures = _plot_profiles(profile_table, output_dir / "figures")

    all_negative = bool((summary_table["distance_hoechst_spearman"] < 0).all())
    all_near_enriched = bool((summary_table["near_to_far_mean_ratio"] > 1).all())
    result: dict[str, Any] = {
        "status": "IN_VIVO_HOECHST_PERFUSION_RELATIVE_VALIDATION_COMPOSITE_IMAGES",
        "images": len(summary_table),
        "downsample": downsample,
        "all_images_negative_distance_correlation": all_negative,
        "all_images_near_vessel_enrichment": all_near_enriched,
        "median_distance_hoechst_spearman": float(
            summary_table["distance_hoechst_spearman"].median()
        ),
        "median_near_to_far_mean_ratio": float(summary_table["near_to_far_mean_ratio"].median()),
        "rcc_transfer": {
            "status": "NOT_COMPUTED",
            "reasons": [
                "The validation images are independent LLC mouse tumours, "
                "not the RCC Xenium section",
                "The TIFF exports are RGB composites rather than raw single-channel acquisitions",
                "Hoechst access validates a perfusion proxy, not administered-antibody delivery",
            ],
        },
        "runtime_seconds": time.time() - started,
        "figures": figures,
    }
    (output_dir / "perfusion_validation.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result
