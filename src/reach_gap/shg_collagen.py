"""Relative SHG collagen texture features without diffusion calibration."""

from __future__ import annotations

import json
import time
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def shg_features(image: np.ndarray) -> dict[str, float]:
    """Extract transparent relative collagen density and orientation features."""

    try:
        from skimage.feature import structure_tensor, structure_tensor_eigenvalues
        from skimage.filters import threshold_otsu
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("SHG validation requires scikit-image") from exc
    values = np.asarray(image, dtype=np.float64)
    if values.ndim > 2:
        values = np.max(values, axis=tuple(range(values.ndim - 2)))
    if values.size == 0:
        raise ValueError("SHG image is empty")
    q99 = float(np.quantile(values, 0.99))
    scaled = np.clip(values / max(q99, 1.0), 0.0, 1.0)
    threshold = (
        float(threshold_otsu(scaled))  # type: ignore[no-untyped-call]
        if np.unique(scaled).size > 1
        else 0.0
    )
    positive = scaled > threshold
    a_rr, a_rc, a_cc = structure_tensor(  # type: ignore[no-untyped-call]
        scaled, sigma=2, order="rc"
    )
    eigenvalues = structure_tensor_eigenvalues(  # type: ignore[no-untyped-call]
        (a_rr, a_rc, a_cc)
    )
    major = eigenvalues[0]
    minor = eigenvalues[1]
    coherence = (major - minor) / np.maximum(major + minor, np.finfo(np.float64).eps)
    active = major > np.quantile(major, 0.50)
    return {
        "intensity_mean_relative": float(np.mean(scaled)),
        "intensity_q90_relative": float(np.quantile(scaled, 0.90)),
        "otsu_threshold_relative": threshold,
        "collagen_positive_fraction": float(np.mean(positive)),
        "orientation_coherence_median": float(np.median(coherence[active])),
        "orientation_coherence_q90": float(np.quantile(coherence[active], 0.90)),
    }


def _read_tiff(path: Path) -> np.ndarray:
    try:
        import tifffile
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("SHG validation requires tifffile") from exc
    return np.asarray(tifffile.imread(path))


def _plot_shg(table: pd.DataFrame, output_dir: Path) -> list[str]:
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover
        return []
    output_dir.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(8, 5))
    groups = [
        group["collagen_positive_fraction"].to_numpy(dtype=np.float64)
        for _, group in table.groupby("tissue_class", sort=True)
    ]
    labels = [str(name) for name, _ in table.groupby("tissue_class", sort=True)]
    plt.boxplot(groups, tick_labels=labels)
    plt.ylabel("Relative SHG-positive fraction")
    plt.title("SHG collagen feature pilot")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    path = output_dir / "shg_relative_collagen_fraction.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return [str(path)]


def benchmark_shg_collagen(
    images: Mapping[str, Sequence[Path]], output_dir: Path
) -> dict[str, Any]:
    """Benchmark the feature extractor across labelled SHG image groups."""

    started = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for tissue_class, paths in images.items():
        for path in paths:
            rows.append(
                {
                    "tissue_class": tissue_class,
                    "image": path.name,
                    **shg_features(_read_tiff(path)),
                }
            )
    table = pd.DataFrame(rows)
    if table.empty:
        raise ValueError("No SHG images were supplied")
    table.to_csv(output_dir / "shg_collagen_features.csv", index=False)
    grouped = (
        table.groupby("tissue_class")
        .agg(
            images=("image", "count"),
            collagen_positive_fraction_median=("collagen_positive_fraction", "median"),
            orientation_coherence_median=("orientation_coherence_median", "median"),
        )
        .reset_index()
    )
    grouped.to_csv(output_dir / "shg_collagen_group_summary.csv", index=False)
    figures = _plot_shg(table, output_dir / "figures")
    result = {
        "status": "SHG_FEATURE_EXTRACTOR_PILOT_NOT_DIFFUSIVITY_CALIBRATION",
        "images": len(table),
        "tissue_classes": sorted(table["tissue_class"].unique().tolist()),
        "transport_coefficient": {
            "status": "NOT_COMPUTED",
            "reasons": [
                "SHG intensity and orientation are not an antibody diffusion coefficient",
                "The pilot images are not registered to the RCC Xenium tissue",
                "No FRAP or tracer transport measurement is supplied",
            ],
        },
        "runtime_seconds": time.time() - started,
        "figures": figures,
    }
    (output_dir / "shg_collagen_benchmark.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result
