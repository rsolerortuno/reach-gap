"""Audit Bordeau et al. administered-trastuzumab supplementary evidence.

The supplementary DOCX contains representative compressed composite figures
rather than raw section-level microscopy.  This module extracts the figures,
quantifies transparent relative features, and records whether the published
penetration direction is reproduced robustly.  It never treats the embedded
figures as raw calibrated microscopy.
"""

from __future__ import annotations

import itertools
import json
import math
import time
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from PIL import Image
from scipy.ndimage import distance_transform_edt

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]

PUBLISHED_ENDPOINTS = {
    "trastuzumab_alone_threshold_positive_penetration_um": 41.30,
    "trastuzumab_alone_threshold_positive_penetration_sd_um": 6.70,
    "trastuzumab_plus_1he_threshold_positive_penetration_um": 58.24,
    "trastuzumab_plus_1he_threshold_positive_penetration_sd_um": 5.40,
    "trastuzumab_alone_penetration_limit_um": 68.60,
    "trastuzumab_plus_1he_penetration_limit_um": 64.54,
    "dose_mg_per_kg": 2.0,
    "time_hours": 24.0,
}


def extract_docx_media(docx_path: Path, output_dir: Path) -> list[Path]:
    """Extract image media from a DOCX archive."""

    output_dir.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []
    with zipfile.ZipFile(docx_path) as archive:
        for name in archive.namelist():
            if not name.startswith("word/media/"):
                continue
            destination = output_dir / Path(name).name
            destination.write_bytes(archive.read(name))
            extracted.append(destination)
    if not extracted:
        raise ValueError("DOCX contains no media images")
    return sorted(extracted)


def _white_runs(values: FloatArray, threshold: float = 0.95) -> list[tuple[int, int]]:
    mask = values > threshold
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for index, selected in enumerate(mask):
        if selected and start is None:
            start = index
        if start is not None and (not selected or index == len(mask) - 1):
            end = index if not selected else index + 1
            runs.append((start, end))
            start = None
    return runs


def split_three_by_three_figure(rgb: NDArray[np.uint8]) -> list[list[NDArray[np.uint8]]]:
    """Split a white-gutter 3x3 supplementary figure into panels."""

    white = np.asarray(np.mean(rgb, axis=2) > 250, dtype=np.float64)
    column_runs = _white_runs(np.mean(white, axis=0))
    row_runs = _white_runs(np.mean(white, axis=1))
    if len(column_runs) < 4 or len(row_runs) < 4:
        raise ValueError("Could not detect the 3x3 white-gutter panel layout")
    x_ranges = [(column_runs[index][1], column_runs[index + 1][0]) for index in range(3)]
    y_ranges = [(row_runs[index][1], row_runs[index + 1][0]) for index in range(3)]
    return [[np.asarray(rgb[y0:y1, x0:x1, :3]) for x0, x1 in x_ranges] for y0, y1 in y_ranges]


def _relative_top_panel_features(panel: NDArray[np.uint8]) -> dict[str, float]:
    image = np.asarray(panel, dtype=np.float64)
    red = image[..., 0]
    green = image[..., 1]
    blue = image[..., 2]
    red_positive = red[red > 0]
    green_positive = green[green > 0]
    if red_positive.size == 0 or green_positive.size == 0:
        raise ValueError("Panel lacks red vascular or green antibody signal")

    red_threshold = max(5.0, float(np.quantile(red_positive, 0.60)))
    vessel = np.asarray(
        (red >= red_threshold) & (red > 1.10 * green) & (red > 1.10 * blue), dtype=np.bool_
    )
    if not np.any(vessel):
        raise ValueError("No vascular pixels were detected")
    green_scale = float(np.quantile(green_positive, 0.99))
    relative_green = np.clip(green / green_scale, 0.0, 1.0)
    antibody = np.asarray(relative_green > 0.20, dtype=np.bool_)
    if not np.any(antibody):
        raise ValueError("No antibody-positive pixels were detected")

    distance_pixels = np.asarray(distance_transform_edt(~vessel), dtype=np.float64)
    diagonal = math.hypot(*panel.shape[:2])
    normalized_distance = distance_pixels[antibody] / diagonal
    return {
        "vascular_fraction": float(np.mean(vessel)),
        "antibody_positive_fraction": float(np.mean(antibody)),
        "antibody_distance_median_fraction_of_diagonal": float(np.median(normalized_distance)),
        "antibody_distance_mean_fraction_of_diagonal": float(np.mean(normalized_distance)),
        "antibody_distance_q90_fraction_of_diagonal": float(np.quantile(normalized_distance, 0.90)),
    }


def _exact_group_difference_pvalue(values: NDArray[np.float64]) -> float:
    observed = float(np.mean(values[3:]) - np.mean(values[:3]))
    extreme = 0
    total = 0
    for selected in itertools.combinations(range(6), 3):
        group_b = np.asarray(selected, dtype=np.int64)
        group_a = np.asarray([index for index in range(6) if index not in selected], dtype=np.int64)
        difference = float(np.mean(values[group_b]) - np.mean(values[group_a]))
        total += 1
        extreme += int(abs(difference) >= abs(observed) - 1e-15)
    return extreme / total


def benchmark_bordeau_supplement(docx_path: Path, output_dir: Path) -> dict[str, Any]:
    """Curate published endpoints and test representative-panel direction robustness."""

    started = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)
    media = extract_docx_media(docx_path, output_dir / "extracted_media")
    figure_candidates: list[tuple[Path, NDArray[np.uint8]]] = []
    for path in media:
        with Image.open(path) as image:
            rgb = np.asarray(image.convert("RGB"))
        try:
            split_three_by_three_figure(rgb)
        except ValueError:
            continue
        figure_candidates.append((path, rgb))
    if len(figure_candidates) < 2:
        raise ValueError("Could not identify both representative tumour-section figures")

    alone_candidates = [item for item in figure_candidates if "image2" in item[0].stem.casefold()]
    combo_candidates = [item for item in figure_candidates if "image3" in item[0].stem.casefold()]
    if alone_candidates and combo_candidates:
        alone_path, alone_rgb = alone_candidates[0]
        combo_path, combo_rgb = combo_candidates[0]
    else:
        # Source Figure S2 is the smaller JPEG and Figure S3 is the larger PNG.
        alone_path, alone_rgb = min(figure_candidates, key=lambda item: item[1].shape[0])
        combo_path, combo_rgb = max(figure_candidates, key=lambda item: item[1].shape[0])
    rows: list[dict[str, Any]] = []
    for group, source_path, rgb in (
        ("trastuzumab_alone", alone_path, alone_rgb),
        ("trastuzumab_plus_1he", combo_path, combo_rgb),
    ):
        panels = split_three_by_three_figure(rgb)
        for section_index, panel in enumerate(panels[0], start=1):
            rows.append(
                {
                    "group": group,
                    "section": section_index,
                    "source_image": source_path.name,
                    **_relative_top_panel_features(panel),
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(output_dir / "bordeau_representative_panel_features.csv", index=False)
    metric = "antibody_distance_mean_fraction_of_diagonal"
    ordered = pd.concat(
        [
            table.loc[table["group"] == "trastuzumab_alone", metric],
            table.loc[table["group"] == "trastuzumab_plus_1he", metric],
        ],
        ignore_index=True,
    ).to_numpy(dtype=np.float64)
    alone_mean = float(np.mean(ordered[:3]))
    combo_mean = float(np.mean(ordered[3:]))
    direction_matches = combo_mean > alone_mean
    pvalue = _exact_group_difference_pvalue(ordered)

    result: dict[str, Any] = {
        "status": "PUBLISHED_ADMINISTERED_TRASTUZUMAB_ENDPOINTS_CURATED_REPRESENTATIVE_FIGURE_ONLY",
        "published_endpoints": PUBLISHED_ENDPOINTS,
        "representative_figure_analysis": {
            "sections_per_group": 3,
            "metric": metric,
            "trastuzumab_alone_mean": alone_mean,
            "trastuzumab_plus_1he_mean": combo_mean,
            "direction_matches_published_threshold_positive_penetration": direction_matches,
            "exact_two_sided_permutation_pvalue": pvalue,
            "interpretation": "DESCRIPTIVE_ONLY_COMPRESSED_REPRESENTATIVE_FIGURES",
        },
        "model_concordance": {
            "status": "NOT_COMPUTED",
            "reasons": [
                "Raw section-level microscopy and animal-level measurements are not deposited",
                "The embedded supplementary figures are compressed and may be rescaled",
                "No reach-gap prediction was prospectively fixed for the SKOV3 experiment",
            ],
        },
        "runtime_seconds": time.time() - started,
    }
    (output_dir / "bordeau_supplement_benchmark.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result
