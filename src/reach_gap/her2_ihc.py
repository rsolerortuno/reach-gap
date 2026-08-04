"""HER2 IHC score benchmark without pretending to calibrate receptor copies."""

from __future__ import annotations

import itertools
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def _her2_image_record(payload: dict[str, Any]) -> dict[str, Any]:
    for record_any in payload["images"]:
        record = dict(record_any)
        if str(record.get("marker_name", "")).casefold() == "her2/neu":
            return record
    raise KeyError("HER2/neu image record is missing")


def _read_rgb_and_mask(
    root: Path, image_record: dict[str, Any], downsample: int
) -> tuple[np.ndarray, np.ndarray]:
    try:
        from PIL import Image
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("HER2 IHC validation requires Pillow") from exc
    image_path = root / str(image_record["file_path"])
    requested_mask = root / str(image_record["mask_path"])
    if requested_mask.exists():
        mask_path = requested_mask
    else:
        candidates = {
            candidate.name.casefold(): candidate
            for candidate in requested_mask.parent.glob("*")
            if candidate.is_file()
        }
        mask_path = candidates.get(requested_mask.name.casefold(), requested_mask)
    with Image.open(image_path) as image_handle:
        rgb = np.asarray(image_handle.convert("RGB"))[::downsample, ::downsample]
    with Image.open(mask_path) as mask_handle:
        mask = np.asarray(mask_handle.convert("L"))[::downsample, ::downsample] > 0
    if rgb.shape[:2] != mask.shape:
        raise ValueError(f"Image/mask shape mismatch: {image_path} versus {mask_path}")
    return rgb, mask


def image_features(rgb: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    """Extract transparent DAB/brown relative intensity features inside a supplied mask."""

    try:
        from skimage.color import rgb2hed
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("HER2 IHC validation requires scikit-image") from exc
    if not np.any(mask):
        raise ValueError("HER2 mask is empty")
    floating = rgb.astype(np.float64) / 255.0
    dab = rgb2hed(floating)[..., 2]
    # The scikit-image HED transform returns positive DAB optical density here.
    dab_strength = np.maximum(dab, 0.0)
    red = floating[..., 0]
    green = floating[..., 1]
    blue = floating[..., 2]
    brown = np.maximum((red - blue) + 0.5 * (green - blue), 0.0)
    return {
        "dab_median": float(np.median(dab_strength[mask])),
        "dab_q90": float(np.quantile(dab_strength[mask], 0.90)),
        "dab_q95": float(np.quantile(dab_strength[mask], 0.95)),
        "brown_median": float(np.median(brown[mask])),
        "brown_q90": float(np.quantile(brown[mask], 0.90)),
        "mask_fraction": float(mask.mean()),
    }


def load_her2_cases(root: Path, downsample: int = 4) -> pd.DataFrame:
    """Parse every case and preserve the dataset's mask semantics explicitly."""

    rows: list[dict[str, Any]] = []
    for json_path in sorted(root.glob("ID*.json")):
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        assessment = dict(payload["diagnostic_assessment"]["HER2/neu"])
        image_record = _her2_image_record(payload)
        rgb, mask = _read_rgb_and_mask(root, image_record, downsample)
        features = image_features(rgb, mask)
        mask_type = str(image_record["mask_type"])
        rows.append(
            {
                "case_id": str(payload["case_id"]),
                "her2_score": int(assessment["her2_score"]),
                "her2_result": str(assessment["result"]),
                "positive_cell_percentage_reported": assessment.get("positive_cell_percentage"),
                "mask_type": mask_type,
                "mask_is_denominator_valid": mask_type == "tumor_cells",
                "magnification": str(image_record["magnification"]),
                **features,
            }
        )
    if not rows:
        raise FileNotFoundError(f"No ID*.json cases found under {root}")
    return pd.DataFrame(rows)


def feature_concordance(table: pd.DataFrame) -> pd.DataFrame:
    """Compare relative image intensity features with pathologist HER2 ordinal scores."""

    rows: list[dict[str, Any]] = []
    score = table["her2_score"].to_numpy(dtype=np.float64)
    for feature in ["dab_median", "dab_q90", "dab_q95", "brown_median", "brown_q90"]:
        values = table[feature].to_numpy(dtype=np.float64)
        statistic = np.nan if np.unique(values).size < 2 else spearmanr(score, values).statistic
        rows.append({"feature": feature, "spearman_vs_her2_score": float(statistic)})
    return pd.DataFrame(rows).sort_values("spearman_vs_her2_score", ascending=False)


def max_statistic_permutation_pvalue(table: pd.DataFrame, features: list[str]) -> tuple[float, int]:
    """Exact multiset-permutation p-value correcting for feature selection.

    The tested statistic is the maximum absolute Spearman correlation over the
    predeclared feature list.  Enumeration is capped to avoid accidental
    combinatorial explosions on future datasets.
    """

    score = table["her2_score"].to_numpy(dtype=np.int64)
    if len(score) < 4 or len(np.unique(score)) < 2:
        return np.nan, 0
    counts = {int(value): int(np.sum(score == value)) for value in np.unique(score)}
    permutations = math.factorial(len(score))
    for count in counts.values():
        permutations //= math.factorial(count)
    if permutations > 100_000:
        return np.nan, permutations

    values = [table[feature].to_numpy(dtype=np.float64) for feature in features]
    observed = max(
        abs(float(spearmanr(score, feature_values).statistic))
        for feature_values in values
        if np.unique(feature_values).size >= 2
    )

    labels = sorted(counts)
    baseline: np.ndarray = np.full(len(score), labels[0], dtype=np.int64)
    assignments: list[np.ndarray] = [baseline]
    # Build unique assignments recursively from combinations of available positions.
    for label in labels[1:]:
        count = counts[label]
        expanded: list[np.ndarray] = []
        for assignment in assignments:
            available = np.flatnonzero(assignment == labels[0])
            for chosen in itertools.combinations(available.tolist(), count):
                candidate = assignment.copy()
                candidate[list(chosen)] = label
                expanded.append(candidate)
        assignments = expanded

    exceed = 0
    evaluated = 0
    for permuted in assignments:
        statistic = max(
            abs(float(spearmanr(permuted, feature_values).statistic))
            for feature_values in values
            if np.unique(feature_values).size >= 2
        )
        exceed += int(statistic >= observed - np.finfo(np.float64).eps)
        evaluated += 1
    return exceed / evaluated, evaluated


def _plot_her2(table: pd.DataFrame, output_dir: Path) -> list[str]:
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover
        return []
    output_dir.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(8, 5))
    plt.scatter(table["her2_score"], table["brown_median"])
    for _, row in table.iterrows():
        plt.annotate(str(row["case_id"]), (row["her2_score"], row["brown_median"]), fontsize=7)
    plt.xlabel("Pathologist HER2 score")
    plt.ylabel("Median relative brown signal")
    plt.title("AHIHCI HER2 image-score benchmark")
    plt.tight_layout()
    path = output_dir / "her2_score_vs_brown_signal.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return [str(path)]


def benchmark_her2_ihc(root: Path, output_dir: Path, downsample: int = 4) -> dict[str, Any]:
    """Run a small ordinal IHC benchmark and enforce the calibration abstention."""

    started = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)
    table = load_her2_cases(root, downsample=downsample)
    concordance = feature_concordance(table)
    table.to_csv(output_dir / "her2_ihc_case_features.csv", index=False)
    concordance.to_csv(output_dir / "her2_ihc_feature_concordance.csv", index=False)
    figures = _plot_her2(table, output_dir / "figures")
    best = concordance.iloc[0]
    candidate_features = ["dab_median", "dab_q90", "dab_q95", "brown_median", "brown_q90"]
    max_statistic_pvalue, permutations_evaluated = max_statistic_permutation_pvalue(
        table, candidate_features
    )
    denominator_valid = table[table["mask_is_denominator_valid"]]
    denominator_valid_spearman = (
        float(
            spearmanr(
                denominator_valid["her2_score"],
                denominator_valid[str(best["feature"])],
            ).statistic
        )
        if len(denominator_valid) >= 3
        else np.nan
    )
    result = {
        "status": "IHC_SCORE_BENCHMARK_ONLY_NOT_MOLECULE_CALIBRATION",
        "cases": len(table),
        "scores_present": sorted(int(value) for value in table["her2_score"].unique()),
        "best_feature": str(best["feature"]),
        "best_spearman_vs_her2_score": float(best["spearman_vs_her2_score"]),
        "max_statistic_exact_permutation_pvalue": max_statistic_pvalue,
        "permutations_evaluated": permutations_evaluated,
        "denominator_valid_only_best_feature_spearman": denominator_valid_spearman,
        "positive_only_masks": int((table["mask_type"] == "only_positive_tumor_cells").sum()),
        "denominator_valid_masks": int(table["mask_is_denominator_valid"].sum()),
        "absolute_antigen_density": {
            "status": "NOT_COMPUTED",
            "reasons": [
                "Ordinal HER2 IHC score is not receptors per cell",
                "No calibration beads or antibodies-bound-per-cell curve is supplied",
                "Positive-only masks cannot define a tumour-cell denominator",
            ],
        },
        "runtime_seconds": time.time() - started,
        "figures": figures,
    }
    (output_dir / "her2_ihc_benchmark.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result
