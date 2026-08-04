"""Relative geometry-expression accessibility analysis for measured RCC targets."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from pydantic import BaseModel, ConfigDict, Field

from reach_gap.visuals import write_horizontal_bar_svg

FloatArray = NDArray[np.float64]


class TargetRankSummary(BaseModel):
    """Monte Carlo relative-rank summary for one target."""

    model_config = ConfigDict(extra="forbid")

    target: str
    tumour_positive_fraction: float = Field(ge=0.0, le=1.0)
    top_rank_probability: float = Field(ge=0.0, le=1.0)
    median_rank: float = Field(ge=1.0)
    rank_q05: float = Field(ge=1.0)
    rank_q95: float = Field(ge=1.0)
    proxy_score_median: float = Field(ge=0.0, le=1.0)
    proxy_score_q05: float = Field(ge=0.0, le=1.0)
    proxy_score_q95: float = Field(ge=0.0, le=1.0)


class RelativeAccessibilityResult(BaseModel):
    """Explicitly relative target ranking with no absolute-drug interpretation."""

    model_config = ConfigDict(extra="forbid")

    version: str
    status: str
    target_count: int
    vessel_definition_count: int
    draws: int
    seed: int
    components: list[str]
    weight_sampling: str
    vessel_definition_sampling: str
    target_summaries: list[TargetRankSummary]
    pairwise_win_probability: dict[str, dict[str, float]]
    leave_one_component_out_top_probability: dict[str, dict[str, float]]
    stable_top_target: str | None
    stable_top_probability: float
    permitted_claims: list[str]
    unsupported_claims: list[str]


def _minmax(values: FloatArray, *, higher_is_better: bool = True) -> FloatArray:
    low = float(np.min(values))
    high = float(np.max(values))
    if high - low <= 1.0e-15:
        scaled = np.full_like(values, 0.5, dtype=np.float64)
    else:
        scaled = np.asarray((values - low) / (high - low), dtype=np.float64)
    return scaled if higher_is_better else 1.0 - scaled


def _load_inputs(repository_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    target_path = repository_root / "results/real_rcc_xenium/target_spatial_summary.csv"
    sensitivity_path = repository_root / "results/real_rcc_xenium/vessel_calling_sensitivity.csv"
    targets = pd.read_csv(target_path)
    sensitivity = pd.read_csv(sensitivity_path)
    required_targets = {
        "target",
        "subset",
        "positive_fraction",
        "measurement",
    }
    required_sensitivity = {
        "definition",
        "target",
        "target_median_um",
        "within_50um",
    }
    if not required_targets.issubset(targets.columns):
        raise ValueError(
            f"Missing target columns: {sorted(required_targets - set(targets.columns))}"
        )
    if not required_sensitivity.issubset(sensitivity.columns):
        missing = sorted(required_sensitivity - set(sensitivity.columns))
        raise ValueError(f"Missing sensitivity columns: {missing}")
    tumour = targets.loc[
        (targets["subset"] == "tumour_region") & (targets["measurement"] == "protein_intensity"),
        ["target", "positive_fraction"],
    ].copy()
    if tumour["target"].duplicated().any():
        raise ValueError("Tumour target rows must be unique")
    merged = sensitivity.merge(tumour, on="target", how="inner", validate="many_to_one")
    if merged.empty:
        raise ValueError("No shared RCC target rows were found")
    counts = merged.groupby("definition")["target"].nunique()
    if counts.nunique() != 1:
        raise ValueError("Every vessel definition must cover the same targets")
    if merged[["target_median_um", "within_50um", "positive_fraction"]].isna().any().any():
        raise ValueError("Relative accessibility inputs contain missing values")
    return tumour.sort_values("target").reset_index(drop=True), merged


def _component_matrix(definition_rows: pd.DataFrame, ordered_targets: list[str]) -> FloatArray:
    rows = definition_rows.set_index("target").loc[ordered_targets]
    coverage = _minmax(rows["positive_fraction"].to_numpy(dtype=np.float64))
    proximity = _minmax(rows["target_median_um"].to_numpy(dtype=np.float64), higher_is_better=False)
    near_source = _minmax(rows["within_50um"].to_numpy(dtype=np.float64))
    return np.column_stack([coverage, proximity, near_source]).astype(np.float64)


def _rank_from_scores(scores: FloatArray) -> NDArray[np.int64]:
    order = np.argsort(-scores, kind="stable")
    ranks = np.empty(len(scores), dtype=np.int64)
    ranks[order] = np.arange(1, len(scores) + 1, dtype=np.int64)
    return ranks


def _simulate(
    matrices: list[FloatArray],
    *,
    draws: int,
    seed: int,
    active_components: tuple[int, ...] = (0, 1, 2),
) -> tuple[FloatArray, NDArray[np.int64]]:
    if draws < 100:
        raise ValueError("At least 100 draws are required")
    if not active_components:
        raise ValueError("At least one component is required")
    rng = np.random.default_rng(seed)
    target_count = matrices[0].shape[0]
    scores = np.empty((draws, target_count), dtype=np.float64)
    ranks = np.empty((draws, target_count), dtype=np.int64)
    for draw in range(draws):
        matrix = matrices[int(rng.integers(0, len(matrices)))][:, active_components]
        weights = rng.dirichlet(np.ones(len(active_components), dtype=np.float64))
        draw_scores = matrix @ weights
        scores[draw] = draw_scores
        ranks[draw] = _rank_from_scores(draw_scores)
    return scores, ranks


def compute_relative_accessibility(
    repository_root: Path,
    *,
    draws: int = 20_000,
    seed: int = 17,
    stable_top_threshold: float = 0.80,
) -> RelativeAccessibilityResult:
    """Rank targets under vessel-definition and objective-weight uncertainty."""

    if not 0.5 <= stable_top_threshold <= 1.0:
        raise ValueError("stable_top_threshold must lie in [0.5, 1]")
    tumour, sensitivity = _load_inputs(repository_root)
    targets = tumour["target"].astype(str).tolist()
    definitions = sorted(sensitivity["definition"].astype(str).unique())
    matrices = [
        _component_matrix(sensitivity.loc[sensitivity["definition"] == definition], targets)
        for definition in definitions
    ]
    scores, ranks = _simulate(matrices, draws=draws, seed=seed)
    top = ranks == 1
    top_probability = top.mean(axis=0)
    target_summaries: list[TargetRankSummary] = []
    fractions = tumour.set_index("target")["positive_fraction"]
    for index, target in enumerate(targets):
        target_summaries.append(
            TargetRankSummary(
                target=target,
                tumour_positive_fraction=float(fractions.loc[target]),
                top_rank_probability=float(top_probability[index]),
                median_rank=float(np.median(ranks[:, index])),
                rank_q05=float(np.quantile(ranks[:, index], 0.05)),
                rank_q95=float(np.quantile(ranks[:, index], 0.95)),
                proxy_score_median=float(np.median(scores[:, index])),
                proxy_score_q05=float(np.quantile(scores[:, index], 0.05)),
                proxy_score_q95=float(np.quantile(scores[:, index], 0.95)),
            )
        )

    pairwise: dict[str, dict[str, float]] = {}
    for left_index, left in enumerate(targets):
        pairwise[left] = {}
        for right_index, right in enumerate(targets):
            if left == right:
                pairwise[left][right] = 0.5
            else:
                pairwise[left][right] = float(
                    np.mean(scores[:, left_index] > scores[:, right_index])
                )

    component_names = ["tumour_positive_fraction", "median_structural_proximity", "within_50um"]
    leave_one_out: dict[str, dict[str, float]] = {}
    for omitted_index, omitted_name in enumerate(component_names):
        active = tuple(index for index in range(3) if index != omitted_index)
        _, ablated_ranks = _simulate(
            matrices,
            draws=max(2_000, draws // 5),
            seed=seed + omitted_index + 1,
            active_components=active,
        )
        leave_one_out[omitted_name] = {
            target: float(np.mean(ablated_ranks[:, index] == 1))
            for index, target in enumerate(targets)
        }

    winner_index = int(np.argmax(top_probability))
    winner_probability = float(top_probability[winner_index])
    stable_target = targets[winner_index] if winner_probability >= stable_top_threshold else None
    return RelativeAccessibilityResult(
        version="0.8.0",
        status="RELATIVE_RCC_GEOMETRY_EXPRESSION_PROXY_ABSOLUTE_REACHABILITY_NOT_COMPUTED",
        target_count=len(targets),
        vessel_definition_count=len(definitions),
        draws=draws,
        seed=seed,
        components=component_names,
        weight_sampling="Uniform over the three-component simplex via Dirichlet(1,1,1)",
        vessel_definition_sampling="Uniform over six preregistered structural vessel definitions",
        target_summaries=target_summaries,
        pairwise_win_probability=pairwise,
        leave_one_component_out_top_probability=leave_one_out,
        stable_top_target=stable_target,
        stable_top_probability=winner_probability,
        permitted_claims=[
            "Targets may be compared on a relative, within-section geometry-expression proxy.",
            "Rank uncertainty includes vessel-definition choice and unknown objective weights.",
            "A stable top target means robustness within this proxy only.",
        ],
        unsupported_claims=[
            "Absolute antibody reachability or engagement.",
            "Clinical target priority or efficacy.",
            "Functional perfusion of structural CD31 definitions.",
            "Equivalence of relative proxy score and drug concentration.",
            "A Bayesian posterior probability that a target is clinically best.",
            "Comparison of proxy scores across datasets or changing target sets.",
        ],
    )


def write_relative_accessibility(
    result: RelativeAccessibilityResult, output_dir: Path
) -> dict[str, str]:
    """Write JSON, target and pairwise tables, ablation table and a rank figure."""

    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "relative_accessibility_v0.8.json"
    target_path = output_dir / "target_rank_summary_v0.8.csv"
    pairwise_path = output_dir / "pairwise_win_probability_v0.8.csv"
    ablation_path = output_dir / "leave_one_component_out_v0.8.csv"
    figure_path = output_dir / "target_top_rank_probability_v0.8.svg"
    json_path.write_text(result.model_dump_json(indent=2), encoding="utf-8")

    ordered = sorted(
        result.target_summaries,
        key=lambda row: (-row.top_rank_probability, row.target),
    )
    with target_path.open("w", encoding="utf-8", newline="") as handle:
        target_writer = csv.DictWriter(handle, fieldnames=list(TargetRankSummary.model_fields))
        target_writer.writeheader()
        for row in ordered:
            target_writer.writerow(row.model_dump(mode="json"))

    targets = [row.target for row in ordered]
    with pairwise_path.open("w", encoding="utf-8", newline="") as handle:
        pairwise_writer = csv.writer(handle)
        pairwise_writer.writerow(["target", *targets])
        for left in targets:
            pairwise_writer.writerow(
                [left] + [result.pairwise_win_probability[left][right] for right in targets]
            )

    with ablation_path.open("w", encoding="utf-8", newline="") as handle:
        ablation_writer = csv.DictWriter(
            handle,
            fieldnames=["omitted_component", "target", "top_rank_probability"],
        )
        ablation_writer.writeheader()
        for omitted, probabilities in result.leave_one_component_out_top_probability.items():
            for target, probability in sorted(probabilities.items()):
                ablation_writer.writerow(
                    {
                        "omitted_component": omitted,
                        "target": target,
                        "top_rank_probability": probability,
                    }
                )

    write_horizontal_bar_svg(
        [row.target for row in ordered],
        [row.top_rank_probability for row in ordered],
        figure_path,
        title="RCC relative geometry-expression top-rank probability",
        x_label="Probability of rank 1 across vessel and objective uncertainty",
        maximum=1.0,
    )
    return {
        "result": str(json_path),
        "targets": str(target_path),
        "pairwise": str(pairwise_path),
        "ablation": str(ablation_path),
        "figure": str(figure_path),
    }


def run_relative_accessibility(
    repository_root: Path,
    output_dir: Path,
    *,
    draws: int = 20_000,
    seed: int = 17,
) -> dict[str, Any]:
    """Compute and persist the relative target analysis."""

    result = compute_relative_accessibility(repository_root, draws=draws, seed=seed)
    paths = write_relative_accessibility(result, output_dir)
    return {"summary": result.model_dump(mode="json"), "paths": paths}
