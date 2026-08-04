from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "reports" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(FIGURES / name, dpi=200, bbox_inches="tight")
    plt.close(fig)


def workflow() -> None:
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.axis("off")

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    top_boxes = [
        (0.6, 4.9, 2.5, 1.45, "Spatial tissue data", "RCC Xenium, H&E,\nprotein imaging"),
        (
            3.8,
            4.9,
            2.5,
            1.45,
            "External validation",
            "Perfusion, HER2 calibration,\nIgG transport priors",
        ),
        (
            7.0,
            4.9,
            2.5,
            1.45,
            "Evidence graph",
            "Same tissue, external,\nmissing and blocked evidence",
        ),
        (
            10.2,
            4.9,
            2.7,
            1.45,
            "Uncertainty analysis",
            "Six vessel definitions and\n20,000 objective-weight draws",
        ),
    ]
    for idx, (x, y, w, h, title, subtitle) in enumerate(top_boxes):
        patch = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.04,rounding_size=0.08",
            linewidth=1.8,
            edgecolor=colors[idx % len(colors)],
            facecolor="white",
        )
        ax.add_patch(patch)
        ax.text(
            x + w / 2, y + 0.98, title, ha="center", va="center", fontsize=13, fontweight="bold"
        )
        ax.text(
            x + w / 2, y + 0.43, subtitle, ha="center", va="center", fontsize=9.5, linespacing=1.3
        )
        if idx < len(top_boxes) - 1:
            nx = top_boxes[idx + 1][0]
            ax.add_patch(
                FancyArrowPatch(
                    (x + w, y + h / 2),
                    (nx, y + h / 2),
                    arrowstyle="-|>",
                    mutation_scale=16,
                    linewidth=1.5,
                )
            )

    ax.text(
        7,
        7.25,
        "reach-gap v0.8.0: from spatial measurements to evidence-aware decisions",
        ha="center",
        fontsize=19,
        fontweight="bold",
    )

    # Decision split
    split_x, split_y = 11.55, 4.9
    ax.add_patch(
        FancyArrowPatch(
            (split_x, split_y), (5.1, 3.25), arrowstyle="-|>", mutation_scale=16, linewidth=1.5
        )
    )
    ax.add_patch(
        FancyArrowPatch(
            (split_x, split_y), (9.9, 3.25), arrowstyle="-|>", mutation_scale=16, linewidth=1.5
        )
    )

    branches = [
        (
            2.5,
            1.0,
            4.8,
            1.9,
            "Relative claim supported",
            (
                "Report target-rank probabilities, pairwise wins\n"
                "and leave-one-component-out sensitivity."
            ),
        ),
        (
            7.5,
            1.0,
            4.8,
            1.9,
            "Absolute claim not supported",
            (
                "Return NOT_COMPUTED, rank the next\n"
                "measurement, and identify what would\n"
                "reduce uncertainty most."
            ),
        ),
    ]
    for idx, (x, y, w, h, title, body) in enumerate(branches):
        patch = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.04,rounding_size=0.08",
            linewidth=1.8,
            edgecolor=colors[(idx + 4) % len(colors)],
            facecolor="white",
        )
        ax.add_patch(patch)
        ax.text(
            x + w / 2, y + 1.32, title, ha="center", va="center", fontsize=13, fontweight="bold"
        )
        ax.text(x + w / 2, y + 0.62, body, ha="center", va="center", fontsize=9.5, linespacing=1.35)

    ax.text(
        7,
        0.35,
        (
            "Research software only — relative geometry-expression evidence is not "
            "therapeutic efficacy or patient guidance."
        ),
        ha="center",
        fontsize=10,
        style="italic",
    )
    save(fig, "reach_gap_workflow.png")


def target_rank() -> None:
    path = (
        ROOT / "results/evidence_synthesis_v0.8/relative_accessibility/target_rank_summary_v0.8.csv"
    )
    df = pd.read_csv(path).sort_values("top_rank_probability")
    fig, ax = plt.subplots(figsize=(10, 5.7))
    bars = ax.barh(df["target"], df["top_rank_probability"] * 100)
    ax.set_xlim(0, 105)
    ax.set_xlabel("Probability of rank 1 across structural and objective uncertainty (%)")
    ax.set_title(
        "Relative RCC target ranking is strongly concentrated on VISTA", fontsize=16, pad=14
    )
    ax.grid(axis="x", alpha=0.25)
    for bar, value in zip(bars, df["top_rank_probability"] * 100, strict=True):
        ax.text(
            min(value + 1.2, 101),
            bar.get_y() + bar.get_height() / 2,
            f"{value:.1f}%",
            va="center",
            fontsize=11,
        )
    fig.subplots_adjust(bottom=0.20)
    fig.text(
        0.5,
        0.035,
        "20,000 draws across six vessel definitions; relative proxy only, not antibody efficacy.",
        ha="center",
        fontsize=9,
        style="italic",
    )
    save(fig, "target_rank_probability.png")


def measurement_priority() -> None:
    path = ROOT / "results/evidence_synthesis_v0.8/evidence/measurement_priorities_v0.8.csv"
    df = pd.read_csv(path).sort_values("uncertainty_share")
    labels = df["requirement"].str.replace("_", " ")
    fig, ax = plt.subplots(figsize=(11, 6))
    bars = ax.barh(labels, df["uncertainty_share"] * 100)
    ax.set_xlabel("Share of unresolved weighted evidence burden (%)")
    ax.set_title("What experiment would reduce uncertainty the most?", fontsize=16, pad=14)
    ax.grid(axis="x", alpha=0.25)
    for bar, value in zip(bars, df["uncertainty_share"] * 100, strict=True):
        ax.text(
            value + 0.5,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.1f}%",
            va="center",
            fontsize=10,
        )
    ax.set_xlim(0, max(df["uncertainty_share"] * 100) * 1.25)
    save(fig, "measurement_priority.png")


def pairwise_heatmap() -> None:
    path = (
        ROOT
        / "results/evidence_synthesis_v0.8/relative_accessibility/pairwise_win_probability_v0.8.csv"
    )
    df = pd.read_csv(path, index_col=0)
    values = df.to_numpy() * 100
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(values, vmin=0, vmax=100, cmap="Blues")
    ax.set_xticks(range(len(df.columns)), df.columns)
    ax.set_yticks(range(len(df.index)), df.index)
    ax.set_xlabel("Opponent")
    ax.set_ylabel("Target")
    ax.set_title("Pairwise probability that the row target ranks higher", fontsize=15, pad=14)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            text_color = "white" if values[i, j] > 55 else "black"
            ax.text(
                j,
                i,
                f"{values[i, j]:.1f}%",
                ha="center",
                va="center",
                color=text_color,
                fontsize=10,
            )
    fig.colorbar(im, ax=ax, label="Win probability (%)", shrink=0.85)
    save(fig, "pairwise_win_probability.png")


def leave_one_out() -> None:
    path = (
        ROOT
        / "results/evidence_synthesis_v0.8/relative_accessibility/leave_one_component_out_v0.8.csv"
    )
    df = pd.read_csv(path)
    pivot = (
        df.pivot(index="omitted_component", columns="target", values="top_rank_probability") * 100
    )
    pivot = pivot[[c for c in ["VISTA", "PD-1", "PD-L1", "LAG-3"] if c in pivot.columns]]
    labels = [s.replace("_", " ") for s in pivot.index]
    x = np.arange(len(labels))
    width = 0.18
    fig, ax = plt.subplots(figsize=(11, 6))
    for idx, target in enumerate(pivot.columns):
        ax.bar(x + (idx - (len(pivot.columns) - 1) / 2) * width, pivot[target], width, label=target)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Probability of rank 1 (%)")
    ax.set_ylim(0, 105)
    ax.set_title(
        "VISTA remains the top relative target after removing any one score component",
        fontsize=15,
        pad=14,
    )
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    ax.grid(axis="y", alpha=0.25)
    save(fig, "leave_one_component_out.png")


def target_geometry() -> None:
    path = ROOT / "results/real_rcc_xenium/target_spatial_summary.csv"
    df = pd.read_csv(path)
    df = df[(df["measurement"] == "protein_intensity") & (df["subset"] == "tumour_region")].copy()
    x = df["median_distance_to_vessel_um_positive"]
    y = df["positive_fraction"] * 100
    sizes = 1000 * df["positive_within_50um"]
    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    ax.scatter(x, y, s=sizes, alpha=0.7)
    offsets = {"VISTA": (8, 7), "PD-1": (-45, 8), "PD-L1": (-55, 8), "LAG-3": (-55, 8)}
    for _, row in df.iterrows():
        ax.annotate(
            row["target"],
            (row["median_distance_to_vessel_um_positive"], row["positive_fraction"] * 100),
            xytext=offsets.get(row["target"], (6, 6)),
            textcoords="offset points",
            fontsize=11,
            fontweight="bold",
        )
    ax.set_xlim(x.min() - 0.18, x.max() + 0.30)
    ax.set_ylim(max(0, y.min() - 2.5), y.max() + 2.0)
    ax.set_xlabel("Median distance of target-positive cells to structural vessel proxy (µm)")
    ax.set_ylabel("Target-positive fraction in pathology-defined tumour region (%)")
    ax.set_title(
        "RCC target geometry combines prevalence and vessel proximity", fontsize=15, pad=14
    )
    ax.grid(alpha=0.25)
    ax.text(
        0.02,
        0.02,
        "Bubble area ∝ fraction of positive cells within 50 µm",
        transform=ax.transAxes,
        fontsize=9,
        style="italic",
    )
    save(fig, "rcc_target_geometry.png")


def validation_summary() -> None:
    build = json.loads((ROOT / "results/build_validation_v0.8.json").read_text())
    coverage = json.loads((ROOT / "results/coverage_v0.8.json").read_text())
    checks = build["checks"]
    rows = [
        ("Tests", f"{checks['pytest']['tests_passed']} passed"),
        (
            "Coverage",
            (
                f"{coverage['coverage_percent']:.2f}% "
                f"(threshold {coverage['threshold_percent']:.0f}%)"
            ),
        ),
        ("Ruff", checks["ruff_check"]["status"]),
        (
            "Mypy strict",
            (
                f"{checks['mypy_full_strict']['status']} · "
                f"{checks['mypy_full_strict']['source_modules']} modules"
            ),
        ),
        (
            "Pyright v0.8",
            (
                f"{checks['pyright_v080_modules']['errors']} errors · "
                f"{checks['pyright_v080_modules']['warnings']} warnings"
            ),
        ),
        (
            "Result invariants",
            (
                f"{checks['result_bundle_validation']['checks']} checks · "
                f"{checks['result_bundle_validation']['issues']} issues"
            ),
        ),
        ("Clean sdist test", f"{checks['sdist_clean_test']['tests_passed']} passed"),
        ("Wheel install", checks["clean_wheel_install"]["status"]),
    ]
    fig, ax = plt.subplots(figsize=(10, 5.8))
    ax.axis("off")
    ax.set_title("v0.8.0 validation summary", fontsize=18, pad=18, fontweight="bold")
    table = ax.table(
        cellText=rows,
        colLabels=["Quality gate", "Result"],
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=[0.35, 0.55],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.7)
    for (r, _c), cell in table.get_celld().items():
        if r == 0:
            cell.set_text_props(weight="bold")
    ax.text(
        0.5,
        0.03,
        (
            "All release checks were executed on the packaged source; absolute RCC "
            "reachability remains intentionally NOT_COMPUTED."
        ),
        ha="center",
        fontsize=9,
        style="italic",
        transform=ax.transAxes,
    )
    save(fig, "validation_summary.png")


def external_validation_panel() -> None:
    panels = [
        (
            ROOT
            / "results/external_validation/her2_receptor_calibration_mckinski2026"
            / "figures/her2_source_protocol_calibration.png",
            "Source-protocol HER2 calibration\nR² = 0.9906",
        ),
        (
            ROOT
            / "results/external_validation/v0.7.1/perfusion_s_biad3159_all4"
            / "Upper_Left_hoechst_distance_profile.png",
            "Independent perfusion-proxy profile\nHoechst signal versus CD31 distance",
        ),
        (
            ROOT
            / "results/external_validation/v0.7.1/breast_xenium_erbb2"
            / "erbb2_tumor_group_sample_comparison.png",
            "Breast Xenium ERBB2 control\n20.623x descriptive tumour-cell mean difference",
        ),
        (
            ROOT / "results/real_rcc_xenium/figures/rna_protein_concordance.png",
            "RCC RNA-protein concordance\nRNA remains distinct from measured protein",
        ),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, (path, title) in zip(axes.flat, panels, strict=True):
        ax.imshow(plt.imread(path))
        ax.set_title(title, fontsize=12, pad=10)
        ax.axis("off")
    fig.suptitle(
        "Independent component validation across imaging and spatial-omics modalities",
        fontsize=18,
        fontweight="bold",
        y=0.99,
    )
    fig.text(
        0.5,
        0.015,
        (
            "External validations constrain methods and priors; they are not silently "
            "transferred into absolute RCC reachability."
        ),
        ha="center",
        fontsize=10,
        style="italic",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    save(fig, "external_validation_panel.png")


def rcc_target_maps_panel() -> None:
    targets = ["VISTA", "PD_1", "PD_L1", "LAG_3"]
    titles = ["VISTA", "PD-1", "PD-L1", "LAG-3"]
    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    for ax, target, title in zip(axes.flat, targets, titles, strict=True):
        path = ROOT / f"results/real_rcc_xenium/figures/target_map_{target}.png"
        ax.imshow(plt.imread(path))
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.axis("off")
    fig.suptitle(
        "Spatial heterogeneity of target-positive fractions across the RCC section",
        fontsize=18,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.01,
        (
            "Within-section protein thresholds; maps are not clinical positivity or "
            "administered-antibody concentration."
        ),
        ha="center",
        fontsize=10,
        style="italic",
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(FIGURES / "rcc_target_maps_panel.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    global FIGURES
    parser = argparse.ArgumentParser(description="Regenerate the portfolio figures.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=FIGURES,
        help="Directory for the nine generated PNG files.",
    )
    args = parser.parse_args()
    FIGURES = args.output_dir.resolve()
    FIGURES.mkdir(parents=True, exist_ok=True)

    workflow()
    target_rank()
    measurement_priority()
    pairwise_heatmap()
    leave_one_out()
    target_geometry()
    validation_summary()
    external_validation_panel()
    rcc_target_maps_panel()

    generated = sorted(FIGURES.glob("*.png"))
    if len(generated) != 9:
        raise RuntimeError(f"Expected 9 portfolio figures, found {len(generated)}")
    print(f"Wrote {len(generated)} figures to {FIGURES}")


if __name__ == "__main__":
    main()
