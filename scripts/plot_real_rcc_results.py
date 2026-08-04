"""Generate compact figures from committed real RCC Xenium summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results" / "real_rcc_xenium",
    )
    args = parser.parse_args()
    out = args.results_dir
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    full_parts = sorted((out / "cell_tables").glob("*.csv.gz"))
    if full_parts:
        scored = pd.concat([pd.read_csv(path) for path in full_parts], ignore_index=True)
    else:
        scored = pd.read_csv(out / "scored_cells_sample.csv.gz")
    targets = ["PD-L1", "VISTA", "PD-1", "LAG-3"]
    for target in targets:
        safe = target.replace("-", "_")
        positive = scored[f"target__{safe}__positive"].astype(float).to_numpy()
        fig, ax = plt.subplots(figsize=(8.5, 4.5))
        hb = ax.hexbin(
            scored["x_um"],
            scored["y_um"],
            C=positive,
            reduce_C_function=np.mean,
            gridsize=260,
            mincnt=3,
        )
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xlabel("X (µm)")
        ax.set_ylabel("Y (µm)")
        ax.set_title(f"{target}: within-section protein-positive fraction")
        cbar = fig.colorbar(hb, ax=ax)
        cbar.set_label("Positive-cell fraction per spatial bin")
        fig.tight_layout()
        fig.savefig(figures / f"target_map_{safe}.png", dpi=180)
        plt.close(fig)

    concordance = pd.read_csv(out / "rna_protein_concordance.csv")
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.bar(concordance["target"], concordance["spearman_raw_signal"])
    ax.axhline(0.0, linewidth=0.8)
    ax.set_ylim(0.0, max(0.30, concordance["spearman_raw_signal"].max() * 1.15))
    ax.set_ylabel("Spearman correlation")
    ax.set_title("Cell-level RNA–protein concordance")  # noqa: RUF001
    for index, value in enumerate(concordance["spearman_raw_signal"]):
        ax.text(index, value + 0.008, f"{value:.3f}", ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(figures / "rna_protein_concordance.png", dpi=180)
    plt.close(fig)

    ranges = pd.read_csv(out / "vessel_robustness_ranges.csv")
    ranges = ranges.set_index("target").loc[targets].reset_index()
    centre = (ranges["target_median_um_min"] + ranges["target_median_um_max"]) / 2.0
    lower = centre - ranges["target_median_um_min"]
    upper = ranges["target_median_um_max"] - centre
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.errorbar(ranges["target"], centre, yerr=np.vstack([lower, upper]), fmt="o", capsize=5)
    ax.set_ylabel("Median distance to vessel proxy (µm)")
    ax.set_title("Sensitivity to six vessel definitions")
    fig.tight_layout()
    fig.savefig(figures / "vessel_definition_sensitivity.png", dpi=180)
    plt.close(fig)

    segmentation = pd.read_csv(out / "segmentation_boundary_metrics_sample.csv.gz")
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for column, label in [
        ("cell_area_relative_error", "Cell"),
        ("nucleus_area_relative_error", "Nucleus"),
    ]:
        values = np.sort(segmentation[column].dropna().to_numpy())
        y = np.arange(1, len(values) + 1) / len(values)
        ax.plot(values, y, label=label)
    ax.set_xlim(0.0, 0.25)
    ax.set_xlabel("Relative polygon-area difference")
    ax.set_ylabel("Empirical cumulative fraction")
    ax.set_title("Segmentation boundary robustness")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / "segmentation_area_robustness.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
