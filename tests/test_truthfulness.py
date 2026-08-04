from __future__ import annotations

import json
from pathlib import Path


def test_readme_quantitative_results_are_artifact_backed() -> None:
    root = Path(__file__).parents[1]
    readme = (root / "README.md").read_text(encoding="utf-8")
    benchmark_path = root / "results" / "simulated" / "benchmark.json"
    index_path = root / "results" / "simulated" / "index.json"
    if "COMMITTED_BENCHMARK_RMSE" in readme:
        benchmark = json.loads(benchmark_path.read_text(encoding="utf-8"))
        value = benchmark["model_comparison"]["mechanistic"]["rmse"]
        assert f"{value:.4f}" in readme
    if "COMMITTED_INDEX_GAP" in readme:
        index = json.loads(index_path.read_text(encoding="utf-8"))
        value = index["expression_reach_gap"]["median"]
        assert f"{value:.3f}" in readme


def test_readme_states_external_validation_not_computed() -> None:
    root = Path(__file__).parents[1]
    readme = (root / "README.md").read_text(encoding="utf-8")
    assert "model concordance**: `NOT_COMPUTED`" in readme
    assert "Real clinical retrospective: `NOT_COMPUTED`" in readme


def test_readme_real_he_numbers_are_artifact_backed() -> None:
    root = Path(__file__).parents[1]
    readme = (root / "README.md").read_text(encoding="utf-8")
    if "COMMITTED_REAL_HE_GEOMETRY" not in readme:
        return
    result = json.loads(
        (root / "results" / "real_rcc_he" / "he_pathology_result.json").read_text(encoding="utf-8")
    )
    tissue = result["tissue"]["estimated_area_mm2"]
    tumour = result["annotations"]["by_name"]["Tumor"]["area_mm2"]
    immune = result["annotations"]["by_name"]["Immune infiltration"]["area_mm2"]
    vessels = result["annotations"]["by_name"]["Blood vessels"]["area_mm2"]
    necrosis = result["annotations"]["by_name"]["Necrosis"]["area_mm2"]
    candidates = result["lumen_candidates"]["count"]
    assert f"{tissue:.3f}" in readme
    assert f"{tumour:.3f}" in readme
    assert f"{immune:.3f}" in readme
    assert f"{vessels:.3f}" in readme
    assert f"{necrosis:.3f}" in readme
    assert f"**{candidates}**" in readme


def test_readme_real_xenium_numbers_are_artifact_backed() -> None:
    root = Path(__file__).parents[1]
    readme = (root / "README.md").read_text(encoding="utf-8")
    if "COMMITTED_REAL_XENIUM" not in readme:
        return
    summary = json.loads(
        (root / "results" / "real_rcc_xenium" / "real_data_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert f"**{summary['cells']:,} cells**" in readme
    for target in ("PD-L1", "VISTA", "PD-1", "LAG-3"):
        result = summary["targets"][target]
        assert f"{result['tumour_positive_fraction'] * 100:.2f}%" in readme
        low, high = result["median_distance_range_across_vessel_definitions_um"]
        assert f"{low:.2f}–{high:.2f} µm" in readme  # noqa: RUF001
    for value in summary["rna_protein_spearman"].values():
        assert f"**{value:.3f}**" in readme


def test_real_xenium_claims_abstain_from_absolute_index() -> None:
    root = Path(__file__).parents[1]
    summary = json.loads(
        (root / "results" / "real_rcc_xenium" / "real_data_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["absolute_index"]["status"] == "NOT_COMPUTED"
    claims = json.loads(
        (root / "results" / "real_rcc_xenium" / "claims.json").read_text(encoding="utf-8")
    )
    assert any("reachable_fraction" in statement for statement in claims["unsupported"])


def test_readme_real_protein_imaging_numbers_are_artifact_backed() -> None:
    root = Path(__file__).parents[1]
    readme = (root / "README.md").read_text(encoding="utf-8")
    if "COMMITTED_REAL_PROTEIN_IMAGING" not in readme:
        return
    import pandas as pd

    result_dir = root / "results" / "real_rcc_imaging"
    summary = json.loads((result_dir / "real_imaging_summary.json").read_text(encoding="utf-8"))
    concordance = pd.read_csv(result_dir / "channel_cell_image_concordance.csv")
    vessels = pd.read_csv(result_dir / "image_vessel_definition_sensitivity.csv")
    assert f"**{summary['cells']:,} cells**" in readme
    for channel in ("CD31", "alphaSMA", "Vimentin", "PanCK", "PD-L1", "VISTA"):
        value = float(
            concordance.loc[concordance["channel"] == channel, "cell_image_spearman"].iloc[0]
        )
        assert f"{value:.3f}" in readme
    balanced = vessels[vessels["definition"] == "balanced_fpr_0_5pct"]
    for target in ("PD-L1", "VISTA"):
        value = float(
            balanced.loc[balanced["target"] == target, "target_median_distance_um"].iloc[0]
        )
        assert f"{value:.2f} µm" in readme


def test_real_protein_imaging_claims_abstain_from_absolute_index() -> None:
    root = Path(__file__).parents[1]
    claims = json.loads(
        (root / "results" / "real_rcc_imaging" / "claims.json").read_text(encoding="utf-8")
    )
    assert any("reachable_fraction" in statement for statement in claims["unsupported"])
    assert "FUNCTIONALLY_PERFUSED_VESSELS_NOT_IDENTIFIED" in claims["abstention_reasons"]


def test_readme_external_validation_numbers_are_artifact_backed() -> None:
    root = Path(__file__).parents[1]
    readme = (root / "README.md").read_text(encoding="utf-8")
    if "COMMITTED_EXTERNAL_COSMX" not in readme:
        return
    cosmx = json.loads(
        (
            root / "results/external_validation/cosmx_gse299786/cosmx_external_validation.json"
        ).read_text(encoding="utf-8")
    )
    sensitivity = cosmx["geometry_sensitivity"]
    assert f"**{cosmx['cells']:,} metadata cells**" in readme
    assert f"**{sensitivity['median_definition_distance_ratio']:.2f}×**" in readme  # noqa: RUF001
    assert f"**{sensitivity['maximum_definition_distance_ratio']:.2f}×**" in readme  # noqa: RUF001
    assert f"**{sensitivity['median_pairwise_balanced_target_rank_spearman']:.2f}**" in readme

    her2 = json.loads(
        (root / "results/external_validation/her2_ihc_ahihci/her2_ihc_benchmark.json").read_text(
            encoding="utf-8"
        )
    )
    assert f"**{her2['best_spearman_vs_her2_score']:.3f}**" in readme
    assert f"**p = {her2['max_statistic_exact_permutation_pvalue']:.6f}**" in readme
    assert f"**{her2['denominator_valid_only_best_feature_spearman']:.3f}**" in readme


def test_external_validation_claims_abstain_from_pharmacological_concordance() -> None:
    root = Path(__file__).parents[1]
    summary = json.loads(
        (root / "results/external_validation/external_validation_summary_v0.6.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["absolute_reachability"]["status"] == "NOT_COMPUTED"
    audit = json.loads(
        (
            root
            / "results/external_validation/external_data_audit/external_data_artifact_audit.json"
        ).read_text(encoding="utf-8")
    )
    assert audit["status"] == "FAIL"
    assert set(audit["invalid_artifacts"]) == {
        "bordeau_pdf_download",
        "bordeau_supplement_download",
    }


def test_readme_v07_numbers_are_artifact_backed() -> None:
    root = Path(__file__).parents[1]
    readme = (root / "README.md").read_text(encoding="utf-8")
    summary = json.loads(
        (root / "results/external_validation/external_validation_summary_v0.7.json").read_text(
            encoding="utf-8"
        )
    )
    perfusion = summary["perfusion_validation"]
    assert f"**{perfusion['median_distance_hoechst_spearman']:.3f}**" in readme
    assert f"**{perfusion['median_near_to_far_mean_ratio']:.2f}×**" in readme  # noqa: RUF001
    her2 = summary["her2_receptor_calibration"]
    assert f"**{her2['tumours_in_workbook']} tumours**" in readme
    assert f"**{her2['uncensored_replicate_pairs']} uncensored replicate pairs**" in readme
    assert f"**R² = {her2['model']['log_log_r2']:.4f}**" in readme
    prior = summary["igg_transport_prior"]["absolute_diffusion_um2_s"]
    assert f"**{prior['broad_low']:.1f}–{prior['broad_high']:.1f} µm²/s**" in readme  # noqa: RUF001
    breast = summary["breast_cohort_audit"]
    assert f"**{breast['cells']:,} unique cells**" in readme
    assert summary["absolute_reachability"]["status"] == "NOT_COMPUTED"


def test_readme_v080_numbers_are_artifact_backed() -> None:
    root = Path(__file__).parents[1]
    readme = (root / "README.md").read_text(encoding="utf-8")
    summary = json.loads(
        (root / "results/evidence_synthesis_v0.8/summary_v0.8.json").read_text(encoding="utf-8")
    )
    evidence = summary["evidence_readiness"]
    relative = summary["relative_accessibility"]
    assert f"**{evidence['absolute_readiness_score']:.1f}/100**" in readme
    assert f"**{relative['stable_top_probability']:.3%}**" in readme
    assert relative["stable_top_target"] in readme
    assert set(summary["absolute_outputs"].values()) == {"NOT_COMPUTED"}
