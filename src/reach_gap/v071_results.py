"""Integrity and scientific-boundary validation for reach-gap v0.7.1 outputs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ValidationIssue:
    code: str
    message: str


@dataclass(frozen=True)
class BundleValidation:
    status: str
    checks: int
    issues: tuple[ValidationIssue, ...]
    summary: dict[str, Any]

    @property
    def ok(self) -> bool:
        return not self.issues


REQUIRED_COMPACT_PATHS = (
    "perfusion_s_biad3159_all4/perfusion_validation_all4.json",
    "perfusion_s_biad3159_all4/perfusion_sensitivity_summary_all4.csv",
    "perfusion_s_biad3159_all4/perfusion_image_robustness_all4.csv",
    "perfusion_s_biad3159_all4/perfusion_distance_profiles_all4_sensitivity.csv",
    "breast_xenium_erbb2/breast_erbb2_validation.json",
    "breast_xenium_erbb2/breast_erbb2_sample_summary.csv",
    "breast_xenium_erbb2/breast_erbb2_group_summary.csv",
    "breast_xenium_erbb2/zarr_array_inventory.json",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _close(left: float, right: float, tolerance: float = 1e-12) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def validate_bundle(root: Path) -> BundleValidation:
    """Validate the compact Colab result bundle and its claims boundary."""

    checks = 0
    issues: list[ValidationIssue] = []

    def check(condition: bool, code: str, message: str) -> None:
        nonlocal checks
        checks += 1
        if not condition:
            issues.append(ValidationIssue(code=code, message=message))

    completion_path = root / "COLAB_COMPLETE.json"
    check(completion_path.is_file(), "missing_completion", str(completion_path))
    if not completion_path.is_file():
        return BundleValidation("FAIL", checks, tuple(issues), {})

    completion = load_json_object(completion_path)
    check(completion.get("version") == "0.7.1", "completion_version", "Expected 0.7.1")
    check(completion.get("status") == "COMPLETE", "completion_status", "Expected COMPLETE")
    check(
        completion.get("perfusion_all_four_fields") is True,
        "perfusion_gate",
        "Four fields not complete",
    )
    check(
        completion.get("breast_erbb2_two_sections") is True,
        "breast_gate",
        "Two sections not complete",
    )

    absolute = completion.get("absolute_rcc_outputs", {})
    for key in ("reachable_fraction", "penetration_depth", "expression_reach_gap"):
        check(
            isinstance(absolute, dict) and absolute.get(key) == "NOT_COMPUTED",
            f"absolute_{key}",
            f"{key} must remain NOT_COMPUTED",
        )

    manifest_rows = completion.get("files", [])
    manifest = {
        str(row["relative_path"]): row
        for row in manifest_rows
        if isinstance(row, dict) and "relative_path" in row
    }
    for relative in REQUIRED_COMPACT_PATHS:
        path = root / relative
        check(path.is_file(), "missing_compact_file", relative)
        row = manifest.get(relative)
        check(row is not None, "missing_manifest_row", relative)
        if path.is_file() and row is not None:
            check(path.stat().st_size == int(row["size_bytes"]), "size_mismatch", relative)
            check(sha256_file(path) == str(row["sha256"]), "hash_mismatch", relative)

    perfusion = load_json_object(root / REQUIRED_COMPACT_PATHS[0])
    check(perfusion.get("version") == "0.7.1", "perfusion_version", "Expected 0.7.1")
    check(perfusion.get("fields") == 4, "perfusion_fields", "Expected four fields")
    check(
        perfusion.get("sensitivity_settings_per_field") == 3,
        "perfusion_sensitivity_count",
        "Expected three settings",
    )
    check(
        perfusion.get("fields_robust_negative_distance_correlation") == 4,
        "perfusion_correlation",
        "Expected 4/4",
    )
    check(
        perfusion.get("fields_robust_near_vessel_enrichment") == 4,
        "perfusion_enrichment",
        "Expected 4/4",
    )
    median_spearman = float(perfusion.get("median_spearman_across_all_sensitivity_runs", math.nan))
    median_ratio = float(
        perfusion.get("median_near_to_far_ratio_across_all_sensitivity_runs", math.nan)
    )
    check(median_spearman < 0.0, "perfusion_direction", "Median correlation must be negative")
    check(median_ratio > 1.0, "perfusion_ratio", "Near/far ratio must exceed one")

    robustness = read_csv_rows(root / REQUIRED_COMPACT_PATHS[2])
    check(len(robustness) == 4, "robustness_rows", "Expected four image rows")
    check(
        all(row.get("all_negative_distance_correlations") == "True" for row in robustness),
        "robustness_negative",
        "Every field must be negative across settings",
    )
    check(
        all(row.get("all_near_to_far_ratios_above_one") == "True" for row in robustness),
        "robustness_ratio",
        "Every field must be enriched across settings",
    )

    sensitivity = read_csv_rows(root / REQUIRED_COMPACT_PATHS[1])
    check(len(sensitivity) == 12, "sensitivity_rows", "Expected 4 x 3 sensitivity rows")
    image_settings: dict[str, set[str]] = {}
    for row in sensitivity:
        image_settings.setdefault(row["image_label"], set()).add(row["red_bleed_correction_alpha"])
    check(
        len(image_settings) == 4
        and all(values == {"0.0", "0.5", "1.0"} for values in image_settings.values()),
        "sensitivity_grid",
        "Each field must contain alpha 0, 0.5 and 1",
    )

    breast = load_json_object(root / REQUIRED_COMPACT_PATHS[4])
    check(breast.get("version") == "0.7.1", "breast_version", "Expected 0.7.1")
    check(breast.get("sections") == 2, "breast_sections", "Expected two sections")
    check(breast.get("cells") == 679197, "breast_cells", "Expected 679,197 cells")
    sample_summary = breast.get("sample_summary", [])
    check(
        isinstance(sample_summary, list) and len(sample_summary) == 2,
        "sample_summary",
        "Expected two sample rows",
    )
    if isinstance(sample_summary, list) and len(sample_summary) == 2:
        total_cells = sum(int(row["cells"]) for row in sample_summary)
        tumor_cells = sum(int(row["provider_tumor_cells"]) for row in sample_summary)
        check(
            total_cells == int(breast["cells"]), "sample_cell_sum", "Sample cell counts do not sum"
        )
        check(
            tumor_cells == int(breast["provider_tumor_cells"]),
            "tumor_cell_sum",
            "Tumour cell counts do not sum",
        )
        rows_by_status = {str(row["provider_status"]): row for row in sample_summary}
        check(
            set(rows_by_status) == {"HER2-2+", "HER2-3+"},
            "provider_statuses",
            "Unexpected HER2 statuses",
        )
        if set(rows_by_status) == {"HER2-2+", "HER2-3+"}:
            fold = float(rows_by_status["HER2-3+"]["ERBB2_RNA_mean_provider_tumor_cells"]) / float(
                rows_by_status["HER2-2+"]["ERBB2_RNA_mean_provider_tumor_cells"]
            )
            reported_fold = float(
                breast["descriptive_tumor_group_mean_fold_HER2_3plus_vs_HER2_2plus"]
            )
            check(
                _close(fold, reported_fold),
                "fold_recalculation",
                "Reported fold does not reproduce",
            )

    sample_rows = read_csv_rows(root / REQUIRED_COMPACT_PATHS[5])
    check(len(sample_rows) == 2, "sample_csv_rows", "Expected two CSV sample rows")
    check(
        sum(int(row["cells"]) for row in sample_rows) == 679197,
        "sample_csv_cells",
        "CSV counts do not sum",
    )

    zarr_inventory = load_json_object(root / REQUIRED_COMPACT_PATHS[7])
    check(
        set(zarr_inventory) == {"S2_Middle_HER2_2plus", "S2_Bottom_HER2_3plus"},
        "zarr_samples",
        "Unexpected inventory samples",
    )
    required_arrays = {
        "cell_features/cell_id",
        "cell_features/data",
        "cell_features/indices",
        "cell_features/indptr",
        "cell_features/csc/data",
        "cell_features/csc/indices",
        "cell_features/csc/indptr",
        "__matrix_encoding__/CSR_FEATURES_BY_CELLS",
    }
    for sample, row in zarr_inventory.items():
        arrays = set(row.get("arrays", [])) if isinstance(row, dict) else set()
        check(required_arrays.issubset(arrays), "zarr_schema", f"Incomplete schema for {sample}")
        check(
            row.get("selected_feature_name") == "ERBB2",
            "zarr_feature",
            f"ERBB2 not selected for {sample}",
        )
        check(
            row.get("selected_feature_index") == 87,
            "zarr_index",
            f"Unexpected ERBB2 index for {sample}",
        )

    summary = {
        "version": "0.7.1",
        "status": (
            "EXTERNAL_PERFUSION_ALL_FOUR_FIELDS_AND_BREAST_XENIUM_ERBB2_RNA_"
            "VALIDATED_MODEL_PHARMACOLOGICAL_CONCORDANCE_NOT_COMPUTED"
        ),
        "perfusion": {
            "fields": int(perfusion["fields"]),
            "sensitivity_runs": len(sensitivity),
            "median_distance_hoechst_spearman": median_spearman,
            "median_near_to_far_mean_ratio": median_ratio,
        },
        "breast_xenium_erbb2": {
            "sections": int(breast["sections"]),
            "cells": int(breast["cells"]),
            "provider_tumor_cells": int(breast["provider_tumor_cells"]),
            "descriptive_tumor_group_mean_fold_HER2_3plus_vs_HER2_2plus": float(
                breast["descriptive_tumor_group_mean_fold_HER2_3plus_vs_HER2_2plus"]
            ),
        },
        "absolute_reachability": {
            "status": "NOT_COMPUTED",
            "reasons": completion.get("reasons_absolute_rcc_remains_not_computed", []),
        },
        "validation": {"checks": checks, "issues": len(issues)},
    }
    return BundleValidation("PASS" if not issues else "FAIL", checks, tuple(issues), summary)


def build_claims(summary: dict[str, Any]) -> dict[str, Any]:
    """Build machine-readable permitted and unsupported v0.7.1 claims."""

    return {
        "version": "0.7.1",
        "status": summary["status"],
        "permitted": [
            {
                "claim": (
                    "All four independent S-BIAD3159 fields show a robust relative Hoechst "
                    "gradient with distance to structural CD31 signal across the locked "
                    "red-channel sensitivity grid."
                ),
                "evidence": summary["perfusion"],
            },
            {
                "claim": (
                    "ERBB2 RNA was extracted from two independent native Xenium "
                    "cell_features Zarr matrices covering 679,197 labelled cells."
                ),
                "evidence": summary["breast_xenium_erbb2"],
            },
            {
                "claim": (
                    "Provider-labelled tumour cells in the HER2-3+ section have a "
                    "descriptively higher mean ERBB2 RNA count than those in the "
                    "HER2-2+ section."
                ),
                "qualifier": (
                    "Two sections only; no cell-level inferential p-value and no clinical "
                    "equivalence claim."
                ),
            },
        ],
        "conditional": [
            "The external IgG interval may be used only as a solver sensitivity prior.",
            (
                "The source-protocol HER2 calibration may be used only within its original "
                "Cy5 assay unless a shared calibrator is measured."
            ),
        ],
        "unsupported": [
            "Functional perfusion of the RCC Xenium section",
            "HER2 surface receptor copies per Xenium cell",
            "Therapeutic-antibody concentration, binding or penetration in either breast section",
            "Absolute RCC reachable_fraction",
            "Absolute RCC penetration_depth",
            "Absolute RCC expression_reach_gap",
            "Model-versus-administered-antibody pharmacological concordance",
        ],
    }
