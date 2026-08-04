"""Build and validate the integrated v0.8 evidence-synthesis result package."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from reach_gap.evidence_synthesis import synthesize_repository_evidence
from reach_gap.relative_accessibility import run_relative_accessibility

REQUIRED_PATHS = (
    "evidence/evidence_graph_v0.8.json",
    "evidence/evidence_requirements_v0.8.csv",
    "evidence/measurement_priorities_v0.8.csv",
    "evidence/evidence_uncertainty_budget_v0.8.svg",
    "relative_accessibility/relative_accessibility_v0.8.json",
    "relative_accessibility/target_rank_summary_v0.8.csv",
    "relative_accessibility/pairwise_win_probability_v0.8.csv",
    "relative_accessibility/leave_one_component_out_v0.8.csv",
    "relative_accessibility/target_top_rank_probability_v0.8.svg",
    "summary_v0.8.json",
    "claims_v0.8.json",
)


@dataclass(frozen=True)
class ValidationIssue:
    """One failed v0.8 package invariant."""

    code: str
    message: str


@dataclass(frozen=True)
class V080Validation:
    """Validation status and machine-readable summary."""

    status: str
    checks: int
    issues: tuple[ValidationIssue, ...]
    summary: dict[str, Any]

    @property
    def ok(self) -> bool:
        """Return whether every invariant passed."""

        return self.status == "PASS"


def sha256_file(path: Path) -> str:
    """Return a lower-case SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected JSON object: {path}")
    return payload


def _write_manifest(output_dir: Path) -> Path:
    rows = []
    for path in sorted(file for file in output_dir.rglob("*") if file.is_file()):
        relative = path.relative_to(output_dir).as_posix()
        if relative == "artifact_manifest_v0.8.json":
            continue
        rows.append(
            {
                "relative_path": relative,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    manifest = {
        "schema_version": "1.0",
        "version": "0.8.0",
        "files": rows,
    }
    path = output_dir / "artifact_manifest_v0.8.json"
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path


def build_v080_package(
    repository_root: Path,
    output_dir: Path,
    *,
    draws: int = 20_000,
    seed: int = 17,
) -> dict[str, Any]:
    """Build evidence graph, relative target analysis, claims and manifest."""

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence = synthesize_repository_evidence(repository_root, output_dir / "evidence")
    relative = run_relative_accessibility(
        repository_root,
        output_dir / "relative_accessibility",
        draws=draws,
        seed=seed,
    )
    evidence_summary = evidence["summary"]
    relative_summary = relative["summary"]
    summary = {
        "version": "0.8.0",
        "status": (
            "EVIDENCE_SYNTHESIS_AND_RELATIVE_RCC_TARGET_RANKING_COMPLETE_"
            "ABSOLUTE_REACHABILITY_NOT_COMPUTED"
        ),
        "evidence_readiness": {
            "absolute_readiness_score": evidence_summary["absolute_readiness_score"],
            "same_tissue_requirements_satisfied": evidence_summary[
                "same_tissue_requirements_satisfied"
            ],
            "requirements_total": evidence_summary["requirements_total"],
            "dominant_uncertainty": evidence_summary["dominant_uncertainty"],
            "measurement_priority": evidence_summary["measurement_priority"],
        },
        "relative_accessibility": {
            "status": relative_summary["status"],
            "targets": relative_summary["target_count"],
            "vessel_definitions": relative_summary["vessel_definition_count"],
            "draws": relative_summary["draws"],
            "stable_top_target": relative_summary["stable_top_target"],
            "stable_top_probability": relative_summary["stable_top_probability"],
        },
        "absolute_outputs": evidence_summary["absolute_outputs"],
        "interpretation": (
            "The v0.8 score measures evidence readiness, and the target ranking is a relative "
            "geometry-expression proxy. Neither is an estimate of therapeutic-"
            "antibody reachability."
        ),
    }
    claims = {
        "version": "0.8.0",
        "status": summary["status"],
        "permitted": [
            (
                "Evidence dependencies and blocked transfers are represented in "
                "a machine-readable graph."
            ),
            "The absolute-readiness score describes evidence completeness under a fixed rubric.",
            (
                "RCC targets are compared using a relative geometry-expression proxy across six "
                "structural vessel definitions and uncertain objective weights."
            ),
            (
                f"{relative_summary['stable_top_target']} is the stable top target "
                "within the relative "
                f"proxy in {relative_summary['stable_top_probability']:.3%} of draws."
            ),
        ],
        "conditional": [
            (
                "The relative ranking is conditional on the measured target set "
                "and within-section thresholds."
            ),
            "Structural CD31 definitions are not functionally perfused-vessel labels.",
            (
                "External calibrations and priors constrain method uncertainty "
                "but do not identify RCC parameters."
            ),
        ],
        "unsupported": [
            "Absolute RCC reachable_fraction.",
            "Absolute RCC penetration_depth.",
            "Absolute RCC expression_reach_gap.",
            "Clinical efficacy or programme-priority claims from the relative target ranking.",
            "Model-versus-administered-antibody concordance.",
        ],
    }
    summary_path = output_dir / "summary_v0.8.json"
    claims_path = output_dir / "claims_v0.8.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    claims_path.write_text(json.dumps(claims, indent=2), encoding="utf-8")
    manifest_path = _write_manifest(output_dir)
    validation = validate_v080_package(output_dir)
    validation_path = output_dir / "validation_v0.8.json"
    validation_path.write_text(
        json.dumps(
            {
                "status": validation.status,
                "checks": validation.checks,
                "issues": [issue.__dict__ for issue in validation.issues],
                "summary": validation.summary,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    if not validation.ok:
        raise RuntimeError(f"v0.8 package validation failed: {validation.issues}")
    _write_manifest(output_dir)
    final_validation = validate_v080_package(output_dir)
    validation_path.write_text(
        json.dumps(
            {
                "status": final_validation.status,
                "checks": final_validation.checks,
                "issues": [issue.__dict__ for issue in final_validation.issues],
                "summary": final_validation.summary,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    if not final_validation.ok:
        raise RuntimeError(f"Final v0.8 package validation failed: {final_validation.issues}")
    _write_manifest(output_dir)
    confirmation = validate_v080_package(output_dir)
    if not confirmation.ok:
        raise RuntimeError(f"v0.8 manifest confirmation failed: {confirmation.issues}")
    return {
        "summary": summary,
        "claims": claims,
        "manifest": str(manifest_path),
        "validation": str(validation_path),
    }


def validate_v080_package(output_dir: Path) -> V080Validation:
    """Validate v0.8 scientific boundaries, deterministic results and file hashes."""

    issues: list[ValidationIssue] = []
    checks = 0

    def check(condition: bool, code: str, message: str) -> None:
        nonlocal checks
        checks += 1
        if not condition:
            issues.append(ValidationIssue(code, message))

    for relative_name in REQUIRED_PATHS:
        check((output_dir / relative_name).is_file(), "missing_file", relative_name)
    manifest_path = output_dir / "artifact_manifest_v0.8.json"
    check(manifest_path.is_file(), "missing_manifest", str(manifest_path))
    if not all((output_dir / relative_name).is_file() for relative_name in REQUIRED_PATHS):
        return V080Validation("FAIL", checks, tuple(issues), {})

    evidence = _load_json(output_dir / "evidence/evidence_graph_v0.8.json")
    relative_payload = _load_json(
        output_dir / "relative_accessibility/relative_accessibility_v0.8.json"
    )
    summary = _load_json(output_dir / "summary_v0.8.json")
    claims = _load_json(output_dir / "claims_v0.8.json")

    check(evidence.get("version") == "0.8.0", "evidence_version", "Expected 0.8.0")
    check(relative_payload.get("version") == "0.8.0", "relative_version", "Expected 0.8.0")
    check(summary.get("version") == "0.8.0", "summary_version", "Expected 0.8.0")
    check(claims.get("version") == "0.8.0", "claims_version", "Expected 0.8.0")
    requirements = evidence.get("requirements", [])
    check(len(requirements) == 8, "requirements", "Expected eight absolute-model requirements")
    if isinstance(requirements, list) and requirements:
        check(
            abs(sum(float(row["weight"]) for row in requirements) - 1.0) < 1.0e-12,
            "weights",
            "Requirement weights must sum to one",
        )
        check(
            abs(sum(float(row["uncertainty_share"]) for row in requirements) - 1.0) < 1.0e-12,
            "uncertainty_shares",
            "Uncertainty shares must sum to one",
        )
    check(
        abs(float(evidence.get("absolute_readiness_score", -1.0)) - 40.5) < 1.0e-9,
        "readiness_score",
        "Unexpected readiness score",
    )
    check(
        evidence.get("same_tissue_requirements_satisfied") == 3,
        "same_tissue_requirements",
        "Expected three satisfied same-tissue requirements",
    )
    check(
        evidence.get("dominant_uncertainty") == "administered_antibody_field",
        "dominant_uncertainty",
        "Unexpected dominant uncertainty",
    )
    absolute = evidence.get("absolute_outputs", {})
    for key in (
        "reachable_fraction",
        "penetration_depth",
        "expression_reach_gap",
        "model_pharmacological_concordance",
    ):
        check(
            isinstance(absolute, dict) and absolute.get(key) == "NOT_COMPUTED",
            f"absolute_{key}",
            f"{key} must remain NOT_COMPUTED",
        )

    check(relative_payload.get("target_count") == 4, "target_count", "Expected four RCC targets")
    check(
        relative_payload.get("vessel_definition_count") == 6,
        "vessel_definitions",
        "Expected six vessel definitions",
    )
    check(relative_payload.get("draws") == 20_000, "draws", "Expected 20,000 draws")
    check(
        relative_payload.get("stable_top_target") == "VISTA",
        "stable_top_target",
        "Expected VISTA as stable relative-proxy top target",
    )
    top_probability = float(relative_payload.get("stable_top_probability", 0.0))
    check(top_probability >= 0.97, "top_probability", "VISTA top probability below 97%")
    pairwise = relative_payload.get("pairwise_win_probability", {})
    if isinstance(pairwise, dict) and "VISTA" in pairwise:
        vista = pairwise["VISTA"]
        check(
            all(float(vista[target]) >= 0.98 for target in ("LAG-3", "PD-1", "PD-L1")),
            "vista_pairwise",
            "VISTA pairwise probability below 98%",
        )
    else:
        check(False, "pairwise_missing", "VISTA pairwise probabilities missing")

    unsupported = claims.get("unsupported", [])
    check(
        any("reachable_fraction" in str(statement) for statement in unsupported),
        "claim_abstention",
        "Absolute reachability abstention missing",
    )

    if manifest_path.is_file():
        manifest = _load_json(manifest_path)
        rows = manifest.get("files", [])
        check(
            isinstance(rows, list) and len(rows) >= len(REQUIRED_PATHS),
            "manifest_rows",
            "Manifest incomplete",
        )
        if isinstance(rows, list):
            for row in rows:
                if not isinstance(row, dict):
                    check(False, "manifest_row", "Manifest row must be an object")
                    continue
                relative_path = str(row.get("relative_path", ""))
                path = output_dir / relative_path
                check(path.is_file(), "manifest_file", relative_path)
                if path.is_file():
                    check(
                        path.stat().st_size == int(row.get("size_bytes", -1)),
                        "manifest_size",
                        relative_path,
                    )
                    check(
                        sha256_file(path) == str(row.get("sha256", "")),
                        "manifest_hash",
                        relative_path,
                    )

    compact_summary = {
        "version": "0.8.0",
        "status": summary["status"],
        "absolute_readiness_score": float(evidence["absolute_readiness_score"]),
        "dominant_uncertainty": evidence["dominant_uncertainty"],
        "stable_top_target": relative_payload["stable_top_target"],
        "stable_top_probability": top_probability,
        "absolute_outputs": absolute,
    }
    return V080Validation("PASS" if not issues else "FAIL", checks, tuple(issues), compact_summary)
