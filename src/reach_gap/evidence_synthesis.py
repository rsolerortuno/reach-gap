"""Evidence graph, readiness audit and measurement-priority synthesis."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

from reach_gap.visuals import write_horizontal_bar_svg

EvidenceClass = Literal[
    "SAME_TISSUE_MEASUREMENT",
    "EXTERNAL_METHOD_VALIDATION",
    "EXTERNAL_ASSAY_CALIBRATION",
    "LITERATURE_PRIOR",
    "PUBLISHED_REFERENCE",
    "MISSING_REQUIREMENT",
    "BLOCKED_OUTPUT",
]
EvidenceState = Literal["SATISFIED", "PARTIAL", "EXTERNAL_ONLY", "MISSING"]
Relation = Literal["SUPPORTS", "CONSTRAINS", "BLOCKED_TRANSFER", "REQUIRED_FOR"]


class EvidenceNode(BaseModel):
    """One auditable evidence object or explicitly missing requirement."""

    model_config = ConfigDict(extra="forbid")

    node_id: str
    title: str
    domain: str
    evidence_class: EvidenceClass
    tissue: str
    organism: str
    assay: str
    status: str
    source_path: str | None = None
    measured_values: dict[str, float | int | str] = Field(default_factory=dict)
    permitted_use: list[str] = Field(default_factory=list)
    prohibited_transfer: list[str] = Field(default_factory=list)


class EvidenceEdge(BaseModel):
    """Directed relation between evidence objects and requirements."""

    model_config = ConfigDict(extra="forbid")

    source: str
    target: str
    relation: Relation
    transferable: bool
    reason: str


class RequirementAssessment(BaseModel):
    """Transparent readiness assessment for one absolute-model requirement."""

    model_config = ConfigDict(extra="forbid")

    requirement: str
    weight: float = Field(gt=0.0, le=1.0)
    evidence_state: EvidenceState
    satisfaction: float = Field(ge=0.0, le=1.0)
    supporting_nodes: list[str]
    reason: str
    next_measurement: str | None
    uncertainty_contribution: float = Field(ge=0.0, le=1.0)
    uncertainty_share: float = Field(ge=0.0, le=1.0)


class EvidenceSynthesis(BaseModel):
    """Machine-readable v0.8 evidence graph and abstention boundary."""

    model_config = ConfigDict(extra="forbid")

    version: str
    status: str
    nodes: list[EvidenceNode]
    edges: list[EvidenceEdge]
    requirements: list[RequirementAssessment]
    absolute_readiness_score: float = Field(ge=0.0, le=100.0)
    satisfied_weight: float = Field(ge=0.0, le=1.0)
    same_tissue_requirements_satisfied: int
    requirements_total: int
    dominant_uncertainty: str
    measurement_priority: list[str]
    absolute_outputs: dict[str, str]
    permitted_claims: list[str]
    unsupported_claims: list[str]
    scoring_rubric: dict[str, float]


_STATE_SATISFACTION: dict[EvidenceState, float] = {
    "SATISFIED": 1.0,
    "PARTIAL": 0.5,
    "EXTERNAL_ONLY": 0.25,
    "MISSING": 0.0,
}

_REQUIREMENT_WEIGHTS: dict[str, float] = {
    "tissue_geometry": 0.08,
    "structural_vasculature": 0.08,
    "spatial_target_localisation": 0.12,
    "functional_perfusion": 0.18,
    "surface_antigen_calibration": 0.18,
    "matrix_transport": 0.14,
    "administered_antibody_field": 0.15,
    "same_tissue_pharmacological_endpoint": 0.07,
}


def _load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected a JSON object: {path}")
    return payload


def _require_status(payload: dict[str, Any], expected_fragment: str, path: Path) -> None:
    status = str(payload.get("status", ""))
    if expected_fragment not in status:
        raise ValueError(f"Unexpected status in {path}: {status!r}")


def _nodes(repository_root: Path) -> list[EvidenceNode]:
    rcc = _load_json(repository_root / "results/real_rcc_xenium/real_data_summary.json")
    imaging = _load_json(repository_root / "results/real_rcc_imaging/real_imaging_summary.json")
    he = _load_json(repository_root / "results/real_rcc_he/he_pathology_result.json")
    perfusion = _load_json(
        repository_root / "results/external_validation/v0.7.1/perfusion_s_biad3159_all4/"
        "perfusion_validation_all4.json"
    )
    breast = _load_json(
        repository_root / "results/external_validation/v0.7.1/breast_xenium_erbb2/"
        "breast_erbb2_validation.json"
    )
    transport = _load_json(
        repository_root / "results/external_validation/igg_transport_prior/igg_transport_prior.json"
    )
    calibration = _load_json(
        repository_root / "results/external_validation/her2_receptor_calibration_mckinski2026/"
        "her2_receptor_calibration.json"
    )
    bordeau = _load_json(
        repository_root / "results/external_validation/bordeau_2021_supplement/"
        "bordeau_supplement_benchmark.json"
    )
    _require_status(rcc, "ABSOLUTE_REACHABILITY_NOT_COMPUTED", Path("real_data_summary.json"))
    _require_status(imaging, "ABSOLUTE_INDEX_NOT_COMPUTED", Path("real_imaging_summary.json"))
    _require_status(he, "TARGET_INDEX_NOT_COMPUTED", Path("he_pathology_result.json"))
    _require_status(perfusion, "ALL_FOUR_FIELDS_COMPLETE", Path("perfusion_validation_all4.json"))
    _require_status(breast, "ERBB2_RNA_SUMMARY_COMPLETE", Path("breast_erbb2_validation.json"))

    return [
        EvidenceNode(
            node_id="rcc_histology_geometry",
            title="RCC histology and pathology geometry",
            domain="tissue_geometry",
            evidence_class="SAME_TISSUE_MEASUREMENT",
            tissue="human RCC section",
            organism="human",
            assay="H&E plus pathology polygons",
            status=str(he["status"]),
            source_path="results/real_rcc_he/he_pathology_result.json",
            measured_values={
                "tissue_area_mm2": float(he["tissue"]["estimated_area_mm2"]),
                "analysis_pixel_size_um": float(he["image"]["analysis_size_um_per_pixel"][0]),
            },
            permitted_use=["Same-section morphology and pathology-region geometry"],
            prohibited_transfer=["H&E lumen candidates as validated perfused vessels"],
        ),
        EvidenceNode(
            node_id="rcc_structural_vasculature",
            title="RCC structural CD31 vessel-proxy geometry",
            domain="structural_vasculature",
            evidence_class="SAME_TISSUE_MEASUREMENT",
            tissue="human RCC section",
            organism="human",
            assay="Xenium protein plus morphology-focus imaging",
            status=str(imaging["status"]),
            source_path="results/real_rcc_imaging/real_imaging_summary.json",
            measured_values={
                "cells": int(imaging["cells"]),
                "vessel_definition_count": int(rcc["vessel_proxy"]["definition_count"]),
                "vessel_fraction_low": float(rcc["vessel_proxy"]["fraction_range"][0]),
                "vessel_fraction_high": float(rcc["vessel_proxy"]["fraction_range"][1]),
            },
            permitted_use=["Relative distance to structural CD31 definitions in the RCC section"],
            prohibited_transfer=["Structural CD31 positivity as proof of functional perfusion"],
        ),
        EvidenceNode(
            node_id="rcc_target_localisation",
            title="RCC target protein and RNA localisation",
            domain="spatial_target_localisation",
            evidence_class="SAME_TISSUE_MEASUREMENT",
            tissue="human RCC section",
            organism="human",
            assay="Xenium gene and protein panel",
            status=str(rcc["status"]),
            source_path="results/real_rcc_xenium/real_data_summary.json",
            measured_values={
                "cells": int(rcc["cells"]),
                "targets": len(rcc["targets"]),
                "pathology_alignment_score": float(rcc["pathology_alignment"]["score"]),
            },
            permitted_use=["Relative within-section target prevalence and spatial geometry"],
            prohibited_transfer=["Fluorescence or RNA counts as surface molecules per cell"],
        ),
        EvidenceNode(
            node_id="external_llc_perfusion",
            title="Four-field in-vivo Hoechst/CD31 perfusion-proxy validation",
            domain="functional_perfusion",
            evidence_class="EXTERNAL_METHOD_VALIDATION",
            tissue="mouse LLC tumours",
            organism="mouse",
            assay="RGB fluorescence microscopy",
            status=str(perfusion["status"]),
            source_path=(
                "results/external_validation/v0.7.1/perfusion_s_biad3159_all4/"
                "perfusion_validation_all4.json"
            ),
            measured_values={
                "fields": int(perfusion["fields"]),
                "sensitivity_runs": int(
                    perfusion["fields"] * perfusion["sensitivity_settings_per_field"]
                ),
                "median_spearman": float(perfusion["median_spearman_across_all_sensitivity_runs"]),
                "median_near_far_ratio": float(
                    perfusion["median_near_to_far_ratio_across_all_sensitivity_runs"]
                ),
            },
            permitted_use=["Validation of a relative perfusion-gradient extraction method"],
            prohibited_transfer=["Functional perfusion labels or source terms in RCC"],
        ),
        EvidenceNode(
            node_id="external_igg_transport_prior",
            title="Published tumour IgG diffusion prior",
            domain="matrix_transport",
            evidence_class="LITERATURE_PRIOR",
            tissue="external tumour systems",
            organism="mixed published models",
            assay="FRAP literature measurements",
            status=str(transport["status"]),
            source_path="results/external_validation/igg_transport_prior/igg_transport_prior.json",
            measured_values={
                "diffusion_low_um2_s": float(transport["absolute_diffusion_um2_s"]["broad_low"]),
                "diffusion_high_um2_s": float(transport["absolute_diffusion_um2_s"]["broad_high"]),
            },
            permitted_use=["Log-uniform solver sensitivity prior"],
            prohibited_transfer=["RCC-specific transport coefficient"],
        ),
        EvidenceNode(
            node_id="external_her2_calibration",
            title="Source-protocol Cy5-to-HER2 receptor calibration",
            domain="surface_antigen_calibration",
            evidence_class="EXTERNAL_ASSAY_CALIBRATION",
            tissue="external HER2-positive tumours",
            organism="human tumour material",
            assay="source Cy5 quantitative-IHC workflow",
            status=str(calibration["status"]),
            source_path=(
                "results/external_validation/her2_receptor_calibration_mckinski2026/"
                "her2_receptor_calibration.json"
            ),
            measured_values={
                "tumours": int(calibration["tumours_in_workbook"]),
                "replicate_pairs": int(calibration["uncensored_replicate_pairs"]),
                "log_log_r2": float(calibration["model"]["log_log_r2"]),
            },
            permitted_use=["Calibration within the original Cy5 protocol"],
            prohibited_transfer=["Conversion of RCC Xenium channels to receptors per cell"],
        ),
        EvidenceNode(
            node_id="external_breast_erbb2_control",
            title="Independent breast Xenium ERBB2 RNA control",
            domain="external_assay_control",
            evidence_class="EXTERNAL_METHOD_VALIDATION",
            tissue="two human breast cancer sections",
            organism="human",
            assay="native Xenium cell_features Zarr",
            status=str(breast["status"]),
            source_path=(
                "results/external_validation/v0.7.1/breast_xenium_erbb2/"
                "breast_erbb2_validation.json"
            ),
            measured_values={
                "cells": int(breast["cells"]),
                "provider_tumor_cells": int(breast["provider_tumor_cells"]),
                "descriptive_fold": float(
                    breast["descriptive_tumor_group_mean_fold_HER2_3plus_vs_HER2_2plus"]
                ),
            },
            permitted_use=[
                "Validation of Xenium sparse-matrix extraction and descriptive RNA contrast"
            ],
            prohibited_transfer=["Surface HER2 calibration or antibody-reach inference"],
        ),
        EvidenceNode(
            node_id="external_administered_antibody_reference",
            title="Published administered-trastuzumab reference",
            domain="administered_antibody_field",
            evidence_class="PUBLISHED_REFERENCE",
            tissue="SKOV3 xenograft",
            organism="mouse xenograft",
            assay="published fluorescence microscopy",
            status=str(bordeau["status"]),
            source_path=(
                "results/external_validation/bordeau_2021_supplement/"
                "bordeau_supplement_benchmark.json"
            ),
            measured_values={
                "dose_mg_per_kg": float(bordeau["published_endpoints"]["dose_mg_per_kg"]),
                "time_hours": float(bordeau["published_endpoints"]["time_hours"]),
                "representative_panel_permutation_p": float(
                    bordeau["representative_figure_analysis"]["exact_two_sided_permutation_pvalue"]
                ),
            },
            permitted_use=["External reference endpoint and negative-result audit"],
            prohibited_transfer=["Same-tissue RCC antibody concentration or model concordance"],
        ),
        EvidenceNode(
            node_id="missing_rcc_functional_perfusion",
            title="RCC same-section functional perfusion",
            domain="functional_perfusion",
            evidence_class="MISSING_REQUIREMENT",
            tissue="human RCC section",
            organism="human",
            assay="not measured",
            status="NOT_COMPUTED",
            permitted_use=[],
            prohibited_transfer=["Do not substitute external LLC perfusion for RCC perfusion"],
        ),
        EvidenceNode(
            node_id="missing_rcc_surface_antigen_calibration",
            title="RCC surface-antigen molecules per cell",
            domain="surface_antigen_calibration",
            evidence_class="MISSING_REQUIREMENT",
            tissue="human RCC section",
            organism="human",
            assay="not measured",
            status="NOT_COMPUTED",
            permitted_use=[],
            prohibited_transfer=["Do not convert RNA or arbitrary fluorescence to receptor copies"],
        ),
        EvidenceNode(
            node_id="missing_rcc_administered_antibody_field",
            title="RCC administered-antibody concentration field",
            domain="administered_antibody_field",
            evidence_class="MISSING_REQUIREMENT",
            tissue="human RCC section",
            organism="human",
            assay="not measured",
            status="NOT_COMPUTED",
            permitted_use=[],
            prohibited_transfer=["Do not infer drug concentration from structural geometry alone"],
        ),
        EvidenceNode(
            node_id="missing_rcc_pharmacological_endpoint",
            title="RCC same-tissue pharmacological endpoint",
            domain="same_tissue_pharmacological_endpoint",
            evidence_class="MISSING_REQUIREMENT",
            tissue="human RCC section",
            organism="human",
            assay="not measured",
            status="NOT_COMPUTED",
            permitted_use=[],
            prohibited_transfer=["Do not claim model-versus-drug concordance without an endpoint"],
        ),
        EvidenceNode(
            node_id="blocked_absolute_rcc_outputs",
            title="Absolute RCC reachability outputs",
            domain="absolute_outputs",
            evidence_class="BLOCKED_OUTPUT",
            tissue="human RCC section",
            organism="human",
            assay="mechanistic synthesis",
            status="NOT_COMPUTED",
            permitted_use=["Explicit abstention and requirement audit"],
            prohibited_transfer=[
                "No absolute reachable_fraction",
                "No absolute penetration_depth",
                "No absolute expression_reach_gap",
            ],
        ),
    ]


def _edges() -> list[EvidenceEdge]:
    return [
        EvidenceEdge(
            source="rcc_histology_geometry",
            target="blocked_absolute_rcc_outputs",
            relation="SUPPORTS",
            transferable=True,
            reason=(
                "Same-section morphology supplies geometry but not source or binding calibration."
            ),
        ),
        EvidenceEdge(
            source="rcc_structural_vasculature",
            target="blocked_absolute_rcc_outputs",
            relation="SUPPORTS",
            transferable=True,
            reason="Structural source geometry is measured under multiple definitions.",
        ),
        EvidenceEdge(
            source="rcc_target_localisation",
            target="blocked_absolute_rcc_outputs",
            relation="SUPPORTS",
            transferable=True,
            reason="Relative target localisation is measured in the RCC section.",
        ),
        EvidenceEdge(
            source="external_llc_perfusion",
            target="missing_rcc_functional_perfusion",
            relation="BLOCKED_TRANSFER",
            transferable=False,
            reason=(
                "Different organism, tumour model and image field; method validation "
                "is transferable, labels are not."
            ),
        ),
        EvidenceEdge(
            source="external_her2_calibration",
            target="missing_rcc_surface_antigen_calibration",
            relation="BLOCKED_TRANSFER",
            transferable=False,
            reason="No shared calibrator between the Cy5 source protocol and RCC Xenium channels.",
        ),
        EvidenceEdge(
            source="external_igg_transport_prior",
            target="blocked_absolute_rcc_outputs",
            relation="CONSTRAINS",
            transferable=True,
            reason="The prior constrains sensitivity ranges but does not identify RCC diffusivity.",
        ),
        EvidenceEdge(
            source="external_administered_antibody_reference",
            target="missing_rcc_administered_antibody_field",
            relation="BLOCKED_TRANSFER",
            transferable=False,
            reason="Published SKOV3 endpoints are not a concentration field in the RCC section.",
        ),
        EvidenceEdge(
            source="external_breast_erbb2_control",
            target="rcc_target_localisation",
            relation="SUPPORTS",
            transferable=True,
            reason="It validates the native Xenium Zarr adapter, not biological calibration.",
        ),
        EvidenceEdge(
            source="missing_rcc_functional_perfusion",
            target="blocked_absolute_rcc_outputs",
            relation="REQUIRED_FOR",
            transferable=False,
            reason=(
                "The vascular boundary condition is unidentified without same-section perfusion."
            ),
        ),
        EvidenceEdge(
            source="missing_rcc_surface_antigen_calibration",
            target="blocked_absolute_rcc_outputs",
            relation="REQUIRED_FOR",
            transferable=False,
            reason="Binding-site density is unidentified without surface molecules per cell.",
        ),
        EvidenceEdge(
            source="missing_rcc_administered_antibody_field",
            target="blocked_absolute_rcc_outputs",
            relation="REQUIRED_FOR",
            transferable=False,
            reason="Absolute source concentration and validation are unavailable.",
        ),
        EvidenceEdge(
            source="missing_rcc_pharmacological_endpoint",
            target="blocked_absolute_rcc_outputs",
            relation="REQUIRED_FOR",
            transferable=False,
            reason="Model concordance cannot be evaluated without a same-tissue endpoint.",
        ),
    ]


def _raw_requirements() -> list[dict[str, Any]]:
    return [
        {
            "requirement": "tissue_geometry",
            "state": "SATISFIED",
            "nodes": ["rcc_histology_geometry"],
            "reason": "Same-section physical geometry and pathology regions are measured.",
            "next": None,
        },
        {
            "requirement": "structural_vasculature",
            "state": "SATISFIED",
            "nodes": ["rcc_structural_vasculature"],
            "reason": "Same-section structural CD31 geometry is available under six definitions.",
            "next": None,
        },
        {
            "requirement": "spatial_target_localisation",
            "state": "SATISFIED",
            "nodes": ["rcc_target_localisation"],
            "reason": "Relative target protein/RNA localisation is available for 465,534 cells.",
            "next": None,
        },
        {
            "requirement": "functional_perfusion",
            "state": "EXTERNAL_ONLY",
            "nodes": ["external_llc_perfusion", "missing_rcc_functional_perfusion"],
            "reason": (
                "The method is externally validated, but RCC vessel perfusion is not measured."
            ),
            "next": "Measure an intravascular perfusion tracer co-registered to RCC CD31.",
        },
        {
            "requirement": "surface_antigen_calibration",
            "state": "EXTERNAL_ONLY",
            "nodes": ["external_her2_calibration", "missing_rcc_surface_antigen_calibration"],
            "reason": (
                "A source-protocol calibration exists, but no shared RCC/Xenium calibrator exists."
            ),
            "next": (
                "Measure target surface molecules per cell with a shared quantitative calibrator."
            ),
        },
        {
            "requirement": "matrix_transport",
            "state": "EXTERNAL_ONLY",
            "nodes": ["external_igg_transport_prior"],
            "reason": "Published FRAP values constrain sensitivity but are not measured in RCC.",
            "next": "Measure IgG diffusivity or a validated transport surrogate in the RCC tissue.",
        },
        {
            "requirement": "administered_antibody_field",
            "state": "MISSING",
            "nodes": [
                "external_administered_antibody_reference",
                "missing_rcc_administered_antibody_field",
            ],
            "reason": "No administered-antibody concentration or engagement field exists in RCC.",
            "next": (
                "Acquire spatial antibody concentration or target-engagement "
                "measurements after dosing."
            ),
        },
        {
            "requirement": "same_tissue_pharmacological_endpoint",
            "state": "MISSING",
            "nodes": ["missing_rcc_pharmacological_endpoint"],
            "reason": (
                "No same-tissue response or engagement endpoint is available for concordance."
            ),
            "next": (
                "Collect a blinded same-sample pharmacological endpoint before model evaluation."
            ),
        },
    ]


def build_evidence_synthesis(repository_root: Path) -> EvidenceSynthesis:
    """Build a fixed-rubric evidence graph from committed v0.7.1 artifacts."""

    if abs(sum(_REQUIREMENT_WEIGHTS.values()) - 1.0) > 1.0e-12:
        raise RuntimeError("Requirement weights must sum to one")
    nodes = _nodes(repository_root)
    edges = _edges()
    raw = _raw_requirements()
    contributions = []
    for row in raw:
        state = row["state"]
        if state not in _STATE_SATISFACTION:
            raise ValueError(f"Unknown evidence state: {state}")
        weight = _REQUIREMENT_WEIGHTS[row["requirement"]]
        satisfaction = _STATE_SATISFACTION[state]
        contributions.append(weight * (1.0 - satisfaction))
    burden_total = sum(contributions)
    if burden_total <= 0.0:
        raise RuntimeError("Uncertainty burden unexpectedly vanished")

    requirements: list[RequirementAssessment] = []
    for row, contribution in zip(raw, contributions, strict=True):
        state = row["state"]
        weight = _REQUIREMENT_WEIGHTS[row["requirement"]]
        satisfaction = _STATE_SATISFACTION[state]
        requirements.append(
            RequirementAssessment(
                requirement=row["requirement"],
                weight=weight,
                evidence_state=state,
                satisfaction=satisfaction,
                supporting_nodes=row["nodes"],
                reason=row["reason"],
                next_measurement=row["next"],
                uncertainty_contribution=contribution,
                uncertainty_share=contribution / burden_total,
            )
        )
    satisfied_weight = sum(row.weight * row.satisfaction for row in requirements)
    priorities = [
        row.requirement
        for row in sorted(
            requirements,
            key=lambda item: (-item.uncertainty_share, item.requirement),
        )
        if row.next_measurement is not None
    ]
    dominant = priorities[0]
    return EvidenceSynthesis(
        version="0.8.0",
        status="EVIDENCE_GRAPH_COMPLETE_ABSOLUTE_RCC_REACHABILITY_NOT_COMPUTED",
        nodes=nodes,
        edges=edges,
        requirements=requirements,
        absolute_readiness_score=100.0 * satisfied_weight,
        satisfied_weight=satisfied_weight,
        same_tissue_requirements_satisfied=sum(
            row.evidence_state == "SATISFIED" for row in requirements
        ),
        requirements_total=len(requirements),
        dominant_uncertainty=dominant,
        measurement_priority=priorities,
        absolute_outputs={
            "reachable_fraction": "NOT_COMPUTED",
            "penetration_depth": "NOT_COMPUTED",
            "expression_reach_gap": "NOT_COMPUTED",
            "model_pharmacological_concordance": "NOT_COMPUTED",
        },
        permitted_claims=[
            (
                "The evidence needed for absolute RCC reachability is represented "
                "as an auditable graph."
            ),
            (
                "Same-tissue measurements, external validations and literature "
                "priors are kept distinct."
            ),
            "The readiness score describes evidence completeness, not biological reachability.",
            (
                "Measurement priorities rank unresolved uncertainty under a fixed "
                "published weighting rubric."
            ),
        ],
        unsupported_claims=[
            "Absolute therapeutic-antibody reachability in RCC.",
            "Transfer of external LLC perfusion labels to RCC vessels.",
            "Conversion of RCC RNA or fluorescence to surface molecules per cell.",
            "Model-versus-drug concordance in the absence of a same-tissue endpoint.",
        ],
        scoring_rubric={
            "SATISFIED": _STATE_SATISFACTION["SATISFIED"],
            "PARTIAL": _STATE_SATISFACTION["PARTIAL"],
            "EXTERNAL_ONLY": _STATE_SATISFACTION["EXTERNAL_ONLY"],
            "MISSING": _STATE_SATISFACTION["MISSING"],
        },
    )


def write_evidence_synthesis(synthesis: EvidenceSynthesis, output_dir: Path) -> dict[str, str]:
    """Write graph, requirement table, priorities and an uncertainty-budget figure."""

    output_dir.mkdir(parents=True, exist_ok=True)
    graph_path = output_dir / "evidence_graph_v0.8.json"
    requirements_path = output_dir / "evidence_requirements_v0.8.csv"
    priorities_path = output_dir / "measurement_priorities_v0.8.csv"
    figure_path = output_dir / "evidence_uncertainty_budget_v0.8.svg"
    graph_path.write_text(synthesis.model_dump_json(indent=2), encoding="utf-8")

    with requirements_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "requirement",
                "weight",
                "evidence_state",
                "satisfaction",
                "uncertainty_contribution",
                "uncertainty_share",
                "reason",
                "next_measurement",
            ],
        )
        writer.writeheader()
        for row in synthesis.requirements:
            writer.writerow(
                {
                    "requirement": row.requirement,
                    "weight": row.weight,
                    "evidence_state": row.evidence_state,
                    "satisfaction": row.satisfaction,
                    "uncertainty_contribution": row.uncertainty_contribution,
                    "uncertainty_share": row.uncertainty_share,
                    "reason": row.reason,
                    "next_measurement": row.next_measurement or "",
                }
            )

    prioritised = sorted(
        (row for row in synthesis.requirements if row.next_measurement is not None),
        key=lambda item: (-item.uncertainty_share, item.requirement),
    )
    with priorities_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "requirement",
                "uncertainty_share",
                "evidence_state",
                "next_measurement",
            ],
        )
        writer.writeheader()
        for rank, row in enumerate(prioritised, start=1):
            writer.writerow(
                {
                    "rank": rank,
                    "requirement": row.requirement,
                    "uncertainty_share": row.uncertainty_share,
                    "evidence_state": row.evidence_state,
                    "next_measurement": row.next_measurement,
                }
            )

    write_horizontal_bar_svg(
        [row.requirement.replace("_", " ") for row in prioritised],
        [row.uncertainty_share for row in prioritised],
        figure_path,
        title="Unresolved evidence contribution to absolute RCC uncertainty",
        x_label="Share of unresolved weighted evidence burden",
        maximum=max(row.uncertainty_share for row in prioritised) * 1.15,
    )
    return {
        "graph": str(graph_path),
        "requirements": str(requirements_path),
        "priorities": str(priorities_path),
        "figure": str(figure_path),
    }


def synthesize_repository_evidence(repository_root: Path, output_dir: Path) -> dict[str, Any]:
    """Build and persist the complete evidence synthesis package."""

    synthesis = build_evidence_synthesis(repository_root)
    paths = write_evidence_synthesis(synthesis, output_dir)
    return {"summary": synthesis.model_dump(mode="json"), "paths": paths}
