from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from reach_gap.evidence_synthesis import (
    build_evidence_synthesis,
    synthesize_repository_evidence,
)


def repository_root() -> Path:
    return Path(__file__).resolve().parents[1]


def test_build_evidence_synthesis_is_truthful() -> None:
    result = build_evidence_synthesis(repository_root())
    assert result.version == "0.8.0"
    assert result.absolute_readiness_score == pytest.approx(40.5)
    assert result.same_tissue_requirements_satisfied == 3
    assert result.requirements_total == 8
    assert result.dominant_uncertainty == "administered_antibody_field"
    assert result.measurement_priority[:3] == [
        "administered_antibody_field",
        "functional_perfusion",
        "surface_antigen_calibration",
    ]
    assert set(result.absolute_outputs.values()) == {"NOT_COMPUTED"}
    assert sum(row.weight for row in result.requirements) == pytest.approx(1.0)
    assert sum(row.uncertainty_share for row in result.requirements) == pytest.approx(1.0)
    blocked = [edge for edge in result.edges if edge.relation == "BLOCKED_TRANSFER"]
    assert len(blocked) == 3
    assert all(not edge.transferable for edge in blocked)


def test_evidence_values_are_artifact_backed() -> None:
    result = build_evidence_synthesis(repository_root())
    nodes = {node.node_id: node for node in result.nodes}
    assert nodes["rcc_target_localisation"].measured_values["cells"] == 465534
    assert nodes["external_llc_perfusion"].measured_values["fields"] == 4
    assert nodes["external_breast_erbb2_control"].measured_values["cells"] == 679197
    assert nodes["external_her2_calibration"].measured_values["log_log_r2"] == pytest.approx(
        0.9906341282693898
    )


def test_synthesize_repository_evidence_writes_package(tmp_path: Path) -> None:
    payload = synthesize_repository_evidence(repository_root(), tmp_path)
    assert payload["summary"]["status"].endswith("NOT_COMPUTED")
    graph = json.loads((tmp_path / "evidence_graph_v0.8.json").read_text(encoding="utf-8"))
    assert graph["absolute_readiness_score"] == pytest.approx(40.5)
    with (tmp_path / "measurement_priorities_v0.8.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["requirement"] == "administered_antibody_field"
    assert (tmp_path / "evidence_uncertainty_budget_v0.8.svg").is_file()


def test_missing_source_artifact_fails_explicitly(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        build_evidence_synthesis(tmp_path)
