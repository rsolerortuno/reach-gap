from __future__ import annotations

import json
from pathlib import Path

from reach_gap.v080_results import (
    REQUIRED_PATHS,
    build_v080_package,
    validate_v080_package,
)


def repository_root() -> Path:
    return Path(__file__).resolve().parents[1]


def test_build_and_validate_v080_package(tmp_path: Path) -> None:
    payload = build_v080_package(repository_root(), tmp_path, draws=20_000, seed=17)
    assert payload["summary"]["version"] == "0.8.0"
    validation = validate_v080_package(tmp_path)
    assert validation.ok, validation.issues
    assert validation.checks >= 60
    stored_validation = json.loads((tmp_path / "validation_v0.8.json").read_text())
    assert stored_validation["checks"] == validation.checks
    assert stored_validation["status"] == validation.status
    assert validation.summary["absolute_readiness_score"] == 40.5
    assert validation.summary["stable_top_target"] == "VISTA"
    assert validation.summary["stable_top_probability"] >= 0.97
    for relative in REQUIRED_PATHS:
        assert (tmp_path / relative).is_file()
    manifest = json.loads((tmp_path / "artifact_manifest_v0.8.json").read_text())
    assert len(manifest["files"]) >= len(REQUIRED_PATHS)


def test_v080_validation_detects_tampering(tmp_path: Path) -> None:
    build_v080_package(repository_root(), tmp_path, draws=20_000, seed=17)
    target = tmp_path / "claims_v0.8.json"
    target.write_text("{}", encoding="utf-8")
    validation = validate_v080_package(tmp_path)
    assert not validation.ok
    assert any(
        issue.code in {"claims_version", "claim_abstention", "manifest_size", "manifest_hash"}
        for issue in validation.issues
    )
