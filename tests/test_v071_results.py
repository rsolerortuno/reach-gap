from __future__ import annotations

import json
import shutil
from pathlib import Path

from reach_gap.v071_results import build_claims, validate_bundle

REPO_ROOT = Path(__file__).resolve().parents[1]
BUNDLE_ROOT = REPO_ROOT / "results/external_validation/v0.7.1"


def test_committed_v071_bundle_passes() -> None:
    result = validate_bundle(BUNDLE_ROOT)
    assert result.ok, result.issues
    assert result.summary["perfusion"]["fields"] == 4
    assert result.summary["breast_xenium_erbb2"]["cells"] == 679197
    assert result.summary["absolute_reachability"]["status"] == "NOT_COMPUTED"


def test_v071_claims_preserve_abstention() -> None:
    result = validate_bundle(BUNDLE_ROOT)
    claims = build_claims(result.summary)
    assert "Absolute RCC reachable_fraction" in claims["unsupported"]


def test_v071_tamper_detection(tmp_path: Path) -> None:
    target = tmp_path / "v0.7.1"
    shutil.copytree(BUNDLE_ROOT, target)
    sample_csv = target / "breast_xenium_erbb2/breast_erbb2_sample_summary.csv"
    sample_csv.write_text(sample_csv.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    result = validate_bundle(target)
    assert not result.ok
    assert any(issue.code in {"size_mismatch", "hash_mismatch"} for issue in result.issues)


def test_v071_summary_is_json_serializable() -> None:
    result = validate_bundle(BUNDLE_ROOT)
    json.dumps(result.summary)
