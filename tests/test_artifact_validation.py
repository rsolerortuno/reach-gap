from __future__ import annotations

from pathlib import Path

from reach_gap.artifact_validation import audit_artifacts, detect_artifact_kind, validate_artifact


def test_content_signature_rejects_html_saved_as_pdf(tmp_path: Path) -> None:
    fake_pdf = tmp_path / "paper.pdf"
    fake_pdf.write_text("<!doctype html><html><body>blocked</body></html>", encoding="utf-8")
    assert detect_artifact_kind(fake_pdf) == "html"
    result = validate_artifact(fake_pdf, "pdf")
    assert result["valid"] is False
    assert result["actual_kind"] == "html"


def test_artifact_audit_accepts_real_signatures_and_writes_report(tmp_path: Path) -> None:
    pdf = tmp_path / "real.pdf"
    pdf.write_bytes(b"%PDF-1.7\nfixture")
    archive = tmp_path / "archive.zip"
    archive.write_bytes(b"PK\x03\x04fixture")
    output = tmp_path / "audit.json"
    result = audit_artifacts(
        {"paper": (pdf, "pdf"), "archive": (archive, "zip")},
        output,
    )
    assert result["status"] == "PASS"
    assert output.exists()
