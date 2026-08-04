"""Content-signature validation for externally downloaded scientific artifacts.

Filename extensions and HTTP status codes are not sufficient evidence that an
artifact is usable.  This module detects common error payloads saved under PDF,
ZIP, GZIP, and TIFF filenames and records an explicit validation status.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def detect_artifact_kind(path: Path) -> str:
    """Detect a small set of scientific file types from content signatures."""

    with path.open("rb") as handle:
        prefix = handle.read(512)
    stripped = prefix.lstrip().lower()
    if prefix.startswith(b"%PDF-"):
        return "pdf"
    if prefix.startswith((b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08")):
        return "zip"
    if prefix.startswith(b"\x1f\x8b"):
        return "gzip"
    if prefix.startswith((b"II*\x00", b"MM\x00*")):
        return "tiff"
    if stripped.startswith((b"<!doctype html", b"<html")):
        return "html"
    if stripped.startswith(b"<?xml") or stripped.startswith(b"<errorbean"):
        return "xml"
    if stripped.startswith((b"{", b"[")):
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            pass
        else:
            return "json"
    return "unknown"


def validate_artifact(path: Path, expected_kind: str) -> dict[str, Any]:
    """Validate one artifact against an expected content type."""

    actual_kind = detect_artifact_kind(path)
    return {
        "path": str(path),
        "expected_kind": expected_kind,
        "actual_kind": actual_kind,
        "valid": actual_kind == expected_kind,
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def audit_artifacts(
    artifacts: Mapping[str, tuple[Path, str]], output_path: Path | None = None
) -> dict[str, Any]:
    """Audit named artifacts and optionally write a JSON report."""

    records = {
        name: validate_artifact(Path(path), expected_kind)
        for name, (path, expected_kind) in artifacts.items()
    }
    result = {
        "status": "PASS" if all(record["valid"] for record in records.values()) else "FAIL",
        "artifacts": records,
        "invalid_artifacts": sorted(
            name for name, record in records.items() if not record["valid"]
        ),
    }
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result
