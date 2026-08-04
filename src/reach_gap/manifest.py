"""Manifest validation utilities."""

from __future__ import annotations

import json
from pathlib import Path

from reach_gap.schemas import ProvenanceManifest


def validate_manifest(path: Path) -> ProvenanceManifest:
    """Validate a JSON manifest and return the typed model."""

    raw = json.loads(path.read_text(encoding="utf-8"))
    return ProvenanceManifest.model_validate(raw)
