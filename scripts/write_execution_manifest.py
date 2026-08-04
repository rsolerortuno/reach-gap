"""Write checksums and execution metadata for committed simulation artefacts."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy

from reach_gap.geometry import GeometryConfig
from reach_gap.schemas import ModelParameters

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "simulated"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    """Record deterministic settings, versions and hashes."""

    names = (
        "features.npz",
        "solution.npz",
        "index.json",
        "claims.json",
        "benchmark.json",
        "benchmark_report.md",
    )
    geometry = GeometryConfig(
        size=42,
        cell_count=900,
        vessel_count=5,
        stroma_level=0.65,
        antigen_level=0.80,
        seed=17,
    )
    payload = {
        "status": "COMPUTED_SIMULATION_ONLY",
        "seed": 17,
        "geometry_config": geometry.__dict__,
        "parameters": ModelParameters().as_dict(),
        "index_draws": 48,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "artefact_sha256": {name: _sha256(RESULTS / name) for name in names},
    }
    (RESULTS / "execution.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
