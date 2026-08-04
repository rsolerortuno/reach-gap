"""Small deterministic report renderer."""

from __future__ import annotations

import json
from pathlib import Path


def render_report(benchmark_path: Path, output_path: Path) -> Path:
    """Render a Markdown summary from a benchmark result document."""

    result = json.loads(benchmark_path.read_text(encoding="utf-8"))
    comparison = result["model_comparison"]
    text = f"""# reach-gap benchmark report

Status: **{result["status"]}**

## Simulation-only model comparison

- Mechanistic RMSE: {comparison["mechanistic"]["rmse"]:.4f}
- Naive weighted-sum RMSE: {comparison["naive_weighted_sum"]["rmse"]:.4f}
- Distance-only RMSE: {comparison["distance_only"]["rmse"]:.4f}

## Interpretation

These values are computed on deterministic simulated tissue generated from the same model
family. They validate implementation and ablations, not real-tissue predictive validity.
The clinical retrospective status is **{result["retrospective"]["status"]}**.
"""
    output_path.write_text(text, encoding="utf-8")
    return output_path
