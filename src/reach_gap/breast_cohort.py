"""Audit independent Xenium breast cell-group cohorts before full spatial extraction."""

from __future__ import annotations

import csv
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd


def _broad_group(value: object) -> str:
    text = str(value)
    return text.split("_", 1)[1] if "_" in text else text


def audit_breast_cell_groups(sample_files: Mapping[str, Path], output_dir: Path) -> dict[str, Any]:
    """Summarise provider cell groups without claiming spatial expression analysis."""

    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    for sample, path in sample_files.items():
        with path.open(newline="", encoding="utf-8") as handle:
            records = [
                {"cell_id": str(row["cell_id"]), "group": str(row["group"])}
                for row in csv.DictReader(handle)
            ]
        table = pd.DataFrame.from_records(records, columns=["cell_id", "group"])
        if table["cell_id"].duplicated().any():
            raise ValueError(f"Duplicate cell identifiers in {sample}")
        table["broad_group"] = table["group"].map(_broad_group)
        for group, count in table["broad_group"].value_counts().items():
            rows.append(
                {
                    "sample": sample,
                    "broad_group": str(group),
                    "cells": int(count),
                    "fraction": float(count / len(table)),
                }
            )
        sample_rows.append(
            {
                "sample": sample,
                "cells": len(table),
                "provider_groups": int(table["group"].nunique()),
                "broad_groups": int(table["broad_group"].nunique()),
            }
        )
    if not rows:
        raise ValueError("No breast cell-group files were supplied")
    group_table = pd.DataFrame(rows)
    sample_table = pd.DataFrame(sample_rows)
    group_table.to_csv(output_dir / "breast_cell_group_composition.csv", index=False)
    sample_table.to_csv(output_dir / "breast_sample_summary.csv", index=False)
    result: dict[str, Any] = {
        "status": "INDEPENDENT_BREAST_COHORT_CELL_GROUP_AUDIT_SPATIAL_EXPRESSION_NOT_COMPUTED",
        "samples": len(sample_table),
        "cells": int(sample_table["cells"].sum()),
        "sample_summaries": sample_rows,
        "spatial_expression": {
            "status": "NOT_COMPUTED",
            "reasons": [
                "The Explorer output bundles were not extracted in this runtime",
                "Cell-group CSVs provide labels but not coordinates or expression values",
                "HER2 state is section-level metadata, not calibrated cell-surface protein",
            ],
        },
    }
    (output_dir / "breast_cell_group_audit.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result
