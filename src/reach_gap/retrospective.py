"""Clinical retrospective schema and deliberately non-default synthetic null checks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

REQUIRED_RETROSPECTIVE_COLUMNS = {
    "programme_id",
    "target",
    "molecule",
    "format",
    "payload",
    "dose",
    "indication",
    "line_of_therapy",
    "target_expression_evidence",
    "outcome",
    "reason_for_discontinuation",
    "source_identifiers",
    "expression_reach_gap",
}


@dataclass(frozen=True)
class RetrospectiveResult:
    """Exploratory difference and permutation result."""

    status: str
    n_programmes: int
    median_difference: float | None
    permutation_p: float | None
    stratified_difference: float | None
    note: str


def validate_retrospective_table(table: pd.DataFrame) -> None:
    """Validate required fields without using labels for model construction."""

    missing = REQUIRED_RETROSPECTIVE_COLUMNS.difference(table.columns)
    if missing:
        raise ValueError(f"Missing retrospective columns: {sorted(missing)}")


def run_retrospective(
    table: pd.DataFrame,
    *,
    permutations: int = 1000,
    seed: int = 17,
    allow_synthetic: bool = False,
) -> RetrospectiveResult:
    """Run the preregistered exploratory comparison on an already-curated table."""

    validate_retrospective_table(table)
    if table["programme_id"].astype(str).str.startswith("SYNTHETIC_").any() and not allow_synthetic:
        return RetrospectiveResult(
            status="NOT_COMPUTED",
            n_programmes=0,
            median_difference=None,
            permutation_p=None,
            stratified_difference=None,
            note="Synthetic schema rows are excluded from scientific analysis defaults.",
        )
    filtered = table[table["outcome"].isin(["success", "failure"])].copy()
    if len(filtered) < 6 or filtered["outcome"].nunique() < 2:
        return RetrospectiveResult(
            status="NOT_COMPUTED",
            n_programmes=len(filtered),
            median_difference=None,
            permutation_p=None,
            stratified_difference=None,
            note="Insufficient verified programmes for the preregistered comparison.",
        )
    values = filtered["expression_reach_gap"].to_numpy(dtype=np.float64)
    labels = (filtered["outcome"] == "success").to_numpy(dtype=np.bool_)
    observed = float(np.median(values[labels]) - np.median(values[~labels]))
    rng = np.random.default_rng(seed)
    null = np.empty(permutations, dtype=np.float64)
    for index in range(permutations):
        shuffled = rng.permutation(labels)
        null[index] = np.median(values[shuffled]) - np.median(values[~shuffled])
    p_value = float((1 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1))

    strata = filtered.groupby(["format", "indication"], dropna=False)
    differences: list[float] = []
    for _, group in strata:
        group_labels = (group["outcome"] == "success").to_numpy(dtype=np.bool_)
        if len(group) >= 4 and group_labels.any() and (~group_labels).any():
            group_values = group["expression_reach_gap"].to_numpy(dtype=np.float64)
            differences.append(
                float(
                    np.median(group_values[group_labels]) - np.median(group_values[~group_labels])
                )
            )
    stratified = float(np.median(differences)) if differences else None
    return RetrospectiveResult(
        status="EXPLORATORY_ONLY",
        n_programmes=len(filtered),
        median_difference=observed,
        permutation_p=p_value,
        stratified_difference=stratified,
        note="Outcome labels were used only in this evaluation function, never for parameters.",
    )
