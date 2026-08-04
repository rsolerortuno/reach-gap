"""Reproducible real-data preparation for the compact public RCC Xenium package.

The workflow computes relative molecular geometry over every cell while refusing to
convert fluorescence to receptor density or endothelial presence to perfusion.
"""

from __future__ import annotations

import json
import platform
import resource
import time
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from reach_gap.segmentation import boundary_polygon_metrics, segmentation_robustness_summary
from reach_gap.xenium import (
    assign_pathology_regions,
    extract_selected_h5_features,
    geojson_vertices,
    infer_annotation_transform,
    load_affine_matrix,
    otsu_threshold,
    read_10x_h5_catalog,
    read_cells,
    resolve_feature_indices,
    robust_scale,
    score_cells,
    sha256_file,
    summarise_h5_features,
    transform_geojson,
    validate_cell_barcode_identity,
    write_partitioned_table,
)

_TARGET_PAIRS: dict[str, tuple[str, str]] = {
    "PD-L1": ("protein__PD_L1", "rna__CD274"),
    "VISTA": ("protein__VISTA", "rna__VSIR"),
    "PD-1": ("protein__PD_1", "rna__PDCD1"),
    "LAG-3": ("protein__LAG_3", "rna__LAG3"),
}


def _find_one(input_dir: Path, names: tuple[str, ...], *, required: bool = True) -> Path | None:
    matches: list[Path] = []
    for name in names:
        matches.extend(path for path in input_dir.rglob(name) if path.is_file())
        matches.extend(path for path in input_dir.rglob(f"{name}.bin") if path.is_file())
    unique = sorted(set(matches))
    if len(unique) > 1:
        raise ValueError(f"Ambiguous input for {names}: {unique}")
    if not unique:
        if required:
            raise FileNotFoundError(f"Missing required input; expected one of {names}")
        return None
    return unique[0]


def target_spatial_summary(scored: pd.DataFrame, diagnostics: dict[str, Any]) -> pd.DataFrame:
    """Summarise relative target positivity and geometry without claiming reachability."""

    rows: list[dict[str, object]] = []
    pathology = scored.get("pathology_region", pd.Series("UNAVAILABLE", index=scored.index))
    tumour_mask = pathology.astype(str).str.contains("Tumor", case=False, regex=False)
    if not tumour_mask.any():
        tumour_mask = scored["in_molecular_tumour_neighbourhood"].astype(bool)
    distance = scored["distance_to_vessel_um"].to_numpy(dtype=np.float64)
    for target, details_any in dict(diagnostics["targets"]).items():
        details = dict(details_any)
        safe = target.replace("-", "_")
        signal = scored[f"target__{safe}__signal"].to_numpy(dtype=np.float64)
        positive = scored[f"target__{safe}__positive"].to_numpy(dtype=np.bool_)
        subsets = {
            "all_cells": np.ones(len(scored), dtype=np.bool_),
            "tumour_region": tumour_mask.to_numpy(dtype=np.bool_),
            "malignant_proxy": scored["cell_is_malignant_proxy"].to_numpy(dtype=np.bool_),
            "immune_high": scored["immune_score"].to_numpy(dtype=np.float64) >= 0.5,
        }
        for subset_name, subset in subsets.items():
            active = subset & positive
            n_subset = int(np.sum(subset))
            n_positive = int(np.sum(active))
            rows.append(
                {
                    "target": target,
                    "measurement": details["measurement"],
                    "subset": subset_name,
                    "cells": n_subset,
                    "positive_cells": n_positive,
                    "positive_fraction": n_positive / max(n_subset, 1),
                    "median_signal": float(np.median(signal[subset])) if n_subset else np.nan,
                    "median_distance_to_vessel_um_positive": (
                        float(np.median(distance[active])) if n_positive else np.nan
                    ),
                    "positive_within_25um": (
                        float(np.mean(distance[active] <= 25.0)) if n_positive else np.nan
                    ),
                    "positive_within_50um": (
                        float(np.mean(distance[active] <= 50.0)) if n_positive else np.nan
                    ),
                    "positive_within_100um": (
                        float(np.mean(distance[active] <= 100.0)) if n_positive else np.nan
                    ),
                    "median_local_caf_positive": (
                        float(
                            np.median(
                                cast(
                                    "pd.Series[Any]", scored.loc[active, "local_caf_score"]
                                ).to_numpy(dtype=np.float64)
                            )
                        )
                        if n_positive
                        else np.nan
                    ),
                    "median_local_ecm_positive": (
                        float(
                            np.median(
                                cast(
                                    "pd.Series[Any]", scored.loc[active, "local_ecm_score"]
                                ).to_numpy(dtype=np.float64)
                            )
                        )
                        if n_positive
                        else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def rna_protein_concordance(scored: pd.DataFrame) -> pd.DataFrame:
    """Measure cell-level RNA/protein agreement for targets measured in both modalities."""

    rows: list[dict[str, object]] = []
    for target, (protein_col, rna_col) in _TARGET_PAIRS.items():
        if protein_col not in scored or rna_col not in scored:
            continue
        protein = scored[protein_col].to_numpy(dtype=np.float64)
        rna = scored[rna_col].to_numpy(dtype=np.float64)
        correlation = spearmanr(protein, rna).statistic
        protein_positive = protein > 0
        rna_positive = rna > 0
        union = int(np.sum(protein_positive | rna_positive))
        rows.append(
            {
                "target": target,
                "spearman_raw_signal": float(correlation),
                "protein_nonzero_fraction": float(np.mean(protein_positive)),
                "rna_nonzero_fraction": float(np.mean(rna_positive)),
                "nonzero_jaccard": int(np.sum(protein_positive & rna_positive)) / max(union, 1),
                "protein_positive_rna_zero_fraction": float(
                    np.mean(protein_positive & ~rna_positive)
                ),
                "rna_positive_protein_zero_fraction": float(
                    np.mean(rna_positive & ~protein_positive)
                ),
            }
        )
    return pd.DataFrame(rows)


def pathology_molecular_summary(scored: pd.DataFrame) -> pd.DataFrame:
    """Summarise marker-derived geometry by supplied pathology region."""

    if "pathology_region" not in scored:
        return pd.DataFrame()
    columns = [
        "endothelial_score",
        "vessel_positive",
        "distance_to_vessel_um",
        "caf_score",
        "local_caf_score",
        "immune_score",
        "epithelial_malignant_score",
    ]
    rows: list[dict[str, object]] = []
    for region, group in scored.groupby("pathology_region", dropna=False):
        row: dict[str, object] = {"pathology_region": str(region), "cells": len(group)}
        for column in columns:
            values = group[column].to_numpy(dtype=np.float64)
            row[f"{column}_median"] = float(np.nanmedian(values))
            row[f"{column}_q90"] = float(np.nanquantile(values, 0.90))
        rows.append(row)
    return pd.DataFrame(rows).sort_values("cells", ascending=False)


def _nearest_distance(coordinates: np.ndarray, mask: np.ndarray) -> np.ndarray:
    if not np.any(mask):
        return np.full(len(coordinates), np.nan, dtype=np.float64)
    tree = cKDTree(coordinates[mask])
    output = np.empty(len(coordinates), dtype=np.float64)
    for start in range(0, len(coordinates), 50_000):
        stop = min(start + 50_000, len(coordinates))
        output[start:stop] = tree.query(coordinates[start:stop], k=1, workers=-1)[0]
    return output


def vessel_definition_sensitivity(
    scored: pd.DataFrame, expression: pd.DataFrame
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Recompute distances under six transparent endothelial definitions."""

    needed = ["cell_id", "protein__CD31", "rna__PECAM1", "rna__VWF", "rna__RAMP2"]
    missing = sorted(set(needed).difference(expression.columns))
    if missing:
        raise KeyError(f"Vessel sensitivity features are missing: {missing}")
    data = scored.merge(expression[needed], on="cell_id", how="left", validate="one_to_one")
    coordinates = data[["x_um", "y_um"]].to_numpy(dtype=np.float64)
    cd31 = robust_scale(data["protein__CD31"].to_numpy(dtype=np.float64))
    endothelial = data["endothelial_score"].to_numpy(dtype=np.float64)
    cd31_otsu = float(np.clip(otsu_threshold(cd31), 0.10, 0.90))
    endothelial_otsu = float(np.clip(otsu_threshold(endothelial), 0.10, 0.90))
    pecam1 = data["rna__PECAM1"].to_numpy(dtype=np.float64) > 0
    vwf = data["rna__VWF"].to_numpy(dtype=np.float64) > 0
    ramp2 = data["rna__RAMP2"].to_numpy(dtype=np.float64) > 0
    rna_count = pecam1.astype(np.int8) + vwf.astype(np.int8) + ramp2.astype(np.int8)
    definitions = {
        "inclusive_endothelial_or_pericyte": data["vessel_positive"].to_numpy(dtype=np.bool_),
        "endothelial_score_otsu": endothelial >= endothelial_otsu,
        "cd31_otsu": cd31 >= cd31_otsu,
        "cd31_high_robust_0_5": cd31 >= 0.5,
        "cd31_plus_any_endothelial_rna": (cd31 >= cd31_otsu) & (rna_count >= 1),
        "cd31_plus_two_endothelial_rna": (cd31 >= cd31_otsu) & (rna_count >= 2),
    }
    tumour = (
        data.get("pathology_region", pd.Series("", index=data.index))
        .astype(str)
        .str.contains("Tumor", case=False, regex=False)
        .to_numpy()
    )
    rows: list[dict[str, object]] = []
    for name, mask in definitions.items():
        distance = _nearest_distance(coordinates, mask)
        base: dict[str, object] = {
            "definition": name,
            "vessel_cells": int(mask.sum()),
            "vessel_fraction": float(mask.mean()),
            "all_median_um": float(np.nanmedian(distance)),
            "all_q90_um": float(np.nanquantile(distance, 0.90)),
            "tumour_median_um": float(np.nanmedian(distance[tumour])),
            "tumour_q90_um": float(np.nanquantile(distance[tumour], 0.90)),
        }
        for target in _TARGET_PAIRS:
            positive = (
                data[f"target__{target.replace('-', '_')}__positive"].to_numpy(dtype=np.bool_)
                & tumour
            )
            row = dict(base)
            row.update(
                {
                    "target": target,
                    "target_positive_tumour_cells": int(positive.sum()),
                    "target_median_um": float(np.nanmedian(distance[positive])),
                    "target_q90_um": float(np.nanquantile(distance[positive], 0.90)),
                    "within_25um": float(np.mean(distance[positive] <= 25.0)),
                    "within_50um": float(np.mean(distance[positive] <= 50.0)),
                    "within_100um": float(np.mean(distance[positive] <= 100.0)),
                }
            )
            rows.append(row)
    frame = pd.DataFrame(rows)
    summary = {
        "status": "COMPUTED_RELATIVE_GEOMETRY_ONLY",
        "cd31_robust_otsu_threshold": cd31_otsu,
        "endothelial_score_otsu_threshold": endothelial_otsu,
        "definitions": {
            name: {
                "vessel_cells": int(mask.sum()),
                "vessel_fraction": float(mask.mean()),
            }
            for name, mask in definitions.items()
        },
        "warning": (
            "None of the definitions establishes functional perfusion; strict definitions "
            "may fragment vessel walls into individual endothelial cells."
        ),
    }
    return frame, summary


def prepare_rcc_xenium_essentials(
    input_dir: Path,
    output_dir: Path,
    *,
    annotation_path: Path | None = None,
    alignment_path: Path | None = None,
    write_full_cell_tables: bool = False,
    seed: int = 17,
) -> dict[str, Any]:
    """Prepare all cells from the compact RCC Xenium package with explicit abstention."""

    start = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)
    inputs: dict[str, Path] = {
        "cells": _find_one(input_dir, ("cells.parquet", "cells.csv.gz")),  # type: ignore[dict-item]
        "h5": _find_one(input_dir, ("cell_feature_matrix.h5",)),  # type: ignore[dict-item]
        "gene_panel": _find_one(input_dir, ("gene_panel.json",)),  # type: ignore[dict-item]
        "protein_panel": _find_one(input_dir, ("protein_panel.json",)),  # type: ignore[dict-item]
        "metrics": _find_one(input_dir, ("metrics_summary.csv",)),  # type: ignore[dict-item]
        "experiment": _find_one(input_dir, ("experiment.xenium",)),  # type: ignore[dict-item]
    }
    cell_boundary_path = _find_one(input_dir, ("cell_boundaries.parquet",), required=False)
    nucleus_boundary_path = _find_one(input_dir, ("nucleus_boundaries.parquet",), required=False)

    cells = read_cells(inputs["cells"])
    catalog = read_10x_h5_catalog(inputs["h5"])
    validate_cell_barcode_identity(cells, catalog)
    summarise_h5_features(inputs["h5"], catalog).to_csv(
        output_dir / "feature_summary.csv", index=False
    )
    selected_indices, resolution = resolve_feature_indices(catalog)
    (output_dir / "marker_resolution.json").write_text(
        json.dumps(resolution, indent=2), encoding="utf-8"
    )
    expression = extract_selected_h5_features(inputs["h5"], catalog, selected_indices)
    scored, diagnostics = score_cells(cells, expression, resolution)

    alignment_record: dict[str, Any] | None = None
    if annotation_path is not None and alignment_path is not None:
        alignment = load_affine_matrix(alignment_path)
        accepted, candidates = infer_annotation_transform(
            geojson_vertices(annotation_path), cells, alignment
        )
        candidate_records = [
            {
                "name": candidate.name,
                "matrix": candidate.matrix.tolist(),
                "scale_x": candidate.scale_x,
                "scale_y": candidate.scale_y,
                "score": candidate.score,
                "fraction_inside": candidate.fraction_inside,
                "accepted": accepted is not None and candidate.name == accepted.name,
            }
            for candidate in candidates
        ]
        (output_dir / "annotation_alignment_candidates.json").write_text(
            json.dumps(candidate_records, indent=2), encoding="utf-8"
        )
        if accepted is not None:
            transformed = output_dir / "pathology_annotations_xenium.geojson"
            transform_geojson(annotation_path, transformed, accepted)
            scored = assign_pathology_regions(scored, transformed)
            alignment_record = next(
                record for record in candidate_records if bool(record["accepted"])
            )
        else:
            diagnostics["warnings"].append(
                "Pathology alignment remained ambiguous; pathology summaries were omitted."
            )

    target_spatial_summary(scored, diagnostics).to_csv(
        output_dir / "target_spatial_summary.csv", index=False
    )
    rna_protein_concordance(scored).to_csv(output_dir / "rna_protein_concordance.csv", index=False)
    pathology_molecular_summary(scored).to_csv(
        output_dir / "pathology_molecular_summary.csv", index=False
    )
    vessel_table, vessel_summary = vessel_definition_sensitivity(scored, expression)
    vessel_table.to_csv(output_dir / "vessel_calling_sensitivity.csv", index=False)
    (output_dir / "vessel_calling_sensitivity.json").write_text(
        json.dumps(vessel_summary, indent=2), encoding="utf-8"
    )
    ranges = (
        vessel_table.groupby("target")
        .agg(
            vessel_definition_count=("definition", "nunique"),
            target_median_um_min=("target_median_um", "min"),
            target_median_um_max=("target_median_um", "max"),
            within_50um_min=("within_50um", "min"),
            within_50um_max=("within_50um", "max"),
            within_100um_min=("within_100um", "min"),
            within_100um_max=("within_100um", "max"),
        )
        .reset_index()
    )
    ranges.to_csv(output_dir / "vessel_robustness_ranges.csv", index=False)

    segmentation_summary: dict[str, object] | None = None
    if cell_boundary_path is not None and nucleus_boundary_path is not None:
        cell_metrics = boundary_polygon_metrics(cell_boundary_path)
        nucleus_metrics = boundary_polygon_metrics(nucleus_boundary_path)
        segmentation, segmentation_summary = segmentation_robustness_summary(
            cells, cell_metrics, nucleus_metrics
        )
        segmentation.sample(n=min(20_000, len(segmentation)), random_state=seed).sort_values(
            "cell_id"
        ).to_csv(
            output_dir / "segmentation_boundary_metrics_sample.csv.gz",
            index=False,
            compression="gzip",
        )
        (output_dir / "segmentation_robustness.json").write_text(
            json.dumps(segmentation_summary, indent=2), encoding="utf-8"
        )

    compact_columns = [
        "cell_id",
        "x_um",
        "y_um",
        "cell_area",
        "nucleus_area",
        "endothelial_score",
        "pericyte_smooth_muscle_score",
        "caf_score",
        "ecm_score",
        "immune_score",
        "epithelial_malignant_score",
        "vessel_positive",
        "distance_to_vessel_um",
        "local_vessel_density",
        "local_caf_score",
        "local_ecm_score",
        "cell_is_malignant_proxy",
        "in_molecular_tumour_neighbourhood",
    ]
    compact_columns.extend(column for column in scored.columns if column.startswith("target__"))
    if "pathology_region" in scored:
        compact_columns.append("pathology_region")
    compact = scored[compact_columns]
    compact.sample(n=min(20_000, len(compact)), random_state=seed).sort_values("cell_id").to_csv(
        output_dir / "scored_cells_sample.csv.gz", index=False, compression="gzip"
    )
    full_parts: list[str] = []
    if write_full_cell_tables:
        full_parts = write_partitioned_table(
            compact, output_dir / "cell_tables", "rcc_cells_scored", rows_per_part=100_000
        )

    diagnostics["pathology_alignment"] = alignment_record
    diagnostics["ecm_measurement_status"] = (
        "NOT_AVAILABLE_IN_GENE_OR_PROTEIN_PANEL"
        if not resolution["groups"].get("ecm")
        else "MARKER_DERIVED"
    )
    (output_dir / "cell_scoring_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2), encoding="utf-8"
    )
    absolute_reasons = [
        "Xenium protein fluorescence is not calibrated to surface antigen density.",
        "CD31-positive structures are not confirmed to be functionally perfused vessels.",
        "The panel provides no direct ECM marker set; ECM was not imputed.",
    ]
    all_inputs = dict(inputs)
    if annotation_path is not None:
        all_inputs["annotation"] = annotation_path
    if alignment_path is not None:
        all_inputs["alignment"] = alignment_path
    if cell_boundary_path is not None:
        all_inputs["cell_boundaries"] = cell_boundary_path
    if nucleus_boundary_path is not None:
        all_inputs["nucleus_boundaries"] = nucleus_boundary_path
    manifest: dict[str, Any] = {
        "schema_version": "1.1",
        "status": "REAL_XENIUM_MOLECULAR_GEOMETRY_PREPARED_ABSOLUTE_INDEX_NOT_COMPUTED",
        "dataset": "Xenium gene and protein expression, FFPE human renal cell carcinoma",
        "cells": len(cells),
        "features": catalog.matrix_shape[0],
        "selected_features": len(selected_indices),
        "matrix_shape": list(catalog.matrix_shape),
        "pathology_alignment": alignment_record,
        "segmentation_robustness": segmentation_summary,
        "input_files": {
            name: {
                "name": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for name, path in all_inputs.items()
        },
        "outputs": {
            "full_cell_parts": [str(Path(path).relative_to(output_dir)) for path in full_parts]
        },
        "vessel_calling_sensitivity": {
            "definitions": int(vessel_table["definition"].nunique()),
            "vessel_fraction_range": [
                float(vessel_table["vessel_fraction"].min()),
                float(vessel_table["vessel_fraction"].max()),
            ],
        },
        "runtime_seconds": time.perf_counter() - start,
        "max_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
        "python": platform.python_version(),
        "absolute_index": {"status": "NOT_COMPUTED", "reasons": absolute_reasons},
    }
    (output_dir / "processing_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    claims = {
        "permitted": [
            "Every Xenium cell was joined exactly to the complete cell-feature matrix.",
            (
                "Relative target signals, pathology regions and distances to marker-derived "
                "endothelial proxies were computed in the same section."
            ),
            "Vessel-calling and segmentation sensitivity were quantified.",
        ],
        "conditional": [
            (
                "Distance results are conditional on marker-derived vessel definitions and "
                "do not establish perfusion."
            ),
            (
                "Protein-positive fractions are within-section thresholds, not clinical assay "
                "positivity rates."
            ),
        ],
        "unsupported": [
            "Absolute receptor density, receptor occupancy or antibody concentration.",
            "A therapeutic-dose reachable_fraction or expression_reach_gap.",
            "Clinical efficacy prediction.",
        ],
        "abstention_reasons": absolute_reasons,
    }
    (output_dir / "claims.json").write_text(json.dumps(claims, indent=2), encoding="utf-8")
    return manifest
