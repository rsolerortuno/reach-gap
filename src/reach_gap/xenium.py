"""Low-memory preparation of Xenium gene-and-protein data for reach-gap.

This module intentionally prepares auditable cell-level inputs. It does not convert
Xenium protein intensity or RNA abundance into absolute surface antigen density.
Any downstream mechanistic index must therefore abstain unless an independent
calibration is supplied.
"""

from __future__ import annotations

import hashlib
import json
import math
import shutil
import zipfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeAlias, cast

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray
from scipy.spatial import cKDTree

FloatArray: TypeAlias = NDArray[np.float64]


@dataclass(frozen=True)
class XeniumFeatureCatalog:
    """Feature metadata and matrix dimensions from a 10x HDF5 matrix."""

    feature_ids: tuple[str, ...]
    feature_names: tuple[str, ...]
    feature_types: tuple[str, ...]
    barcodes: tuple[str, ...]
    matrix_shape: tuple[int, int]


@dataclass(frozen=True)
class AlignmentCandidate:
    """One candidate mapping of annotation coordinates to Xenium microns."""

    name: str
    matrix: NDArray[np.float64]
    scale_x: float
    scale_y: float
    score: float
    fraction_inside: float


DEFAULT_MARKERS: dict[str, tuple[str, ...]] = {
    "endothelial": (
        "CD31",
        "PECAM1",
        "VWF",
        "EMCN",
        "KDR",
        "ENG",
        "RAMP2",
        "PLVAP",
    ),
    "pericyte_smooth_muscle": (
        "alphaSMA",
        "ACTA2",
        "RGS5",
        "CSPG4",
        "MCAM",
        "PDGFRB",
        "DES",
    ),
    "caf": (
        "alphaSMA",
        "Vimentin",
        "VIM",
        "FAP",
        "PDGFRA",
        "PDGFRB",
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "DCN",
        "LUM",
        "SPARC",
        "FN1",
    ),
    "ecm": (
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "COL4A1",
        "COL4A2",
        "COL6A1",
        "COL6A2",
        "FN1",
        "LAMA4",
        "LAMB1",
        "DCN",
        "LUM",
        "SPARC",
    ),
    "epithelial_malignant": (
        "PanCK",
        "E-cadherin",
        "Beta-catenin",
        "EPCAM",
        "KRT7",
        "KRT8",
        "KRT18",
        "KRT19",
        "CA9",
        "PAX8",
        "KIM1",
        "HAVCR1",
    ),
    "immune": (
        "CD45",
        "PTPRC",
        "CD3E",
        "CD4",
        "CD8A",
        "CD20",
        "MS4A1",
        "CD68",
        "CD163",
        "CD11c",
        "ITGAX",
    ),
}

DEFAULT_TARGET_ALIASES: dict[str, tuple[str, ...]] = {
    "PD-L1": ("PD-L1", "CD274"),
    "VISTA": ("VISTA", "VSIR"),
    "PD-1": ("PD-1", "PDCD1"),
    "LAG-3": ("LAG-3", "LAG3"),
}

EXPECTED_RCC_FILES: dict[str, dict[str, Any]] = {
    "Xenium_V1_Human_Kidney_FFPE_Protein_updated_outs.zip": {
        "size": 36_149_509_228,
        "md5": "76d46bac8060f8bc3ecb450e03b4f3f6",
    },
    "Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_image.ome.tif": {
        "size": 3_720_697_771,
        "md5": "96ad5f699c7d6280cdf6af1c13f39515",
    },
    "Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_imagealignment.csv": {
        "size": 129,
        "md5": "e78bbee6561b9c037cd3eb839f63272",
    },
    "Xenium_V1_Human_Kidney_FFPE_Protein_updated_annotation.geojson": {
        "size": 65_584,
        "md5": "b5e848d7147f25817568d5871592eb56",
    },
}


def _decode(values: NDArray[Any]) -> tuple[str, ...]:
    decoded: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        else:
            decoded.append(str(value))
    return tuple(decoded)


def md5_file(path: Path, *, chunk_size: int = 16 * 1024 * 1024) -> str:
    """Compute a streaming MD5 for comparison with provider checksums."""

    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_file(path: Path, *, chunk_size: int = 16 * 1024 * 1024) -> str:
    """Compute a streaming SHA-256 for the run provenance manifest."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def verify_expected_files(
    raw_dir: Path,
    *,
    verify_large_md5: bool = True,
    expected: Mapping[str, Mapping[str, Any]] = EXPECTED_RCC_FILES,
) -> dict[str, Any]:
    """Validate names, sizes and provider MD5 values for the RCC download."""

    report: dict[str, Any] = {
        "files": {},
        "all_present": True,
        "all_verified": True,
        "all_md5_computed": True,
        "verification_policy": (
            "size_and_provider_md5_for_all_files"
            if verify_large_md5
            else "size_for_large_files_and_provider_md5_below_1GB"
        ),
    }
    for name, specification in expected.items():
        path = raw_dir / name
        entry: dict[str, Any] = {
            "path": str(path),
            "present": path.exists(),
            "expected_size": int(specification["size"]),
            "expected_md5": str(specification["md5"]),
        }
        if not path.exists():
            report["all_present"] = False
            report["all_verified"] = False
            report["files"][name] = entry
            continue
        size = path.stat().st_size
        entry["observed_size"] = size
        entry["size_matches"] = size == int(specification["size"])
        should_hash = verify_large_md5 or size < 1_000_000_000
        if should_hash:
            observed_md5 = md5_file(path)
            entry["observed_md5"] = observed_md5
            entry["md5_matches"] = observed_md5 == str(specification["md5"])
        else:
            entry["observed_md5"] = None
            entry["md5_matches"] = None
            entry["md5_status"] = "NOT_COMPUTED_BY_CONFIGURATION"
            report["all_md5_computed"] = False
        if not entry["size_matches"] or entry.get("md5_matches") is False:
            report["all_verified"] = False
        report["files"][name] = entry
    return report


def inspect_zip(zip_path: Path) -> pd.DataFrame:
    """Return an inventory of every member without extracting the archive."""

    rows: list[dict[str, Any]] = []
    with zipfile.ZipFile(zip_path) as archive:
        for info in archive.infolist():
            rows.append(
                {
                    "member": info.filename,
                    "basename": Path(info.filename).name,
                    "uncompressed_bytes": info.file_size,
                    "compressed_bytes": info.compress_size,
                    "compression": info.compress_type,
                    "crc32": f"{info.CRC:08x}",
                    "is_dir": info.is_dir(),
                }
            )
    return pd.DataFrame(rows)


def _members_by_basename(inventory: pd.DataFrame) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    for basename, group in inventory.groupby("basename", sort=False):
        mapping[str(basename)] = [str(value) for value in group["member"].tolist()]
    return mapping


def select_essential_members(inventory: pd.DataFrame) -> list[str]:
    """Select the minimum full-resolution files needed for cell-level preparation."""

    by_name = _members_by_basename(inventory)
    selected: list[str] = []
    preferred_groups = (
        ("cells.parquet", "cells.csv.gz"),
        ("cell_feature_matrix.h5",),
        ("metrics_summary.csv",),
        ("gene_panel.json",),
        ("protein_panel.json",),
        ("experiment.xenium",),
        ("analysis_summary.html",),
        ("overview_scan.png",),
    )
    for alternatives in preferred_groups:
        for basename in alternatives:
            candidates = by_name.get(basename, [])
            if candidates:
                selected.append(candidates[0])
                break
    # Secondary analysis CSVs are small and useful for QC, but never pull image payloads.
    for row in inventory.itertuples(index=False):
        member = str(row.member)
        if (
            not bool(row.is_dir)
            and member.lower().endswith(".csv")
            and "/analysis/" in f"/{member.lower()}"
            and int(cast(Any, row.uncompressed_bytes)) <= 100_000_000
        ):
            selected.append(member)
    return sorted(set(selected))


def extract_members(zip_path: Path, members: Sequence[str], output_dir: Path) -> dict[str, str]:
    """Extract selected members while preserving their relative paths."""

    output_dir.mkdir(parents=True, exist_ok=True)
    extracted: dict[str, str] = {}
    with zipfile.ZipFile(zip_path) as archive:
        archive_names = set(archive.namelist())
        for member in members:
            if member not in archive_names:
                raise KeyError(f"Archive member not found: {member}")
            destination = output_dir / member
            destination.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(member) as source, destination.open("wb") as target:
                shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
            extracted[member] = str(destination)
    return extracted


def find_extracted(output_dir: Path, basename: str) -> Path | None:
    """Find one extracted file by basename, rejecting ambiguous duplicates."""

    matches = [path for path in output_dir.rglob(basename) if path.is_file()]
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(f"Multiple extracted files named {basename}: {matches}")
    return matches[0]


def read_cells(path: Path) -> pd.DataFrame:
    """Read and validate the Xenium cell summary."""

    if path.suffix == ".parquet" or path.name.endswith(".parquet.bin"):
        try:
            table = pd.read_parquet(path)
        except (ImportError, ModuleNotFoundError):
            from reach_gap.parquet_lite import read_flat_parquet

            table = read_flat_parquet(
                path,
                columns=[
                    "cell_id",
                    "x_centroid",
                    "y_centroid",
                    "cell_area",
                    "nucleus_area",
                ],
            )
    elif path.name.endswith(".csv.gz"):
        table = pd.read_csv(path, compression="gzip")
    else:
        table = pd.read_csv(path)
    required = {"cell_id", "x_centroid", "y_centroid", "cell_area", "nucleus_area"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Xenium cells table is missing: {sorted(missing)}")
    if table["cell_id"].duplicated().any():
        raise ValueError("Xenium cells table contains duplicate cell_id values")
    table = table.copy()
    table["cell_id"] = table["cell_id"].astype(str)
    table["x_um"] = pd.to_numeric(table["x_centroid"], errors="raise")
    table["y_um"] = pd.to_numeric(table["y_centroid"], errors="raise")
    if not np.isfinite(table[["x_um", "y_um"]].to_numpy(dtype=np.float64)).all():
        raise ValueError("Cell coordinates contain non-finite values")
    return table


def read_10x_h5_catalog(path: Path) -> XeniumFeatureCatalog:
    """Read feature metadata from a standard 10x cell-feature HDF5 matrix."""

    import h5py

    with h5py.File(str(path), "r") as handle:
        if "matrix" not in handle:
            raise ValueError("HDF5 file has no /matrix group")
        matrix = cast(h5py.Group, handle["matrix"])
        features = cast(h5py.Group, matrix["features"])
        names_key = "name" if "name" in features else "gene_names"
        ids_key = "id" if "id" in features else "genes"
        feature_ids_ds = cast(h5py.Dataset, features[ids_key])
        feature_names_ds = cast(h5py.Dataset, features[names_key])
        feature_ids = _decode(np.asarray(feature_ids_ds[...]))
        feature_names = _decode(np.asarray(feature_names_ds[...]))
        if "feature_type" in features:
            feature_type_ds = cast(h5py.Dataset, features["feature_type"])
            feature_types = _decode(np.asarray(feature_type_ds[...]))
        else:
            feature_types = tuple("Gene Expression" for _ in feature_names)
        barcodes_ds = cast(h5py.Dataset, matrix["barcodes"])
        shape_ds = cast(h5py.Dataset, matrix["shape"])
        barcodes = _decode(np.asarray(barcodes_ds[...]))
        raw_shape = tuple(int(value) for value in np.asarray(shape_ds[...]))
    if len(raw_shape) != 2:
        raise ValueError(f"Unexpected matrix shape: {raw_shape}")
    if raw_shape != (len(feature_names), len(barcodes)):
        raise ValueError(
            "HDF5 shape does not match feature/barcode arrays: "
            f"{raw_shape}, {len(feature_names)}, {len(barcodes)}"
        )
    return XeniumFeatureCatalog(
        feature_ids=feature_ids,
        feature_names=feature_names,
        feature_types=feature_types,
        barcodes=barcodes,
        matrix_shape=raw_shape,
    )


def validate_cell_barcode_identity(cells: pd.DataFrame, catalog: XeniumFeatureCatalog) -> None:
    """Require exact cell-identifier identity before joining matrix and geometry."""

    if len(cells) != len(catalog.barcodes):
        raise ValueError(
            f"Cell table/HDF5 barcode count mismatch: {len(cells)} vs {len(catalog.barcodes)}"
        )
    cell_ids = set(cells["cell_id"].astype(str))
    matrix_barcodes = set(catalog.barcodes)
    if cell_ids != matrix_barcodes:
        missing_from_cells = sorted(matrix_barcodes.difference(cell_ids))[:10]
        missing_from_matrix = sorted(cell_ids.difference(matrix_barcodes))[:10]
        raise ValueError(
            "Cell table/HDF5 barcode identities differ; refusing to impute unmatched cells. "
            f"Examples missing from cells: {missing_from_cells}; "
            f"examples missing from matrix: {missing_from_matrix}"
        )


def _is_protein_feature_type(feature_type: str) -> bool:
    """Recognize Xenium/10x protein feature labels across software versions."""

    normalized = _normalise_token(feature_type)
    return "PROTEIN" in normalized or "ANTIBODY" in normalized


def canonical_feature_name(name: str, feature_type: str) -> str:
    """Create a stable, collision-resistant column name."""

    prefix = "protein" if _is_protein_feature_type(feature_type) else "rna"
    cleaned = "".join(character if character.isalnum() else "_" for character in name)
    cleaned = "_".join(part for part in cleaned.split("_") if part)
    return f"{prefix}__{cleaned}"


def _normalise_token(value: str) -> str:
    return "".join(character for character in value.upper() if character.isalnum())


def resolve_feature_indices(
    catalog: XeniumFeatureCatalog,
    *,
    markers: Mapping[str, Sequence[str]] = DEFAULT_MARKERS,
    targets: Mapping[str, Sequence[str]] = DEFAULT_TARGET_ALIASES,
    include_all_proteins: bool = True,
) -> tuple[list[int], dict[str, Any]]:
    """Resolve requested markers against exact normalized feature names."""

    normalized_to_indices: dict[str, list[int]] = {}
    for index, name in enumerate(catalog.feature_names):
        normalized_to_indices.setdefault(_normalise_token(name), []).append(index)

    requested_groups: dict[str, Sequence[str]] = {**markers, **targets}
    resolution: dict[str, Any] = {"groups": {}, "unresolved": {}}
    selected: set[int] = set()
    for group, aliases in requested_groups.items():
        resolved: list[dict[str, Any]] = []
        unresolved: list[str] = []
        group_indices: set[int] = set()
        for alias in aliases:
            candidates = normalized_to_indices.get(_normalise_token(alias), [])
            if not candidates:
                unresolved.append(alias)
                continue
            for index in candidates:
                if index in group_indices:
                    continue
                group_indices.add(index)
                selected.add(index)
                resolved.append(
                    {
                        "alias": alias,
                        "index": index,
                        "feature_id": catalog.feature_ids[index],
                        "feature_name": catalog.feature_names[index],
                        "feature_type": catalog.feature_types[index],
                        "column": canonical_feature_name(
                            catalog.feature_names[index], catalog.feature_types[index]
                        ),
                    }
                )
        resolution["groups"][group] = resolved
        if unresolved:
            resolution["unresolved"][group] = unresolved
    if include_all_proteins:
        for index, feature_type in enumerate(catalog.feature_types):
            if _is_protein_feature_type(feature_type):
                selected.add(index)
    ordered = sorted(selected)
    resolution["selected_indices"] = ordered
    resolution["selected_count"] = len(ordered)
    return ordered, resolution


def summarise_h5_features(
    path: Path,
    catalog: XeniumFeatureCatalog,
    *,
    chunk_nnz: int = 5_000_000,
) -> pd.DataFrame:
    """Compute totals and positive-cell fractions for all features by streaming NNZ arrays."""

    import h5py

    feature_count, cell_count = catalog.matrix_shape
    sums = np.zeros(feature_count, dtype=np.float64)
    positive_cells = np.zeros(feature_count, dtype=np.int64)
    with h5py.File(str(path), "r") as handle:
        matrix = cast(h5py.Group, handle["matrix"])
        indices_ds = cast(h5py.Dataset, matrix["indices"])
        data_ds = cast(h5py.Dataset, matrix["data"])
        nnz = int(data_ds.shape[0])
        for start in range(0, nnz, chunk_nnz):
            stop = min(start + chunk_nnz, nnz)
            indices = np.asarray(indices_ds[start:stop], dtype=np.int64)
            values = np.asarray(data_ds[start:stop], dtype=np.float64)
            sums += np.bincount(indices, weights=values, minlength=feature_count)
            positive_cells += np.bincount(indices, minlength=feature_count)
    return pd.DataFrame(
        {
            "feature_index": np.arange(feature_count, dtype=np.int64),
            "feature_id": catalog.feature_ids,
            "feature_name": catalog.feature_names,
            "feature_type": catalog.feature_types,
            "total_signal": sums,
            "mean_signal_per_cell": sums / max(cell_count, 1),
            "positive_cells": positive_cells,
            "positive_fraction": positive_cells / max(cell_count, 1),
        }
    )


def extract_selected_h5_features(
    path: Path,
    catalog: XeniumFeatureCatalog,
    selected_indices: Sequence[int],
    *,
    chunk_cells: int = 20_000,
) -> pd.DataFrame:
    """Read selected features from CSC HDF5 without materializing the full matrix."""

    import h5py

    feature_count, cell_count = catalog.matrix_shape
    selected = np.asarray(sorted(set(int(value) for value in selected_indices)), dtype=np.int64)
    if selected.size == 0:
        raise ValueError("No selected features")
    if selected.min() < 0 or selected.max() >= feature_count:
        raise IndexError("Selected feature index is outside matrix bounds")
    lookup = np.full(feature_count, -1, dtype=np.int64)
    lookup[selected] = np.arange(selected.size, dtype=np.int64)
    output = np.zeros((cell_count, selected.size), dtype=np.float32)
    with h5py.File(str(path), "r") as handle:
        matrix = cast(h5py.Group, handle["matrix"])
        indptr_ds = cast(h5py.Dataset, matrix["indptr"])
        indices_ds = cast(h5py.Dataset, matrix["indices"])
        data_ds = cast(h5py.Dataset, matrix["data"])
        for cell_start in range(0, cell_count, chunk_cells):
            cell_stop = min(cell_start + chunk_cells, cell_count)
            indptr = np.asarray(indptr_ds[cell_start : cell_stop + 1], dtype=np.int64)
            nnz_start = int(indptr[0])
            nnz_stop = int(indptr[-1])
            local_indptr = indptr - nnz_start
            indices = np.asarray(indices_ds[nnz_start:nnz_stop], dtype=np.int64)
            values = np.asarray(data_ds[nnz_start:nnz_stop], dtype=np.float32)
            local_cells = np.repeat(
                np.arange(cell_stop - cell_start, dtype=np.int64), np.diff(local_indptr)
            )
            mapped = lookup[indices]
            keep = mapped >= 0
            block = output[cell_start:cell_stop]
            np.add.at(block, (local_cells[keep], mapped[keep]), values[keep])
    columns = [
        canonical_feature_name(catalog.feature_names[index], catalog.feature_types[index])
        for index in selected
    ]
    if len(columns) != len(set(columns)):
        raise ValueError("Canonical selected-feature columns are not unique")
    table = pd.DataFrame(output, columns=columns)
    table.insert(0, "cell_id", list(catalog.barcodes))
    return table


def robust_scale(
    values: ArrayLike,
    *,
    lower_quantile: float = 0.05,
    upper_quantile: float = 0.99,
) -> FloatArray:
    """Map non-negative signal to [0, 1] with log1p and robust quantiles."""

    array = np.maximum(np.asarray(values, dtype=np.float64), 0.0)
    transformed = np.log1p(array)
    low = float(np.quantile(transformed, lower_quantile))
    high = float(np.quantile(transformed, upper_quantile))
    if not math.isfinite(low) or not math.isfinite(high) or high <= low:
        return np.zeros_like(transformed)
    return np.asarray(np.clip((transformed - low) / (high - low), 0.0, 1.0), dtype=np.float64)


def otsu_threshold(values: ArrayLike, *, bins: int = 256) -> float:
    """Compute Otsu's one-dimensional threshold without an image dependency."""

    array = np.asarray(values, dtype=np.float64)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return 0.5
    low, high = float(array.min()), float(array.max())
    if high <= low:
        return high
    hist, edges = np.histogram(array, bins=bins, range=(low, high))
    hist = hist.astype(np.float64)
    probabilities = hist / max(hist.sum(), 1.0)
    centres = (edges[:-1] + edges[1:]) / 2.0
    cumulative_weight = np.cumsum(probabilities)
    cumulative_mean = np.cumsum(probabilities * centres)
    total_mean = cumulative_mean[-1]
    denominator = cumulative_weight * (1.0 - cumulative_weight)
    between = np.zeros_like(denominator)
    valid = denominator > 0
    between[valid] = (
        total_mean * cumulative_weight[valid] - cumulative_mean[valid]
    ) ** 2 / denominator[valid]
    return float(centres[int(np.argmax(between))])


def _columns_for_resolution(
    resolution: Mapping[str, Any], group: str, available_columns: Iterable[str]
) -> list[str]:
    available = set(available_columns)
    return [
        str(item["column"])
        for item in resolution["groups"].get(group, [])
        if str(item["column"]) in available
    ]


def group_score(table: pd.DataFrame, columns: Sequence[str]) -> FloatArray:
    """Compute a conservative group score as the mean of per-feature robust scales."""

    if not columns:
        return np.zeros(len(table), dtype=np.float64)
    scaled = np.column_stack([robust_scale(table[column].to_numpy()) for column in columns])
    return np.asarray(np.mean(scaled, axis=1), dtype=np.float64)


def threshold_aligned_scale(values: ArrayLike, threshold: float) -> FloatArray:
    """Scale values so the declared positivity threshold maps exactly to 0.5."""

    array = np.clip(np.asarray(values, dtype=np.float64), 0.0, 1.0)
    threshold = float(np.clip(threshold, 1.0e-6, 1.0 - 1.0e-6))
    below = 0.5 * array / threshold
    above = 0.5 + 0.5 * (array - threshold) / (1.0 - threshold)
    return np.where(array <= threshold, below, above)


def local_mean_scores(
    coordinates: NDArray[np.float64],
    values: FloatArray,
    *,
    neighbours: int = 24,
    maximum_distance_um: float = 100.0,
    chunk_size: int = 50_000,
) -> FloatArray:
    """Average a cell score over nearby cells using a bounded nearest-neighbour query."""

    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("coordinates must have shape (n_cells, 2)")
    if len(values) != len(coordinates):
        raise ValueError("values and coordinates differ in length")
    count = len(coordinates)
    if count == 0:
        return np.array([], dtype=np.float64)
    k = min(max(neighbours, 1), count)
    tree = cKDTree(coordinates)
    output = np.empty(count, dtype=np.float64)
    for start in range(0, count, chunk_size):
        stop = min(start + chunk_size, count)
        distances, indices = tree.query(coordinates[start:stop], k=k, workers=-1)
        distances_array = np.asarray(distances, dtype=np.float64)
        indices_array = np.asarray(indices, dtype=np.int64)
        if distances_array.ndim == 1:
            distances_array = distances_array[:, None]
            indices_array = indices_array[:, None]
        valid = distances_array <= maximum_distance_um
        weights = valid.astype(np.float64)
        denominator = np.maximum(weights.sum(axis=1), 1.0)
        output[start:stop] = (values[indices_array] * weights).sum(axis=1) / denominator
    return output


def nearest_vessel_geometry(
    coordinates: NDArray[np.float64],
    vessel_positive: NDArray[np.bool_],
    *,
    density_radius_um: float = 100.0,
    density_neighbours: int = 32,
    chunk_size: int = 50_000,
) -> tuple[FloatArray, FloatArray]:
    """Compute nearest endothelial-cell distance and a bounded local vessel-density proxy."""

    vessel_coordinates = coordinates[vessel_positive]
    if len(vessel_coordinates) == 0:
        return (
            np.full(len(coordinates), np.nan, dtype=np.float64),
            np.zeros(len(coordinates), dtype=np.float64),
        )
    tree = cKDTree(vessel_coordinates)
    distances = np.empty(len(coordinates), dtype=np.float64)
    local_density = np.empty(len(coordinates), dtype=np.float64)
    k = min(max(density_neighbours, 1), len(vessel_coordinates))
    for start in range(0, len(coordinates), chunk_size):
        stop = min(start + chunk_size, len(coordinates))
        nearest, _ = tree.query(coordinates[start:stop], k=1, workers=-1)
        distances[start:stop] = np.asarray(nearest, dtype=np.float64)
        neighbourhood, _ = tree.query(coordinates[start:stop], k=k, workers=-1)
        neighbourhood_array = np.asarray(neighbourhood, dtype=np.float64)
        if neighbourhood_array.ndim == 1:
            neighbourhood_array = neighbourhood_array[:, None]
        local_density[start:stop] = np.mean(neighbourhood_array <= density_radius_um, axis=1)
    return distances, local_density


def score_cells(
    cells: pd.DataFrame,
    expression: pd.DataFrame,
    resolution: Mapping[str, Any],
    *,
    local_neighbours: int = 24,
    local_radius_um: float = 100.0,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Create transparent marker scores, vessel geometry and target signals."""

    merged = cells.merge(expression, on="cell_id", how="left", validate="one_to_one")
    feature_columns = [column for column in expression.columns if column != "cell_id"]
    merged[feature_columns] = merged[feature_columns].fillna(0.0)
    coordinates = merged[["x_um", "y_um"]].to_numpy(dtype=np.float64)

    score_groups = (
        "endothelial",
        "pericyte_smooth_muscle",
        "caf",
        "ecm",
        "epithelial_malignant",
        "immune",
    )
    group_columns: dict[str, list[str]] = {}
    for group in score_groups:
        columns = _columns_for_resolution(resolution, group, merged.columns)
        group_columns[group] = columns
        merged[f"{group}_score"] = group_score(merged, columns)

    vessel_raw = np.maximum(
        merged["endothelial_score"].to_numpy(dtype=np.float64),
        0.5 * merged["pericyte_smooth_muscle_score"].to_numpy(dtype=np.float64),
    )
    vessel_threshold = otsu_threshold(vessel_raw)
    # Guard against pathological Otsu splits by requiring a nontrivial high-confidence tail.
    vessel_threshold = float(np.clip(vessel_threshold, 0.25, 0.85))
    vessel_signal = threshold_aligned_scale(vessel_raw, vessel_threshold)
    vessel_positive = vessel_signal >= 0.5
    if vessel_positive.mean() < 0.001:
        fallback = float(np.quantile(vessel_raw, 0.99))
        vessel_signal = threshold_aligned_scale(vessel_raw, max(fallback, 1.0e-6))
        vessel_positive = vessel_signal >= 0.5
        vessel_threshold = fallback

    local_caf = local_mean_scores(
        coordinates,
        merged["caf_score"].to_numpy(dtype=np.float64),
        neighbours=local_neighbours,
        maximum_distance_um=local_radius_um,
    )
    local_ecm = local_mean_scores(
        coordinates,
        merged["ecm_score"].to_numpy(dtype=np.float64),
        neighbours=local_neighbours,
        maximum_distance_um=local_radius_um,
    )
    vessel_distance, local_vessel_density = nearest_vessel_geometry(coordinates, vessel_positive)
    merged["vessel_signal"] = vessel_signal
    merged["vessel_positive"] = vessel_positive
    merged["distance_to_vessel_um"] = vessel_distance
    merged["local_vessel_density"] = local_vessel_density
    merged["local_caf_score"] = np.clip(local_caf, 0.0, 1.0)
    merged["local_ecm_score"] = np.clip(local_ecm, 0.0, 1.0)

    malignant_raw = merged["epithelial_malignant_score"].to_numpy(dtype=np.float64)
    malignant_threshold = float(np.clip(otsu_threshold(malignant_raw), 0.20, 0.85))
    malignant = malignant_raw >= malignant_threshold
    merged["cell_is_malignant_proxy"] = malignant
    if np.any(malignant):
        malignant_tree = cKDTree(coordinates[malignant])
        distance_to_malignant, _ = malignant_tree.query(coordinates, k=1, workers=-1)
        merged["distance_to_malignant_proxy_um"] = distance_to_malignant
        merged["in_molecular_tumour_neighbourhood"] = distance_to_malignant <= 150.0
    else:
        merged["distance_to_malignant_proxy_um"] = np.nan
        merged["in_molecular_tumour_neighbourhood"] = False

    target_diagnostics: dict[str, Any] = {}
    for target in DEFAULT_TARGET_ALIASES:
        columns = _columns_for_resolution(resolution, target, merged.columns)
        if not columns:
            continue
        # Prefer direct protein intensity; use RNA only as a separately identified fallback.
        protein_columns = [column for column in columns if column.startswith("protein__")]
        active_columns = protein_columns or columns
        target_raw = group_score(merged, active_columns)
        threshold = float(np.clip(otsu_threshold(target_raw), 0.10, 0.90))
        signal = threshold_aligned_scale(target_raw, threshold)
        safe_name = target.replace("-", "_")
        merged[f"target__{safe_name}__signal"] = signal
        merged[f"target__{safe_name}__positive"] = signal >= 0.5
        target_diagnostics[target] = {
            "columns": active_columns,
            "measurement": "protein_intensity" if protein_columns else "rna_proxy",
            "raw_threshold": threshold,
            "positive_cells": int(np.sum(signal >= 0.5)),
            "positive_fraction": float(np.mean(signal >= 0.5)),
        }

    diagnostics = {
        "group_columns": group_columns,
        "vessel_threshold": vessel_threshold,
        "vessel_positive_cells": int(vessel_positive.sum()),
        "vessel_positive_fraction": float(vessel_positive.mean()),
        "malignant_proxy_threshold": malignant_threshold,
        "malignant_proxy_cells": int(malignant.sum()),
        "targets": target_diagnostics,
        "warnings": [
            "Endothelial-cell presence is not a measurement of vessel perfusion.",
            "Xenium protein signal is scaled mean fluorescence intensity, "
            "not antigen molecules per cell.",
            "The molecular tumour neighbourhood is a fallback proxy until "
            "pathology alignment is verified.",
        ],
    }
    return merged, diagnostics


def load_affine_matrix(path: Path) -> NDArray[np.float64]:
    """Read the 3x3 affine transformation supplied by Xenium Explorer."""

    matrix = np.loadtxt(path, delimiter=",")
    if matrix.shape != (3, 3):
        raise ValueError(f"Expected a 3x3 alignment matrix, observed {matrix.shape}")
    if not np.allclose(matrix[2], np.array([0.0, 0.0, 1.0]), atol=1.0e-8):
        raise ValueError("Alignment matrix final row is not [0, 0, 1]")
    return np.asarray(matrix, dtype=np.float64)


def geojson_vertices(path: Path) -> NDArray[np.float64]:
    """Flatten Polygon and MultiPolygon vertices from a GeoJSON file."""

    document = json.loads(path.read_text(encoding="utf-8"))
    vertices: list[tuple[float, float]] = []
    for feature in document.get("features", []):
        geometry = feature.get("geometry", {})
        geometry_type = geometry.get("type")
        coordinates = geometry.get("coordinates", [])
        polygons = coordinates if geometry_type == "MultiPolygon" else [coordinates]
        if geometry_type not in {"Polygon", "MultiPolygon"}:
            continue
        for polygon in polygons:
            for ring in polygon:
                for coordinate in ring:
                    vertices.append((float(coordinate[0]), float(coordinate[1])))
    if not vertices:
        raise ValueError("GeoJSON contains no polygon vertices")
    return np.asarray(vertices, dtype=np.float64)


def _bbox_score(
    transformed: NDArray[np.float64],
    cell_bbox: tuple[float, float, float, float],
) -> tuple[float, float]:
    x_min, x_max, y_min, y_max = cell_bbox
    margin_x = max((x_max - x_min) * 0.10, 100.0)
    margin_y = max((y_max - y_min) * 0.10, 100.0)
    inside = (
        (transformed[:, 0] >= x_min - margin_x)
        & (transformed[:, 0] <= x_max + margin_x)
        & (transformed[:, 1] >= y_min - margin_y)
        & (transformed[:, 1] <= y_max + margin_y)
    )
    fraction_inside = float(np.mean(inside))
    tx_min, ty_min = transformed.min(axis=0)
    tx_max, ty_max = transformed.max(axis=0)
    cell_width = max(x_max - x_min, 1.0)
    cell_height = max(y_max - y_min, 1.0)
    width_ratio = max((tx_max - tx_min) / cell_width, 1.0e-9)
    height_ratio = max((ty_max - ty_min) / cell_height, 1.0e-9)
    shape_penalty = abs(math.log(width_ratio)) + abs(math.log(height_ratio))
    centre_distance = math.hypot(
        ((tx_min + tx_max) - (x_min + x_max)) / (2.0 * cell_width),
        ((ty_min + ty_max) - (y_min + y_max)) / (2.0 * cell_height),
    )
    score = fraction_inside - 0.15 * shape_penalty - 0.10 * centre_distance
    return score, fraction_inside


def infer_annotation_transform(
    vertices: NDArray[np.float64],
    cells: pd.DataFrame,
    affine: NDArray[np.float64],
    *,
    candidate_pixel_sizes_um: Sequence[float] = (1.0, 0.2125, 0.425, 0.5),
) -> tuple[AlignmentCandidate | None, list[AlignmentCandidate]]:
    """Evaluate transform direction and pixel scaling; abstain when ambiguous."""

    homogeneous = np.column_stack([vertices, np.ones(len(vertices), dtype=np.float64)])
    cell_bbox = (
        float(cells["x_um"].min()),
        float(cells["x_um"].max()),
        float(cells["y_um"].min()),
        float(cells["y_um"].max()),
    )
    transformations = {
        "identity": np.eye(3, dtype=np.float64),
        "affine": affine,
        "inverse_affine": np.linalg.inv(affine),
    }
    candidates: list[AlignmentCandidate] = []
    seen: set[tuple[float, ...]] = set()
    for transform_name, matrix in transformations.items():
        transformed_pixels = (matrix @ homogeneous.T).T[:, :2]
        for scale in candidate_pixel_sizes_um:
            signature = tuple(np.round((matrix * float(scale)).ravel(), 10).tolist())
            if signature in seen:
                continue
            seen.add(signature)
            transformed = transformed_pixels * np.array([scale, scale])
            score, fraction_inside = _bbox_score(transformed, cell_bbox)
            candidates.append(
                AlignmentCandidate(
                    name=f"{transform_name}_scale_{scale:g}",
                    matrix=np.asarray(matrix, dtype=np.float64),
                    scale_x=float(scale),
                    scale_y=float(scale),
                    score=score,
                    fraction_inside=fraction_inside,
                )
            )
    candidates.sort(key=lambda candidate: candidate.score, reverse=True)
    best = candidates[0]
    runner_up = candidates[1]
    if best.fraction_inside < 0.70 or best.score - runner_up.score < 0.05:
        return None, candidates
    return best, candidates


def transform_geojson(
    input_path: Path,
    output_path: Path,
    candidate: AlignmentCandidate,
) -> None:
    """Transform a GeoJSON annotation layer into inferred Xenium micron coordinates."""

    document = json.loads(input_path.read_text(encoding="utf-8"))

    def transform_coordinate(coordinate: Sequence[float]) -> list[float]:
        vector = np.array([float(coordinate[0]), float(coordinate[1]), 1.0])
        result = candidate.matrix @ vector
        return [float(result[0] * candidate.scale_x), float(result[1] * candidate.scale_y)]

    for feature in document.get("features", []):
        geometry = feature.get("geometry", {})
        geometry_type = geometry.get("type")
        coordinates = geometry.get("coordinates", [])
        if geometry_type == "Polygon":
            geometry["coordinates"] = [
                [transform_coordinate(coordinate) for coordinate in ring] for ring in coordinates
            ]
        elif geometry_type == "MultiPolygon":
            geometry["coordinates"] = [
                [[transform_coordinate(coordinate) for coordinate in ring] for ring in polygon]
                for polygon in coordinates
            ]
    document.setdefault("reach_gap", {})["inferred_transform"] = {
        "name": candidate.name,
        "matrix": candidate.matrix.tolist(),
        "scale_x": candidate.scale_x,
        "scale_y": candidate.scale_y,
        "score": candidate.score,
        "fraction_inside": candidate.fraction_inside,
    }
    output_path.write_text(json.dumps(document), encoding="utf-8")


def _annotation_names(document: Mapping[str, Any]) -> list[str]:
    names: list[str] = []
    for feature in document.get("features", []):
        properties = feature.get("properties", {})
        classification = properties.get("classification") or {}
        name = properties.get("name") or classification.get("name")
        names.append(str(name or "UNLABELLED"))
    return names


def assign_pathology_regions(cells: pd.DataFrame, geojson_path: Path) -> pd.DataFrame:
    """Assign polygon labels to cell centroids with explicit overlap handling."""

    from shapely.geometry import Point, shape
    from shapely.strtree import STRtree

    document = json.loads(geojson_path.read_text(encoding="utf-8"))
    geometries = [shape(feature["geometry"]) for feature in document.get("features", [])]
    names = _annotation_names(document)
    if not geometries:
        output = cells.copy()
        output["pathology_region"] = "UNANNOTATED"
        return output
    tree = STRtree(geometries)
    try:
        from shapely import points

        point_array = points(
            cells["x_um"].to_numpy(dtype=np.float64),
            cells["y_um"].to_numpy(dtype=np.float64),
        )
        pairs = tree.query(point_array, predicate="intersects")
        matched_by_cell: dict[int, set[str]] = {}
        if pairs.size:
            for cell_index, geometry_index in pairs.T:
                matched_by_cell.setdefault(int(cell_index), set()).add(names[int(geometry_index)])
        labels = [
            "|".join(sorted(matched_by_cell[index])) if index in matched_by_cell else "UNANNOTATED"
            for index in range(len(cells))
        ]
    except (ImportError, AttributeError):
        labels = []
        for x, y in cells[["x_um", "y_um"]].itertuples(index=False, name=None):
            point = Point(float(x), float(y))
            candidate_indices = tree.query(point, predicate="intersects")
            if len(candidate_indices) == 0:
                labels.append("UNANNOTATED")
            else:
                matched = sorted({names[int(index)] for index in candidate_indices})
                labels.append("|".join(matched))
    output = cells.copy()
    output["pathology_region"] = labels
    return output


def write_partitioned_table(
    table: pd.DataFrame,
    output_dir: Path,
    stem: str,
    *,
    rows_per_part: int = 100_000,
) -> list[str]:
    """Write bounded-size table parts, preferring Zstandard-compressed Parquet."""

    if rows_per_part < 1:
        raise ValueError("rows_per_part must be positive")
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []
    try:
        import pyarrow  # noqa: F401

        suffix = "parquet"

        def writer(frame: pd.DataFrame, path: Path) -> None:
            frame.to_parquet(path, index=False, compression="zstd")

    except ImportError:
        suffix = "csv.gz"

        def writer(frame: pd.DataFrame, path: Path) -> None:
            frame.to_csv(path, index=False, compression="gzip")

    for part, start in enumerate(range(0, len(table), rows_per_part)):
        stop = min(start + rows_per_part, len(table))
        path = output_dir / f"{stem}.part{part:04d}.{suffix}"
        writer(table.iloc[start:stop], path)
        paths.append(str(path))
    if not paths:
        path = output_dir / f"{stem}.part0000.{suffix}"
        writer(table, path)
        paths.append(str(path))
    return paths


def read_partitioned_table(paths: Sequence[Path]) -> pd.DataFrame:
    """Read table parts produced by :func:`write_partitioned_table`."""

    frames: list[pd.DataFrame] = []
    for path in sorted(paths):
        if path.suffix == ".parquet":
            frames.append(pd.read_parquet(path))
        elif path.name.endswith(".csv.gz"):
            frames.append(pd.read_csv(path, compression="gzip"))
        else:
            frames.append(pd.read_csv(path))
    if not frames:
        raise FileNotFoundError("No table parts were supplied")
    return pd.concat(frames, ignore_index=True)


def build_target_tables(
    scored: pd.DataFrame,
    output_dir: Path,
    diagnostics: Mapping[str, Any],
) -> dict[str, list[str]]:
    """Write generic reach-gap input tables for every resolved target."""

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, list[str]] = {}
    if "pathology_region" in scored:
        pathology_tumour = (
            scored["pathology_region"].astype(str).str.contains("Tumor", case=False, regex=False)
        )
        in_tumour_region = pathology_tumour.to_numpy(dtype=np.bool_)
        region_definition = "pathology_tumour_polygon"
    else:
        in_tumour_region = scored["in_molecular_tumour_neighbourhood"].to_numpy(dtype=np.bool_)
        region_definition = "molecular_tumour_neighbourhood_150um"

    for target, target_info in diagnostics.get("targets", {}).items():
        safe_name = target.replace("-", "_")
        signal_column = f"target__{safe_name}__signal"
        if signal_column not in scored:
            continue
        table = pd.DataFrame(
            {
                "cell_id": scored["cell_id"].astype(str),
                "x_um": scored["x_um"].astype(float),
                "y_um": scored["y_um"].astype(float),
                # Historical schema name: this means 'cell lies in tumour region', not malignancy.
                "is_tumour": in_tumour_region,
                "target_signal": scored[signal_column].astype(float),
                "vessel_signal": scored["vessel_signal"].astype(float),
                "ecm_score": scored["local_ecm_score"].astype(float),
                "caf_score": scored["local_caf_score"].astype(float),
                "cell_is_malignant_proxy": scored["cell_is_malignant_proxy"].astype(bool),
                "distance_to_vessel_um": scored["distance_to_vessel_um"].astype(float),
                "target_measurement": str(target_info["measurement"]),
                "tumour_region_definition": region_definition,
            }
        )
        outputs[target] = write_partitioned_table(table, output_dir, f"reach_gap_cells_{safe_name}")
    return outputs


def write_manifest(
    output_path: Path,
    *,
    source_files: Sequence[Path],
    zip_inventory: pd.DataFrame,
    extracted: Mapping[str, str],
    diagnostics: Mapping[str, Any],
    alignment: AlignmentCandidate | None,
) -> None:
    """Write an auditable preparation manifest with explicit unsupported claims."""

    manifest = {
        "schema_version": "1.0",
        "dataset": (
            "Xenium In Situ Gene and Protein Expression data for FFPE Human Renal Cell Carcinoma"
        ),
        "platform": "Xenium Onboard Analysis 4.0",
        "licence": "CC BY 4.0",
        "source_files": [
            {
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": (
                    sha256_file(path)
                    if path.stat().st_size < 1_000_000_000
                    else "NOT_COMPUTED_LARGE_FILE"
                ),
            }
            for path in source_files
            if path.exists()
        ],
        "zip_members": len(zip_inventory),
        "zip_uncompressed_bytes": int(zip_inventory["uncompressed_bytes"].sum()),
        "extracted": dict(extracted),
        "diagnostics": diagnostics,
        "annotation_alignment": None
        if alignment is None
        else {
            "name": alignment.name,
            "matrix": alignment.matrix.tolist(),
            "scale_x": alignment.scale_x,
            "scale_y": alignment.scale_y,
            "score": alignment.score,
            "fraction_inside": alignment.fraction_inside,
        },
        "calibration": {
            "antigen_molecules_per_cell": "NOT_COMPUTED",
            "antigen_nM_per_signal": None,
            "reason": (
                "Xenium protein values are scaled mean fluorescence intensity, "
                "not absolute surface density."
            ),
        },
        "perfusion": {
            "measured": False,
            "reason": "CD31/endothelial presence does not establish functional perfusion.",
        },
        "permitted_claims": [
            "Cell-level RNA and protein signals were prepared from the complete "
            "Xenium cell-feature matrix.",
            "Distances are measured to high-confidence endothelial-cell proxies "
            "in the same section.",
            "Target and barrier scores are relative, transparent marker-derived quantities.",
        ],
        "unsupported_claims": [
            "Absolute antibody penetration or receptor occupancy.",
            "Clinical efficacy prediction.",
            "Perfused-vessel identification.",
            "Absolute surface-antigen density.",
        ],
    }
    output_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def prepare_rcc_xenium(
    *,
    raw_dir: Path,
    output_dir: Path,
    verify_large_md5: bool = True,
    force: bool = False,
) -> dict[str, Any]:
    """Run the complete in-Drive preparation pipeline with resumable checkpoints."""

    output_dir.mkdir(parents=True, exist_ok=True)
    verification_path = output_dir / "download_verification.json"
    if verification_path.exists() and not force:
        verify_report = json.loads(verification_path.read_text(encoding="utf-8"))
        cached_sizes_valid = all(
            (raw_dir / name).exists()
            and (raw_dir / name).stat().st_size == int(specification["size"])
            for name, specification in EXPECTED_RCC_FILES.items()
        )
        if not (
            verify_report.get("all_present")
            and verify_report.get("all_verified")
            and cached_sizes_valid
        ):
            verify_report = verify_expected_files(raw_dir, verify_large_md5=verify_large_md5)
    else:
        verify_report = verify_expected_files(raw_dir, verify_large_md5=verify_large_md5)
    verification_path.write_text(json.dumps(verify_report, indent=2), encoding="utf-8")
    if not verify_report["all_present"]:
        missing = [name for name, entry in verify_report["files"].items() if not entry["present"]]
        raise FileNotFoundError(f"Missing downloaded files: {missing}")
    if not verify_report["all_verified"]:
        raise ValueError("At least one downloaded file failed size or MD5 verification")

    zip_path = raw_dir / "Xenium_V1_Human_Kidney_FFPE_Protein_updated_outs.zip"
    he_path = raw_dir / "Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_image.ome.tif"
    alignment_path = raw_dir / "Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_imagealignment.csv"
    annotation_path = raw_dir / "Xenium_V1_Human_Kidney_FFPE_Protein_updated_annotation.geojson"

    inventory_path = output_dir / "zip_inventory.csv"
    if inventory_path.exists() and not force:
        inventory = pd.read_csv(inventory_path)
    else:
        inventory = inspect_zip(zip_path)
        inventory.to_csv(inventory_path, index=False)

    extracted_manifest_path = output_dir / "extracted_members.json"
    extracted_dir = output_dir / "extracted"
    extracted: dict[str, str]
    if extracted_manifest_path.exists() and not force:
        cached = json.loads(extracted_manifest_path.read_text(encoding="utf-8"))
        cached_paths_valid = cached and all(
            Path(path).exists() and Path(path).stat().st_size > 0 for path in cached.values()
        )
        if cached_paths_valid:
            extracted = {str(key): str(value) for key, value in cached.items()}
        else:
            members = select_essential_members(inventory)
            extracted = extract_members(zip_path, members, extracted_dir)
    else:
        members = select_essential_members(inventory)
        extracted = extract_members(zip_path, members, extracted_dir)
    extracted_manifest_path.write_text(json.dumps(extracted, indent=2), encoding="utf-8")

    cells_path = find_extracted(extracted_dir, "cells.parquet") or find_extracted(
        extracted_dir, "cells.csv.gz"
    )
    h5_path = find_extracted(extracted_dir, "cell_feature_matrix.h5")
    if cells_path is None or h5_path is None:
        raise FileNotFoundError("Essential cells or HDF5 matrix file was not extracted")
    cells = read_cells(cells_path)
    catalog = read_10x_h5_catalog(h5_path)
    validate_cell_barcode_identity(cells, catalog)

    feature_summary_path = output_dir / "feature_summary.csv"
    resolution_path = output_dir / "marker_resolution.json"
    selected_parts_manifest = output_dir / "selected_expression_parts.json"
    selected_expression: pd.DataFrame
    if (
        feature_summary_path.exists()
        and resolution_path.exists()
        and selected_parts_manifest.exists()
        and not force
    ):
        selected_part_paths = [
            Path(value) for value in json.loads(selected_parts_manifest.read_text(encoding="utf-8"))
        ]
        if selected_part_paths and all(path.exists() for path in selected_part_paths):
            resolution = json.loads(resolution_path.read_text(encoding="utf-8"))
            selected_indices = [int(value) for value in resolution["selected_indices"]]
            selected_expression = read_partitioned_table(selected_part_paths)
        else:
            selected_indices, resolution = resolve_feature_indices(catalog)
            selected_expression = extract_selected_h5_features(h5_path, catalog, selected_indices)
            selected_part_paths = [
                Path(value)
                for value in write_partitioned_table(
                    selected_expression, output_dir / "tables", "selected_expression"
                )
            ]
    else:
        feature_summary = summarise_h5_features(h5_path, catalog)
        feature_summary.to_csv(feature_summary_path, index=False)
        selected_indices, resolution = resolve_feature_indices(catalog)
        resolution_path.write_text(json.dumps(resolution, indent=2), encoding="utf-8")
        selected_expression = extract_selected_h5_features(h5_path, catalog, selected_indices)
        selected_part_paths = [
            Path(value)
            for value in write_partitioned_table(
                selected_expression, output_dir / "tables", "selected_expression"
            )
        ]
    if not feature_summary_path.exists():
        summarise_h5_features(h5_path, catalog).to_csv(feature_summary_path, index=False)
    if not resolution_path.exists():
        resolution_path.write_text(json.dumps(resolution, indent=2), encoding="utf-8")
    selected_expression_parts = [str(path) for path in selected_part_paths]
    selected_parts_manifest.write_text(
        json.dumps(selected_expression_parts, indent=2), encoding="utf-8"
    )

    diagnostics_path = output_dir / "cell_scoring_diagnostics.json"
    scored_parts_manifest = output_dir / "scored_cell_parts.json"
    if diagnostics_path.exists() and scored_parts_manifest.exists() and not force:
        scored_part_paths = [
            Path(value) for value in json.loads(scored_parts_manifest.read_text(encoding="utf-8"))
        ]
        if scored_part_paths and all(path.exists() for path in scored_part_paths):
            diagnostics = json.loads(diagnostics_path.read_text(encoding="utf-8"))
            scored = read_partitioned_table(scored_part_paths)
        else:
            scored, diagnostics = score_cells(cells, selected_expression, resolution)
            scored_part_paths = []
    else:
        scored, diagnostics = score_cells(cells, selected_expression, resolution)
        scored_part_paths = []

    alignment_candidates_path = output_dir / "annotation_alignment_candidates.json"
    transformed_path = output_dir / "pathology_annotations_xenium.geojson"
    best_alignment: AlignmentCandidate | None = None
    if transformed_path.exists() and alignment_candidates_path.exists() and not force:
        candidate_records = json.loads(alignment_candidates_path.read_text(encoding="utf-8"))
        accepted = [record for record in candidate_records if record.get("accepted")]
        if accepted:
            record = accepted[0]
            best_alignment = AlignmentCandidate(
                name=str(record["name"]),
                matrix=np.asarray(record["matrix"], dtype=np.float64),
                scale_x=float(record["scale_x"]),
                scale_y=float(record["scale_y"]),
                score=float(record["score"]),
                fraction_inside=float(record["fraction_inside"]),
            )
            if "pathology_region" not in scored.columns:
                scored = assign_pathology_regions(scored, transformed_path)
    else:
        alignment = load_affine_matrix(alignment_path)
        vertices = geojson_vertices(annotation_path)
        best_alignment, candidates = infer_annotation_transform(vertices, cells, alignment)
        candidate_records = [
            {
                "name": candidate.name,
                "matrix": candidate.matrix.tolist(),
                "scale_x": candidate.scale_x,
                "scale_y": candidate.scale_y,
                "score": candidate.score,
                "fraction_inside": candidate.fraction_inside,
                "accepted": best_alignment is not None and candidate.name == best_alignment.name,
            }
            for candidate in candidates
        ]
        alignment_candidates_path.write_text(
            json.dumps(candidate_records, indent=2), encoding="utf-8"
        )
        if best_alignment is not None:
            transform_geojson(annotation_path, transformed_path, best_alignment)
            scored = assign_pathology_regions(scored, transformed_path)
        else:
            diagnostics.setdefault("warnings", []).append(
                "Pathology alignment was ambiguous; tumour regions use the molecular "
                "150 µm neighbourhood proxy."
            )

    # Re-write scored parts whenever pathology annotations were newly added or a
    # checkpoint was incomplete.
    if not scored_part_paths or (
        best_alignment is not None and "pathology_region" in scored.columns
    ):
        scored_part_paths = [
            Path(value)
            for value in write_partitioned_table(scored, output_dir / "tables", "cells_reach_gap")
        ]
    scored_parts = [str(path) for path in scored_part_paths]
    scored_parts_manifest.write_text(json.dumps(scored_parts, indent=2), encoding="utf-8")
    diagnostics_path.write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")

    target_outputs = build_target_tables(scored, output_dir / "targets", diagnostics)
    write_manifest(
        output_dir / "processing_manifest.json",
        source_files=[zip_path, he_path, alignment_path, annotation_path],
        zip_inventory=inventory,
        extracted=extracted,
        diagnostics=diagnostics,
        alignment=best_alignment,
    )
    result = {
        "status": "PREPARED_WITH_ABSOLUTE_INDEX_NOT_COMPUTED",
        "cells": len(cells),
        "features": catalog.matrix_shape[0],
        "selected_features": len(selected_indices),
        "selected_expression_parts": selected_expression_parts,
        "scored_cell_parts": scored_parts,
        "targets": target_outputs,
        "pathology_alignment": None if best_alignment is None else best_alignment.name,
        "output_dir": str(output_dir),
        "absolute_index": {
            "status": "NOT_COMPUTED",
            "reason": (
                "No independent conversion from Xenium protein intensity to "
                "surface antigen density is available."
            ),
        },
    }
    (output_dir / "run_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result
