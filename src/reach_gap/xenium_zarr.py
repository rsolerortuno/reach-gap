"""Native Xenium ``cell_features`` Zarr decoding for reach-gap v0.7.1.

The provider schema stores a feature-by-cell sparse matrix under ``cell_features``.
Feature names and dimensions are group attributes, while cell identifiers are packed
integer pairs.  This module keeps schema detection explicit and rejects ambiguous
basename matches—the ambiguity that caused the v0.7.1 Colab failure.
"""

from __future__ import annotations

import contextlib
import json
import re
from collections.abc import Mapping, Sequence
from typing import Any, TypeAlias

import numpy as np
from numpy.typing import NDArray
from scipy import sparse

FloatArray: TypeAlias = NDArray[np.float64]
IntArray: TypeAlias = NDArray[np.int64]


def decode_string_array(values: NDArray[Any]) -> list[str]:
    """Decode an arbitrary one-dimensional string-like array."""

    decoded: list[str] = []
    for value in np.asarray(values).ravel():
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        else:
            decoded.append(str(value))
    return decoded


def decode_zarr_string_attribute(value: Any) -> list[str]:
    """Decode list-like Zarr attributes, including JSON-encoded variants."""

    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.startswith("["):
            with contextlib.suppress(json.JSONDecodeError):
                value = json.loads(stripped)
        else:
            return [value]
    return decode_string_array(np.asarray(value, dtype=object))


def xenium_packed_cell_ids_to_strings(packed_cell_ids: NDArray[Any]) -> list[str]:
    """Convert 10x uint32 ``[prefix, dataset_suffix]`` IDs to Xenium IDs.

    Each hexadecimal nibble of the prefix is shifted to the alphabet ``a``-``p``.
    The second integer is retained after a hyphen, matching provider cell-group CSVs.
    """

    values = np.asarray(packed_cell_ids)
    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError(
            f"Expected cell_features/cell_id with shape (cells, 2); observed {values.shape}"
        )

    shifted_hex_alphabet = "abcdefghijklmnop"
    cell_ids: list[str] = []
    for prefix, suffix in values[:, :2]:
        prefix_int = int(prefix)
        suffix_int = int(suffix)
        if prefix_int < 0 or prefix_int > 0xFFFFFFFF:
            raise ValueError(f"Packed cell prefix outside uint32 range: {prefix_int}")
        if suffix_int < 0:
            raise ValueError(f"Packed cell suffix must be non-negative: {suffix_int}")
        hexadecimal = f"{prefix_int:08x}"
        shifted = "".join(shifted_hex_alphabet[int(character, 16)] for character in hexadecimal)
        cell_ids.append(f"{shifted}-{suffix_int}")
    return cell_ids


def resolve_array_name(array_names: Sequence[str], candidates: Sequence[str]) -> str:
    """Resolve exactly one array path without ambiguous basename selection."""

    normalized_names = [name.strip("/") for name in array_names]
    normalized_candidates = [candidate.strip("/") for candidate in candidates]

    for candidate in normalized_candidates:
        if candidate in normalized_names:
            return candidate

    suffix_matches = sorted(
        {
            name
            for name in normalized_names
            if any(name.endswith("/" + candidate) for candidate in normalized_candidates)
        }
    )
    if len(suffix_matches) != 1:
        raise KeyError(
            "Could not resolve one Zarr array for candidates "
            f"{list(candidates)}. Matches: {suffix_matches}. "
            f"Available arrays: {sorted(normalized_names)}"
        )
    return suffix_matches[0]


def build_native_cell_features_matrix(
    *,
    attributes: Mapping[str, Any],
    arrays: Mapping[str, NDArray[Any]],
) -> tuple[sparse.spmatrix, list[str], list[str], str]:
    """Build a native Xenium sparse matrix from materialised components.

    This pure function is the regression-test seam used without requiring Zarr.
    CSR and CSC encodings are both accepted; the resulting orientation is always
    ``features x cells``.
    """

    required_attributes = {"number_cells", "number_features", "feature_keys"}
    missing_attributes = sorted(required_attributes - set(attributes))
    if missing_attributes:
        raise KeyError(
            "The cell_features group is missing required attributes: "
            f"{missing_attributes}. Available attributes: {sorted(attributes)}"
        )

    number_cells = int(attributes["number_cells"])
    number_features = int(attributes["number_features"])
    if number_cells <= 0 or number_features <= 0:
        raise ValueError("number_cells and number_features must be positive")

    feature_names = decode_zarr_string_attribute(attributes["feature_keys"])
    if len(feature_names) != number_features:
        raise ValueError(
            "feature_keys length does not equal number_features: "
            f"{len(feature_names)} != {number_features}"
        )

    array_names = tuple(sorted(name.strip("/") for name in arrays))
    cell_id_name = resolve_array_name(array_names, ("cell_id",))
    cell_ids = xenium_packed_cell_ids_to_strings(np.asarray(arrays[cell_id_name]))
    if len(cell_ids) != number_cells:
        raise ValueError(
            f"cell_id length does not equal number_cells: {len(cell_ids)} != {number_cells}"
        )
    if len(set(cell_ids)) != len(cell_ids):
        raise ValueError("Converted Xenium cell IDs are not unique")

    if {"data", "indices", "indptr"}.issubset(array_names):
        data = np.asarray(arrays["data"])
        indices: IntArray = np.asarray(arrays["indices"], dtype=np.int64)
        indptr: IntArray = np.asarray(arrays["indptr"], dtype=np.int64)
        matrix: sparse.spmatrix = sparse.csr_matrix(
            (data, indices, indptr), shape=(number_features, number_cells)
        )
        encoding = "CSR_FEATURES_BY_CELLS"
    elif {"csc/data", "csc/indices", "csc/indptr"}.issubset(array_names):
        data = np.asarray(arrays["csc/data"])
        indices = np.asarray(arrays["csc/indices"], dtype=np.int64)
        indptr = np.asarray(arrays["csc/indptr"], dtype=np.int64)
        matrix = sparse.csc_matrix((data, indices, indptr), shape=(number_features, number_cells))
        encoding = "CSC_FEATURES_BY_CELLS"
    else:
        raise KeyError(
            "Neither the native CSR nor CSC sparse representation is complete "
            f"under cell_features. Available arrays: {list(array_names)}"
        )

    if matrix.shape != (number_features, number_cells):  # pragma: no cover - SciPy enforces shape
        raise ValueError(
            f"Unexpected sparse matrix shape {matrix.shape}; expected "
            f"({number_features}, {number_cells})"
        )
    matrix_nnz = matrix.getnnz()
    if matrix_nnz != len(data):  # pragma: no cover - constructed directly from data
        raise ValueError(f"Sparse matrix nnz mismatch: {matrix_nnz} != {len(data)}")

    return matrix, feature_names, cell_ids, encoding


def find_feature_index(feature_names: Sequence[str], aliases: Sequence[str]) -> int:
    """Find exactly one normalized feature matching the supplied aliases."""

    normalized_aliases = {re.sub(r"[^A-Z0-9]", "", alias.upper()) for alias in aliases}
    normalized_names = [
        re.sub(r"[^A-Z0-9]", "", feature_name.upper()) for feature_name in feature_names
    ]
    matches = [index for index, name in enumerate(normalized_names) if name in normalized_aliases]
    if len(matches) != 1:
        matching_names = [feature_names[index] for index in matches]
        raise ValueError(
            f"Expected exactly one feature for aliases {list(aliases)}, "
            f"found indices {matches} and names {matching_names}"
        )
    return matches[0]


def find_erbb2_feature_index(feature_names: Sequence[str]) -> int:
    """Find the unique ERBB2/HER2 row."""

    return find_feature_index(feature_names, ("ERBB2", "HER2"))


def extract_feature_vector(
    matrix: sparse.spmatrix,
    feature_index: int,
    feature_count: int,
    barcode_count: int,
) -> FloatArray:
    """Extract one feature while accepting either explicit matrix orientation."""

    if feature_index < 0 or feature_index >= feature_count:
        raise IndexError(feature_index)
    if matrix.shape == (feature_count, barcode_count):
        values = matrix.getrow(feature_index).toarray().ravel()
    elif matrix.shape == (barcode_count, feature_count):
        values = matrix.getcol(feature_index).toarray().ravel()
    else:
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match "
            f"{feature_count} features and {barcode_count} barcodes"
        )
    return np.asarray(values, dtype=np.float64)
