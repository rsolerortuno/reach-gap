"""Zarr file I/O adapter for :mod:`reach_gap.xenium_zarr`."""

from __future__ import annotations

import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol, cast

import numpy as np
from numpy.typing import NDArray
from scipy import sparse

from .xenium_zarr import build_native_cell_features_matrix


class ArrayLike(Protocol):
    def __getitem__(self, key: object) -> Any: ...


@dataclass(frozen=True)
class XeniumSparseMatrix:
    matrix: sparse.spmatrix
    feature_names: tuple[str, ...]
    cell_ids: tuple[str, ...]
    array_names: tuple[str, ...]
    root_prefix: str
    encoding: str


def detect_zarr_root_prefix(path: Path) -> str:
    """Return the shallowest root containing a Zarr v2 group marker."""

    with zipfile.ZipFile(path, mode="r") as archive:
        candidates: list[str] = []
        for name in archive.namelist():
            normalized = name.rstrip("/")
            if normalized.endswith("/.zgroup"):
                candidates.append(normalized[: -len("/.zgroup")])
            elif normalized == ".zgroup":
                candidates.append("")
        if not candidates:
            return ""
        return min(candidates, key=lambda value: (value.count("/"), len(value)))


def _list_zarr_arrays(group: Any) -> list[str]:
    names: list[str] = []

    def visitor(name: str, obj: Any) -> None:
        if hasattr(obj, "shape") and hasattr(obj, "dtype"):
            names.append(name.strip("/"))

    group.visititems(visitor)
    return sorted(names)


def read_xenium_sparse_matrix(path: Path) -> XeniumSparseMatrix:
    """Read a zipped native Xenium matrix using an installed Zarr implementation."""

    try:
        import zarr  # type: ignore[import-not-found]
        from zarr.storage import ZipStore  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError(
            "Reading a real Xenium Zarr requires the project's pinned zarr dependency"
        ) from exc

    root_prefix = detect_zarr_root_prefix(path)
    store = ZipStore(str(path), mode="r")
    try:
        root = zarr.open_group(store=store, mode="r", path=root_prefix or None)
        arrays = _list_zarr_arrays(root)
        try:
            cell_features = root["cell_features"]
        except KeyError as exc:
            raise KeyError(
                f"Expected the native Xenium 'cell_features' group. Available arrays: {arrays}"
            ) from exc

        attributes = dict(cell_features.attrs)
        local_names = _list_zarr_arrays(cell_features)
        materialized: dict[str, NDArray[Any]] = {
            name: np.asarray(cast(ArrayLike, cell_features[name])[:]) for name in local_names
        }
        matrix, feature_names, cell_ids, encoding = build_native_cell_features_matrix(
            attributes=attributes,
            arrays=materialized,
        )
        return XeniumSparseMatrix(
            matrix=matrix,
            feature_names=tuple(feature_names),
            cell_ids=tuple(cell_ids),
            array_names=tuple([*arrays, f"__matrix_encoding__/{encoding}"]),
            root_prefix=root_prefix,
            encoding=encoding,
        )
    finally:
        store.close()
