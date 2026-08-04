from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from reach_gap.xenium_zarr import (
    build_native_cell_features_matrix,
    decode_string_array,
    decode_zarr_string_attribute,
    extract_feature_vector,
    find_erbb2_feature_index,
    find_feature_index,
    resolve_array_name,
    xenium_packed_cell_ids_to_strings,
)


def test_packed_cell_ids_match_xenium_encoding() -> None:
    values = np.asarray([[0x01234567, 1], [0x89ABCDEF, 7]], dtype=np.uint32)
    assert xenium_packed_cell_ids_to_strings(values) == ["abcdefgh-1", "ijklmnop-7"]


def test_packed_cell_ids_reject_wrong_shape() -> None:
    with pytest.raises(ValueError, match="shape"):
        xenium_packed_cell_ids_to_strings(np.asarray([1, 2], dtype=np.uint32))


def test_packed_ids_reject_invalid_values() -> None:
    with pytest.raises(ValueError, match="uint32"):
        xenium_packed_cell_ids_to_strings(np.asarray([[-1, 1]], dtype=np.int64))
    with pytest.raises(ValueError, match="non-negative"):
        xenium_packed_cell_ids_to_strings(np.asarray([[1, -1]], dtype=np.int64))


def test_string_decoding_variants() -> None:
    assert decode_string_array(np.asarray([b"ERBB2", "ACTB"], dtype=object)) == ["ERBB2", "ACTB"]
    assert decode_zarr_string_attribute(b'["ERBB2", "ACTB"]') == ["ERBB2", "ACTB"]
    assert decode_zarr_string_attribute("ERBB2") == ["ERBB2"]
    assert decode_zarr_string_attribute("[not-json") == ["[not-json"]


def test_array_resolution_rejects_ambiguous_legacy_basename() -> None:
    names = ["cell_features/data", "cell_features/csc/data"]
    with pytest.raises(KeyError, match="Matches"):
        resolve_array_name(names, ("matrix/data", "data"))


def test_array_resolution_supports_unique_suffix_and_exact_path() -> None:
    assert resolve_array_name(["cell_features/cell_id"], ("cell_id",)) == "cell_features/cell_id"
    assert resolve_array_name(["cell_id"], ("cell_id",)) == "cell_id"


def _csr_fixture() -> tuple[dict[str, object], dict[str, np.ndarray]]:
    attributes: dict[str, object] = {
        "number_cells": 3,
        "number_features": 2,
        "feature_keys": ["ACTB", "ERBB2"],
    }
    arrays = {
        "cell_id": np.asarray([[0x01234567, 1], [0x89ABCDEF, 1], [0x11111111, 1]], dtype=np.uint32),
        "data": np.asarray([4, 5, 6], dtype=np.int32),
        "indices": np.asarray([0, 1, 2], dtype=np.int32),
        "indptr": np.asarray([0, 1, 3], dtype=np.int32),
    }
    return attributes, arrays


def test_native_csr_feature_extraction() -> None:
    attributes, arrays = _csr_fixture()
    matrix, features, cell_ids, encoding = build_native_cell_features_matrix(
        attributes=attributes, arrays=arrays
    )
    assert encoding == "CSR_FEATURES_BY_CELLS"
    assert features == ["ACTB", "ERBB2"]
    assert len(cell_ids) == 3
    index = find_erbb2_feature_index(features)
    np.testing.assert_array_equal(extract_feature_vector(matrix, index, 2, 3), [0.0, 5.0, 6.0])


def test_native_csc_fallback() -> None:
    attributes = {
        "number_cells": 3,
        "number_features": 2,
        "feature_keys": '["ACTB", "ERBB2"]',
    }
    arrays = {
        "cell_id": np.asarray([[1, 1], [2, 1], [3, 1]], dtype=np.uint32),
        "csc/data": np.asarray([4, 5, 6], dtype=np.int32),
        "csc/indices": np.asarray([0, 1, 1], dtype=np.int32),
        "csc/indptr": np.asarray([0, 1, 2, 3], dtype=np.int32),
    }
    matrix, features, _, encoding = build_native_cell_features_matrix(
        attributes=attributes, arrays=arrays
    )
    assert encoding == "CSC_FEATURES_BY_CELLS"
    np.testing.assert_array_equal(extract_feature_vector(matrix, 1, 2, 3), [0.0, 5.0, 6.0])
    assert features == ["ACTB", "ERBB2"]


def test_native_matrix_rejects_malformed_metadata() -> None:
    cell_ids = np.asarray([[1, 1]], dtype=np.uint32)
    with pytest.raises(KeyError, match="required attributes"):
        build_native_cell_features_matrix(attributes={}, arrays={"cell_id": cell_ids})
    with pytest.raises(ValueError, match="positive"):
        build_native_cell_features_matrix(
            attributes={"number_cells": 0, "number_features": 1, "feature_keys": ["ERBB2"]},
            arrays={"cell_id": np.empty((0, 2), dtype=np.uint32)},
        )
    with pytest.raises(ValueError, match="feature_keys length"):
        build_native_cell_features_matrix(
            attributes={"number_cells": 1, "number_features": 2, "feature_keys": ["ERBB2"]},
            arrays={"cell_id": cell_ids},
        )
    with pytest.raises(ValueError, match="cell_id length"):
        build_native_cell_features_matrix(
            attributes={"number_cells": 2, "number_features": 1, "feature_keys": ["ERBB2"]},
            arrays={"cell_id": cell_ids},
        )
    with pytest.raises(ValueError, match="not unique"):
        build_native_cell_features_matrix(
            attributes={"number_cells": 2, "number_features": 1, "feature_keys": ["ERBB2"]},
            arrays={"cell_id": np.asarray([[1, 1], [1, 1]], dtype=np.uint32)},
        )


def test_native_matrix_rejects_incomplete_sparse_arrays() -> None:
    with pytest.raises(KeyError, match="Neither"):
        build_native_cell_features_matrix(
            attributes={"number_cells": 1, "number_features": 1, "feature_keys": ["ERBB2"]},
            arrays={"cell_id": np.asarray([[1, 1]], dtype=np.uint32), "data": np.asarray([1])},
        )


def test_feature_matching_rejects_zero_or_multiple_hits() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        find_feature_index(["ACTB"], ["ERBB2", "HER2"])
    with pytest.raises(ValueError, match="exactly one"):
        find_feature_index(["ERBB2", "HER-2"], ["ERBB2", "HER2"])


def test_feature_vector_validation_and_transposed_orientation() -> None:
    transposed = sparse.csr_matrix(np.asarray([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]))
    np.testing.assert_array_equal(extract_feature_vector(transposed, 1, 2, 3), [2.0, 4.0, 6.0])
    with pytest.raises(IndexError):
        extract_feature_vector(transposed, 2, 2, 3)
    with pytest.raises(ValueError, match="does not match"):
        extract_feature_vector(sparse.csr_matrix(np.ones((4, 4))), 1, 2, 3)


def test_detect_zarr_root_prefix(tmp_path: Path) -> None:
    import zipfile

    from reach_gap.xenium_zarr_io import detect_zarr_root_prefix

    archive = tmp_path / "matrix.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("nested/root/.zgroup", "{}")
        handle.writestr("nested/root/cell_features/.zgroup", "{}")
    assert detect_zarr_root_prefix(archive) == "nested/root"


def test_detect_zarr_root_prefix_returns_empty_when_absent(tmp_path: Path) -> None:
    import zipfile

    from reach_gap.xenium_zarr_io import detect_zarr_root_prefix

    archive = tmp_path / "not_zarr.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("README.txt", "not a zarr")
    assert detect_zarr_root_prefix(archive) == ""
