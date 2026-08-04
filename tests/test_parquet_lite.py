from __future__ import annotations

import struct

import numpy as np

from reach_gap.parquet_lite import _CompactReader, _decode_plain


def test_compact_reader_decodes_simple_struct() -> None:
    # Field 1, compact I32, value 7 encoded as zig-zag 14; then STOP.
    reader = _CompactReader(bytes([0x15, 0x0E, 0x00]))
    assert reader.read_struct() == {1: 7}
    assert reader.position == 3


def test_plain_numeric_decoding() -> None:
    payload = struct.pack("<ddd", 1.5, -2.0, 9.25)
    observed = _decode_plain(payload, 5, 3)
    assert isinstance(observed, np.ndarray)
    np.testing.assert_allclose(observed, [1.5, -2.0, 9.25])


def test_plain_byte_array_decoding() -> None:
    payload = struct.pack("<I", 3) + b"abc" + struct.pack("<I", 1) + b"z"
    assert _decode_plain(payload, 6, 2) == ["abc", "z"]


def test_plain_decoding_accepts_small_zero_padding() -> None:
    numeric = struct.pack("<dd", 1.0, 2.0) + (b"\x00" * 8)
    observed = _decode_plain(numeric, 5, 2)
    assert isinstance(observed, np.ndarray)
    np.testing.assert_allclose(observed, [1.0, 2.0])

    text = struct.pack("<I", 3) + b"abc" + (b"\x00" * 8)
    assert _decode_plain(text, 6, 1) == ["abc"]


def test_plain_decoding_rejects_nonzero_padding() -> None:
    payload = struct.pack("<d", 1.0) + b"\x01"
    try:
        _decode_plain(payload, 5, 1)
    except ValueError as error:
        assert "trailing bytes" in str(error)
    else:
        raise AssertionError("non-zero padding must be rejected")
