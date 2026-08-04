"""Minimal, audited Parquet reader for flat Xenium tables.

The project normally uses ``pyarrow`` through pandas. Some execution sandboxes do not
provide an Arrow wheel, even though the Xenium cell and boundary tables are otherwise
small enough to process. This module implements only the narrow Parquet subset emitted
by the public Xenium RCC bundle:

* flat optional columns;
* PLAIN value encoding;
* Zstandard compression;
* DataPage V1 or V2;
* BYTE_ARRAY, INT32, INT64, FLOAT and DOUBLE physical types;
* no null values in requested columns.

Anything outside that subset raises an explicit error. It is not a general Parquet
implementation and must never silently reinterpret unsupported encodings.
"""

from __future__ import annotations

import ctypes
import ctypes.util
import struct
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

import numpy as np
import pandas as pd

# Thrift compact protocol type tags.
_STOP = 0
_BOOLEAN_TRUE = 1
_BOOLEAN_FALSE = 2
_BYTE = 3
_I16 = 4
_I32 = 5
_I64 = 6
_DOUBLE = 7
_BINARY = 8
_LIST = 9
_SET = 10
_MAP = 11
_STRUCT = 12

# Parquet enums used by the supported Xenium files.
_PAGE_DATA = 0
_PAGE_DICTIONARY = 2
_PAGE_DATA_V2 = 3
_ENCODING_PLAIN = 0
_CODEC_ZSTD = 6

# Parquet physical types.
_TYPE_INT32 = 1
_TYPE_INT64 = 2
_TYPE_FLOAT = 4
_TYPE_DOUBLE = 5
_TYPE_BYTE_ARRAY = 6


class UnsupportedParquetError(ValueError):
    """Raised when an input requires Parquet features outside the audited subset."""


class _CompactReader:
    """Read anonymous Thrift compact-protocol values without generated classes."""

    def __init__(self, data: bytes) -> None:
        self.data = data
        self.position = 0

    def _read_byte(self) -> int:
        if self.position >= len(self.data):
            raise EOFError("Unexpected end of Thrift compact payload")
        value = self.data[self.position]
        self.position += 1
        return value

    def _read_varint(self) -> int:
        shift = 0
        value = 0
        while True:
            byte = self._read_byte()
            value |= (byte & 0x7F) << shift
            if byte < 0x80:
                return value
            shift += 7
            if shift > 70:
                raise ValueError("Invalid oversized Thrift varint")

    def _read_integer(self) -> int:
        raw = self._read_varint()
        return (raw >> 1) ^ -(raw & 1)

    def _read_binary(self) -> bytes:
        length = self._read_varint()
        stop = self.position + length
        if stop > len(self.data):
            raise EOFError("Truncated Thrift binary value")
        value = self.data[self.position : stop]
        self.position = stop
        return value

    def _read_value(self, type_tag: int) -> Any:
        if type_tag == _BOOLEAN_TRUE:
            return True
        if type_tag == _BOOLEAN_FALSE:
            return False
        if type_tag == _BYTE:
            raw = self._read_byte()
            return raw - 256 if raw >= 128 else raw
        if type_tag in {_I16, _I32, _I64}:
            return self._read_integer()
        if type_tag == _DOUBLE:
            stop = self.position + 8
            if stop > len(self.data):
                raise EOFError("Truncated Thrift double")
            value = struct.unpack_from("<d", self.data, self.position)[0]
            self.position = stop
            return value
        if type_tag == _BINARY:
            return self._read_binary()
        if type_tag in {_LIST, _SET}:
            header = self._read_byte()
            size = header >> 4
            element_type = header & 0x0F
            if size == 15:
                size = self._read_varint()
            return [self._read_value(element_type) for _ in range(size)]
        if type_tag == _MAP:
            size = self._read_varint()
            if size == 0:
                return {}
            types = self._read_byte()
            key_type = types >> 4
            value_type = types & 0x0F
            return {self._read_value(key_type): self._read_value(value_type) for _ in range(size)}
        if type_tag == _STRUCT:
            return self.read_struct()
        raise ValueError(f"Unknown Thrift compact type tag: {type_tag}")

    def read_struct(self) -> dict[int, Any]:
        fields: dict[int, Any] = {}
        previous_field_id = 0
        while True:
            header = self._read_byte()
            if header == _STOP:
                return fields
            delta = header >> 4
            type_tag = header & 0x0F
            field_id = previous_field_id + delta if delta else self._read_integer()
            previous_field_id = field_id
            fields[field_id] = self._read_value(type_tag)


@dataclass(frozen=True)
class _SchemaColumn:
    name: str
    physical_type: int
    repetition_type: int | None


@dataclass(frozen=True)
class _ColumnChunk:
    name: str
    physical_type: int
    codec: int
    num_values: int
    total_compressed_size: int
    data_page_offset: int
    dictionary_page_offset: int | None
    repetition_type: int | None


@dataclass(frozen=True)
class ParquetFileInfo:
    """Compact metadata summary used for validation and provenance."""

    num_rows: int
    columns: tuple[str, ...]
    row_groups: int
    created_by: str | None


class _ZstdDecoder:
    """Small ctypes binding to the system libzstd shared library."""

    def __init__(self) -> None:
        library_name = ctypes.util.find_library("zstd") or "libzstd.so.1"
        self._library = ctypes.CDLL(library_name)
        self._library.ZSTD_decompress.argtypes = [
            ctypes.c_void_p,
            ctypes.c_size_t,
            ctypes.c_void_p,
            ctypes.c_size_t,
        ]
        self._library.ZSTD_decompress.restype = ctypes.c_size_t
        self._library.ZSTD_isError.argtypes = [ctypes.c_size_t]
        self._library.ZSTD_isError.restype = ctypes.c_uint
        self._library.ZSTD_getErrorName.argtypes = [ctypes.c_size_t]
        self._library.ZSTD_getErrorName.restype = ctypes.c_char_p

    def decompress(self, payload: bytes, expected_size: int) -> bytes:
        if expected_size < 0:
            raise ValueError("Negative expected Zstandard output size")
        destination = ctypes.create_string_buffer(expected_size)
        source = ctypes.create_string_buffer(payload)
        observed_size = int(
            self._library.ZSTD_decompress(
                destination,
                expected_size,
                source,
                len(payload),
            )
        )
        if self._library.ZSTD_isError(observed_size):
            message = self._library.ZSTD_getErrorName(observed_size).decode("utf-8")
            raise ValueError(f"Zstandard decompression failed: {message}")
        if observed_size != expected_size:
            raise ValueError(
                "Zstandard output length differs from Parquet page metadata: "
                f"{observed_size} != {expected_size}"
            )
        return destination.raw[:observed_size]


_ZSTD = _ZstdDecoder()


def _read_footer(path: Path) -> dict[int, Any]:
    with path.open("rb") as handle:
        if handle.read(4) != b"PAR1":
            raise ValueError(f"Not an Apache Parquet file: {path}")
        handle.seek(-8, 2)
        metadata_length = struct.unpack("<I", handle.read(4))[0]
        if handle.read(4) != b"PAR1":
            raise ValueError(f"Parquet footer magic is missing: {path}")
        handle.seek(-8 - metadata_length, 2)
        payload = handle.read(metadata_length)
    reader = _CompactReader(payload)
    metadata = reader.read_struct()
    if reader.position != len(payload):
        raise ValueError("Parquet footer contains unread bytes")
    return metadata


def _decode_text(value: bytes | str) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _schema(metadata: dict[int, Any]) -> list[_SchemaColumn]:
    elements = metadata.get(2)
    if not isinstance(elements, list) or len(elements) < 2:
        raise UnsupportedParquetError("Parquet schema is missing or nested unexpectedly")
    root = elements[0]
    if int(root.get(5, -1)) != len(elements) - 1:
        raise UnsupportedParquetError("Only flat Parquet schemas are supported")
    columns: list[_SchemaColumn] = []
    for element in elements[1:]:
        if 1 not in element or 4 not in element:
            raise UnsupportedParquetError("Nested or untyped schema element")
        columns.append(
            _SchemaColumn(
                name=_decode_text(element[4]),
                physical_type=int(element[1]),
                repetition_type=(int(element[3]) if 3 in element else None),
            )
        )
    return columns


def inspect_flat_parquet(path: Path) -> ParquetFileInfo:
    """Return metadata without loading any row data."""

    metadata = _read_footer(path)
    schema = _schema(metadata)
    created_by = _decode_text(metadata[6]) if 6 in metadata else None
    row_groups = metadata.get(4, [])
    return ParquetFileInfo(
        num_rows=int(metadata.get(3, 0)),
        columns=tuple(column.name for column in schema),
        row_groups=len(row_groups),
        created_by=created_by,
    )


def _chunks_for_row_group(
    row_group: dict[int, Any], schema: Sequence[_SchemaColumn]
) -> dict[str, _ColumnChunk]:
    raw_chunks = row_group.get(1)
    if not isinstance(raw_chunks, list) or len(raw_chunks) != len(schema):
        raise UnsupportedParquetError("Row-group columns do not match flat schema")
    chunks: dict[str, _ColumnChunk] = {}
    for schema_column, raw_chunk in zip(schema, raw_chunks, strict=True):
        column_metadata = raw_chunk.get(3)
        if not isinstance(column_metadata, dict):
            raise UnsupportedParquetError("Encrypted or external Parquet columns are unsupported")
        path_in_schema = column_metadata.get(3, [])
        if [_decode_text(value) for value in path_in_schema] != [schema_column.name]:
            raise UnsupportedParquetError("Nested Parquet column paths are unsupported")
        chunks[schema_column.name] = _ColumnChunk(
            name=schema_column.name,
            physical_type=int(column_metadata[1]),
            codec=int(column_metadata[4]),
            num_values=int(column_metadata[5]),
            total_compressed_size=int(column_metadata[7]),
            data_page_offset=int(column_metadata[9]),
            dictionary_page_offset=(int(column_metadata[11]) if 11 in column_metadata else None),
            repetition_type=schema_column.repetition_type,
        )
    return chunks


def _read_page_header(handle: BinaryIO, offset: int) -> tuple[dict[int, Any], int]:
    handle.seek(offset)
    # Xenium page headers are tens of bytes. The generous cap avoids byte-wise I/O
    # while preserving a hard upper bound for corrupt inputs.
    payload = handle.read(65_536)
    reader = _CompactReader(payload)
    header = reader.read_struct()
    if reader.position > 16_384:
        raise UnsupportedParquetError("Unexpectedly large Parquet page header")
    return header, reader.position


def _strip_v1_levels(payload: bytes, *, optional: bool) -> bytes:
    if not optional:
        return payload
    if len(payload) < 4:
        raise ValueError("Truncated DataPage V1 definition-level prefix")
    level_length = struct.unpack_from("<I", payload, 0)[0]
    stop = 4 + level_length
    if stop > len(payload):
        raise ValueError("Truncated DataPage V1 definition levels")
    return payload[stop:]


def _trim_zero_padding(payload: bytes, expected: int, label: str) -> bytes:
    """Accept only the small all-zero page padding emitted by the Xenium writer."""

    if len(payload) < expected:
        raise ValueError(f"{label} page size mismatch: {len(payload)} < {expected}")
    trailing = payload[expected:]
    if trailing and (len(trailing) > 16 or any(trailing)):
        raise ValueError(f"{label} page contains unsupported trailing bytes: {len(trailing)}")
    return payload[:expected]


def _decode_plain(
    payload: bytes,
    physical_type: int,
    count: int,
) -> list[Any] | np.ndarray:
    if physical_type == _TYPE_DOUBLE:
        expected = count * 8
        payload = _trim_zero_padding(payload, expected, "DOUBLE")
        return np.frombuffer(payload, dtype="<f8").copy()
    if physical_type == _TYPE_FLOAT:
        expected = count * 4
        payload = _trim_zero_padding(payload, expected, "FLOAT")
        return np.frombuffer(payload, dtype="<f4").copy()
    if physical_type == _TYPE_INT32:
        expected = count * 4
        payload = _trim_zero_padding(payload, expected, "INT32")
        return np.frombuffer(payload, dtype="<i4").copy()
    if physical_type == _TYPE_INT64:
        expected = count * 8
        payload = _trim_zero_padding(payload, expected, "INT64")
        return np.frombuffer(payload, dtype="<i8").copy()
    if physical_type == _TYPE_BYTE_ARRAY:
        values: list[str] = []
        position = 0
        for _ in range(count):
            if position + 4 > len(payload):
                raise ValueError("Truncated BYTE_ARRAY length")
            length = struct.unpack_from("<I", payload, position)[0]
            position += 4
            stop = position + length
            if stop > len(payload):
                raise ValueError("Truncated BYTE_ARRAY value")
            values.append(payload[position:stop].decode("utf-8"))
            position = stop
        trailing = payload[position:]
        if trailing and (len(trailing) > 16 or any(trailing)):
            raise ValueError("BYTE_ARRAY page contains unsupported trailing bytes")
        return values
    raise UnsupportedParquetError(f"Unsupported Parquet physical type: {physical_type}")


def _read_chunk(handle: BinaryIO, chunk: _ColumnChunk) -> list[Any] | np.ndarray:
    if chunk.codec != _CODEC_ZSTD:
        raise UnsupportedParquetError(
            f"Only Zstandard-compressed Xenium columns are supported, got codec {chunk.codec}"
        )
    if chunk.dictionary_page_offset is not None:
        raise UnsupportedParquetError("Dictionary-encoded pages are unsupported")
    offset = chunk.data_page_offset
    values: list[Any] = []
    numeric_parts: list[np.ndarray] = []
    values_read = 0
    while values_read < chunk.num_values:
        header, header_size = _read_page_header(handle, offset)
        page_type = int(header[1])
        uncompressed_size = int(header[2])
        compressed_size = int(header[3])
        handle.seek(offset + header_size)
        compressed_payload = handle.read(compressed_size)
        if len(compressed_payload) != compressed_size:
            raise EOFError("Truncated Parquet page payload")

        if page_type == _PAGE_DICTIONARY:
            raise UnsupportedParquetError("Dictionary pages are unsupported")
        if page_type == _PAGE_DATA_V2:
            page = header.get(8)
            if not isinstance(page, dict):
                raise ValueError("DataPage V2 header is missing")
            count = int(page[1])
            null_count = int(page[2])
            encoding = int(page[4])
            repetition_length = int(page[5])
            definition_length = int(page[6])
            if null_count != 0:
                raise UnsupportedParquetError("Null values are unsupported")
            if encoding != _ENCODING_PLAIN:
                raise UnsupportedParquetError(
                    f"Only PLAIN data pages are supported, got encoding {encoding}"
                )
            levels_length = repetition_length + definition_length
            value_payload = compressed_payload[levels_length:]
            expected_value_size = uncompressed_size - levels_length
            if bool(page.get(7, True)):
                value_payload = _ZSTD.decompress(value_payload, expected_value_size)
            elif len(value_payload) != expected_value_size:
                raise ValueError("Uncompressed DataPage V2 size mismatch")
        elif page_type == _PAGE_DATA:
            page = header.get(5)
            if not isinstance(page, dict):
                raise ValueError("DataPage V1 header is missing")
            count = int(page[1])
            encoding = int(page[2])
            if encoding != _ENCODING_PLAIN:
                raise UnsupportedParquetError(
                    f"Only PLAIN data pages are supported, got encoding {encoding}"
                )
            uncompressed_payload = _ZSTD.decompress(compressed_payload, uncompressed_size)
            value_payload = _strip_v1_levels(
                uncompressed_payload,
                optional=chunk.repetition_type == 1,
            )
        else:
            raise UnsupportedParquetError(f"Unsupported Parquet page type: {page_type}")

        decoded = _decode_plain(value_payload, chunk.physical_type, count)
        if isinstance(decoded, np.ndarray):
            numeric_parts.append(decoded)
        else:
            values.extend(decoded)
        values_read += count
        offset += header_size + compressed_size

    if values_read != chunk.num_values:
        raise ValueError(
            f"Column {chunk.name} value count mismatch: {values_read} != {chunk.num_values}"
        )
    if numeric_parts:
        return np.asarray(np.concatenate(numeric_parts))
    return values


def read_flat_parquet(
    path: Path,
    columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Read selected columns from the audited flat Xenium Parquet subset."""

    metadata = _read_footer(path)
    schema = _schema(metadata)
    available = [column.name for column in schema]
    requested = available if columns is None else list(columns)
    missing = sorted(set(requested).difference(available))
    if missing:
        raise KeyError(f"Parquet columns are missing: {missing}")

    parts: dict[str, list[list[Any] | np.ndarray[Any, Any]]] = {name: [] for name in requested}
    row_groups = metadata.get(4, [])
    with path.open("rb") as handle:
        for row_group in row_groups:
            chunks = _chunks_for_row_group(row_group, schema)
            for name in requested:
                decoded = _read_chunk(handle, chunks[name])
                parts[name].append(decoded)

    data: dict[str, Any] = {}
    for name in requested:
        column_parts = parts[name]
        if not column_parts:
            data[name] = []
        elif isinstance(column_parts[0], np.ndarray):
            data[name] = np.concatenate(
                [part for part in column_parts if isinstance(part, np.ndarray)]
            )
        else:
            flattened: list[Any] = []
            for part in column_parts:
                if isinstance(part, np.ndarray):
                    flattened.extend(part.tolist())
                else:
                    flattened.extend(part)
            data[name] = flattened
    frame = pd.DataFrame(data)
    expected_rows = int(metadata.get(3, 0))
    if len(frame) != expected_rows:
        raise ValueError(f"Parquet row count mismatch: {len(frame)} != {expected_rows}")
    return frame


def iter_flat_parquet_row_groups(
    path: Path,
    columns: Sequence[str],
) -> Iterable[pd.DataFrame]:
    """Yield selected flat columns one row group at a time for bounded memory."""

    metadata = _read_footer(path)
    schema = _schema(metadata)
    available = [column.name for column in schema]
    missing = sorted(set(columns).difference(available))
    if missing:
        raise KeyError(f"Parquet columns are missing: {missing}")
    with path.open("rb") as handle:
        for row_group in metadata.get(4, []):
            chunks = _chunks_for_row_group(row_group, schema)
            data = {name: _read_chunk(handle, chunks[name]) for name in columns}
            yield pd.DataFrame(data)
