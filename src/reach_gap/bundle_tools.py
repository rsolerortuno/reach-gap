"""Streaming tools for extracting only essential Xenium bundle members.

These utilities avoid loading a multi-gigabyte ZIP or HDF5 matrix into memory. Large
members can be written directly as independently verifiable parts small enough for
Drive/chat transfer, without first materialising a second full copy on disk.
"""

from __future__ import annotations

import hashlib
import json
import math
import zipfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

from reach_gap.xenium import inspect_zip, select_essential_members


def _safe_member_destination(root: Path, member: str) -> Path:
    destination = (root / member).resolve()
    root_resolved = root.resolve()
    if root_resolved not in destination.parents and destination != root_resolved:
        raise ValueError(f"Unsafe ZIP member path: {member}")
    return destination


def _stream_member_to_file(
    archive: zipfile.ZipFile,
    member: str,
    destination: Path,
    *,
    chunk_size: int,
) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    sha256 = hashlib.sha256()
    md5 = hashlib.md5()
    total = 0
    with archive.open(member) as source, destination.open("wb") as target:
        while chunk := source.read(chunk_size):
            target.write(chunk)
            sha256.update(chunk)
            md5.update(chunk)
            total += len(chunk)
    return {
        "path": str(destination),
        "bytes": total,
        "sha256": sha256.hexdigest(),
        "md5": md5.hexdigest(),
    }


def _stream_member_to_parts(
    archive: zipfile.ZipFile,
    member: str,
    output_dir: Path,
    *,
    expected_size: int,
    part_size_bytes: int,
    chunk_size: int,
) -> dict[str, Any]:
    basename = Path(member).name
    part_count = math.ceil(expected_size / part_size_bytes)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_sha256 = hashlib.sha256()
    source_md5 = hashlib.md5()
    parts: list[dict[str, Any]] = []
    with archive.open(member) as source:
        total = 0
        for index in range(1, part_count + 1):
            part_name = f"{basename}.part{index:03d}-of-{part_count:03d}"
            part_path = output_dir / part_name
            part_sha256 = hashlib.sha256()
            part_md5 = hashlib.md5()
            written = 0
            with part_path.open("wb") as target:
                while written < part_size_bytes:
                    request = min(chunk_size, part_size_bytes - written)
                    chunk = source.read(request)
                    if not chunk:
                        break
                    target.write(chunk)
                    part_sha256.update(chunk)
                    part_md5.update(chunk)
                    source_sha256.update(chunk)
                    source_md5.update(chunk)
                    written += len(chunk)
                    total += len(chunk)
            if written == 0:
                part_path.unlink(missing_ok=True)
                break
            parts.append(
                {
                    "index": index,
                    "name": part_name,
                    "size": written,
                    "sha256": part_sha256.hexdigest(),
                    "md5": part_md5.hexdigest(),
                }
            )
        if source.read(1):
            raise ValueError(f"ZIP member {member} exceeds its declared size")
    if total != expected_size:
        raise ValueError(f"Extracted byte count mismatch for {member}: {total} != {expected_size}")
    split_manifest = {
        "format_version": 2,
        "source": {
            "name": basename,
            "member": member,
            "size": expected_size,
            "md5Checksum": source_md5.hexdigest(),
            "sha256": source_sha256.hexdigest(),
        },
        "part_size_bytes": part_size_bytes,
        "parts": parts,
        "reassembly": {
            "unix": f"cat {basename}.part* > {basename}",
            "verification": "Verify every part SHA-256, then source MD5/SHA-256.",
        },
    }
    manifest_path = output_dir / f"{basename}.split-manifest.json"
    manifest_path.write_text(json.dumps(split_manifest, indent=2), encoding="utf-8")
    return {
        "mode": "split_parts",
        "member": member,
        "bytes": total,
        "sha256": source_sha256.hexdigest(),
        "md5": source_md5.hexdigest(),
        "part_count": len(parts),
        "manifest": str(manifest_path),
        "parts_dir": str(output_dir),
    }


def extract_xenium_essentials_low_memory(
    zip_path: Path,
    output_dir: Path,
    *,
    part_size_bytes: int = 95_000_000,
    chunk_size: int = 8 * 1024 * 1024,
    selected_members: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Extract essential Xenium members, splitting large members during streaming."""

    if part_size_bytes <= 0 or chunk_size <= 0:
        raise ValueError("part_size_bytes and chunk_size must be positive")
    output_dir.mkdir(parents=True, exist_ok=True)
    inventory = inspect_zip(zip_path)
    members = (
        list(selected_members)
        if selected_members is not None
        else select_essential_members(inventory)
    )
    if not members:
        raise ValueError("No essential Xenium members were found in the archive")
    inventory_by_member = {
        str(row.member): {
            "uncompressed_bytes": _coerce_int(row.uncompressed_bytes),
            "compressed_bytes": _coerce_int(row.compressed_bytes),
            "crc32": str(row.crc32),
        }
        for row in inventory.itertuples(index=False)
    }
    results: list[dict[str, Any]] = []
    with zipfile.ZipFile(zip_path) as archive:
        archive_names = set(archive.namelist())
        for member in members:
            if member not in archive_names:
                raise KeyError(f"Selected archive member not found: {member}")
            size = _coerce_int(inventory_by_member[member]["uncompressed_bytes"])
            if size > part_size_bytes:
                result = _stream_member_to_parts(
                    archive,
                    member,
                    output_dir / "split" / Path(member).name,
                    expected_size=size,
                    part_size_bytes=part_size_bytes,
                    chunk_size=chunk_size,
                )
            else:
                destination = _safe_member_destination(output_dir / "files", member)
                result = {
                    "mode": "single_file",
                    "member": member,
                    **_stream_member_to_file(
                        archive,
                        member,
                        destination,
                        chunk_size=chunk_size,
                    ),
                }
            result["zip_crc32"] = inventory_by_member[member]["crc32"]
            result["compressed_bytes"] = inventory_by_member[member]["compressed_bytes"]
            results.append(result)
    package_manifest = {
        "status": "ESSENTIAL_XENIUM_PACKAGE_EXTRACTED",
        "source_zip": {
            "path": str(zip_path),
            "bytes": zip_path.stat().st_size,
        },
        "part_size_bytes": part_size_bytes,
        "selected_member_count": len(members),
        "members": results,
        "memory_model": (
            "Streaming extraction; memory bounded by chunk_size plus ZIP decompressor buffers. "
            "No full archive member is loaded into RAM."
        ),
    }
    manifest_path = output_dir / "essential_package_manifest.json"
    manifest_path.write_text(json.dumps(package_manifest, indent=2), encoding="utf-8")
    package_manifest["manifest_path"] = str(manifest_path)
    return package_manifest


def _coerce_int(value: object) -> int:
    """Convert an inventory scalar to an integer after explicit validation."""

    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (str, bytes, bytearray)):
        return int(value)
    raise TypeError(f"Expected integer-compatible inventory value, got {type(value).__name__}")
