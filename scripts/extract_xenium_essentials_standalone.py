#!/usr/bin/env python3
"""Extract transferable essential files from a Xenium output ZIP with bounded memory.

This script uses only the Python standard library. It reads the ZIP central directory,
selects the minimum cell-level files needed by reach-gap, and streams large members
directly into independently checksummed parts. It never loads the source ZIP or a full
archive member into RAM and never creates a second unsplit copy of a large member.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
import zipfile
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path, PurePosixPath
from typing import BinaryIO

DEFAULT_PART_SIZE = 95_000_000
DEFAULT_CHUNK_SIZE = 8 * 1024 * 1024
SMALL_ANALYSIS_LIMIT = 100_000_000


@dataclass(frozen=True)
class MemberInfo:
    member: str
    basename: str
    uncompressed_bytes: int
    compressed_bytes: int
    crc32: str
    is_dir: bool


def _digest_stream(source: BinaryIO, target: BinaryIO, chunk_size: int) -> dict[str, object]:
    sha256 = hashlib.sha256()
    md5 = hashlib.md5()
    total = 0
    while True:
        chunk = source.read(chunk_size)
        if not chunk:
            break
        target.write(chunk)
        sha256.update(chunk)
        md5.update(chunk)
        total += len(chunk)
    return {"bytes": total, "sha256": sha256.hexdigest(), "md5": md5.hexdigest()}


def inspect_archive(zip_path: Path) -> list[MemberInfo]:
    with zipfile.ZipFile(zip_path) as archive:
        return [
            MemberInfo(
                member=info.filename,
                basename=PurePosixPath(info.filename).name,
                uncompressed_bytes=int(info.file_size),
                compressed_bytes=int(info.compress_size),
                crc32=f"{info.CRC:08x}",
                is_dir=info.is_dir(),
            )
            for info in archive.infolist()
        ]


def select_essential_members(inventory: Iterable[MemberInfo]) -> list[MemberInfo]:
    rows = list(inventory)
    by_basename: dict[str, list[MemberInfo]] = {}
    for row in rows:
        by_basename.setdefault(row.basename, []).append(row)

    selected: dict[str, MemberInfo] = {}
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
            candidates = sorted(by_basename.get(basename, []), key=lambda row: row.member)
            if candidates:
                selected[candidates[0].member] = candidates[0]
                break

    for row in rows:
        lowered = f"/{row.member.lower()}"
        if (
            not row.is_dir
            and lowered.endswith(".csv")
            and "/analysis/" in lowered
            and row.uncompressed_bytes <= SMALL_ANALYSIS_LIMIT
        ):
            selected[row.member] = row

    return [selected[key] for key in sorted(selected)]


def _safe_destination(root: Path, member: str) -> Path:
    destination = (root / PurePosixPath(member)).resolve()
    resolved_root = root.resolve()
    if destination != resolved_root and resolved_root not in destination.parents:
        raise ValueError(f"Unsafe ZIP member path: {member}")
    return destination


def _write_single(
    archive: zipfile.ZipFile,
    member: MemberInfo,
    output_root: Path,
    chunk_size: int,
) -> dict[str, object]:
    destination = _safe_destination(output_root / "files", member.member)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with archive.open(member.member) as source, destination.open("wb") as target:
        digests = _digest_stream(source, target, chunk_size)
    if int(digests["bytes"]) != member.uncompressed_bytes:
        raise ValueError(f"Byte-count mismatch for {member.member}")
    return {
        "mode": "single_file",
        "member": member.member,
        "path": str(destination),
        **digests,
        "zip_crc32": member.crc32,
        "compressed_bytes": member.compressed_bytes,
    }


def _write_parts(
    archive: zipfile.ZipFile,
    member: MemberInfo,
    output_root: Path,
    part_size: int,
    chunk_size: int,
) -> dict[str, object]:
    basename = member.basename
    target_dir = output_root / "split" / basename
    target_dir.mkdir(parents=True, exist_ok=True)
    part_count = math.ceil(member.uncompressed_bytes / part_size)
    source_sha256 = hashlib.sha256()
    source_md5 = hashlib.md5()
    parts: list[dict[str, object]] = []
    total = 0

    with archive.open(member.member) as source:
        for index in range(1, part_count + 1):
            part_name = f"{basename}.part{index:03d}-of-{part_count:03d}"
            part_path = target_dir / part_name
            part_sha256 = hashlib.sha256()
            part_md5 = hashlib.md5()
            written = 0
            with part_path.open("wb") as target:
                while written < part_size:
                    chunk = source.read(min(chunk_size, part_size - written))
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
            raise ValueError(f"ZIP member exceeds declared size: {member.member}")

    if total != member.uncompressed_bytes:
        raise ValueError(
            f"Byte-count mismatch for {member.member}: {total} != {member.uncompressed_bytes}"
        )

    split_manifest = {
        "format_version": 2,
        "source": {
            "name": basename,
            "member": member.member,
            "size": member.uncompressed_bytes,
            "md5Checksum": source_md5.hexdigest(),
            "sha256": source_sha256.hexdigest(),
            "zip_crc32": member.crc32,
        },
        "part_size_bytes": part_size,
        "parts": parts,
        "reassembly": {
            "unix": f"cat {basename}.part* > {basename}",
            "powershell": (
                "$parts = Get-ChildItem '"
                + basename
                + ".part*' | Sort-Object Name; "
                + "$out = [System.IO.File]::Create('"
                + basename
                + "'); try { foreach ($p in $parts) { "
                + "$in = [System.IO.File]::OpenRead($p.FullName); "
                + "try { $in.CopyTo($out) } finally { $in.Dispose() } } } "
                + "finally { $out.Dispose() }"
            ),
            "verification": "Verify each part SHA-256, then reconstructed source MD5/SHA-256.",
        },
    }
    manifest_path = target_dir / f"{basename}.split-manifest.json"
    manifest_path.write_text(json.dumps(split_manifest, indent=2) + "\n", encoding="utf-8")
    return {
        "mode": "split_parts",
        "member": member.member,
        "bytes": total,
        "sha256": source_sha256.hexdigest(),
        "md5": source_md5.hexdigest(),
        "part_count": len(parts),
        "manifest": str(manifest_path),
        "parts_dir": str(target_dir),
        "zip_crc32": member.crc32,
        "compressed_bytes": member.compressed_bytes,
    }


def extract(
    zip_path: Path,
    output_dir: Path,
    *,
    part_size: int,
    chunk_size: int,
    list_only: bool,
) -> dict[str, object]:
    started = time.time()
    inventory = inspect_archive(zip_path)
    selected = select_essential_members(inventory)
    if not selected:
        raise RuntimeError("No expected Xenium essential members were found")

    output_dir.mkdir(parents=True, exist_ok=True)
    inventory_path = output_dir / "xenium_bundle_inventory.json"
    inventory_path.write_text(
        json.dumps([asdict(row) for row in inventory], indent=2) + "\n", encoding="utf-8"
    )
    selection_path = output_dir / "selected_essential_members.json"
    selection_path.write_text(
        json.dumps([asdict(row) for row in selected], indent=2) + "\n", encoding="utf-8"
    )

    if list_only:
        result: dict[str, object] = {
            "status": "ESSENTIAL_MEMBERS_LISTED_NOT_EXTRACTED",
            "source_zip": {"path": str(zip_path), "bytes": zip_path.stat().st_size},
            "archive_member_count": len(inventory),
            "selected_member_count": len(selected),
            "selected_uncompressed_bytes": sum(row.uncompressed_bytes for row in selected),
            "inventory": str(inventory_path),
            "selection": str(selection_path),
        }
    else:
        outputs: list[dict[str, object]] = []
        with zipfile.ZipFile(zip_path) as archive:
            for position, member in enumerate(selected, start=1):
                print(
                    f"[{position}/{len(selected)}] {member.member} "
                    f"({member.uncompressed_bytes / 1e6:.1f} MB)",
                    flush=True,
                )
                if member.uncompressed_bytes > part_size:
                    outputs.append(_write_parts(archive, member, output_dir, part_size, chunk_size))
                else:
                    outputs.append(_write_single(archive, member, output_dir, chunk_size))
        result = {
            "status": "ESSENTIAL_XENIUM_PACKAGE_EXTRACTED",
            "source_zip": {"path": str(zip_path), "bytes": zip_path.stat().st_size},
            "archive_member_count": len(inventory),
            "selected_member_count": len(selected),
            "part_size_bytes": part_size,
            "chunk_size_bytes": chunk_size,
            "members": outputs,
            "memory_model": (
                "Bounded by one chunk plus ZIP decompressor buffers; no full archive "
                "member is loaded and no unsplit copy is created for a large selected member."
            ),
        }

    result["elapsed_wall_seconds"] = round(time.time() - started, 3)
    result["python"] = sys.version
    result["platform"] = sys.platform
    result["pid"] = os.getpid()
    manifest_path = output_dir / "essential_package_manifest.json"
    manifest_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    result["manifest_path"] = str(manifest_path)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract only transferable reach-gap inputs from a Xenium output ZIP."
    )
    parser.add_argument("zip_path", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("xenium-essential-package"))
    parser.add_argument("--part-size-mb", type=int, default=95)
    parser.add_argument("--chunk-size-mb", type=int, default=8)
    parser.add_argument(
        "--list-only",
        action="store_true",
        help="Inspect and record selected members without extracting their bytes.",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    if arguments.part_size_mb <= 0 or arguments.chunk_size_mb <= 0:
        raise SystemExit("part and chunk sizes must be positive")
    result = extract(
        arguments.zip_path,
        arguments.output_dir,
        part_size=arguments.part_size_mb * 1_000_000,
        chunk_size=arguments.chunk_size_mb * 1024 * 1024,
        list_only=arguments.list_only,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
