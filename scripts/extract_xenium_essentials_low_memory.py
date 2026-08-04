"""Run low-memory Xenium essential-member extraction without installing the CLI."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from reach_gap.bundle_tools import extract_xenium_essentials_low_memory


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("zip_path", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("xenium-essential-package"))
    parser.add_argument("--part-size-mb", type=int, default=95)
    arguments = parser.parse_args()
    result = extract_xenium_essentials_low_memory(
        arguments.zip_path,
        arguments.output_dir,
        part_size_bytes=arguments.part_size_mb * 1_000_000,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
