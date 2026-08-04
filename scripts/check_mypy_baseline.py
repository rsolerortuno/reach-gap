"""Fail when legacy strict-Mypy debt increases above a recorded ceiling."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum-errors", type=int, required=True)
    parser.add_argument("--output", type=Path, default=Path("results/mypy_legacy_current.txt"))
    args = parser.parse_args()

    command = [
        sys.executable,
        "-m",
        "mypy",
        "--strict",
        "--python-version",
        "3.11",
        "--follow-imports=skip",
        "src/reach_gap",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    text = completed.stdout + completed.stderr
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    match = re.search(r"Found (\d+) errors? in", text)
    errors = int(match.group(1)) if match else (0 if completed.returncode == 0 else -1)
    if errors < 0:
        print(text)
        print("Could not parse Mypy error count.", file=sys.stderr)
        return 2
    print(f"legacy_mypy_errors={errors}; maximum={args.maximum_errors}")
    if errors > args.maximum_errors:
        print("Legacy strict-Mypy debt increased.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
