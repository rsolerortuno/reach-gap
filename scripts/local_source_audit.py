"""Dependency-free local source audit used when Ruff is unavailable."""

from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATHS = [ROOT / "src", ROOT / "tests", ROOT / "scripts"]


def main() -> None:
    """Check syntax, line length and clinical-label isolation in parameter code."""

    failures: list[str] = []
    for directory in PATHS:
        for path in sorted(directory.rglob("*.py")):
            source = path.read_text(encoding="utf-8")
            try:
                ast.parse(source, filename=str(path))
            except SyntaxError as error:
                failures.append(f"{path}: syntax error: {error}")
            for line_number, line in enumerate(source.splitlines(), start=1):
                if len(line) > 100:
                    failures.append(f"{path}:{line_number}: line length {len(line)}")
    parameter_source = (
        (ROOT / "src" / "reach_gap" / "config.py").read_text(encoding="utf-8").lower()
    )
    for forbidden in ("outcome", "retrospective", "success", "failure"):
        if forbidden in parameter_source:
            failures.append(f"config.py contains forbidden clinical label token: {forbidden}")
    if failures:
        raise SystemExit("\n".join(failures))
    print("local source audit: PASS")


if __name__ == "__main__":
    main()
