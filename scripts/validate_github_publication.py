from __future__ import annotations

import argparse
import hashlib
import json
import re
import tomllib
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlsplit

import yaml

ROOT = Path(__file__).resolve().parents[1]
REQUIRED_README_HEADINGS = {
    "The problem",
    "The idea",
    "Headline results",
    "How the pipeline works",
    "Quick start",
    "Tests and reproducibility",
    "Explicitly uncomputed claims",
    "Limitations",
    "Future plans",
    "Conclusion",
}
REQUIRED_FIGURES = {
    "reach_gap_workflow.png",
    "rcc_target_geometry.png",
    "rcc_target_maps_panel.png",
    "target_rank_probability.png",
    "pairwise_win_probability.png",
    "leave_one_component_out.png",
    "measurement_priority.png",
    "external_validation_panel.png",
    "validation_summary.png",
}
REQUIRED_WORKFLOWS = {"ci.yml", "release.yml"}
MARKDOWN_LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def local_markdown_targets(path: Path) -> list[Path]:
    targets: list[Path] = []
    for match in MARKDOWN_LINK.finditer(path.read_text(encoding="utf-8")):
        raw = match.group(1).strip().split()[0].strip("<>")
        parsed = urlsplit(raw)
        if parsed.scheme or parsed.netloc or raw.startswith("#"):
            continue
        relative = unquote(parsed.path)
        if not relative:
            continue
        targets.append((path.parent / relative).resolve())
    return targets


def validate() -> dict[str, Any]:
    issues: list[str] = []
    with (ROOT / "pyproject.toml").open("rb") as handle:
        project = tomllib.load(handle)["project"]
    version = str(project["version"])
    if version != "0.8.0":
        issues.append(f"Unexpected package version: {version}")

    readme = ROOT / "README.md"
    readme_text = readme.read_text(encoding="utf-8")
    headings = {line[3:].strip() for line in readme_text.splitlines() if line.startswith("## ")}
    missing_headings = sorted(REQUIRED_README_HEADINGS - headings)
    if missing_headings:
        issues.append(f"README missing headings: {missing_headings}")

    markdown_files = [
        ROOT / "README.md",
        ROOT / "CHANGELOG.md",
        ROOT / "CONTRIBUTING.md",
        ROOT / "RELEASE_NOTES_v0.8.0.md",
        ROOT / "GITHUB_RELEASE_BODY.md",
        ROOT / "reports/README.md",
        ROOT / "reports/model_card.md",
        ROOT / "docs/GITHUB_RELEASE.md",
    ]
    broken_links: list[str] = []
    for markdown in markdown_files:
        if not markdown.is_file():
            broken_links.append(f"missing markdown file: {markdown.relative_to(ROOT)}")
            continue
        for target in local_markdown_targets(markdown):
            try:
                target.relative_to(ROOT.resolve())
            except ValueError:
                broken_links.append(
                    f"link escapes repository: {markdown.relative_to(ROOT)} -> {target}"
                )
                continue
            if not target.exists():
                broken_links.append(
                    f"broken link: {markdown.relative_to(ROOT)} -> {target.relative_to(ROOT)}"
                )
    if broken_links:
        issues.extend(broken_links)

    figure_dir = ROOT / "reports" / "figures"
    figure_files = {path.name for path in figure_dir.glob("*.png")}
    missing_figures = sorted(REQUIRED_FIGURES - figure_files)
    if missing_figures:
        issues.append(f"Missing portfolio figures: {missing_figures}")
    undersized_figures = sorted(
        path.name for path in figure_dir.glob("*.png") if path.stat().st_size < 25_000
    )
    if undersized_figures:
        issues.append(f"Suspiciously small figures: {undersized_figures}")

    yaml_paths = sorted((ROOT / ".github").rglob("*.yml"))
    parsed_yaml: dict[str, Any] = {}
    for path in yaml_paths:
        try:
            parsed_yaml[str(path.relative_to(ROOT))] = yaml.safe_load(
                path.read_text(encoding="utf-8")
            )
        except yaml.YAMLError as exc:
            issues.append(f"Invalid YAML {path.relative_to(ROOT)}: {exc}")
    workflow_names = {path.name for path in (ROOT / ".github/workflows").glob("*.yml")}
    if not REQUIRED_WORKFLOWS.issubset(workflow_names):
        issues.append(f"Missing workflows: {sorted(REQUIRED_WORKFLOWS - workflow_names)}")

    ci_text = (ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")
    release_text = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    for command in ["pytest -q", "ruff check", "mypy src/reach_gap", "python -m build"]:
        if command not in ci_text and command not in release_text:
            issues.append(f"Workflows do not contain required command: {command}")
    for command in [
        "validate-v080-results",
        "generate_portfolio_figures.py",
        "sha256sum",
        "softprops/action-gh-release@v2",
    ]:
        if command not in release_text:
            issues.append(f"Release workflow missing: {command}")

    build_validation_path = ROOT / "results/build_validation_v0.8.json"
    build_validation = json.loads(build_validation_path.read_text(encoding="utf-8"))
    coverage_validation = json.loads(
        (ROOT / "results/coverage_v0.8.json").read_text(encoding="utf-8")
    )
    checks = build_validation["checks"]
    if checks["pytest"]["status"] != "PASS" or checks["pytest"]["tests_failed"] != 0:
        issues.append("Bundled pytest validation is not PASS")
    if checks["coverage"]["percent"] < checks["coverage"]["threshold_percent"]:
        issues.append("Bundled coverage is below its threshold")
    if checks["scientific_status"]["status"] != "PASS_WITH_ABSTENTION":
        issues.append("Scientific status does not preserve the required abstention")

    files_for_manifest = [
        ROOT / "README.md",
        ROOT / ".github/workflows/ci.yml",
        ROOT / ".github/workflows/release.yml",
        ROOT / "scripts/generate_portfolio_figures.py",
        ROOT / "scripts/validate_github_publication.py",
        *sorted(figure_dir.glob("*.png")),
    ]
    manifest = [
        {
            "path": str(path.relative_to(ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in files_for_manifest
    ]
    return {
        "schema_version": "1.0",
        "repository_version": version,
        "status": "PASS" if not issues else "FAIL",
        "checks": {
            "readme_required_headings": len(REQUIRED_README_HEADINGS),
            "markdown_files_checked": len(markdown_files),
            "yaml_files_checked": len(yaml_paths),
            "portfolio_figures": len(figure_files),
            "bundled_tests_passed": checks["pytest"]["tests_passed"],
            "bundled_coverage_percent": coverage_validation["coverage_percent"],
        },
        "issues": issues,
        "publication_manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate the GitHub publication structure for reach-gap."
    )
    parser.add_argument("--output", type=Path, help="Optional JSON output path")
    args = parser.parse_args()
    payload = validate()
    rendered = json.dumps(payload, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    if payload["issues"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
