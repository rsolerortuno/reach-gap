from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from reach_gap.benchmark import run_benchmark, write_benchmark
from reach_gap.cli import app
from reach_gap.manifest import validate_manifest
from reach_gap.reporting import render_report

runner = CliRunner()


def test_manifest_and_cli(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "source_id": "test",
                "platform": "simulation",
                "segmentation_version": "v1",
                "vessel_definition": "synthetic",
            }
        ),
        encoding="utf-8",
    )
    assert validate_manifest(manifest).source_id == "test"
    result = runner.invoke(app, ["validate-manifest", str(manifest)])
    assert result.exit_code == 0


def test_simulate_solve_index_cli(tmp_path: Path) -> None:
    features = tmp_path / "features.npz"
    result = runner.invoke(app, ["simulate", "--output", str(features), "--size", "18"])
    assert result.exit_code == 0
    solution = tmp_path / "solution.npz"
    result = runner.invoke(app, ["solve", str(features), "--output", str(solution)])
    assert result.exit_code == 0
    output_dir = tmp_path / "index"
    result = runner.invoke(
        app,
        ["index", str(features), "--output-dir", str(output_dir), "--draws", "4"],
    )
    assert result.exit_code == 0
    assert (output_dir / "claims.json").exists()


def test_benchmark_and_report(tmp_path: Path) -> None:
    result = run_benchmark(quick=True)
    assert result["status"] == "PASS"
    output = write_benchmark(result, tmp_path / "benchmark.json")
    report = render_report(output, tmp_path / "report.md")
    assert "simulation" in report.read_text(encoding="utf-8").lower()


def test_validate_artifact_cli_rejects_wrong_content(tmp_path: Path) -> None:
    fake_pdf = tmp_path / "paper.pdf"
    fake_pdf.write_text("<html>not a PDF</html>", encoding="utf-8")
    result = runner.invoke(
        app,
        ["validate-artifact", str(fake_pdf), "--expected-kind", "pdf"],
    )
    assert result.exit_code == 1
    assert '"actual_kind": "html"' in result.output


def test_validate_v071_results_cli(tmp_path: Path) -> None:
    bundle = Path(__file__).resolve().parents[1] / "results/external_validation/v0.7.1"
    output_dir = tmp_path / "validated"
    result = runner.invoke(
        app,
        ["validate-v071-results", str(bundle), "--output-dir", str(output_dir)],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["status"] == "PASS"
    assert payload["summary"]["absolute_reachability"]["status"] == "NOT_COMPUTED"
    assert (output_dir / "validation_v0.7.1.json").is_file()
    assert (output_dir / "claims_v0.7.1.json").is_file()


def test_audit_xenium_zarr_erbb2_cli(monkeypatch: object, tmp_path: Path) -> None:
    import numpy as np
    from scipy import sparse

    from reach_gap import cli
    from reach_gap.xenium_zarr_io import XeniumSparseMatrix

    matrix = XeniumSparseMatrix(
        matrix=sparse.csr_matrix(np.asarray([[1, 0, 2], [0, 5, 6]], dtype=float)),
        feature_names=("ACTB", "ERBB2"),
        cell_ids=("abcdefgh-1", "ijklmnop-1", "aaaaaaaa-1"),
        array_names=("cell_features/data", "__matrix_encoding__/CSR_FEATURES_BY_CELLS"),
        root_prefix="",
        encoding="CSR_FEATURES_BY_CELLS",
    )
    monkeypatch.setattr(
        cli,
        "read_xenium_sparse_matrix",
        lambda _: matrix,
    )
    output = tmp_path / "audit.json"
    result = runner.invoke(
        app,
        ["audit-xenium-zarr-erbb2", str(tmp_path / "matrix.zarr.zip"), "--output", str(output)],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["selected_feature_name"] == "ERBB2"
    assert payload["ERBB2_RNA_total_counts"] == 11.0
    assert output.is_file()


def test_v080_cli_commands(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    evidence_dir = tmp_path / "evidence"
    result = runner.invoke(
        app,
        [
            "synthesize-evidence",
            "--repository-root",
            str(root),
            "--output-dir",
            str(evidence_dir),
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["absolute_readiness_score"] == 40.5

    relative_dir = tmp_path / "relative"
    result = runner.invoke(
        app,
        [
            "rank-relative-accessibility",
            "--repository-root",
            str(root),
            "--output-dir",
            str(relative_dir),
            "--draws",
            "2000",
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["stable_top_target"] == "VISTA"


def test_validate_v080_results_cli(tmp_path: Path) -> None:
    from reach_gap.v080_results import build_v080_package

    root = Path(__file__).resolve().parents[1]
    package = tmp_path / "v080"
    build_v080_package(root, package)
    result = runner.invoke(app, ["validate-v080-results", str(package)])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["status"] == "PASS"
    assert payload["summary"]["stable_top_target"] == "VISTA"
