"""Typer command-line interface."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import typer

from reach_gap.artifact_validation import validate_artifact
from reach_gap.benchmark import run_benchmark, write_benchmark
from reach_gap.bordeau_validation import benchmark_bordeau_supplement
from reach_gap.breast_cohort import audit_breast_cell_groups
from reach_gap.bundle_tools import extract_xenium_essentials_low_memory
from reach_gap.cosmx import prepare_cosmx_external_validation
from reach_gap.evidence_synthesis import synthesize_repository_evidence
from reach_gap.features import (
    extract_features,
    ingest_cell_table,
    load_features,
    save_features,
)
from reach_gap.geometry import GeometryConfig, simulate_geometry
from reach_gap.he_pathology import prepare_he_pathology_rcc
from reach_gap.her2_calibration import benchmark_her2_receptor_calibration
from reach_gap.her2_ihc import benchmark_her2_ihc
from reach_gap.indexing import compute_index, write_index_outputs
from reach_gap.manifest import validate_manifest
from reach_gap.perfusion import benchmark_perfusion_tiffs
from reach_gap.real_rcc import prepare_rcc_xenium_essentials
from reach_gap.real_rcc_imaging import (
    discover_rcc_protein_image_inputs,
    prepare_rcc_protein_imaging,
)
from reach_gap.relative_accessibility import run_relative_accessibility
from reach_gap.reporting import render_report
from reach_gap.retrospective import run_retrospective
from reach_gap.schemas import ModelParameters
from reach_gap.sensitivity import sobol_sensitivity
from reach_gap.shg_collagen import benchmark_shg_collagen
from reach_gap.solver import save_solution, solve_transport_robust
from reach_gap.transport_priors import build_igg_transport_prior
from reach_gap.v071_results import build_claims, validate_bundle
from reach_gap.v080_results import build_v080_package, validate_v080_package
from reach_gap.xenium import prepare_rcc_xenium
from reach_gap.xenium_zarr import extract_feature_vector, find_erbb2_feature_index
from reach_gap.xenium_zarr_io import read_xenium_sparse_matrix

app = typer.Typer(help="Mechanistic spatial antibody reachability with uncertainty and abstention.")


@app.command()
def simulate(
    output: Path = typer.Option(Path("simulation_features.npz"), "--output"),
    size: int = typer.Option(48, "--size"),
    seed: int = typer.Option(17, "--seed"),
) -> None:
    """Generate a deterministic synthetic tumour section and extracted features."""

    geometry = simulate_geometry(GeometryConfig(size=size, seed=seed))
    features = extract_features(geometry, antigen_calibration_nM_per_signal=250.0)
    save_features(features, output)
    typer.echo(str(output))


@app.command()
def ingest(
    cells: Path,
    manifest: Path,
    output_dir: Path = typer.Option(Path("run"), "--output-dir"),
    grid_size: int = typer.Option(64, "--grid-size", min=8),
) -> None:
    """Ingest a generic spatial cell table and provenance manifest."""

    features_path, manifest_path = ingest_cell_table(
        cells, manifest, output_dir, grid_size=grid_size
    )
    typer.echo(f"features={features_path}\nmanifest={manifest_path}")


@app.command()
def features(
    output: Path = typer.Option(Path("features.npz"), "--output"),
    seed: int = typer.Option(17, "--seed"),
) -> None:
    """Create model-ready features from a deterministic simulation."""

    geometry = simulate_geometry(GeometryConfig(seed=seed))
    save_features(extract_features(geometry), output)
    typer.echo(str(output))


@app.command()
def solve(
    features_path: Path,
    output: Path = typer.Option(Path("solution.npz"), "--output"),
) -> None:
    """Solve the steady reaction-diffusion system."""

    result = solve_transport_robust(load_features(features_path), ModelParameters())
    if not result.converged:
        raise typer.Exit(code=2)
    save_solution(result, output)
    typer.echo(str(output))


@app.command(name="index")
def index_command(
    features_path: Path,
    output_dir: Path = typer.Option(Path("index_run"), "--output-dir"),
    draws: int = typer.Option(32, "--draws"),
    seed: int = typer.Option(17, "--seed"),
) -> None:
    """Compute uncertainty-aware reachability metrics and claims."""

    output, claims = compute_index(
        load_features(features_path), ModelParameters(), draws=draws, seed=seed
    )
    index_path, claims_path = write_index_outputs(output, claims, output_dir)
    typer.echo(f"index={index_path}\nclaims={claims_path}")


@app.command()
def sensitivity(
    features_path: Path,
    output: Path = typer.Option(Path("sensitivity.json"), "--output"),
    sample_power: int = typer.Option(3, "--sample-power"),
) -> None:
    """Run global Sobol sensitivity analysis."""

    result = sobol_sensitivity(
        load_features(features_path), ModelParameters(), sample_power=sample_power
    )
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    typer.echo(str(output))


@app.command(name="extract-xenium-essentials")
def extract_xenium_essentials_command(
    zip_path: Path,
    output_dir: Path = typer.Option(Path("xenium-essential-package"), "--output-dir"),
    part_size_mb: int = typer.Option(95, "--part-size-mb", min=1, max=500),
) -> None:
    """Stream essential Xenium members into transferable files or verified parts."""

    result = extract_xenium_essentials_low_memory(
        zip_path,
        output_dir,
        part_size_bytes=part_size_mb * 1_000_000,
    )
    typer.echo(json.dumps(result, indent=2))


@app.command(name="prepare-he-rcc")
def prepare_he_rcc_command(
    he_path: Path,
    annotation_path: Path,
    alignment_path: Path | None = typer.Option(None, "--alignment-path"),
    output_dir: Path = typer.Option(Path("reach-gap-he-analysis"), "--output-dir"),
    analysis_level: int = typer.Option(4, "--analysis-level", min=0),
) -> None:
    """Prepare real RCC H&E geometry without asserting target reachability."""

    result = prepare_he_pathology_rcc(
        he_path=he_path,
        annotation_path=annotation_path,
        alignment_path=alignment_path,
        output_dir=output_dir,
        analysis_level=analysis_level,
    )
    typer.echo(json.dumps(result, indent=2))


@app.command(name="prepare-xenium-rcc")
def prepare_xenium_rcc_command(
    raw_dir: Path,
    output_dir: Path = typer.Option(Path("reach-gap-analysis"), "--output-dir"),
    verify_large_md5: bool = typer.Option(True, "--verify-large-md5/--skip-large-md5"),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Prepare the public Xenium RCC gene-and-protein dataset by streaming its bundle."""

    result = prepare_rcc_xenium(
        raw_dir=raw_dir,
        output_dir=output_dir,
        verify_large_md5=verify_large_md5,
        force=force,
    )
    typer.echo(json.dumps(result, indent=2))


@app.command(name="prepare-xenium-essential-rcc")
def prepare_xenium_essential_rcc_command(
    input_dir: Path,
    output_dir: Path = typer.Option(Path("reach-gap-rcc-essential"), "--output-dir"),
    annotation_path: Path | None = typer.Option(None, "--annotation-path"),
    alignment_path: Path | None = typer.Option(None, "--alignment-path"),
    write_full_cell_tables: bool = typer.Option(False, "--write-full-cell-tables/--sample-only"),
) -> None:
    """Prepare the RCC Xenium dataset from compact essential files, without the large ZIP."""

    result = prepare_rcc_xenium_essentials(
        input_dir=input_dir,
        output_dir=output_dir,
        annotation_path=annotation_path,
        alignment_path=alignment_path,
        write_full_cell_tables=write_full_cell_tables,
    )
    typer.echo(json.dumps(result, indent=2))


@app.command(name="prepare-rcc-protein-imaging")
def prepare_rcc_protein_imaging_command(
    scored_cells_dir: Path,
    h5_path: Path,
    morphology_dir: Path,
    output_dir: Path = typer.Option(Path("reach-gap-rcc-protein-imaging"), "--output-dir"),
    qc_mask_dir: Path | None = typer.Option(None, "--qc-mask-dir"),
    background_dir: Path | None = typer.Option(None, "--background-dir"),
    level: int = typer.Option(4, "--level", min=0),
    write_cell_tables: bool = typer.Option(False, "--write-cell-tables/--summary-only"),
) -> None:
    """Add morphology-focus protein images to the real RCC geometry analysis."""

    channel_paths, qc_paths, background_paths = discover_rcc_protein_image_inputs(
        morphology_dir,
        qc_mask_dir=qc_mask_dir,
        background_dir=background_dir,
    )
    result = prepare_rcc_protein_imaging(
        scored_cells_dir=scored_cells_dir,
        h5_path=h5_path,
        channel_paths=channel_paths,
        output_dir=output_dir,
        qc_mask_paths=qc_paths or None,
        background_paths=background_paths or None,
        level=level,
        write_cell_tables=write_cell_tables,
    )
    typer.echo(json.dumps(result, indent=2))


@app.command(name="prepare-cosmx-external-validation")
def prepare_cosmx_external_validation_command(
    input_dir: Path,
    output_dir: Path = typer.Option(Path("cosmx-external-validation"), "--output-dir"),
) -> None:
    """Validate CosMx flat-file adapters and relative RNA geometry in pixel units."""

    samples: dict[str, dict[str, Path]] = {}
    for metadata in sorted(input_dir.glob("*_metadata_file.csv.gz")):
        prefix = metadata.name.replace("_metadata_file.csv.gz", "")
        accession = prefix.split("_")[0]
        expression = input_dir / f"{prefix}_exprMat_file.csv.gz"
        polygons = input_dir / f"{prefix}_polygons.csv.gz"
        if expression.exists() and polygons.exists():
            samples[accession] = {
                "metadata": metadata,
                "expression": expression,
                "polygons": polygons,
            }
    if not samples:
        raise typer.BadParameter("No complete CosMx sample triplets were found")
    result = prepare_cosmx_external_validation(samples, output_dir)
    typer.echo(json.dumps(result, indent=2))


@app.command(name="benchmark-her2-ihc")
def benchmark_her2_ihc_command(
    input_dir: Path,
    output_dir: Path = typer.Option(Path("her2-ihc-benchmark"), "--output-dir"),
    downsample: int = typer.Option(4, "--downsample", min=1, max=32),
) -> None:
    """Benchmark relative HER2 IHC intensity against ordinal pathology scores."""

    result = benchmark_her2_ihc(input_dir, output_dir, downsample=downsample)
    typer.echo(json.dumps(result, indent=2))


@app.command(name="benchmark-shg-collagen")
def benchmark_shg_collagen_command(
    input_dir: Path,
    output_dir: Path = typer.Option(Path("shg-collagen-benchmark"), "--output-dir"),
) -> None:
    """Benchmark relative SHG collagen features without deriving diffusivity."""

    images: dict[str, list[Path]] = {}
    for tissue_dir in sorted(path for path in input_dir.iterdir() if path.is_dir()):
        paths = sorted(tissue_dir.glob("*.tif")) + sorted(tissue_dir.glob("*.tiff"))
        if paths:
            images[tissue_dir.name] = paths
    if not images:
        paths = sorted(input_dir.glob("*.tif")) + sorted(input_dir.glob("*.tiff"))
        if paths:
            images["unlabelled"] = paths
    result = benchmark_shg_collagen(images, output_dir)
    typer.echo(json.dumps(result, indent=2))


@app.command(name="audit-breast-cell-groups")
def audit_breast_cell_groups_command(
    input_dir: Path,
    output_dir: Path = typer.Option(Path("breast-cell-groups"), "--output-dir"),
) -> None:
    """Audit independent breast-section cell-group composition."""

    samples = {path.stem: path for path in sorted(input_dir.glob("*_cell_groups.csv"))}
    if not samples:
        raise typer.BadParameter("No *_cell_groups.csv files were found")
    result = audit_breast_cell_groups(samples, output_dir)
    typer.echo(json.dumps(result, indent=2))


@app.command(name="benchmark-perfusion")
def benchmark_perfusion_command(
    input_dir: Path,
    output_dir: Path = typer.Option(Path("perfusion-benchmark"), "--output-dir"),
    downsample: int = typer.Option(4, "--downsample", min=1, max=32),
) -> None:
    """Benchmark relative Hoechst perfusion against CD31 distance in RGB TIFFs."""

    images = {path.stem: path for path in sorted(input_dir.glob("Fig_5C_*.tif"))}
    images.update({path.stem: path for path in sorted(input_dir.glob("Fig_5C_*.tiff"))})
    if not images:
        raise typer.BadParameter("No Fig_5C RGB TIFFs were found")
    corrections = {label: 1.0 for label in images if "Bottom" in label}
    result = benchmark_perfusion_tiffs(
        images,
        output_dir,
        downsample=downsample,
        red_correction_by_label=corrections,
    )
    typer.echo(json.dumps(result, indent=2))


@app.command(name="benchmark-her2-receptor-calibration")
def benchmark_her2_receptor_calibration_command(
    workbook_path: Path,
    output_dir: Path = typer.Option(Path("her2-receptor-calibration"), "--output-dir"),
) -> None:
    """Fit the published source-protocol Cy5-to-HER2 receptor calibration."""

    result = benchmark_her2_receptor_calibration(workbook_path, output_dir)
    typer.echo(json.dumps(result, indent=2))


@app.command(name="build-igg-transport-prior")
def build_igg_transport_prior_command(
    output_dir: Path = typer.Option(Path("igg-transport-prior"), "--output-dir"),
) -> None:
    """Build a literature-derived IgG transport sensitivity prior."""

    result = build_igg_transport_prior(output_dir)
    typer.echo(json.dumps(result, indent=2))


@app.command(name="benchmark-bordeau-supplement")
def benchmark_bordeau_supplement_command(
    docx_path: Path,
    output_dir: Path = typer.Option(Path("bordeau-supplement"), "--output-dir"),
) -> None:
    """Curate administered-trastuzumab endpoints from the valid supplementary DOCX."""

    result = benchmark_bordeau_supplement(docx_path, output_dir)
    typer.echo(json.dumps(result, indent=2))


@app.command(name="validate-artifact")
def validate_artifact_command(
    path: Path,
    expected_kind: str = typer.Option(..., "--expected-kind"),
) -> None:
    """Validate PDF/ZIP/GZIP/TIFF content instead of trusting its extension."""

    supported = {"pdf", "zip", "gzip", "tiff", "json", "html", "xml"}
    if expected_kind not in supported:
        raise typer.BadParameter(f"expected kind must be one of {sorted(supported)}")
    result = validate_artifact(path, expected_kind)
    typer.echo(json.dumps(result, indent=2))
    if not result["valid"]:
        raise typer.Exit(code=1)


@app.command(name="validate-v071-results")
def validate_v071_results_command(
    bundle_dir: Path,
    output_dir: Path | None = typer.Option(None, "--output-dir"),
) -> None:
    """Validate the completed v0.7.1 Colab bundle and its abstention boundary."""

    validation = validate_bundle(bundle_dir)
    payload = {
        "status": validation.status,
        "checks": validation.checks,
        "issues": [{"code": issue.code, "message": issue.message} for issue in validation.issues],
        "summary": validation.summary,
    }
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        validation_path = output_dir / "validation_v0.7.1.json"
        claims_path = output_dir / "claims_v0.7.1.json"
        validation_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        if validation.ok:
            claims_path.write_text(
                json.dumps(build_claims(validation.summary), indent=2),
                encoding="utf-8",
            )
    typer.echo(json.dumps(payload, indent=2))
    if not validation.ok:
        raise typer.Exit(code=1)


@app.command(name="audit-xenium-zarr-erbb2")
def audit_xenium_zarr_erbb2_command(
    zarr_zip: Path,
    output: Path = typer.Option(Path("xenium_erbb2_audit.json"), "--output"),
) -> None:
    """Read a native Xenium Zarr ZIP and summarize its ERBB2 feature vector."""

    loaded = read_xenium_sparse_matrix(zarr_zip)
    feature_index = find_erbb2_feature_index(loaded.feature_names)
    values = extract_feature_vector(
        loaded.matrix,
        feature_index,
        len(loaded.feature_names),
        len(loaded.cell_ids),
    )
    result = {
        "version": "0.7.1",
        "source": str(zarr_zip),
        "root_prefix": loaded.root_prefix,
        "encoding": loaded.encoding,
        "array_names": list(loaded.array_names),
        "features": len(loaded.feature_names),
        "cells": len(loaded.cell_ids),
        "selected_feature_index": feature_index,
        "selected_feature_name": loaded.feature_names[feature_index],
        "ERBB2_RNA_mean_all_matrix_cells": float(np.mean(values)),
        "ERBB2_RNA_median_all_matrix_cells": float(np.median(values)),
        "ERBB2_RNA_positive_fraction_all_matrix_cells": float(np.mean(values > 0)),
        "ERBB2_RNA_total_counts": float(np.sum(values)),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    typer.echo(json.dumps(result, indent=2))


@app.command(name="synthesize-evidence")
def synthesize_evidence_command(
    repository_root: Path = typer.Option(Path("."), "--repository-root"),
    output_dir: Path = typer.Option(Path("evidence-synthesis-v0.8"), "--output-dir"),
) -> None:
    """Build the v0.8 evidence graph, readiness audit and measurement priorities."""

    result = synthesize_repository_evidence(repository_root, output_dir)
    typer.echo(json.dumps(result["summary"], indent=2))


@app.command(name="rank-relative-accessibility")
def rank_relative_accessibility_command(
    repository_root: Path = typer.Option(Path("."), "--repository-root"),
    output_dir: Path = typer.Option(Path("relative-accessibility-v0.8"), "--output-dir"),
    draws: int = typer.Option(20_000, "--draws", min=100),
    seed: int = typer.Option(17, "--seed"),
) -> None:
    """Rank RCC targets on a relative geometry-expression proxy only."""

    result = run_relative_accessibility(
        repository_root,
        output_dir,
        draws=draws,
        seed=seed,
    )
    typer.echo(json.dumps(result["summary"], indent=2))


@app.command(name="build-v080-package")
def build_v080_package_command(
    repository_root: Path = typer.Option(Path("."), "--repository-root"),
    output_dir: Path = typer.Option(Path("results/evidence_synthesis_v0.8"), "--output-dir"),
    draws: int = typer.Option(20_000, "--draws", min=100),
    seed: int = typer.Option(17, "--seed"),
) -> None:
    """Build and validate the complete v0.8 evidence-synthesis package."""

    result = build_v080_package(
        repository_root,
        output_dir,
        draws=draws,
        seed=seed,
    )
    typer.echo(json.dumps(result["summary"], indent=2))


@app.command(name="validate-v080-results")
def validate_v080_results_command(
    package_dir: Path,
) -> None:
    """Validate v0.8 evidence, ranking and abstention invariants."""

    validation = validate_v080_package(package_dir)
    payload = {
        "status": validation.status,
        "checks": validation.checks,
        "issues": [issue.__dict__ for issue in validation.issues],
        "summary": validation.summary,
    }
    typer.echo(json.dumps(payload, indent=2))
    if not validation.ok:
        raise typer.Exit(code=1)


@app.command()
def retrospective(
    table_path: Path,
    output: Path = typer.Option(Path("retrospective.json"), "--output"),
) -> None:
    """Run the preregistered evaluation-only retrospective."""

    result = run_retrospective(pd.read_csv(table_path))
    output.write_text(json.dumps(result.__dict__, indent=2), encoding="utf-8")
    typer.echo(str(output))


@app.command()
def benchmark(
    output: Path = typer.Option(Path("results/simulated/benchmark.json"), "--output"),
    quick: bool = typer.Option(False, "--quick"),
    seed: int = typer.Option(17, "--seed"),
) -> None:
    """Run the complete simulation benchmark and all required ablations."""

    result = run_benchmark(quick=quick, seed=seed)
    write_benchmark(result, output)
    typer.echo(str(output))


@app.command()
def report(
    benchmark_path: Path,
    output: Path = typer.Option(Path("benchmark_report.md"), "--output"),
) -> None:
    """Render a Markdown report from a benchmark result."""

    render_report(benchmark_path, output)
    typer.echo(str(output))


@app.command(name="validate-manifest")
def validate_manifest_command(path: Path) -> None:
    """Validate a provenance manifest."""

    manifest = validate_manifest(path)
    typer.echo(manifest.model_dump_json(indent=2))


if __name__ == "__main__":
    app()
