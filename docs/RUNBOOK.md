# Runbook

## Install and verify

```bash
uv sync --extra dev
uv run ruff check .
uv run ruff format --check .
uv run mypy --strict src/reach_gap
uv run pytest
uv run reach benchmark --quick --output /tmp/benchmark.json
uv run reach validate-manifest results/manifest.json
```

## Committed artefacts

Run these as separate commands:

```bash
PYTHONPATH=src python -m reach_gap.cli benchmark \
  --output results/simulated/benchmark.json --seed 17
PYTHONPATH=src python -m reach_gap.cli report \
  results/simulated/benchmark.json \
  --output results/simulated/benchmark_report.md
PYTHONPATH=src python scripts/generate_reference_run.py
PYTHONPATH=src python scripts/write_execution_manifest.py
```

`results/simulated/execution.json` records versions, settings and SHA-256 hashes.

## Real spatial input

Export a cell table with columns
`cell_id,x_um,y_um,is_tumour,target_signal,vessel_signal,ecm_score,caf_score`.
For Visium spots, use the spot identifier as `cell_id` and retain the spot-level resolution in downstream
claims. For Xenium, CosMx or MERFISH, provide segmented-cell centroids. The adapter deliberately does not
infer markers from gene names and does not convert RNA counts to surface density.

Create a JSON manifest matching `ProvenanceManifest`, including source identifiers, licence,
segmentation version, vessel definition and, where available,
`antigen_calibration_nM_per_signal`. If calibration is absent, ingestion succeeds but indexing returns
`INSUFFICIENT_EVIDENCE`.

```bash
reach ingest cells.csv manifest.json --output-dir run/
reach solve run/features.npz --output run/solution.npz
reach index run/features.npz --output-dir run/index/
```

## Platform notes

- Xenium/CosMx/MERFISH: use cell centroids and explicit validated marker channels where available.
- Visium: spot-level output is not cell-level engagement.
- Vessel detection must be supplied or independently validated; an empty mask forces abstention.
- Endothelial presence is not perfusion and must be recorded separately.

## Retrospective curation

Populate `examples/retrospective_template.csv` from public primary sources and registries. Record source
identifiers for every field. Do not alter model parameters after inspecting outcomes. The bundled
synthetic rows are schema tests and are excluded from scientific analysis by default.

## Public Xenium RCC workflow

The Colab notebook is the recommended entry point:

```text
notebooks/prepare_reach_gap_xenium_rcc.ipynb
```

For a local or cloud machine with the four provider files in one directory:

```bash
uv sync --extra xenium
uv run reach prepare-xenium-rcc /path/to/raw-folder \
  --output-dir /path/to/raw-folder/reach-gap-analysis \
  --verify-large-md5
```

The command is resumable by default. Use `--force` only to discard reusable stage outputs. The first MD5
pass over the 36.1 GB bundle is I/O-bound. Subsequent runs reuse `download_verification.json` when file
sizes remain unchanged.

Prepared cell tables are partitioned under `tables/` and `targets/`. The absolute index remains
`NOT_COMPUTED` because this dataset has no conversion from scaled protein fluorescence to surface-antigen
density. See `docs/XENIUM_RCC.md`.

## Real RCC H&E/pathology run

```bash
uv sync --extra dev --extra xenium
uv run reach prepare-he-rcc \
  /path/to/Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_image.ome.tif \
  /path/to/Xenium_V1_Human_Kidney_FFPE_Protein_updated_annotation.geojson \
  --alignment-path /path/to/Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_imagealignment.csv \
  --output-dir results/real_rcc_he \
  --analysis-level 4
```

The command reads one OME-TIFF pyramid level and does not load the native full-resolution image. The
output must abstain from target-specific indexing.

## Reduce the Xenium bundle without high RAM

Run this beside the original ZIP, preferably on the local machine that already stores it:

```bash
uv run reach extract-xenium-essentials \
  Xenium_V1_Human_Kidney_FFPE_Protein_updated_outs.zip \
  --output-dir xenium-essential-package \
  --part-size-mb 95
```

The extractor streams data with an 8 MB buffer. Large members are split while being decompressed; no
full member is loaded into RAM and no second unsplit copy is required. Upload the resulting
`essential_package_manifest.json`, split manifests and part files to Drive.

## Compact RCC Xenium package — version 0.4

Required files:

```text
cells.parquet
cell_feature_matrix.h5
gene_panel.json
protein_panel.json
metrics_summary.csv
experiment.xenium
cell_boundaries.parquet          # optional but required for committed segmentation QC
nucleus_boundaries.parquet       # optional but required for committed segmentation QC
```

Run:

```bash
uv sync --extra xenium --extra dev
uv run reach prepare-xenium-essential-rcc /data/xenium-essential-package \
  --annotation-path /data/pathology_annotations.geojson \
  --alignment-path /data/he_imagealignment.csv \
  --output-dir results/real_rcc_xenium
```

Add `--write-full-cell-tables` only when the complete per-cell derivative package is required. The default
keeps a deterministic 20,000-cell sample and all aggregate results.

Expected terminal status:

```text
REAL_XENIUM_MOLECULAR_GEOMETRY_PREPARED_ABSOLUTE_INDEX_NOT_COMPUTED
```

Validation:

```bash
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv run mypy --strict src/reach_gap
python -m compileall -q src scripts tests
```

The real-data run must not create `reachable_fraction` or `expression_reach_gap` unless the manifest
contains both a surface-antigen calibration and perfusion evidence.

## RCC morphology-focus protein imaging — version 0.5

Required exact files in `morphology_focus/`:

```text
ch0028_cd31.ome.tif
ch0032_alphasma.ome.tif
ch0031_vimentin.ome.tif
ch0030_panck.ome.tif
ch0021_pd-l1.ome.tif
ch0020_vista.ome.tif
```

Run after the compact RCC cell workflow has produced its scored cell-table parts:

```bash
uv sync --extra xenium --extra dev
uv run reach prepare-rcc-protein-imaging \
  /data/v0.4/cell_tables \
  /data/cell_feature_matrix.h5 \
  /data/morphology_focus \
  --qc-mask-dir /data/aux_outputs/morphology_focus_qc_masks \
  --background-dir /data/aux_outputs/background_qc_images \
  --output-dir results/real_rcc_imaging \
  --level 4
```

Use `--write-cell-tables` to export all per-cell image-derived measurements. The default commits only
aggregate summaries and figures. Expected status:

```text
REAL_PROTEIN_IMAGE_GEOMETRY_PREPARED_ABSOLUTE_INDEX_NOT_COMPUTED
```

## Version 0.7 external component benchmarks

```bash
reach benchmark-perfusion /path/to/S-BIAD3159_LLC_CD31_Hoechst \
  --output-dir results/external_validation/perfusion_s_biad3159

reach benchmark-her2-receptor-calibration \
  /path/to/41598_2026_42898_MOESM1_ESM.xlsx \
  --output-dir results/external_validation/her2_receptor_calibration_mckinski2026

reach build-igg-transport-prior \
  --output-dir results/external_validation/igg_transport_prior

reach benchmark-bordeau-supplement \
  /path/to/NIHMS1684776-supplement-1.docx \
  --output-dir results/external_validation/bordeau_2021_supplement

reach audit-breast-cell-groups /path/to/10x_breast_12sample_Xenium \
  --output-dir results/external_validation/breast_xenium_cell_groups
```

The perfusion command subtracts red from blue in labels containing `Bottom` because the source CD8
composite is magenta. Review the source-channel legend before applying this rule to another dataset.

### Version 0.7.1 completion and audit

Validate the compact Colab completion bundle and emit locked claims:

```bash
reach validate-v071-results /path/to/v0.7.1_colab_completion \
  --output-dir results/external_validation
```

Audit one native provider Zarr matrix directly:

```bash
reach audit-xenium-zarr-erbb2 /path/to/cell_feature_matrix.zarr.zip \
  --output /tmp/erbb2-zarr-summary.json
```

The Zarr audit accepts the provider `cell_features` feature-by-cell CSR layout and a complete CSC
fallback. Ambiguous arrays, duplicate packed cell IDs, dimension mismatches and non-unique ERBB2 feature
matches fail explicitly.
