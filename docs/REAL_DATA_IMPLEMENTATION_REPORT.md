# Real-data implementation report — v0.4.0

## What changed

Version 0.4 replaces the previous “molecular execution pending” state with a complete run over the public
RCC Xenium compact package. The workflow no longer requires the 36.1 GB archive once the essential files
have been extracted.

Implemented components:

- Exact barcode identity validation across 465,534 cells.
- Sparse HDF5 extraction of 51 selected features from a 543 × 465,534 matrix.
- Distinct RNA and protein channels with protein preferred for target positivity.
- Pathology-coordinate transform selection with an acceptance margin.
- Endothelial/pericyte, CAF, immune and malignant-proxy scores.
- Cell-centroid distance and local density relative to marker-derived vessel proxies.
- PD-L1, VISTA, PD-1 and LAG-3 relative target summaries.
- RNA–protein concordance analysis.
- Six-definition vascular sensitivity analysis.
- Streaming cell/nucleus polygon reconstruction from 23.2 million boundary vertices.
- Compact Parquet fallback for the narrow flat/ZSTD/PLAIN Xenium schema.
- A new `prepare-xenium-essential-rcc` CLI command.
- Separate compact source repository and full per-cell analysis package.

## What broke and how it was fixed

### Arrow was unavailable

The execution environment could not install `pyarrow`. Three alternatives were attempted: pandas' Arrow
engine, package installation and direct wheel retrieval. None was available. A narrow `parquet_lite`
reader was therefore implemented for the exact provider schema. It rejects dictionary encoding, nested
schemas, nulls and unsupported codecs rather than silently guessing.

The boundary files contained up to 16 all-zero bytes of page padding. Initial strict decoding rejected
this. The reader was changed to accept only short all-zero padding; non-zero trailing bytes remain an
error. The complete 11.6 million cell-boundary and 11.6 million nucleus-boundary rows then decoded.

### Initial boundary aggregation was too slow

A Python loop over approximately 925,000 polygons exceeded the command window. It was replaced with
vectorised shoelace and perimeter reductions per row group. Both boundary tables then completed in
11.40 seconds with a peak resident set of approximately 555 MiB.

### The first all-cell command exceeded the orchestration timeout

The expensive scoring stages completed and wrote validated checkpoints before the command was stopped.
Final tables were resumed from five non-overlapping, duplicate-free cell partitions. The end-to-end wall
clock is therefore recorded as `NOT_COMPUTED_RESUMED_FROM_CHECKPOINTS`; no runtime was reconstructed or
invented.

## Scientific decisions

- Protein intensity is used for relative target positivity where measured; RNA remains a separate
  comparison, never a hidden substitute.
- No ECM score is imputed because the measured panel did not resolve a direct ECM marker set.
- Vascular geometry is reported across six definitions because the choice materially changes distances.
- Pathology polygons are used only after the transform passes a quantitative support and separation rule.
- Absolute mechanistic indexing abstains because antigen calibration and perfusion are absent.

## Validation

- 48 tests passed.
- 87.93% coverage over non-I/O modules.
- Complete cell/HDF5 identifier equality passed.
- Cell boundary coverage: 100%.
- Nucleus boundary coverage: 98.66%.
- Local source audit and `compileall` passed.
- Ruff and strict mypy were not executable in this sandbox because the required packages were absent;
  both remain configured in GitHub Actions and are not reported as locally passed.

## Highest-risk remaining assumptions

1. Marker-derived endothelial cells trace a therapeutically relevant vessel network.
2. A 2D centroid-to-cell proxy is informative about 3D vessel-wall distance.
3. Relative fluorescence thresholds preserve enough biology to compare spatial target patterns, despite
   lacking receptor-density calibration.
