# Implementation report — v0.7.1

## Integrated modules

- `xenium_zarr.py`: native `cell_features` schema detection, unambiguous array resolution, packed cell-ID
  conversion, feature-by-cell CSR reconstruction, CSC fallback, feature lookup and extraction.
- `xenium_zarr_io.py`: bounded provider-Zarr I/O adapter kept separate from schema and scientific logic.
- `v071_results.py`: manifest/hash validation, cross-file scientific invariants, summary generation and
  machine-readable claim boundaries.

## CLI additions

```text
reach validate-v071-results
reach audit-xenium-zarr-erbb2
```

The first command validates the compact Colab bundle and emits claims. The second reads a native Xenium
Zarr ZIP, resolves ERBB2 exactly once and writes a descriptive extraction summary.

## Regression coverage

The tests reproduce the exact ambiguity that caused the failed Colab run (`data` matched both direct CSR
and CSC paths), packed-ID conversion, root-prefix handling, CSR extraction, CSC fallback, feature
ambiguity, compact-result invariants, hash reproduction and tamper detection. CLI wiring is tested without
requiring the multi-gigabyte source bundles.

## Full-package integration

The patch was merged into the verified v0.7.0 source distribution. Version metadata, optional Xenium
requirements, README, methods, runbook, decisions, limitations, tests and release reports were updated in
one source tree. The complete original suite and the new tests were then rerun before rebuilding the wheel
and source distribution.
