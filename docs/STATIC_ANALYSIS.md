# Static-analysis report — v0.6.1

Version 0.6.1 removes the disclosed legacy typing baseline rather than hiding it.

## Verified gates

- `ruff check .`: PASS.
- `ruff format --check .`: PASS.
- `mypy src/reach_gap` under strict project configuration: PASS across 25 source modules.
- Pyright 1.1.411 on `src/reach_gap`: 0 errors, 0 warnings.
- Full tests: 68 passed.
- Configured coverage: 86.49% in the clean Python 3.13 verification run, threshold 85%.
- `python -m compileall src tests`: PASS.

## What changed

The release replaces ambiguous NumPy aliases with explicit `TypeAlias` declarations, narrows HDF5
groups and datasets, types array-like public inputs, normalises SciPy union returns, removes untyped
writer lambdas, and validates scalar conversions at I/O boundaries. Third-party `tifffile` source is
skipped by a narrowly scoped Mypy override because the installed release contains Python 3.12 type
syntax while the package still targets Python 3.11+. Calls into untyped scikit-image functions carry
three line-level `no-untyped-call` suppressions; no reach-gap module is ignored.

## Interpretation

This is an engineering and maintainability improvement. It does not validate perfusion, antigen copy
number, tissue diffusion, antibody concentration or therapeutic efficacy.
