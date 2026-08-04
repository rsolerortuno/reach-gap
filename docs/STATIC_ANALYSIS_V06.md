# Static-analysis policy — v0.6.0

Version 0.6 separates new-code quality gates from disclosed legacy typing debt.

## Passing gates

- `ruff check .`: PASS.
- `ruff format --check .`: PASS.
- Strict Mypy on `cosmx.py`, `her2_ihc.py`, `shg_collagen.py` and `artifact_validation.py`: PASS.
- Pyright 1.1.411 on the same four modules: 0 errors, 0 warnings.
- Full tests: 68 passed.
- Configured coverage: 86.83%, threshold 85%.

## Legacy Mypy baseline

Strict Mypy across the complete package still reports **224 errors in 17 files** when targeted to Python
3.11 with external imports skipped. The current Python 3.13/stub environment reports 254 errors, showing
that absolute counts depend on interpreter and stub resolution.

The CI policy now:

1. requires strict Mypy and Pyright to pass for every module added in v0.6;
2. runs the complete legacy package and fails if the Python 3.11 error count increases above 224;
3. keeps the full log in `results/mypy_legacy_baseline_v0.6.txt`.

This prevents new debt while avoiding a false claim that the legacy package is fully typed. Reducing the
baseline remains desirable, but it is not represented as scientific validation.
