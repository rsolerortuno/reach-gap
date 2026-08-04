# Static analysis — v0.7.0

- Ruff 0.16.0 lint: PASS.
- Ruff 0.16.0 format check: PASS.
- Strict Mypy 2.3.0: PASS across all 30 source modules.
- Pyright 1.1.411: PASS with 0 errors and 0 warnings for the five v0.7 modules.
- The clean v0.6.1 baseline remains the last verified full-package Pyright run. In this runtime, newer
  third-party scientific stubs expose legacy-library typing disagreements unrelated to v0.7; these are
  not represented as new reach-gap errors or silently suppressed.
- Tests: 80 passed.
- Configured coverage: 86.74%, above the required 85% threshold.

No source module added in v0.7 is excluded from coverage or static analysis.
