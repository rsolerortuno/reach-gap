# Static analysis and build validation — v0.7.1

## Full merged package

- Tests: **100 passed**, 0 failed.
- Configured coverage: **87.12%**, above the required 85% threshold.
- Ruff 0.16.0 lint: **PASS**.
- Ruff 0.16.0 format check: **PASS** across 101 files.
- Strict Mypy 2.3.0: **PASS** across all 33 source modules.
- Pyright 1.1.411 on the four modified v0.7.1 files, including CLI wiring: **0 errors, 0 warnings**.
- Python byte compilation: **PASS**.
- Result-bundle validation: **69 checks, 0 issues**.

## Full-package Pyright scope

A full 33-module Pyright run was attempted and reported 49 errors, all in unchanged legacy modules and
arising from third-party NumPy, SciPy and pandas stub disagreements already disclosed for v0.7.0. The four
v0.7.1-modified files are clean. The last verified all-package Pyright baseline remains v0.6.1; no Pyright
error is suppressed or represented as a new reach-gap defect.

## Scientific status

`EXTERNAL_PERFUSION_ALL_FOUR_FIELDS_AND_BREAST_XENIUM_ERBB2_RNA_VALIDATED_MODEL_PHARMACOLOGICAL_CONCORDANCE_NOT_COMPUTED`

Absolute RCC outputs remain `NOT_COMPUTED`.
