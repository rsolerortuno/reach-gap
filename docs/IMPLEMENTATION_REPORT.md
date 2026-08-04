# Implementation report — v0.6.1

Version 0.6.1 is a typing-hardening release built on the unchanged v0.6 external-validation results.

## Delivered

- Full-package strict Mypy pass across 25 source modules.
- Full-package Pyright 1.1.411 pass with 0 errors and 0 warnings.
- Explicit NumPy array aliases and array-like interfaces.
- Typed HDF5 group/dataset narrowing.
- Typed SciPy return normalisation and I/O scalar conversion.
- CI upgraded from four-module gates plus a 224-error baseline to full-package gates.
- 68 tests remain passing; configured coverage remains above 85%.

## Scientific status

No real absolute reachability output was added. The status remains:

`EXTERNAL_ADAPTERS_AND_RELATIVE_FEATURES_VALIDATED_PHARMACOLOGICAL_CONCORDANCE_NOT_COMPUTED`

The required evidence for an absolute index remains measured administered-antibody distribution,
functional perfusion, calibrated surface antigen capacity and tissue-specific transport calibration.
