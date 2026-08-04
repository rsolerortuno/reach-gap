# reach-gap v0.8.0

## Evidence-aware relative accessibility

Version 0.8.0 adds two major capabilities:

1. a machine-readable evidence graph that separates same-tissue measurements, external validations, priors, missing inputs and blocked outputs;
2. an uncertainty-aware relative RCC target ranking across six vessel definitions and 20,000 objective-weight draws.

## Headline results

- evidence readiness: 40.5/100 under the fixed eight-layer audit rubric;
- same-tissue requirements satisfied: 3 of 8;
- stable relative top target: VISTA;
- VISTA rank-1 probability: 98.525%;
- 465,534 RCC cells analysed;
- 679,197 breast Xenium cells included in the independent ERBB2 control;
- absolute RCC reachability remains `NOT_COMPUTED`.

## Validation

- 120 tests passed;
- 88.30% coverage;
- Ruff and strict Mypy passed;
- Pyright reported 0 errors on v0.8 modules;
- 67 result invariants passed with 0 issues;
- clean wheel installation and clean sdist test passed.

## Release assets

- `reach_gap-0.8.0-py3-none-any.whl`
- `reach_gap-0.8.0.tar.gz`
- `SHA256SUMS`
- `build_validation_v0.8.json`
- `artifact_manifest_v0.8.json`

## Scientific scope

The v0.8 ranking is a relative same-section geometry-expression proxy. It is not an estimate of therapeutic efficacy or absolute antibody penetration.
