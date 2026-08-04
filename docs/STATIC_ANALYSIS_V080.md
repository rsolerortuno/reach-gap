# Static analysis and build validation — v0.8.0

## Full package

- Tests: **120 passed**, 0 failed.
- Configured line coverage: **88.30%**, above the required 85% threshold.
- Ruff 0.16.0 lint: **PASS**.
- Ruff 0.16.0 format check: **PASS** across 116 files.
- Strict Mypy 2.3.0: **PASS** across all 37 source modules.
- Pyright 1.1.411 on the four v0.8 modules: **0 errors, 0 warnings**.
- Full-package Pyright: **0 errors, 1 warning**.
- Python byte compilation and `tabnanny`: **PASS**.
- Local source audit and quick simulation benchmark: **PASS**.
- v0.8 compact-result validation: **70 checks, 0 issues**.

## Full-package Pyright warning

The sole warning is `reportMissingModuleSource` for the optional `pyarrow` import in the unchanged
`xenium.py` module. Pyright found the import interface but could not resolve its implementation source in
this runtime. It is not a type error, does not affect the four v0.8 modules, and no warning is suppressed.

## Clean distributions

The wheel is installed into an isolated environment and its version, CLI and v0.8 validator are run.
The source distribution is extracted independently and the complete 120-test suite is repeated from that
clean tree.

## Scientific status

`EVIDENCE_SYNTHESIS_AND_RELATIVE_RCC_TARGET_RANKING_COMPLETE_ABSOLUTE_REACHABILITY_NOT_COMPUTED`

The evidence-readiness score and Monte Carlo rank frequencies are audit quantities, not posterior
probabilities of therapeutic reachability. Absolute RCC `reachable_fraction`, `penetration_depth`,
`expression_reach_gap` and pharmacological concordance remain `NOT_COMPUTED`.
