# reach-gap v0.8.0 — final release report

## Release status

**COMPLETE WITH EXPLICIT ABSOLUTE-REACHABILITY ABSTENTION**

Version 0.8.0 converts the v0.7.1 evidence layers into a machine-readable dependency graph and adds a
same-section RCC relative geometry-expression ranking. It does not present either output as a posterior
probability of antibody delivery or clinical success.

## Scientific results

- Evidence-readiness audit: **40.5/100**.
- Same-tissue requirements fully satisfied: **3/8**.
- Dominant unresolved requirement: **administered-antibody field**
  (**25.21%** of the unresolved weighted burden).
- RCC relative target analysis: **20,000 draws** across
  **6 structural-vessel definitions**.
- Stable top target in this proxy: **VISTA**, rank 1 in
  **98.525%** of draws.
- VISTA pairwise win frequencies: **99.105%**
  versus PD-1, **99.420%** versus PD-L1 and
  **100.000%** versus LAG-3.
- VISTA remains rank 1 in at least
  **91.175%**
  of every leave-one-component-out analysis.

The ranking is conditional on one RCC section, four measured targets, within-section positivity definitions
and structural rather than functionally perfused vessels. It is not a clinical target-priority claim.

## Experimental-design result

The next measurements ranked by unresolved weighted contribution are:

1. administered-antibody concentration or target-engagement field in the RCC tissue;
2. functional perfusion in the same section;
3. surface-antigen calibration in molecules per cell;
4. RCC matrix-transport measurement; and
5. a same-tissue pharmacological endpoint.

This is the most useful new capability in v0.8: the tool identifies which measurement would reduce the
largest identifiable evidence gap instead of collapsing incompatible evidence into a confidence number.

## Verification

- **120 tests passed** from the working tree and again from a clean sdist extraction.
- **88.30% line coverage**, threshold 85%.
- Ruff 0.16.0 lint and format: **PASS**.
- Strict Mypy 2.3.0: **PASS across 37 modules**.
- Pyright 1.1.411 on v0.8 modules: **0 errors, 0 warnings**.
- Full-package Pyright: **0 errors, 1 optional-dependency source warning**.
- Compact-result validator: **67 checks, 0 issues**.
- Wheel clean installation, CLI, quick benchmark, `compileall`, `tabnanny` and local source audit: **PASS**.

## Preserved abstentions

The following RCC outputs remain `NOT_COMPUTED`:

- absolute `reachable_fraction`;
- absolute `penetration_depth`;
- absolute `expression_reach_gap`; and
- model-versus-administered-antibody pharmacological concordance.

## Critical assessment

Overall assessment: **8.8/10**. The project is strong in transparent scientific software, evidence
separation, reproducibility and experimental-design reasoning. A higher score requires a matched dataset
containing administered antibody, functional perfusion, calibrated surface antigen and a pharmacological
endpoint in the same tissue, rather than additional software-only complexity.
