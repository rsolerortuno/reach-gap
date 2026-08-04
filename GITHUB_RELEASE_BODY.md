# reach-gap v0.8.0 — evidence-aware relative antibody accessibility

This release adds an auditable evidence graph and uncertainty-aware relative RCC target ranking.

## Highlights

- 465,534-cell RCC spatial analysis.
- Six structural-vessel definitions.
- 20,000 uncertainty draws.
- VISTA ranked first in 98.525% of draws under the declared relative objective.
- Evidence readiness scored 40.5/100; 3 of 8 absolute-reachability requirements are supported by same-tissue data.
- Missing measurements are ranked by expected uncertainty reduction.
- Absolute `reachable_fraction`, `penetration_depth`, `expression_reach_gap` and pharmacological concordance remain explicitly `NOT_COMPUTED`.

## Quality gates

- 120 tests passed.
- 88.30% coverage.
- Ruff PASS.
- Strict Mypy PASS across 37 modules.
- Pyright v0.8: 0 errors, 0 warnings.
- 67 result invariants, 0 issues.
- Clean wheel and sdist validation PASS.

## Assets

Download the wheel for installation or the source distribution for a reproducible source build. Verify downloads with `SHA256SUMS`.

> Research software only. Relative target ranking is not therapeutic efficacy or patient guidance.

## Automated publication gate

The release is created only after the tag matches the package version and the workflow repeats compilation, Ruff, strict Mypy, all 120 tests, coverage enforcement, scientific result validation, figure regeneration, package inspection and clean-wheel smoke testing.
