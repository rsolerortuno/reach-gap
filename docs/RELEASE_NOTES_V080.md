# Release notes — reach-gap v0.8.0

## Added

- Machine-readable evidence graph separating same-tissue measurements, external validations, assay
  calibrations, literature priors, missing requirements and blocked outputs.
- Fixed evidence-readiness and unresolved-uncertainty audit.
- Ranked measurement priorities for completing absolute RCC validation.
- Relative RCC geometry-expression target analysis across six vessel definitions and 20,000 objective
  uncertainty draws.
- Pairwise target win probabilities and leave-one-component-out stability analysis.
- Four new CLI commands and compact SVG figures.
- Integrated v0.8 manifest and 70-invariant validator.
- 120-test validation suite with 88.30% configured line coverage.

## Headline results

- Absolute evidence-readiness score: **40.5/100**, explicitly not a biological reachability score.
- Three of eight same-tissue requirements fully satisfied.
- Largest unresolved weighted contribution: administered-antibody field (**25.21%**).
- VISTA ranks first in **98.525%** of relative geometry-expression draws.
- VISTA remains first in more than **91%** of every leave-one-component-out analysis.

## Preserved abstentions

Absolute RCC `reachable_fraction`, `penetration_depth`, `expression_reach_gap` and pharmacological
concordance remain `NOT_COMPUTED`.
