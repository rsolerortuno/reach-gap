# reach-gap v0.8.0 model card

## Summary

`reach-gap` is research software for combining spatial tumour measurements with a simplified antibody-transport model and an explicit evidence graph. Version 0.8.0 supports mechanistic simulation, real RCC spatial preparation, external component validation, evidence-readiness auditing and relative target ranking.

## Intended use

- exploratory antibody-target accessibility research;
- spatial geometry and target-localisation analysis;
- sensitivity analysis across vessel definitions and uncertain objectives;
- evidence-gap auditing and experimental prioritisation;
- reproducible portfolio and methods demonstrations.

## Not intended for

- patient treatment decisions;
- clinical response prediction;
- dosing recommendations;
- claims of absolute tumour penetration without matched measurements;
- ranking targets as universally superior across diseases or cohorts.

## Evidence in v0.8.0

- 465,534-cell RCC spatial dataset with RNA, protein and geometry outputs;
- pathology-defined tumour regions and structural vessel proxies;
- six vessel definitions used in sensitivity analysis;
- four-field independent Hoechst–CD31 validation;
- source-protocol HER2 quantitative calibration;
- literature-curated IgG transport priors;
- breast Xenium ERBB2 RNA control across 679,197 cells;
- negative-result retention for the Bordeau representative-image reanalysis.

## Main outputs

- evidence graph and readiness score;
- unresolved uncertainty budget and measurement priorities;
- relative target rank distributions;
- pairwise win probabilities;
- leave-one-component-out robustness;
- explicit `NOT_COMPUTED` absolute outputs.

## Main result

VISTA ranked first in 98.525% of 20,000 relative proxy draws across six structural vessel definitions. This reflects same-section protein prevalence and geometry under the declared objective. It is not an efficacy claim.

## Validation

- 120 tests passed;
- 88.30% coverage;
- strict Mypy passed on 37 source modules;
- Ruff passed;
- Pyright returned 0 errors on the v0.8 modules;
- 67 result invariants passed with 0 issues;
- wheel and source distribution were tested after clean installation/extraction.

## Known limitations

The absolute RCC outputs remain unidentified because the same tissue lacks functional perfusion, calibrated surface-antigen copies, administered-antibody concentration or engagement, RCC-specific transport and a pharmacological endpoint. The real analysis currently represents one section rather than a cohort.

## Ethical and safety considerations

The software must not be used to recommend patient treatment. Relative target scores can be misinterpreted if presented without their evidence scope. Every downstream use should preserve the release labels and abstentions.
