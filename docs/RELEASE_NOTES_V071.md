# reach-gap v0.7.1 release notes

- Completed S-BIAD3159 perfusion-proxy analysis on 4/4 fields and 12 locked sensitivity runs.
- Added native 10x Xenium `cell_features` Zarr support with feature-by-cell CSR, CSC fallback and packed
  cell IDs.
- Extracted ERBB2 RNA from 679,197 labelled cells across independent HER2-2+ and HER2-3+ sections.
- Recorded a descriptive 20.623-fold tumour-group mean difference without treating cells as independent
  biological replicates.
- Added `validate-v071-results` and `audit-xenium-zarr-erbb2` CLI commands.
- Added hash-chain, cross-file invariant, schema-ambiguity and tamper-detection regression tests.
- Preserved all absolute RCC and pharmacological-concordance abstentions.
- Merged into the full v0.7.0 source tree and rebuilt as installable v0.7.1 wheel and source distribution.
