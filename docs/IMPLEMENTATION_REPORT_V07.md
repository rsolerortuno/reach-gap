# Implementation report — v0.7.0

## Added modules

- `perfusion.py`: RGB TIFF metadata, CD31 segmentation, Hoechst-to-CD31 distance profiles, red-channel
  contamination sensitivity and result plots.
- `her2_calibration.py`: dependency-light XLSX reader, uncensored replicate extraction, source-protocol
  log-log calibration, censor-aware prediction boundaries and score summaries.
- `transport_priors.py`: unit-checked FRAP observation registry and conservative IgG sensitivity prior.
- `bordeau_validation.py`: DOCX media extraction, 3×3 panel parsing, published endpoint curation and
  representative-panel negative-result audit.
- `breast_cohort.py`: uniqueness and composition audit for independent Xenium breast cell groups.

## New CLI commands

```text
reach benchmark-perfusion
reach benchmark-her2-receptor-calibration
reach build-igg-transport-prior
reach benchmark-bordeau-supplement
reach audit-breast-cell-groups
```

## Real-data outputs

The release includes compact results under `results/external_validation/` for:

- two S-BIAD3159 perfusion TIFFs;
- the McKinski quantitative HER2 workbook;
- Netti/Davies/Ramanujan transport priors;
- the valid Bordeau supplementary DOCX;
- two Xenium breast cell-group tables; and
- a content-signature audit of eight source artifacts.

## Abstention changes

Version 0.7 narrows several unknowns without erasing them:

- perfusion extraction is independently validated, but not measured in RCC;
- receptor-copy calibration exists for another assay, but is not transferred to Xenium;
- a defensible IgG diffusion prior exists, but is not an RCC measurement;
- administered-antibody endpoints are valid and curated, but raw data/model concordance are absent.

The RCC absolute solver is therefore still gated.
