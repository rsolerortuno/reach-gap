# External validation v0.7

## Scope

Version 0.7 adds four evidence layers that were absent or invalid in v0.6:

1. an independent **in-vivo perfusion-proxy** image benchmark using Hoechst 33342 and CD31;
2. a published **quantitative HER2 receptor calibration** workbook;
3. a literature-derived **IgG tumour diffusion prior** from FRAP measurements; and
4. a valid administered-trastuzumab supplementary DOCX with published penetration endpoints.

These layers constrain uncertainty and validate components. They do not jointly identify absolute
antibody reachability in the RCC Xenium section.

## S-BIAD3159 perfusion-proxy benchmark

The full Drive package contains four paired CZI/TIFF views. The release run processed the two TIFFs that
could be transferred through the 100 MB connector limit: `Fig_5C_Upper_Left.tif` and
`Fig_5C_Bottom_Left.tif`. Both are 5,568 × 5,568 RGB ImageJ TIFFs. The TIFF X-resolution corresponds to
approximately 0.1036 µm/pixel; the analysis used a deterministic fourfold stride.

The green channel was segmented as CD31 after robust q99.5 scaling, Otsu thresholding and component
filtering. The blue channel was treated as the Hoechst perfusion proxy. Because CD8 is displayed in
magenta in the bottom panel, that panel was also analysed after subtracting the red channel from blue.

Across the two transferred views:

- distance–Hoechst correlations were negative in both images;
- the median Spearman correlation was **-0.160**;
- mean relative Hoechst within 10 µm of CD31 was a median **2.38-fold** greater than at 50–100 µm.

This validates a relative perfusion-gradient extraction method in independent LLC mouse tumours. It does
not prove which RCC CD31 structures were perfused and does not measure antibody delivery.

## Quantitative HER2 source-protocol calibration

The McKinski 2026 workbook contains 50 tumour rows. Sixty-five uncensored replicate pairs linked raw Cy5
MFI to HER2 receptors per cell. The empirical log-log relation was:

```text
log10(receptors/cell) = 0.786114 × log10(raw Cy5 MFI) + 3.775379
```

The raw-scale R² of the log-log prediction was **0.9906** and rank correlation was **1.000**. The source
assay reported an LLOQ of 10,375 and ULOQ of 178,649 receptors per cell.

The fitted relation is deliberately labelled
`SOURCE_PROTOCOL_HER2_RECEPTOR_CALIBRATION_NOT_XENIUM_TRANSFER`. No shared calibrator was measured in
both the source Cy5 assay and Xenium, so it cannot convert the RCC protein channels to receptor copies.

## IgG transport prior

Four tumour FRAP measurements from Netti et al. yield central IgG diffusion values of **8.7–19.7 µm²/s**
and a broad reported-interval envelope of **5.4–31.2 µm²/s**. Davies et al. reported tumour/free
diffusion ratios of approximately **0.30–0.50** across rhabdomyosarcoma clones. Ramanujan et al. provides
supporting collagen-gel and xenograft comparisons.

Version 0.7 exposes 5.4–31.2 µm²/s as a log-uniform **sensitivity prior only**. It is not an RCC-specific
transport coefficient, and no universal mapping from SHG brightness to diffusion is emitted.

## Administered-trastuzumab reference

The valid 1.44 MB Bordeau supplementary DOCX contains representative SKOV3 sections after 2 mg/kg
trastuzumab, with red CD31 and green trastuzumab, plus the competitive-inhibitor condition. The published
aggregate threshold-positive penetration endpoints are curated as **41.30 ± 6.70 µm** for trastuzumab
alone and **58.24 ± 5.40 µm** with 1HE at 24 h.

A transparent analysis of the three compressed representative panels per group did **not** reproduce the
published direction using a scale-free mean antibody-distance metric: group means were 0.0340 and
0.0279 of the panel diagonal, exact two-sided permutation p = 0.60. This negative result is retained.
The embedded figures may be rescaled and are not the raw animal-level data used for the publication.
Consequently, model concordance remains `NOT_COMPUTED`.

## Independent breast cohort audit

The two newly supplied Xenium cell-group files contain **679,197 unique cells**:

- S2-Middlle HER2-2+: 345,556 cells and 11 provider groups;
- S2-Bottom HER2-3+: 333,641 cells and 15 provider groups.

The output bundles and H&E images are present in Drive, but were too large for direct transfer into this
runtime. The committed result is therefore a cohort/cell-label audit, not spatial expression analysis.

## Aggregate status

```text
EXTERNAL_PERFUSION_AND_CALIBRATION_PRIORS_VALIDATED_
MODEL_PHARMACOLOGICAL_CONCORDANCE_NOT_COMPUTED
```

Absolute real-tissue `reachable_fraction`, `penetration_depth` and `expression_reach_gap` remain
`NOT_COMPUTED`.
