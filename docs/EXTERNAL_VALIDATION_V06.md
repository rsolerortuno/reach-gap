# External validation report — v0.6.0

## Scope

Version 0.6 tests whether the software can preserve geometry, modality semantics and abstention across
independent public data formats. It does **not** validate therapeutic antibody penetration. The release
adds four CosMx samples, an 11-case HER2 IHC benchmark, a 12-image SHG collagen pilot and content-signature
validation for downloaded artifacts.

The aggregate status is:

```text
EXTERNAL_ADAPTERS_AND_RELATIVE_FEATURES_VALIDATED_PHARMACOLOGICAL_CONCORDANCE_NOT_COMPUTED
```

## CosMx cross-platform validation

Four GSE299786 samples were processed with explicit `(fov, cell_ID)` keys. The workflow reconstructed
polygon area, perimeter and centroid, joined selected RNA features, called three RNA-endothelial proxy
definitions and measured target-to-proxy distance **within each FOV**. Distances remain in pixels because
no physical scale was supplied in the flat files used for this run.

| Sample | Metadata cells | Polygon/expression match | Median area error | q95 area error | Median centroid error |
|---|---:|---:|---:|---:|---:|
| Lung adenocarcinoma TMA1 | 43,565 | 100.00% | 2.29% | 12.62% | 3.35 px |
| Lung adenocarcinoma TMA2 | 48,602 | 80.50% | 2.39% | 13.50% | 3.41 px |
| Mesothelioma TMA1 | 60,395 | 100.00% | 3.83% | 16.61% | 3.94 px |
| Mesothelioma TMA2 | 46,955 | 100.00% | 3.15% | 14.38% | 3.81 px |

The 80.50% result in Lung TMA2 is preserved as a source-table discrepancy: 39,123 cells are shared by
metadata, polygons and expression, while metadata contains 48,602 rows. Missing cells are not fabricated
or silently dropped from the audit denominator.

### Definition sensitivity

The median ratio between the strict and inclusive target-distance estimates across sample–target pairs
was **9.01×**, with a maximum of **10.75×**. Therefore an RNA-endothelial threshold is a dominant
structural assumption, not a minor tuning parameter.

Under the balanced definition, pairwise Spearman correlations of the six-target distance ranking ranged
from **−0.60 to 0.71**, with a median of **−0.23**. Target proximity rankings therefore did not
generalise across the four samples. The release explicitly records:

```text
DEFINITION_SENSITIVE_AND_TARGET_RANKING_NOT_STABLE
```

This negative result prevents a visually attractive but unsupported cross-tumour target ranking.
RNA-positive endothelial cells remain `RNA_ENDOTHELIAL_PROXY_NOT_PERFUSION`.

## HER2 IHC ordinal benchmark

The AHIHCI benchmark contains 11 cases spanning HER2 scores 0, 1, 2 and 3. Relative brown-signal features
were extracted inside the supplied masks while retaining whether each mask represented all tumour cells
or only positive tumour cells.

- Best prespecified feature: `brown_median`.
- Spearman versus ordinal HER2 score: **0.921**.
- Exact max-statistic permutation p-value across the five prespecified image features: **0.000758** from
  3,960 unique score assignments.
- Spearman after restricting to the nine denominator-valid `tumor_cells` masks: **0.868**.
- Two masks contain only positive tumour cells and cannot define a positivity denominator.

The exact permutation test corrects for selecting the best result among the five candidate features. It
does not turn an 11-case ordinal benchmark into a quantitative receptor calibration. HER2 score remains
an ordered pathology category, not receptors per cell or binding capacity.

## SHG collagen pilot

The SHG feature extractor was run on 12 images: four tumour-adjacent colon, four cancerous colon and four
normal colon images. It reports relative Otsu-positive fraction and structure-tensor orientation
coherence.

| Pilot class | Images | Median SHG-positive fraction | Median orientation coherence |
|---|---:|---:|---:|
| Adjacent colon | 4 | 0.0730 | 0.3715 |
| Cancerous colon | 4 | 0.0654 | 0.2883 |
| Normal colon | 4 | 0.0995 | 0.3987 |

These 12 images are a feature-extractor smoke test, not a powered biological comparison. They are not
registered to the RCC section, and neither intensity nor coherence is an antibody diffusion coefficient.
No transport parameter is derived.

## External artifact integrity

File extensions were no longer trusted. Content-signature validation found:

- The downloaded AHIHCI ZIP fragment is a valid ZIP payload.
- The CosMx expression file is valid GZIP.
- The representative SHG file is valid TIFF.
- The nominal Bordeau PDF is actually an HTML response.
- The nominal Bordeau supplementary ZIP is an XML error response stating that the article is not in the
  Europe PMC open-access subset.

Consequently, the two Bordeau fallback files are excluded from pharmacological validation even though the
Colab downloader recorded successful HTTP downloads. This is a deliberate fail-closed decision.

## 10x breast data

The Colab run downloaded the `S1_Top` Xenium section, its H&E and associated files. The prespecified
contrasting `S2_Mid` and `S2_Bot` sections were not downloaded because the original token names did not
match the provider filenames. A single section was not substituted for the planned multi-section
validation. Status:

```text
NOT_ANALYZED_INCOMPLETE_PRIORITY_SET
```

## Claims boundary

Supported in v0.6:

- cross-platform parsing and geometry reconstruction;
- explicit source-table match rates;
- relative RNA proxy sensitivity in pixels;
- ordinal HER2 image-score association;
- relative SHG feature extraction;
- content-signature validation and fail-closed artifact handling.

Not supported:

- functionally perfused vessels from RNA or CD31 alone;
- conversion of CosMx pixels to micrometres without scale metadata;
- receptors per cell from IHC score or image intensity;
- antibody diffusivity from SHG;
- measured trastuzumab-distribution concordance;
- real-tissue `reachable_fraction`, `penetration_depth`, `expression_reach_gap`, efficacy or clinical
  success prediction.

## Files

- `results/external_validation/external_validation_summary_v0.6.json`
- `results/external_validation/claims_v0.6.json`
- `results/external_validation/cosmx_gse299786/`
- `results/external_validation/her2_ihc_ahihci/`
- `results/external_validation/shg_collagen_pilot/`
- `results/external_validation/external_data_audit/`
