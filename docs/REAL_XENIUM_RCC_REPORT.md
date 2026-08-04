# Real RCC Xenium molecular-geometry report

## Status

`REAL_XENIUM_MOLECULAR_GEOMETRY_PREPARED_ABSOLUTE_INDEX_NOT_COMPUTED`

This analysis uses the complete compact cell-level output from the public FFPE human renal cell
carcinoma Xenium gene-and-protein dataset. It is an application of the spatial preparation layer, not a
validation against measured therapeutic-antibody concentration.

## Inputs and integrity

The analysis used:

- `cells.parquet`: 465,534 cell centroids and provider summary areas.
- `cell_feature_matrix.h5`: 543 features × 465,534 cells.
- 405 gene-expression and 27 protein-expression features, plus controls.
- Complete cell and nucleus boundary Parquet tables.
- Protein and gene panel metadata.
- Provider H&E alignment matrix and pathology polygons.

Cell IDs were required to be exactly identical as sets between the cell table and HDF5 matrix. All input
files have SHA-256 values in `results/real_rcc_xenium/processing_manifest.json`.

The accepted annotation mapping was `affine_scale_0.2125`, with score 0.8511 and 91.01% of sampled
polygon vertices inside the cell-coordinate support. The runner-up separation criterion passed; the
alignment was not selected by visual preference.

## Relative target results

Within the pathology-defined tumour region (including tumour-associated immune and stromal cells):

| Target | Cells | Protein-positive | Fraction | Median distance to primary vessel proxy | Within 50 µm |
|---|---:|---:|---:|---:|---:|
| PD-L1 | 335,789 | 15,592 | 4.64% | 16.63 µm | 89.72% |
| VISTA | 335,789 | 147,547 | 43.94% | 13.39 µm | 94.52% |
| PD-1 | 335,789 | 59,132 | 17.61% | 16.78 µm | 94.43% |
| LAG-3 | 335,789 | 77,540 | 23.09% | 16.64 µm | 93.28% |

The positivity cutoffs are Otsu-derived within-section thresholds after deterministic robust scaling.
They do not correspond to clinical IHC cutoffs, receptor copies per cell or an externally calibrated
assay threshold.

## Vessel-definition robustness

The primary proxy combined endothelial and pericyte/smooth-muscle marker scores and classified 77,406
cells (16.63%) as vessel-associated. Six transparent definitions were evaluated, from a strict
CD31-plus-two-endothelial-RNA rule to a broad CD31 threshold. Vessel-positive fractions ranged from
8.84% to 29.56%.

Target-positive tumour-cell median distances varied as follows:

| Target | Minimum | Maximum | Within-50-µm range |
|---|---:|---:|---:|
| PD-L1 | 7.00 µm | 26.68 µm | 77.61–97.90% |
| VISTA | 7.35 µm | 19.35 µm | 87.59–98.97% |
| PD-1 | 8.70 µm | 22.71 µm | 86.88–99.15% |
| LAG-3 | 8.62 µm | 23.12 µm | 84.71–98.77% |

The qualitative conclusion that most target-positive cells lie within 100 µm of at least one
marker-derived endothelial proxy was robust. The exact median distance was not robust enough to be
reported without the range. None of the definitions identifies functional perfusion or a vessel-wall
boundary.

## RNA versus protein

Cell-level raw-signal Spearman correlations were:

| Target | Spearman correlation |
|---|---:|
| PD-L1 | 0.022 |
| VISTA | 0.110 |
| PD-1 | 0.168 |
| LAG-3 | 0.247 |

This provides direct evidence, in this section, against silently treating target RNA as quantitative
surface protein. It does not establish general correlations across tissues, platforms or antibodies.

## Segmentation robustness

Boundary metrics were reconstructed for all 465,534 cell polygons and 459,288 nucleus polygons. Median
relative differences between shoelace polygon area and the provider summary were 4.90% for cells and
6.86% for nuclei. Cell-area differences were below 10% for 97.00% of cells; nucleus-area differences
were below 10% for 77.84% of cells. The median polygon had 25 vertices.

The differences are consistent with provider rounding, rasterisation or boundary-processing conventions,
but this interpretation is an inference. The analysis reports them as sensitivity, not proof that one
representation is correct.

## Why the absolute index abstains

Three required physical inputs remain unavailable:

1. Xenium protein fluorescence has no matched conversion to surface-antigen density.
2. CD31/endothelial presence does not establish which vessels were perfused at dosing time.
3. The measured panel contains no direct ECM-density or collagen-architecture measurement; ECM was not
   imputed from unrelated markers.

Consequently, real-tissue `reachable_fraction`, `expression_reach_gap`, `penetration_depth` and barrier
attribution remain `NOT_COMPUTED`. Filling those fields with nominal assumptions would make parameter
uncertainty larger than the claimed signal and violate the project contract.

## Artefacts

- `real_data_summary.json`: compact machine-readable headline results.
- `target_spatial_summary.csv`: target × subset summaries.
- `rna_protein_concordance.csv`: paired modality comparisons.
- `vessel_calling_sensitivity.csv`: six vessel definitions.
- `segmentation_robustness.json`: polygon robustness statistics.
- `figures/`: target maps and sensitivity plots.
- `claims.json`: permitted, conditional and unsupported statements.
- `full_cell_tables_manifest.json`: checksums for separately distributed per-cell derivatives.
