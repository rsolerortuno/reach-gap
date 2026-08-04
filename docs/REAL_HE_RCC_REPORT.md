# Real RCC H&E and pathology report

## Evidence level

This report covers **real histology and supplied pathology annotations only**. It does not contain
Xenium cell coordinates, RNA/protein measurements, calibrated surface-antigen density, perfusion, drug
concentration or receptor occupancy. Consequently, every target-specific reachability metric remains
`NOT_COMPUTED` with `INSUFFICIENT_EVIDENCE`.

## Source integrity

The post-Xenium H&E image was supplied as 40 bytewise parts. Every part matched the SHA-256 recorded in
the split manifest. The reconstructed file contained 3,720,697,771 bytes and matched the provider MD5
`96ad5f699c7d6280cdf6af1c13f39515`.

The image is a tiled BigTIFF/OME-TIFF with seven pyramid levels. Its native dimensions are 60,680 ×
24,096 RGB pixels, with OME physical pixel sizes 0.273770745 µm × 0.273773107 µm. This corresponds to a
field approximately 6.597 mm × 16.613 mm.

## Computed pathology geometry

Exact polygon areas were calculated in native H&E pixel space using the physical pixel sizes declared
in OME-XML.

| Supplied annotation | Area (mm²) |
|---|---:|
| Tumor | 53.324 |
| Immune infiltration | 8.583 |
| Adipose tissue | 4.066 |
| Blood vessels | 2.126 |
| Hemorrhage | 1.365 |
| Necrosis | 0.205 |

The labels are pathologist-provided analysis regions. In particular, `Blood vessels` is one supplied
polygonal region; it is not a complete segmentation of all vessels and does not establish perfusion.
Annotations may overlap, so their areas must not be summed as disjoint tissue fractions.

## Bounded-memory H&E analysis

The committed run used pyramid level 4 (3,793 × 1,506 RGB pixels; approximately 4.38 µm/pixel). A
background-adaptive tissue mask estimated 76.794 mm² of tissue. Relative hematoxylin and eosin channels
were calculated for morphology QC and summarized within each supplied annotation.

The analysis completed in 18.29 seconds after reconstruction, with maximum resident memory of
approximately 1.00 GiB. The native 4.38-billion-value RGB array was never materialized in memory.

## Lumen-candidate result

The conservative candidate generator marked 772 bright intratissue spaces, occupying 2.383% of the
estimated tissue mask. These candidates are explicitly classified as
`UNVALIDATED_LUMEN_CANDIDATE`. Visual review shows that the set includes tissue clefts, tubular or ductal
spaces, adipose-associated spaces and likely processing artefacts in addition to possible lumina.

They are therefore **excluded from the mechanistic solver and from all vascular-distance metrics**.
Using them as vessels would create a plausible-looking but scientifically indefensible accessibility
map.

## Outputs

- `results/real_rcc_he/he_pathology_result.json`
- `results/real_rcc_he/claims.json`
- `results/real_rcc_he/processing_manifest.json`
- `results/real_rcc_he/unvalidated_lumen_candidates.json`
- `results/real_rcc_he/qc/pathology_overlay.png`
- `results/real_rcc_he/qc/tissue_mask.png`
- `results/real_rcc_he/qc/hematoxylin_proxy.png`
- `results/real_rcc_he/qc/eosin_proxy.png`
- `results/real_rcc_he/qc/unvalidated_lumen_candidates.png`

## What is required for the next gate

The full 36.1 GB bundle itself is not required in the chat runtime. The low-memory extractor needs to be
run once beside the bundle to produce only:

1. `cells.parquet` or `cells.csv.gz`;
2. `cell_feature_matrix.h5`;
3. `gene_panel.json` and `protein_panel.json` when present;
4. `experiment.xenium`, metrics and small analysis CSVs.

Large essential members can be emitted directly as 95 MB verified parts. Once those parts are in Drive,
the molecular layer can be reconstructed and processed here without transferring irrelevant image and
transcript payloads.
