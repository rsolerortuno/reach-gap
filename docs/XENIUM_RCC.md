# Xenium RCC real-data preparation

## Dataset

The supported real-data workflow targets the public 10x Genomics dataset **Xenium In Situ Gene and
Protein Expression data for FFPE Human Renal Cell Carcinoma**, analyzed with Xenium Onboard Analysis
4.0 and licensed under CC BY 4.0.

The workflow expects the original provider filenames:

- `Xenium_V1_Human_Kidney_FFPE_Protein_updated_outs.zip`
- `Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_image.ome.tif`
- `Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_imagealignment.csv`
- `Xenium_V1_Human_Kidney_FFPE_Protein_updated_annotation.geojson`

Provider sizes and MD5 values are encoded in `reach_gap.xenium.EXPECTED_RCC_FILES`. The CLI checks all
provider MD5 values by default. The Colab notebook uses a practical first-pass policy: exact size checks
for every file and MD5 checks for files below 1 GB, with skipped large hashes recorded explicitly; full
large-file MD5 verification is one configuration switch. MD5 is used only for provider transfer
integrity, while generated artefacts use stronger provenance hashes where practical.

## What is processed

The adapter inventories the complete ZIP and selectively extracts the full cell table, full cell-feature
HDF5 matrix, gene and protein panels, metrics, experiment manifest, analysis summary and small secondary
analysis CSV files. It deliberately does not extract the complete transcript table or all full-resolution
protein images during the first pass.

Every cell in the HDF5 matrix is processed. Selected expression outputs retain all protein features and
RNA markers required for transparent endothelial, pericyte, CAF/ECM, immune and epithelial/malignant
proxy scores. Feature-level summaries are computed across the full matrix.

Distances are measured to high-confidence endothelial-cell proxies using the CD31 protein signal and
available endothelial RNA markers. This is **not** a perfusion measurement.

## Targets

The initial target adapters resolve the four measured immune-checkpoint proteins:

- PD-L1 / CD274
- VISTA / VSIR
- PD-1 / PDCD1
- LAG-3 / LAG3

Protein expression is preferred where present. RNA is retained only as an explicitly labelled fallback.
The target positivity split is a deterministic Otsu threshold on robust-scaled signal; the threshold and
resolved features are recorded in `cell_scoring_diagnostics.json` and `marker_resolution.json`.

## Tumour-region semantics

The historical generic input column `is_tumour` means **the cell lies in the tumour analysis region**. It
does not mean the cell is malignant. This distinction is necessary for immune and stromal targets.

The adapter attempts to transform the supplied pathology polygons into Xenium coordinates. It evaluates
both affine directions and several documented/likely pixel scales against the observed cell-coordinate
bounds. The polygons are used only when one candidate is both well supported and separated from the
runner-up. Otherwise the workflow abstains from pathology assignment and uses a clearly labelled
150 µm molecular tumour-neighbourhood proxy.

## Absolute index status

The real-data preparation does not produce an absolute mechanistic reachability number. Xenium protein
features are scaled mean fluorescence intensities, not surface molecules per cell or extracellular molar
concentrations. The manifest therefore records:

```text
absolute_index.status = NOT_COMPUTED
antigen_calibration_nM_per_signal = null
```

The prepared target tables can be ingested by `reach`, but the index must return
`INSUFFICIENT_EVIDENCE` until an independent calibration is supplied. Relative maps, marker-derived
barrier features and vascular-distance geometry remain valid descriptive outputs within their stated
limitations.

## Colab

Open `notebooks/prepare_reach_gap_xenium_rcc.ipynb` in Google Colab. The notebook mounts Drive, locates
the four original files, installs optional dependencies, writes the adapter locally and runs the resumable
pipeline. Outputs are partitioned into bounded files so individual parts can be inspected or transferred
without moving the 36.1 GB source archive. The QC stage also attempts to read only the smallest safe H&E
pyramid level; if no safe level exists, it records `NOT_COMPUTED` rather than loading the 3.72 GB image
into memory.

## Version 0.3 status

The H&E/pathology layer is now computed and committed. It verifies image integrity, physical scale,
annotation geometry, tissue masking and morphology summaries. The target-specific molecular layer is
still pending the reduced essential package.

Do not treat the 772 H&E bright-space candidates as vessels. They are deliberately excluded from the
solver. The next molecular gate requires the cell table and cell-feature HDF5 matrix, not the complete
image/transcript payload of the 36.1 GB bundle.

## Version 0.4 execution outcome

The compact package was sufficient to execute the molecular layer over all 465,534 cells without Colab
or the 36.1 GB ZIP. Direct protein channels were available for PD-L1, VISTA, PD-1, LAG-3 and CD31. Cell
and nucleus boundaries were processed with bounded memory.

The output is deliberately split into two levels:

1. **Committed summaries:** target fractions, pathology summaries, RNA–protein concordance, vascular
   sensitivity, segmentation robustness, plots and claims.
2. **Separate per-cell package:** full scored cell and target tables with a SHA-256 manifest.

Absolute reachability remains unavailable for the calibration and perfusion reasons recorded in
`results/real_rcc_xenium/claims.json`.
