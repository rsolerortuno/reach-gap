# Real RCC morphology-focus protein imaging report — v0.5.0

## Scope

Version 0.5 adds six full-section Xenium morphology-focus protein images to the existing RCC cell-level
RNA/protein analysis: CD31, alphaSMA, Vimentin, PanCK, PD-L1 and VISTA. The images, their provider QC
masks and protein-cycle background maps were processed without loading the native planes into memory.
The result is an image-derived vascular/stromal geometry layer and a matched image-versus-cell
measurement audit. It is not an antibody-distribution measurement.

## Inputs and execution

- Cells: **465,534**.
- Native morphology plane: **54,002 × 27,328 pixels** at **0.2125 µm/pixel**.
- Analysis pyramid level: **4**, corresponding to **3.4 µm/pixel**.
- Channels: CD31, alphaSMA, Vimentin, PanCK, PD-L1 and VISTA.
- QC: matched morphology-focus QC masks and `background_02_*` acquisition images.
- Runtime: **32.33 seconds**.
- Peak resident memory: **917.0 MiB**.

The OME-TIFF tiles use JPEG 2000. Because the active registry did not provide `imagecodecs`, the pipeline
reads compressed tile byte ranges directly and decodes one tile at a time with Pillow. All reconstructed
source channels matched the split-manifest checksums.

## Image–cell measurement concordance

The table compares local morphology-focus intensity at the cell centroid with the provider HDF5
cell-aggregate protein value. Spearman correlation is descriptive; neither side is an absolute surface
copy-number measurement.

| Channel | Spearman | Image TPR at fixed ~0.5% negative-tail FPR |
|---|---:|---:|
| CD31 | **0.761** | 0.261 |
| alphaSMA | **0.836** | 0.199 |
| Vimentin | **0.853** | 0.189 |
| PanCK | **0.476** | 0.441 |
| PD-L1 | **0.334** | 0.485 |
| VISTA | **0.812** | 0.158 |

The threshold was not optimized for accuracy. It was fixed to the 99.5th percentile of local image
signal among cells negative under the deterministic HDF5 within-section threshold. This makes the false
positive target transparent and prevents clinical or outcome labels from entering calibration.

The low PD-L1 concordance shows that cell-aggregate intensity and local pixel signal are not
interchangeable. VISTA aggregate positivity is broad, whereas strong local VISTA pixel signal is much
more restricted. Neither should be interpreted as membrane copy number.

## Image-defined vascular geometry

CD31 structures were filtered by a minimum connected area of 25 µm². Three negative-tail definitions
were preregistered for robustness:

| Definition | CD31 pixel fraction | Components | Median tumour-cell distance | PD-L1-positive median | VISTA-positive median |
|---|---:|---:|---:|---:|---:|
| Inclusive, 1% tail | 6.52% | 16,906 | 14.42 µm | 19.23 µm | 13.60 µm |
| Balanced, 0.5% tail | 5.10% | 15,763 | 17.00 µm | **21.77 µm** | **14.42 µm** |
| Strict, 0.1% tail | 2.65% | 12,349 | 24.04 µm | 30.41 µm | 20.68 µm |

Absolute distances depend materially on the CD31 definition. The qualitative ordering is stable across
all three definitions: VISTA-positive tumour cells are closer to image-CD31 structures than PD-L1-positive
tumour cells. This is a relative spatial result, not evidence that those structures were perfused.

Under the balanced definition, 56.25% of PD-L1-positive tumour cells and 70.87% of VISTA-positive tumour
cells were within 25 µm of an image-CD31 structure. At 100 µm, the values were 89.76% and 95.89%,
respectively.

## Perivascular stroma

The median excess alphaSMA and Vimentin signals decreased with distance from balanced image-CD31
structures:

| Distance | alphaSMA median / q90 | Vimentin median / q90 |
|---|---:|---:|
| 0–10 µm | 2 / 12 | 19 / 56 |
| 10–25 µm | 1 / 5 | 15 / 52 |
| 25–50 µm | 1 / 4 | 14 / 50 |
| 50–100 µm | 0 / 3 | 11 / 48 |
| ≥100 µm | 0 / 2 | 1 / 35 |

This supports a measurable perivascular stromal gradient. It does not calibrate a matrix diffusion
coefficient: Vimentin is not ECM, alphaSMA is not collagen, and both can be expressed outside CAF or
perivascular compartments.

## Target-image assignment

Within pathology-defined tumour regions:

- PD-L1: HDF5-positive **4.64%**, image-positive **2.91%**, positive by both **2.35%**.
- VISTA: HDF5-positive **43.94%**, image-positive **6.09%**, positive by both **5.94%**.

The difference is especially large for VISTA. The provider HDF5 value is a cell-level aggregate, while
the image rule selects strong local signal near the centroid. The two measurements answer different
questions and are retained separately.

## Scientific status

```text
REAL_PROTEIN_IMAGE_GEOMETRY_PREPARED_ABSOLUTE_INDEX_NOT_COMPUTED
```

The real mechanistic index remains uncomputed because:

- CD31 presence does not identify functionally perfused vessels.
- Fluorescence intensity is not calibrated to surface-antigen molecules.
- alphaSMA/Vimentin do not provide a calibrated matrix transport coefficient.
- No administered-antibody distribution, receptor occupancy or engagement is measured.

The committed outputs are in `results/real_rcc_imaging/`. Claims and abstention reasons are machine
readable in `claims.json`.
