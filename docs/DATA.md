# Data

The committed repository contains only deterministic simulated tissue fixtures and small tabular schema
examples. No human spatial dataset, immunoPET dataset or clinical-trial database is bundled.

Real-data adapters accept flat cell tables and optional raster arrays. Users are responsible for
checking licences, consent restrictions and platform export terms. Every input manifest must include a
source URL or accession, licence, checksum, platform, coordinate unit, segmentation version, vessel
marker definition, antigen calibration method and whether perfusion was measured.

## Public RCC H&E input committed in v0.3

The repository does not commit the 3.72 GB source image. It commits derived QC images, JSON metrics,
claims and provenance hashes under `results/real_rcc_he/`.

Source names:

- `Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_image.ome.tif`
- `Xenium_V1_Human_Kidney_FFPE_Protein_updated_he_imagealignment.csv`
- `Xenium_V1_Human_Kidney_FFPE_Protein_updated_annotation.geojson`

The full molecular ZIP remains external. Use `extract-xenium-essentials` to create a reduced,
checksummed package.

## Public RCC Xenium compact package

Version 0.4 consumes cell-level derivatives from the public FFPE human renal cell carcinoma Xenium
RNA/protein dataset. The source archive and image are not redistributed. The committed processing
manifest records names, byte sizes and SHA-256 values for the exact cell matrix, boundaries, panels,
pathology annotation and alignment used.

Full per-cell outputs are generated data, not source data. They are distributed in a separate checksummed
analysis package; aggregate summaries and deterministic samples remain in the repository.

## RCC morphology-focus images — version 0.5

The source repository does not redistribute the multi-gigabyte provider OME-TIFF files. Committed outputs
record checksums for CD31, alphaSMA, Vimentin, PanCK, PD-L1 and VISTA, their QC masks and the five
`background_02_*` images. Aggregate results and figures are under `results/real_rcc_imaging/`.

Required provider names are resolved exactly by `prepare-rcc-protein-imaging`; no channel is inferred from
position alone. The images are processed at a declared OME pyramid level with physical scale preserved.

## Version 0.7 sources

- S-BIAD3159: four paired CZI/TIFF confocal fields with Hoechst 33342, CD31 and CD4/CD8. The committed
  release run uses the two TIFFs transferable under the connector limit; all eight source files remain in
  Drive.
- McKinski 2026 quantitative HER2: article PDF and `41598_2026_42898_MOESM1_ESM.xlsx`.
- Netti 2000, Ramanujan 2002 and Davies 2002: valid source PDFs for FRAP/transport priors.
- Bordeau 2021: valid `NIHMS1684776-supplement-1.docx`; invalid historical fallback files remain excluded.
- 10x breast S2-Middlle and S2-Bottom: Explorer ZIP, post-Xenium H&E, alignment and cell-group CSVs in
  external storage. Only the cell-group CSVs are committed as derived summaries.
