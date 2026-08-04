# Data for non-experts

## RCC Xenium

Xenium measures many RNA molecules and selected protein signals while preserving the location of individual cells in tissue. In this project it provides the main RCC cell map, target measurements and spatial geometry.

## H&E imaging

H&E is a standard tissue stain used by pathologists. It helps define tumour, immune, necrotic and other tissue regions. The project uses these annotations to restrict and interpret the spatial analysis.

## Morphology-focus protein imaging

These images measure markers such as CD31, alphaSMA, Vimentin, PD-L1 and VISTA. They are useful for checking local assignment and spatial organisation, but fluorescence values are not automatically molecules per cell.

## Perfusion-proxy microscopy

Hoechst injected into living animals can act as a relative proxy for vascular access. The independent images used here validate the extraction method, but they come from LLC tumours rather than the RCC section.

## HER2 calibration

A published workbook links source-protocol Cy5 intensity to HER2 receptors per cell. The relationship is strong within that assay, but it cannot be applied to Xenium without a shared calibrator.

## IgG transport studies

Published FRAP experiments measure how antibodies move through tumour tissue. These values define a sensitivity range, not an RCC-specific measurement.

## Why raw data are not included

Several source datasets are very large and remain subject to the original providers' terms. The repository therefore includes compact derived tables, figures, manifests and checksums rather than redistributing the complete raw files.
