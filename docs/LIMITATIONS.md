# Limitations

This file was written before results were generated.

- Transcript abundance is not surface protein density. Without an assay-specific calibration,
  antigen capacity is unknown and the tool must abstain.
- A two-dimensional section is not a three-dimensional tumour. Out-of-plane transport can make
  nearest-vessel distances and apparent barriers wrong.
- Vascular perfusion is not vascular presence. A vessel marker does not prove that the vessel
  delivered the molecule at the analysed time.
- Parameter ranges are wide and may still omit structural uncertainty.
- The clinical retrospective is confounded and underpowered; dose, format, indication, line of
  therapy, payload, safety, development strategy and commercial decisions can dominate outcome.
- Nothing in the committed release is validated against measured drug concentration in the same
  tissue used to derive the spatial features.
- The model omits convection, interstitial pressure gradients, FcRn recycling, target shedding,
  nonspecific binding, bivalent avidity, heterogeneous vascular permeability and systemic PK.
- ECM and CAF effects are represented as scalar diffusion attenuation, not explicit fibres, pores or
  anisotropic transport.
- The steady-state approximation can be wrong for short exposures or rapidly clearing formats.
- Simulation performance and relative real-tissue geometry cannot support a clinical or target-selection claim without calibrated delivery measurements.

## Xenium RCC-specific limitations

- The sample represents one stage III renal cell carcinoma donor and cannot establish generality.
- Xenium protein values are scaled per-cell mean fluorescence intensity, not receptor copy number.
- CD31-positive/endothelial cells do not identify functional perfusion or vessel lumen topology.
- Cell-centroid distances are not exact distances to the nearest vessel wall.
- The panel resolves CAF-associated signals but no direct ECM marker set used by this workflow; collagen
  density, fibre orientation, pore size and interstitial pressure remain unmeasured.
- Pathology polygons are supplied in the post-Xenium H&E coordinate system; use requires an inferred
  transform that may abstain.
- The first pass does not use all individual transcript coordinates or all protein-image pixels, although
  it processes every cell and inventories the full archive.
- No measured therapeutic-antibody concentration, receptor occupancy or matched pharmacology is present.

## H&E-only geometry layer

The committed real RCC run contains histology and supplied pathology polygons, but no cell-level target
measurements. H&E colour and morphology are not substitutes for RNA, protein or antigen density.

Bright intratissue spaces detected from H&E are not validated vessels. They can include tissue tears,
tubules, ducts, adipose spaces and processing artefacts. They are retained only for QC and are excluded
from the mechanistic model.

The supplied `Blood vessels` annotation is a broad pathology region rather than a complete vascular
network. Vascular presence also does not establish vascular perfusion.


## Version 0.4 molecular layer

- Target positivity is based on deterministic within-section intensity thresholds. These are not clinical
  assay cutoffs and cannot be compared directly across sections without normalization controls.
- RNA–protein correlations are descriptive for one donor, one section and one panel.
- Six vessel definitions expose threshold sensitivity but do not span all plausible vascular models.
- Boundary-area agreement tests representation consistency, not biological segmentation accuracy.
- The absolute mechanistic index remains uncomputed; the real results are molecular geometry, not drug
  penetration.

## Version 0.5 morphology-focus image layer

- CD31 image structures improve continuous geometry but still do not identify functional perfusion.
- The 3.4 µm/pixel analysis level resolves section-wide structure, not submicron membrane localization.
- Local centroid-neighbourhood signal and HDF5 cell aggregates are different measurements; agreement is
  neither expected nor forced.
- Negative-tail thresholds are transparent within-section rules, not clinical cutoffs or cross-study
  normalizations.
- alphaSMA and Vimentin provide stromal context but no calibrated ECM transport coefficient.
- Background correlations are relative QC and cannot prove that residual signal is or is not artefactual.
- The image data remain from one section and one donor.

## Version 0.7–0.7.1 external calibration layer

- S-BIAD3159 validates a relative Hoechst/CD31 gradient in four independent LLC mouse-tumour fields, not
  perfusion of the RCC vessels. The four fields are not independent animals.
- The quantitative HER2 curve is specific to a Cy5 IHC protocol and cannot calibrate Xenium intensity.
- The 5.4–31.2 µm²/s IgG interval is a literature sensitivity prior, not an RCC-specific measurement.
- Bordeau representative images are compressed and not the animal-level raw data. Their scale-free panel
  analysis did not reproduce the published group direction.
- ERBB2 was extracted from two breast Xenium sections, but the sections are not donor-level replication.
  Cell-level counts are descriptive and are not clinical HER2 scores or surface-receptor measurements.
- The external perfusion, breast expression, Cy5 calibration, IgG prior and RCC geometry are not
  co-registered and cannot be multiplied into an absolute reachability estimate.

## v0.8-specific limitations

The evidence-readiness score is sensitive to its explicit requirement weights and categorical
satisfaction rubric. It is intended for audit and experimental planning, not cross-project benchmarking.
The relative target scores are normalised within one RCC section and four measured targets, so they are
not comparable across datasets or target sets. Two spatial components—median distance and fraction
within 50 µm—are correlated summaries of structural proximity; leave-one-component-out analyses are
reported to expose this dependence. Rank frequencies integrate design-choice uncertainty, not biological
sampling uncertainty. Functional perfusion, receptor copies, drug exposure and pharmacological outcome
remain absent.
