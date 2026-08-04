# Assumption ledger

| ID | Assumption | Why needed | Failure mode | Coverage |
|---|---|---|---|---|
| A1 | A 2D section can approximate local transport topology. | Spatial assays are usually sectional. | Out-of-plane vessels or barriers reverse inferred distances. | Limitation and geometry-regime abstention; not experimentally resolved. |
| A2 | Vessel masks represent accessible concentration boundaries. | Provides source locations. | Non-perfused vessels are treated as sources. | Explicit limitation; real adapter requires perfusion evidence field but cannot verify it. |
| A3 | Binding is locally quasi-equilibrated relative to the reported steady state. | Reduces two-state kinetics to a saturable sink. | Slow association/dissociation or transient dosing invalidates occupancy. | Analytical and monotonicity tests only. |
| A4 | Effective antigen capacity can be represented in nM-equivalent tissue volume. | Couples spatial target density to consumption. | RNA/protein signal cannot be calibrated, causing arbitrary barrier strength. | Mandatory calibration flag and abstention test. |
| A5 | ECM and CAF scores modify scalar diffusivity. | Encodes spatial hindrance. | Anisotropy, convection, pressure and pore topology are missed. | Robustness tests; unresolved for real tissue. |
| A6 | A steady state is meaningful for the chosen exposure. | Makes the solver tractable. | Plasma concentration falls before penetration equilibrates. | Fixed boundary scenario and explicit limitation; transient PK is not modelled. |
| A7 | Parameter ranges cover relevant biology. | Enables uncertainty propagation. | Ranges can still be too narrow or structurally wrong. | Sobol analysis; no empirical guarantee. |
| A8 | Cell-centre sampling represents cell engagement. | Connects field to segmented cells. | Large cells or membrane polarity make centre values misleading. | Segmentation-jitter band. |
| A9 | Barrier penalties are interpretable enough for attribution. | Provides a qualitative dominant barrier. | Correlated penalties are not identifiable. | Winner-frequency and margin abstention. |
| A10 | Simulation generated with the same model validates implementation, not biology. | Provides numerical ground truth. | Tautological biological validation claim. | README truthfulness test and explicit label. |

## A-011 — Endothelial-cell proxies approximate vascular geometry

- **Needed for:** distance-to-vessel and local vascular-density features.
- **Assumption:** high CD31/endothelial-marker cells trace vessel locations closely enough at cell-level
  resolution.
- **Failure mode:** isolated endothelial cells, segmentation errors or non-perfused vascular structures
  are interpreted as sources.
- **Coverage:** marker resolution, threshold and positive fraction are recorded; perfusion remains false
  in the manifest. Visual QC is required.

## A-012 — Local CAF/ECM scores summarize a transport-relevant neighbourhood

- **Needed for:** diffusion attenuation fields.
- **Assumption:** the mean of marker-derived scores among up to 24 cells within 100 µm reflects the local
  matrix/fibroblast environment.
- **Failure mode:** acellular matrix, anisotropic collagen or barriers visible only in morphology are
  missed.
- **Coverage:** neighbour count, radius and resolved markers are deterministic and reported. H&E/image
  validation remains outstanding.

## A-013 — Molecular tumour neighbourhood is an acceptable fallback only

- **Needed for:** region membership when pathology alignment abstains.
- **Assumption:** cells within 150 µm of high epithelial/malignant-proxy cells lie in a tumour-associated
  region.
- **Failure mode:** clear-cell RCC with weak epithelial markers, normal tubules or tumour-adjacent immune
  aggregates are misclassified.
- **Coverage:** the fallback is labelled in every target table and is never presented as pathologist ground
  truth.

## A-014 — HDF5 protein feature-type naming varies by 10x software

- **Needed because:** protein features must not be omitted when a matrix uses `Antibody Capture` rather
  than `Protein Expression`.
- **Failure mode:** an unrelated feature type containing these words could be misclassified.
- **Mitigation:** feature names and resolved columns are emitted in `marker_resolution.json`; tests cover
  both expected labels.
- **Status:** tested.

## A-015 — A reduced transfer-integrity check is acceptable for the first Colab pass

- **Needed because:** a mandatory approximately 40 GB MD5 read can exhaust practical session time.
- **Failure mode:** a same-size corruption in a large source file could escape the initial size check.
- **Mitigation:** the verification policy and non-computed hashes are explicit; ZIP extraction checks
  member CRCs; full provider MD5 remains one configuration switch away.
- **Status:** partially tested; full source checksums require the real-data runtime.

## A-016 — H&E bright spaces are not vessels

- **Needed because:** morphology can reveal empty or weakly stained spaces, but several renal and
  processing structures have the same appearance.
- **Failure mode:** tubules, ducts, adipose-associated spaces, tears and processing clefts are treated as
  vascular sources, producing spuriously short vessel distances.
- **Mitigation:** candidates are emitted only as `UNVALIDATED_LUMEN_CANDIDATE`, visually auditable and
  excluded from the solver and all vascular-distance metrics.
- **Status:** enforced by tests and by the committed real-data claims artefact.

## A-017 — Supplied pathology polygons are overlapping analysis regions

- **Needed because:** physical areas are useful for QC and later alignment.
- **Failure mode:** summing polygon areas as mutually exclusive tissue fractions double-counts overlap;
  the broad `Blood vessels` polygon is mistaken for vessel segmentation.
- **Mitigation:** areas are reported per label, overlap is retained, and no compositional percentage is
  derived from the sum. The blood-vessel label is explicitly described as a pathologist region rather
  than perfusion evidence.
- **Status:** documented and represented in the real-data output schema.

## A-018 — Pyramid-level morphology preserves only coarse geometry

- **Needed because:** the native H&E contains billions of RGB values and should not be materialized in
  memory.
- **Failure mode:** capillaries and small lumina disappear at the selected approximately 4.38 µm/pixel
  level; stain proxies are relative rather than quantitative.
- **Mitigation:** the selected level and physical scale are recorded, target indexing abstains, and the
  native image remains available for later region-specific tile reads.
- **Status:** bounded-memory execution measured; biological validation remains unresolved.


## A-019 — Within-section protein thresholds are descriptive

- **Needed for:** defining relative target-positive populations without an external clinical cutoff.
- **Failure mode:** Otsu separates technical intensity modes rather than biologically meaningful target
  abundance.
- **Mitigation:** thresholds are labelled within-section, raw modality and positive fractions are emitted,
  and no clinical or pharmacological meaning is assigned.
- **Status:** computed, not externally validated.

## A-020 — Cell-centroid distance approximates vessel-wall proximity

- **Needed for:** a scalable geometry statistic across 465,534 cells.
- **Failure mode:** large cells, elongated endothelial cells and vessel-lumen topology make centroid
  distance differ from membrane-to-wall distance.
- **Mitigation:** six vessel definitions and cell-boundary robustness are reported; the result is never
  described as exact penetration distance.
- **Status:** sensitivity quantified; physical validation outstanding.

## A-021 — Provider boundary order defines valid polygons

- **Needed for:** streaming shoelace-area and perimeter reconstruction.
- **Failure mode:** unordered or self-intersecting vertices yield misleading areas.
- **Mitigation:** polygons cover every cell, provider-summary area differences are reported, and unsupported
  Parquet structure fails explicitly.
- **Status:** computationally checked; no independent segmentation ground truth.

## A-027 — Image-CD31 is a relative source-geometry proxy

**Assumption:** connected CD31 morphology-focus structures are informative about potential vessel-wall
geometry.

**Failure mode:** CD31-positive structures may be nonperfused, collapsed, discontinuous or unrelated to
drug delivery at the dosing time.

**Coverage:** three threshold definitions, component filtering and explicit `FUNCTIONALLY_PERFUSED_VESSELS_NOT_IDENTIFIED` abstention.

## A-028 — Centroid-local image signal can audit cell aggregates

**Assumption:** a 3×3 neighbourhood at 3.4 µm/pixel provides a useful local comparison with HDF5
cell-aggregate protein.

**Failure mode:** membrane signal, irregular segmentation, neighbouring cells and subcellular localization
can make the centroid neighbourhood unrepresentative.

**Coverage:** image and HDF5 measurements remain separate; concordance and assignment disagreements are
reported rather than reconciled.

## A-029 — alphaSMA/Vimentin gradients describe stroma but not transport

**Assumption:** spatial gradients in alphaSMA and Vimentin are useful relative descriptors of
perivascular stroma.

**Failure mode:** either marker can occur in tumour, immune, endothelial or smooth-muscle compartments,
and neither measures collagen density, pore size or anisotropy.

**Coverage:** the gradients are descriptive only and cannot parameterize real `D(x)`.

## A-030 — CosMx local pixel coordinates support within-FOV relative geometry

- **Needed for:** cross-platform adapter and neighbourhood validation.
- **Failure mode:** pixel distances are interpreted as micrometres or compared across differently scaled
  acquisitions.
- **Mitigation:** all outputs retain `distance_unit=px`; distances are calculated within FOV; the absolute
  solver abstains.
- **Status:** computationally validated, physically uncalibrated.

## A-031 — RNA marker combinations are only endothelial proxies

- **Needed for:** testing geometry when no vascular protein/perfusion channel is available.
- **Failure mode:** ambient RNA, dropout and marker choice alter the source set and create a false vessel
  map.
- **Mitigation:** three definitions are reported; approximately ninefold distance sensitivity and unstable
  target ranks are explicit.
- **Status:** sensitivity quantified; perfusion not measured.

## A-032 — Ordinal HER2 score can benchmark relative image features

- **Needed for:** testing whether image features preserve an ordered pathology signal.
- **Failure mode:** score categories, stain protocol and positive-only masks are mistaken for receptor-copy
  calibration.
- **Mitigation:** exact permutation testing, denominator-valid sensitivity and explicit
  `NOT_MOLECULE_CALIBRATION` status.
- **Status:** ordinal association supported; absolute antigen capacity unavailable.

## A-033 — SHG texture is descriptive rather than a transport coefficient

- **Needed for:** validating collagen-feature extraction.
- **Failure mode:** relative brightness or orientation is inserted directly into `D(x)`.
- **Mitigation:** no diffusivity is emitted; the pilot is not registered to RCC and lacks FRAP/tracer data.
- **Status:** feature extraction tested; transport calibration absent.

## A-034 — Successful download status does not establish valid scientific content

- **Needed for:** robust automated data acquisition.
- **Failure mode:** access-denied HTML/XML is saved under `.pdf` or `.zip` and treated as evidence.
- **Mitigation:** content signatures and hashes are audited before analysis.
- **Status:** enforced; two invalid Bordeau payloads detected.

## A-035 — Composite Hoechst intensity is a relative perfusion proxy

- **Needed for:** independent validation that vascular proximity and in-vivo tracer access can diverge.
- **Failure mode:** RGB compositing, CD4/CD8 channel overlap and display scaling alter the blue channel.
- **Mitigation:** analyse a clean CD4 panel and a red-corrected CD8 panel; retain relative results only.
- **Status:** supported in two transferred views; RCC transfer unsupported.

## A-036 — Source Cy5 calibration is internally quantitative but externally nonexchangeable

- **Needed for:** testing receptor-copy calibration mechanics and censoring behavior.
- **Failure mode:** applying the curve to Xenium fluorescence produces false copy numbers.
- **Mitigation:** hard source-MFI bounds and explicit `NOT_XENIUM_TRANSFER` status.

## A-037 — Published tumour FRAP spans a plausible IgG sensitivity interval

- **Needed for:** replacing an uncited generic diffusion range with measured scenario bounds.
- **Failure mode:** the RCC tissue lies outside the external-study envelope.
- **Mitigation:** broad log-uniform prior and no claim that it is an RCC measurement.

## A-038 — Supplementary representative sections are descriptive, not raw validation data

- **Needed for:** auditing administered-antibody evidence when raw images are unavailable.
- **Failure mode:** rescaling/compression changes pixel distances and thresholds.
- **Mitigation:** scale-free metrics, exact descriptive test and no model-concordance claim.

## v0.8 audit assumptions

- Evidence readiness is defined by the fixed eight-layer rubric in `EVIDENCE_SYNTHESIS_V080.md`.
- External-only evidence receives partial audit credit but never satisfies a same-tissue requirement.
- Each of the six structural-vessel definitions is considered equally plausible for the relative proxy.
- All convex combinations of the three relative-proxy components are explored uniformly over the
  simplex; no component weights are learned from outcomes.
- The relative analysis is conditional on the four targets measured in the RCC panel.
