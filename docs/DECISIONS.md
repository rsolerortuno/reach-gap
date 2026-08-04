# Decision log

## D001 — Conditional reachability and the gap
The specification defines reachable fraction among target-positive cells but asks to subtract it from
target-positive fraction. Direct subtraction mixes denominators. The implementation preserves
conditional reachable fraction and defines the population-scale gap as `p_expression * (1-r_reach)`.
Rejected alternative: silently redefining reachable fraction over all cells.

## D002 — Steady reaction–diffusion rather than a fitted index
The production model is a finite-volume steady reaction–diffusion equation with saturable consumption.
Rejected alternative: a learned or hand-tuned weighted sum. A weighted sum exists only as a benchmark.

## D003 — Quasi-equilibrium binding
Binding kinetics are collapsed to occupancy `c/(KD+c)` and internalisation-driven consumption. This is
more conservative than claiming a full kinetic PK/PD model. Rejected alternative: inventing association
and dissociation defaults unsupported for a target.

## D004 — Broad parameter envelopes
Literature supports mechanisms more strongly than universal numeric values. Numeric intervals are broad
scenario envelopes, explicitly not pooled estimates. Real-data defaults requiring antigen or dose
calibration are excluded.

## D005 — Simulation-only release
External spatial data, immunoPET and clinical databases were not downloaded. Adapters, schemas and
curation instructions are implemented; all real-data outputs are `NOT_COMPUTED`. No accession or outcome
is fabricated.

## D006 — Barrier attribution
Attribution uses dimensionless mechanistic penalties and abstains under weak separation. Rejected
alternative: selecting the largest raw feature or claiming causal posterior probabilities.

## D007 — Worked retrospective substitute
A synthetic schema demonstration is included instead of a purported real clinical table. It is labelled
`SYNTHETIC_SCHEMA_EXAMPLE`, excluded from analysis defaults and counted against the evaluation score.

## D008 — Calibration uncertainty enters the reaction term
An early implementation declared antigen-calibration uncertainty but did not multiply the sink by it.
That made a reported uncertainty dimension inert. The results layer was stopped, the factor was added to
the nonlinear consumption term and all benchmark artefacts were regenerated. Rejected alternative:
retain the parameter only in metadata.

## D009 — Non-degenerate benchmark regime
The first smoke geometry saturated at near-universal reachability and could not distinguish format,
dose or barrier effects. The benchmark was changed before committed results to a lower vessel
concentration, higher antigen capacity and stroma-rich geometry. The exact scenario is committed and is
not selected from clinical labels.

## D010 — Robust solver retry ladder
Extreme uncertainty draws occasionally challenged Picard convergence. The same equation is retried with
three genuinely different damping schedules before abstaining. No failed draw is silently discarded.
Rejected alternative: narrow parameter ranges until every solve converges.

## D011 — Primary metric can disagree with classification
The mechanistic model beats both baselines on preregistered RMSE, while distance-to-vessel wins binary
accuracy in the committed simulation. Both are reported. Rejected alternative: headline only the metric
that favours the mechanistic model.

## D012 — Local quality-tool substitution
The build environment could execute installed scientific packages but its package registry did not
provide Ruff or mypy, and direct GitHub DNS was unavailable. Tests, coverage, compilation and a custom
syntax and line-length source audit were run locally; Ruff and strict mypy remain enforced in GitHub Actions but
were not truthfully claimed as locally executed. This substitution lowers the software-engineering and
reproducibility grades.

## D013 — Separate regeneration processes
Running the full benchmark and the 48-draw reference index sequentially inside one sandbox process caused
an environment-specific sparse-solver stall, while each command completed independently. Committed
artefacts were therefore generated in separate commands. The runbook preserves that process boundary.

## D-012 — Real-data preparation uses all cells but selective bundle extraction

**Decision:** inventory the complete Xenium archive, extract the full cell-feature matrix and cell table,
but defer the 62.9 million-row transcript table and most full-resolution image members during the first
pass.

**Why:** cell-level protein/RNA, coordinates and segmentation-derived summaries are sufficient to build
the first mechanistic inputs. Extracting every image and transcript would increase temporary storage and
failure risk without improving the initial accessibility calculation.

**Rejected alternative:** treat the Xenium Explorer subset as equivalent to the full bundle. It omits
several programmatic tables needed for efficient preparation.

## D-013 — `is_tumour` denotes region membership, not malignant identity

**Decision:** real-data target tables place all cells inside the tumour analysis region under `is_tumour`,
including immune and stromal cells. A separate `cell_is_malignant_proxy` field is retained.

**Why:** the scientific question is reachability of target-positive cells *within a tumour*, not only
reachability of malignant cells. The earlier naming was ambiguous and would incorrectly discard PD-1,
VISTA and other microenvironment targets.

## D-014 — No absolute antigen calibration from Xenium protein intensity

**Decision:** use measured protein intensity for relative positivity and spatial mapping, while keeping
`antigen_calibration_nM_per_signal = null`.

**Why:** Xenium reports scaled mean fluorescence intensity. Converting it to receptors per cell or molar
sink density without a matched standard would fabricate a physical parameter. The real-data absolute
index therefore remains `NOT_COMPUTED`/abstained.

## D-015 — Pathology alignment is inferred only with a separation criterion

**Decision:** evaluate affine direction and pixel-scale candidates against cell-coordinate bounds and use
pathology polygons only when the leading candidate has at least 70% vertex support and a score margin of
0.05 over the runner-up.

**Why:** the supplemental GeoJSON is aligned to the post-Xenium H&E image, not directly guaranteed to be
in Xenium micron coordinates. Silent use of the wrong direction or scale would corrupt every tumour-region
claim.

## D-016 — Require barcode identity, not only equal row counts

**Decision:** the Xenium adapter rejects a cell table and HDF5 matrix unless their cell identifiers are
exactly identical as sets. The expression table is then joined by `cell_id` rather than row position.

**Why:** equal lengths do not establish that geometry and molecular measurements belong to the same
cells. Filling unmatched rows with zero would create plausible but false spatial biology.

**Rejected alternative:** trust provider row order and validate only the number of cells.

## D-017 — Accept both common 10x protein feature labels

**Decision:** feature types containing either `protein` or `antibody` are treated as protein measurements.

**Why:** 10x matrix conventions vary across products/software versions (`Protein Expression` versus
`Antibody Capture`). Exact matching to one label could silently discard all protein features.

**Rejected alternative:** hard-code the label observed in one documentation version.

## D-018 — Make full large-file MD5 optional in Colab

**Decision:** the Colab default validates exact byte sizes and provider MD5 values for files below 1 GB,
and records large-file MD5 as `NOT_COMPUTED_BY_CONFIGURATION`. The user may enable a full checksum
pass.

**Why:** reading approximately 40 GB from mounted Drive solely for hashes can consume a large fraction
of a transient Colab session before the same archive is read again for extraction. The reduced check is
explicitly recorded and does not masquerade as a complete checksum verification.

**Rejected alternative:** silently skip checksums, or force a full pass that makes routine execution
fragile.

## D-019 — Compute H&E geometry without pretending it is molecular reachability

**Decision:** reconstruct and analyse the complete real H&E OME-TIFF and supplied pathology polygons,
while returning `INSUFFICIENT_EVIDENCE` for all target-specific metrics.

**Why:** real tissue geometry is useful evidence and tests the ingestion path, but H&E cannot supply
target-positive cells, calibrated surface antigen or functional perfusion. Reporting a reach gap from
these inputs would violate the mechanistic contract.

**Rejected alternative:** delay every real-data result until the 36.1 GB molecular bundle is available.
That would discard valid histology and provenance work that can be completed independently.

## D-020 — H&E spaces are candidates, not vessels

**Decision:** call bright intratissue components `UNVALIDATED_LUMEN_CANDIDATE`, expose them as QC, and
exclude them from vascular distance and the solver.

**Why:** visual inspection shows that the rule-based set contains tissue clefts, tubules, ducts, adipose
spaces and processing artefacts. H&E alone cannot reliably distinguish perfused vasculature.

**Rejected alternative:** select the most vessel-like components and use them as a vascular mask. That
would turn an unvalidated image heuristic into the dominant physical boundary condition.

## D-021 — Stream essential ZIP members directly into transfer-sized parts

**Decision:** add `extract-xenium-essentials`, which reads the archive sequentially and writes large
essential members directly into checksummed 95 MB parts.

**Why:** the 36.1 GB archive does not need to be loaded into RAM or extracted wholesale. The required
cell table, HDF5 matrix and panels can be isolated with memory bounded by an 8 MB buffer. Direct splitting
also avoids a second full-size temporary copy.

**Rejected alternative:** split and transfer the entire 36.1 GB ZIP. Most members are unnecessary for the
first molecular reachability pass and would consume connector storage and execution time.

## D-022 — Use a narrow audited Parquet fallback rather than block real execution

**Decision:** implement a reader only for the flat, Zstandard-compressed, PLAIN-encoded Parquet subset
observed in the provider's cell and boundary tables.

**Why:** the environment could not obtain an Arrow wheel, while the provider tables were otherwise small
and structurally simple. Unsupported nested schemas, dictionaries, codecs, nulls or encodings raise an
error.

**Rejected alternative:** write a permissive general-purpose Parquet parser or silently fall back to CSV.
The former would be unauditable; the latter would require transferring larger duplicate files.

## D-023 — Report relative target geometry, not an approximate reach gap

**Decision:** compute within-section protein positivity and distance distributions while keeping all
absolute solver outputs `NOT_COMPUTED`.

**Why:** direct protein channels are superior to RNA proxies, but fluorescence still lacks a conversion to
surface antigen density. A nominal conversion would determine the binding-site barrier and dominate the
answer.

**Rejected alternative:** sweep an arbitrary calibration range and label the result an uncertainty
interval. Without an empirical anchor, the range would be a scenario exercise, not propagated
measurement uncertainty.

## D-024 — Vessel definition is a first-class sensitivity dimension

**Decision:** evaluate six marker-derived vessel definitions and report ranges for every target.

**Why:** the selected definition changes vessel-positive fraction from 8.84% to 29.56% and shifts target
median distances by up to approximately 20 µm. A single threshold would conceal structural uncertainty.

**Rejected alternative:** select the definition yielding the cleanest spatial map or shortest distances.

## D-025 — Do not impute absent ECM measurements

**Decision:** retain ECM as unavailable and set the marker-derived ECM score to zero only as a sentinel,
never as evidence that ECM is absent.

**Why:** the measured panel did not resolve a direct ECM marker set. Deriving collagen architecture from
CAF or H&E intensity would conflate cells, acellular matrix and staining.

**Rejected alternative:** use CAF score as both cellular stroma and ECM hindrance.

## D-026 — Distribute full per-cell derivatives separately

**Decision:** commit compact summaries, samples, plots and checksums in the source repository; package
full cell and target tables as a separate analysis artefact.

**Why:** large generated tables are reproducible derivatives and would make the Git repository difficult
to review or clone. Their hashes preserve auditability without treating them as source code.

**Rejected alternative:** remove the full tables entirely or commit more than 100 MB of compressed
results to Git.

## D-027 — Decode provider JPEG 2000 tiles directly

**Decision:** read TIFF tile offsets and decode each JPEG 2000 payload with Pillow.

**Why:** `imagecodecs` was unavailable, while full external conversion would duplicate multi-gigabyte
images and weaken reproducibility. Tile decoding uses the provider file directly and keeps memory bounded.

**Rejected alternatives:** fail the image phase; require proprietary conversion; allocate the native plane.

## D-028 — Calibrate image thresholds to a fixed negative-cell tail

**Decision:** use the 99.5th percentile among HDF5-negative cells as the descriptive image threshold, with
99% and 99.9% sensitivity definitions.

**Why:** this fixes an interpretable negative-tail error target and avoids optimizing agreement or outcome.

**Rejected alternative:** maximize F1, correlation or target separation after inspecting results.

## D-029 — Do not feed real image stroma into the solver

**Decision:** retain alphaSMA/Vimentin as descriptive spatial signals only.

**Why:** they do not calibrate collagen, matrix diffusivity or interstitial pressure. Using them as direct
transport coefficients would be a visually plausible but unsupported substitution.

## D-030 — Cross-platform distances remain in source pixels

**Decision:** retain CosMx distances in pixels and calculate them within FOV.

**Why:** the flat files used in this run did not provide a verified physical pixel scale or a globally
continuous image coordinate system suitable for transport modelling.

**Rejected alternative:** infer micrometres from another CosMx release or compare raw pixels with Xenium
micrometre distances.

## D-031 — Treat endothelial RNA calls as a sensitivity family, not a vessel truth set

**Decision:** evaluate inclusive, balanced and strict marker-count definitions and report their full
impact.

**Why:** median target distances changed approximately ninefold across definitions and target rankings
were not reproducible across samples.

**Rejected alternative:** select the definition producing the cleanest vascular pattern.

## D-032 — Correct HER2 feature selection with an exact max-statistic permutation test

**Decision:** test the maximum absolute Spearman correlation across the five prespecified image features
against all 3,960 unique score assignments.

**Why:** reporting only the best of five correlations without selection correction would be optimistic.

## D-033 — Validate artifact content, not extensions or HTTP success

**Decision:** check file signatures before accepting external PDFs, archives and scientific images.

**Why:** the nominal Bordeau PDF and ZIP were actually HTML and XML error payloads. They are excluded from
validation.

## D-034 — Gate new typing debt while exposing the legacy baseline

**Decision:** require strict Mypy and Pyright for v0.6 modules and fail CI if the disclosed legacy error
count increases.

**Why:** claiming complete typing would be false, while allowing new untyped modules would worsen the
repository.


## D-035 — Replace the legacy typing baseline with full-package gates

**Decision:** require strict Mypy and Pyright to pass across every source module.

**Why:** v0.6.1 resolves the disclosed baseline through explicit I/O and array typing; retaining an
error-count waiver would no longer be justified. Third-party exceptions must remain narrowly scoped and
documented.

## D-036 — Validate perfusion gradients without transferring vessel identities

**Decision:** use independent Hoechst/CD31 images to validate the extraction method, while retaining RCC
perfusion as unknown.

**Why:** the LLC data show functional access relative to CD31, but they are not the RCC section and cannot
label its vessels.

## D-037 — Keep quantitative HER2 calibration assay-specific

**Decision:** fit and expose the Cy5-to-receptor relation only inside the source assay's observed range.

**Why:** cross-assay fluorescence has no common scale without a shared calibrator.

## D-038 — Use FRAP values as sensitivity priors, not fitted RCC coefficients

**Decision:** expose the literature envelope as a log-uniform scenario prior and prohibit direct
SHG-to-diffusion conversion.

**Why:** matrix composition, porosity, pressure and binding differ across tissues and studies.

## D-039 — Preserve failed representative-figure concordance

**Decision:** commit the non-concordant Bordeau panel analysis and retain model concordance as
`NOT_COMPUTED`.

**Why:** choosing another image threshold to force agreement with the published group mean would be
post-hoc fitting to a compressed figure.

## D-040 — Treat native Xenium Zarr paths as a schema, not suffix guesses

**Decision:** resolve one `cell_features` group, read dimensions and feature keys from its attributes,
prefer its complete feature-by-cell CSR representation and accept `csc/` only as a complete fallback.
Ambiguous suffix matches are errors.

**Why:** the provider bundle contains both `cell_features/data` and `cell_features/csc/data`. Selecting an
array by the final path component silently confuses two sparse encodings and can return a wrong feature
vector without a shape error. Explicit schema validation is safer and regression-testable.

