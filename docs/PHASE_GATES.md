# Phase gates and adversarial review

## Phase 0 — Design

**Passed:** transport derivation, exact metrics, preregistered predictions, retrospective protocol,
parameter registry, assumption ledger and limitations were written before committed results.

**Adversarial finding:** transcript signal cannot identify binding capacity and vessel markers cannot
identify perfusion. **Decision:** real runs may ingest these inputs but must abstain without calibration;
no default conversion was invented.

## Phase 1 — Skeleton

**Passed:** Python package, `src/` layout, Typer commands, typed schemas, CI, licence and empty-to-green
suite.

**Failure/substitution:** Ruff and mypy were unavailable from the sandbox registry. CI retains both;
local validation used compilation, pytest, coverage and a source audit.

## Phase 2 — Geometry simulator

**Passed:** deterministic vessels, tumour mask, ECM, CAF, antigen and cell centroids with tunable barrier
regimes.

**Adversarial finding:** the initial geometry saturated reachability. **Fix:** move the committed
benchmark to a stroma-rich, lower-source, higher-capacity regime before recording results.

## Phase 3 — Solver

**Passed:** finite-volume reaction–diffusion solver, harmonic face diffusivity, analytical slab limit,
zero-sink limit and monotonicity tests.

**Adversarial finding:** extreme draws can slow Picard convergence. **Fix:** three damping schedules,
then explicit failure; parameter ranges were not narrowed to hide the problem.

## Phase 4 — Spatial features

**Passed:** vessel distance, bounded ECM/CAF fields, calibrated antigen capacity, generic cell-table
adapter and manifest validation.

**Adversarial finding:** a 2D centroid and a vessel marker are weaker evidence than membrane exposure and
perfused 3D vasculature. **Unfixed:** requires experimental inputs outside this release.

## Phase 5 — Index and uncertainty

**Passed:** reachable fraction, population reach, expression–reach gap, penetration depth, interval
propagation, abstention and `claims.json`.

**Failure fixed:** antigen-calibration uncertainty was initially metadata-only. It now multiplies the
reaction sink and all committed artefacts were regenerated.

**Adversarial finding:** barrier penalties were correlated. **Decision:** emit
`BARRIER_INDISTINGUISHABLE` unless winner frequency and margin both pass.

## Phase 6 — Sensitivity

**Passed:** variance-based first- and total-order estimates over all seven uncertain parameters.

**Result:** internalisation had the largest total-order estimate in the committed simulation; antigen
calibration was second. The small sample is a screening analysis, not a precise universal ranking.

## Phase 7 — Benchmark

**Passed:** mechanistic, weighted-sum and distance-only ablations; size, antigen, dose and internalisation
checks; uncertainty niches; segmentation/cell-calling perturbations.

**Adversarial result:** the mechanistic model won RMSE but distance-only won binary accuracy. Both remain
headline results.

## Phase 8 — Retrospective scaffold

**Passed:** schema, curation protocol, confounders, permutation analysis and synthetic schema test.

**Not computed:** no real programme table was fabricated. Distance-matched and clinical nulls remain
`NOT_COMPUTED_REAL_DATA`.

## Phase 9 — Documentation and packaging

**Passed:** README truthfulness tests, results hashes, runbook, reports and distributable ZIP.

**Unresolved:** local Ruff and strict-mypy execution could not be confirmed in the sandbox. This is
explicitly counted against the final grade.

## Phase 10 — Real RCC H&E/pathology geometry (v0.3 extension)

**Passed:** all 40 H&E parts were verified against the supplied SHA-256 manifest, reconstructed
byte-for-byte, and the complete file matched the provider MD5. The OME pyramid, physical scale,
pathology polygons, tissue mask and stain-proxy QC were processed with bounded memory.

**Computed real result:** 53.324 mm² of supplied tumour annotation, 8.583 mm² of immune-infiltration
annotation, 76.794 mm² estimated tissue at the analysis level, and 772 bright intratissue candidates.
The candidates are deliberately labelled `UNVALIDATED_LUMEN_CANDIDATE`.

**Adversarial finding:** visual inspection showed that morphology-only bright spaces include renal
lumina, tissue clefts, adipose-associated spaces and artefacts. **Decision:** none enters the vascular
boundary condition, distance-to-vessel calculation or mechanistic index.

**Scientific gate:** target-specific real-tissue reachability remains `NOT_COMPUTED`. The run abstains
because cell coordinates, molecular target signal, surface-antigen calibration and functionally
perfused vessels are not available from H&E alone.

**Next substitution:** extract only essential cell/molecular members beside the 36.1 GB ZIP using
streaming and direct 95 MB part output. Transferring the full archive is neither necessary nor desirable.

