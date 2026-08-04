# Methods

## Scope and preregistration

`reach-gap` is a hypothesis-generating accessibility instrument. It is not a whole-body
pharmacokinetic model, a clinical outcome predictor, or evidence that a programme would have
succeeded or failed. The simulation benchmark and all predictions below were written before
benchmark outputs were generated.

## Transport model

Let `c(x,t)` be free interstitial antibody concentration and `b(x,t)` the fraction of locally
available target occupied. On a two-dimensional tissue section the mechanistic core is

\[
\frac{\partial c}{\partial t}=\nabla\cdot(D(x)\nabla c)
-k_{int}\alpha_B B_{max}(x)\frac{c}{K_D+c}-k_{clear}c.
\]

At vessel pixels, `c = c_v`, where `c_v` is the stated local vascular concentration scenario.
At the outer tissue boundary the normal flux is zero. The quasi-equilibrium bound fraction is

\[
b(x)=\frac{c(x)}{K_D+c(x)}.
\]

The effective diffusivity is

\[
D(x)=D_0\exp[-\beta_E E(x)-\beta_F F(x)],
\]

where `E` and `F` are bounded spatial ECM and fibroblast scores. This exponential form guarantees
positive diffusivity and treats the two image-derived scores as modifiers of a cited free-tissue
diffusion range, not as independently fitted clinical weights.

The steady state is solved by Picard linearisation. At iteration `m`, the saturable sink is
represented by

\[
q^{(m)}(x)=k_{clear}+\frac{k_{int}\alpha_B B_{max}(x)}{K_D+c^{(m)}(x)},
\]

and the finite-volume system

\[
-\nabla\cdot(D\nabla c^{(m+1)})+q^{(m)}c^{(m+1)}=0
\]

is solved with harmonic face diffusivities. Iteration stops only when the relative infinity-norm
change is below tolerance; otherwise the run is invalid.

### Dimensionless interpretation

For a characteristic vessel spacing `L`, the local linearised Thiele modulus is

\[
\phi=L\sqrt{q/D},
\]

and the corresponding penetration length is `ell = sqrt(D/q)`. The saturation ratio is
`Sigma = c_v/K_D`. A binding-consumption group is

\[
Da_B=\frac{k_{int}\alpha_B B_{max}L^2}{D(K_D+c_v)}.
\]

`alpha_B` is a dimensionless calibration multiplier sampled through the full pipeline. Large `phi` or `Da_B` predicts shallow penetration; large `Sigma` partially relieves the
binding-site barrier. These quantities are reported for interpretation but are not fitted scores.

## Metrics

For tumour cells `i`, target-positive indicator `z_i`, and engagement threshold `tau`:

- `target_positive_fraction = sum(z_i) / N`.
- `reachable_fraction = sum(z_i I[b_i >= tau]) / sum(z_i)`, conditional on target positivity.
- `population_reachable_fraction = sum(z_i I[b_i >= tau]) / N`.
- `expression_reach_gap = target_positive_fraction * (1 - reachable_fraction)`. This resolves the
  unit mismatch that would arise from subtracting a conditional fraction directly from a
  population fraction.
- `penetration_depth` is the largest vessel distance among cells with `b_i >= tau`, reported only
  when vessels are detected.
- `decision_stability` is true only when the full reported interval lies on one side of the stated
  reachable-fraction decision threshold.

Uncertainty is propagated by sampling every uncertain physical or calibration parameter from its
registered range and rerunning the solver. Point estimates are medians and intervals are empirical
5th and 95th percentiles. A point estimate without its interval is not a valid index output.

## Barrier attribution

At each target-positive cell and uncertainty draw, three non-negative mechanistic penalties are
computed: squared distance relative to the draw-specific penetration length, logarithmic matrix
resistance `log(D0/D)`, and the local binding-consumption group. Normalised penalties are interpreted
as attribution weights, not causal posterior probabilities. A barrier is named only when one barrier
wins at least 60% of draws and its mean weight exceeds the runner-up by at least 0.15. Otherwise the
output is `BARRIER_INDISTINGUISHABLE`.

## Abstention

`INSUFFICIENT_EVIDENCE` is emitted when no vascular structures are present, antigen calibration is
missing, solver convergence fails, the geometry exceeds the registered maximum vessel distance, or
the uncertainty interval spans the decision threshold. Geometry and calibration failures are not
silently repaired.

## Global sensitivity

First-order and total-order Sobol indices are estimated with two independent quasi-random matrices
and Saltelli-style hybrid matrices. The output of interest is expression–reach gap. Monte Carlo error
is expected at the small default sample used in the committed smoke benchmark; production use should
increase the sample count.

## Benchmark preregistration

Primary simulated comparison: cell-level bound-fraction RMSE and reachable-cell classification
accuracy for the mechanistic solver, an explicitly naive weighted sum, and distance-to-vessel alone.
The benchmark is simulation-only and cannot establish real-tissue validity.

Predictions, fixed before execution:

1. The expression–reach gap will be largest for high-density, rapidly internalised targets in
   stroma-rich tumours.
2. Parameter uncertainty will span the decision threshold for a substantial share of niches,
   forcing abstention.
3. A clinical retrospective, when curated, will most likely be null or weakly suggestive with a
   wide interval.
4. Smaller formats will show greater penetration depth than IgG under matched local vascular
   concentration; increasing antigen density or internalisation will reduce penetration; increasing
   dose will partially rescue it.

## Clinical retrospective preregistration

The primary exploratory statistic is the difference in median expression–reach gap between publicly
documented successful and unsuccessful programmes, with a permutation interval. Required fields are
target, molecule, format, payload, dose, indication, line of therapy, target-expression evidence,
outcome, reason for discontinuation, and source identifiers. Dose, format, indication, line of
therapy, commercial strategy, and incomplete reporting are declared confounders. No model parameter,
threshold, or weight may be learned from the outcome column. The success criterion is a directionally
consistent difference whose permutation p-value is below 0.05 and whose interval excludes zero.
No post-hoc subgroup is permitted. The real retrospective is `NOT_COMPUTED` in this release.

## Morphology-focus image geometry — version 0.5

The RCC protein images are externally stored, tiled pyramidal OME-TIFF planes. The workflow decodes one
declared pyramid level by reading compressed tile byte ranges and decoding each JPEG 2000 payload
independently. Native-resolution planes are never allocated. Physical pixel size is read from OME metadata
and multiplied by `2^level`.

For each channel, the local image measurement for cell `i` is the mean of the 3×3 image neighbourhood
centred on the cell centroid. It is compared with the matched HDF5 cell-aggregate protein value using
Spearman correlation. This is a measurement-assignment audit, not calibration to receptor copies.

The HDF5 reference-positive rule uses robust scaling followed by deterministic Otsu thresholding. Image
thresholds are fixed from the negative-cell tail:

```text
threshold_q = quantile(local image signal among HDF5-negative cells, q)
```

The primary definition uses `q = 0.995`; sensitivity uses 0.99 and 0.999. Thresholds are not selected to
maximize concordance or any clinical endpoint.

Image-CD31 structural masks retain connected components with area at least 25 µm². Euclidean distance
transforms are converted to microns. Target-positive tumour-cell distances are reported at each CD31
definition. The balanced definition is the descriptive reference, while the full range drives the
structural uncertainty statement.

alphaSMA and Vimentin are baseline-subtracted using the image mode and summarized in prespecified
0–10, 10–25, 25–50, 50–100 and ≥100 µm bands from image-CD31 structures. They remain relative stromal
signals and do not set `D(x)` in a real mechanistic run.

Provider QC masks exclude invalid pixels. Protein-cycle background images are used only for relative
residual spatial-correlation QC. Neither operation creates an antigen or transport calibration.

## External four-field perfusion and native Xenium Zarr completion — version 0.7.1

The S-BIAD3159 completion applies the v0.7 locked segmentation and distance-profile workflow to all four
fields. Each field is evaluated at red-channel subtraction coefficients 0, 0.5 and 1. A field is called
robust only when every sensitivity setting retains a negative Spearman association between relative
Hoechst signal and CD31 distance and enrichment within 10 µm versus 50–100 µm. Fields, rather than pixels,
are the independent image units; no RCC vessel identity is transferred.

The breast reader resolves a single `cell_features` group and obtains `feature_keys`, `number_features`
and `number_cells` from group attributes. Packed two-integer cell identifiers are converted to canonical
Xenium strings. The preferred representation is feature-by-cell CSR (`data`, `indices`, `indptr`), with
`csc/` accepted only as a complete fallback. Ambiguous suffix matches fail rather than selecting the first
array. ERBB2 is required to match exactly one normalized feature name.

Counts are joined to provider cell groups by canonical cell ID. Section and provider-labelled tumour-group
means, medians and positive fractions are descriptive. Cells are not treated as biological replicates and
no cell-level inferential p-value is emitted. ERBB2 RNA is not converted to HER2 receptors per cell.


## v0.8 evidence graph and relative target proxy

The v0.8 evidence-readiness audit uses eight declared requirement weights that sum to one and categorical
satisfaction values of 1.00 (`SATISFIED`), 0.50 (`PARTIAL`), 0.25 (`EXTERNAL_ONLY`) and 0.00
(`MISSING`). The readiness score is 100 times the weighted satisfaction sum. Unresolved uncertainty
contribution is `weight × (1 - satisfaction)` and is normalised across unresolved requirements. This is
an audit rubric, not a probabilistic model.

The relative RCC target proxy uses tumour-region protein-positive fraction, target-positive median
structural-vessel distance and target-positive fraction within 50 µm. For each draw, one of six vessel
definitions is sampled uniformly and component weights are sampled from `Dirichlet(1,1,1)`. Components
are min-max normalised within the four-target set for that draw. Rank frequencies describe robustness to
these declared choices; they are not Bayesian posterior probabilities or absolute reachability estimates.
