# Relative RCC geometry-expression accessibility — v0.8.0

## Question

Given the measurements already available in the same RCC section, which of the four measured targets is
most consistently favoured by **relative target prevalence and structural proximity**? This question is
narrower than antibody reachability and can be answered without inventing perfusion, receptor-copy or
drug-concentration measurements.

## Inputs

For PD-L1, VISTA, PD-1 and LAG-3, the analysis uses:

1. tumour-region protein-positive fraction;
2. target-positive median distance to a structural vessel definition; and
3. the fraction of target-positive tumour cells within 50 µm of that definition.

All six preregistered RCC structural-vessel definitions are retained. No external perfusion,
calibration, diffusion or drug-distribution measurement enters the target ranking.

## Uncertainty propagation

Each of 20,000 deterministic-seed draws:

1. samples one of the six vessel definitions uniformly;
2. normalises the three components within the four-target set;
3. samples objective weights uniformly over the three-component simplex with
   `Dirichlet(1, 1, 1)`; and
4. ranks the four targets by the resulting convex score.

The resulting top-rank frequencies are a **design-space robustness measure**, not a Bayesian posterior
probability and not a clinical probability of success. Scores are relative to this target set and this
section; adding targets or changing the dataset changes the normalisation.

## Result

| Target | Tumour protein-positive fraction | Probability of rank 1 | Median rank | Proxy score, 5th–95th percentile |
|---|---:|---:|---:|---:|
| VISTA | 43.94% | **98.525%** | 1 | 0.853–1.000 |
| PD-1 | 17.61% | 0.895% | 2 | 0.141–0.817 |
| PD-L1 | 4.64% | 0.580% | 4 | 0.000–0.606 |
| LAG-3 | 23.09% | 0.000% | 3 | 0.195–0.642 |

VISTA beats PD-1 in **99.105%**, PD-L1 in **99.420%** and LAG-3 in **100%** of draws.

The result is not driven by a single component. VISTA remains top in:

- **91.175%** of draws when tumour-positive fraction is omitted;
- **92.225%** when median structural proximity is omitted; and
- **97.225%** when the within-50-µm component is omitted.

The stable result is therefore:

```text
VISTA_IS_STABLE_TOP_TARGET_WITHIN_RELATIVE_RCC_GEOMETRY_EXPRESSION_PROXY
```

## Interpretation

This says that VISTA combines broad measured tumour-region positivity with favourable structural
proximity more consistently than the other three targets in this section. It does **not** show that a
VISTA antibody reaches more cells, binds more molecules, has better pharmacology or should be advanced
clinically. Functional perfusion, surface density, affinity, internalisation, exposure and toxicity are
outside the proxy.

Machine-readable outputs and the leave-one-component-out audit are under
`results/evidence_synthesis_v0.8/relative_accessibility/`.
