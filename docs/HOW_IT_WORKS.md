# How reach-gap works

## A non-technical explanation

An antibody enters a tumour through blood vessels, moves through tissue and binds to its target. A target can therefore be abundant but still difficult to reach if the relevant cells are far from perfused vessels or if transport is limited.

`reach-gap` separates this problem into four questions:

1. **Where are the target-positive cells?**
   Spatial RNA and protein measurements locate cells and measure target signal.
2. **Where are the possible vascular sources?**
   Structural vessel definitions are built from spatial and image-derived evidence.
3. **Which transport assumptions are actually supported?**
   Same-tissue measurements are separated from external validation and literature priors.
4. **Is there enough evidence to report an absolute result?**
   When required measurements are absent, the tool returns `NOT_COMPUTED` instead of inventing a number.

## The v0.8 evidence graph

Every required evidence layer is assigned a state:

- `SAME_TISSUE` — directly measured in the analysed RCC section;
- `EXTERNAL_ONLY` — measured or validated elsewhere;
- `PRIOR_ONLY` — literature-derived sensitivity information;
- `MISSING` — not available;
- `BLOCKED` — a requested output cannot be computed from the available evidence.

The readiness score summarises completeness under a fixed audit rubric. It is not a probability of biological success.

## Relative target ranking

When absolute reachability is unavailable, v0.8 can still compare targets on a restricted same-section proxy using:

- target-positive fraction in tumour regions;
- median structural-vessel proximity;
- fraction of target-positive cells within 50 µm.

The analysis repeats the ranking across six vessel definitions and 20,000 randomly sampled objective weights. It reports rank distributions rather than one deterministic score.

## Why abstention matters

RNA is not automatically surface protein. Protein intensity is not automatically receptor copies. Structural CD31 is not automatically functional perfusion. A literature diffusion coefficient is not automatically the coefficient of the analysed tumour. `reach-gap` preserves these distinctions and blocks unsupported transfers.
