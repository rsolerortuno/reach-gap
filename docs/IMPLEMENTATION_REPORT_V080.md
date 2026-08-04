# Implementation report — v0.8.0

## Added modules

- `evidence_synthesis.py`: typed evidence nodes and edges, fixed readiness rubric, uncertainty budget,
  blocked-transfer audit and measurement priorities.
- `relative_accessibility.py`: deterministic Monte Carlo target ranking across six structural-vessel
  definitions and uniformly uncertain objective weights, including pairwise and component-ablation
  analyses.
- `v080_results.py`: integrated package builder, SHA-256 manifest and scientific-boundary validator.
- `visuals.py`: dependency-light, accessible SVG bar charts for compact release artifacts.

## Added CLI commands

```text
reach synthesize-evidence
reach rank-relative-accessibility
reach build-v080-package
reach validate-v080-results
```

## Scientific design controls

- Same-tissue measurements and external evidence have different node classes.
- Transfer-blocking edges are first-class outputs rather than documentation-only caveats.
- Absolute outputs are validated to remain `NOT_COMPUTED`.
- The target proxy uses only same-section RCC target and structural-geometry measurements.
- Objective weights are not fitted to a desired ranking.
- Six vessel definitions are sampled rather than selecting the most favourable one.
- Leave-one-component-out analyses test whether the stable top target depends on one feature.
- Manifest hashes detect modification of every compact result artifact.

## Integrated outputs

`results/evidence_synthesis_v0.8/` contains:

- evidence graph JSON;
- requirement and measurement-priority CSV files;
- uncertainty-budget SVG;
- relative target result JSON;
- target, pairwise and ablation CSV files;
- target-rank SVG;
- aggregate summary and claims;
- validation report; and
- SHA-256 artifact manifest.
