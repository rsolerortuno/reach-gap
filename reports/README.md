# Reports

This directory contains portfolio-oriented summaries derived from the machine-readable outputs in `results/`.

## Figures

- `reach_gap_workflow.png` — architecture and abstention logic.
- `rcc_target_geometry.png` — target prevalence versus structural vessel proximity.
- `rcc_target_maps_panel.png` — four-target spatial heterogeneity across the RCC section.
- `target_rank_probability.png` — probability of rank 1 across uncertainty draws.
- `pairwise_win_probability.png` — pairwise target-ranking probabilities.
- `leave_one_component_out.png` — robustness after removing each score component.
- `measurement_priority.png` — experiments ranked by unresolved evidence burden.
- `external_validation_panel.png` — perfusion, calibration, breast Xenium and RNA–protein checks.
- `validation_summary.png` — release quality gates.

## Metrics

The `metrics/` folder contains selected compact CSV and JSON outputs copied from the canonical `results/` tree for easier review. The canonical results remain the source of truth.

To regenerate the portfolio figures, install the optional visualization dependencies with `python -m pip install -e ".[viz]"` and run `python scripts/generate_portfolio_figures.py`.

The GitHub workflows regenerate these figures into a temporary directory before packaging a release. The command fails unless exactly nine PNG files are produced.
