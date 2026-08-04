# Critical evaluation — v0.7.0

## Overall assessment: 8.0/10

Version 0.7 is a meaningful scientific improvement rather than another software-only increment. It adds
independent in-vivo perfusion-proxy evidence, a quantitative receptor calibration source, measured IgG
transport priors and a valid administered-antibody reference. It also reports a failed attempt to recover
the publication-level penetration direction from compressed representative figures.

### Biologics discovery hiring manager — 8.0/10

The project now demonstrates awareness of the entire evidence chain: vessel presence, perfusion,
transport, binding capacity and administered-antibody distribution. The strongest signal is the refusal
to transfer external calibrations silently.

### Methods reviewer — 7.7/10

Strengths include fixed thresholds, explicit channel-contamination sensitivity, source-protocol bounds,
unit conversions and preservation of the Bordeau negative result. The main weakness is that each
calibration layer comes from a different tissue or assay.

### Skeptical statistician — 7.7/10

The exact six-panel permutation calculation is appropriate but underpowered and descriptive. The
perfusion analysis currently covers two transferred views; all four TIFFs should be processed in the
same locked workflow. The HER2 fit has strong within-source association but no external calibration
holdout.

### Software engineer — 8.9/10

The new modules are typed, tested and CLI-accessible. Inputs fail explicitly. Results and claim scope are
machine-readable. Full-package strict Mypy remains clean; the new modules also pass Pyright 1.1.411.

### Reproducibility auditor — 8.5/10

All locally used external artifacts pass content-signature validation and are hashed. Large Drive files
remain a practical transfer limitation. A Colab runner is supplied for the complete four-image and breast
bundle execution.

## What prevents a higher score

A 9/10 release would require a single matched dataset containing administered antibody, structural and
functional vasculature, target protein/copy calibration, tissue geometry and animal/sample-level raw
measurements, evaluated without fitting to the endpoint. The current evidence layers are scientifically
useful but are not co-registered.
