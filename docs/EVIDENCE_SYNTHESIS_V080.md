# Evidence synthesis — v0.8.0

## Purpose

Version 0.8 turns the accumulated reach-gap evidence into an explicit dependency graph. The graph keeps
same-tissue RCC measurements, external method validations, assay-specific calibrations, literature
priors, missing requirements and blocked outputs as separate node classes. Edges record which evidence
supports a method, constrains a sensitivity range, is blocked from transfer, or is required before an
absolute output can be emitted.

The graph is an audit of **evidence completeness**. It is not a biological reachability model and the
readiness score is not a probability that an antibody will reach its target.

## Required evidence layers

The fixed v0.8 rubric assigns the following weights to the eight requirements for absolute RCC
reachability:

| Requirement | Weight | v0.8 evidence state | Satisfaction |
|---|---:|---|---:|
| Tissue geometry | 0.08 | SATISFIED | 1.00 |
| Structural vasculature | 0.08 | SATISFIED | 1.00 |
| Spatial target localisation | 0.12 | SATISFIED | 1.00 |
| Functional perfusion | 0.18 | EXTERNAL_ONLY | 0.25 |
| Surface-antigen calibration | 0.18 | EXTERNAL_ONLY | 0.25 |
| Matrix transport | 0.14 | EXTERNAL_ONLY | 0.25 |
| Administered-antibody field | 0.15 | MISSING | 0.00 |
| Same-tissue pharmacological endpoint | 0.07 | MISSING | 0.00 |

The categorical satisfaction rubric is fixed at 1.00 for `SATISFIED`, 0.50 for `PARTIAL`, 0.25 for
`EXTERNAL_ONLY` and 0.00 for `MISSING`. The weighted result is an **absolute-evidence readiness score of
40.5/100**. Three of eight requirements are fully satisfied in the same RCC tissue.

These weights are an explicit engineering and experimental-planning convention. They have not been
validated as a clinical or regulatory score. Their role is to make assumptions reviewable and to prevent
an external prior from being silently treated as a same-tissue measurement.

## Uncertainty budget and next measurements

Unresolved weighted burden is normalised to one. The largest contributions are:

1. administered-antibody concentration or engagement field: **25.21%**;
2. functional perfusion in the RCC section: **22.69%**;
3. RCC surface-antigen calibration: **22.69%**;
4. RCC matrix transport measurement: **17.65%**; and
5. a blinded same-tissue pharmacological endpoint: **11.76%**.

This is a measurement-priority ranking under the fixed rubric, not a cost-effectiveness analysis. The
first requested experiment is a spatial administered-antibody or target-engagement measurement after
dosing, because it both supplies the missing drug field and creates a direct validation target. A
co-registered intravascular perfusion tracer and a shared quantitative surface-antigen calibrator remain
nearly equal priorities.

## Blocked transfers retained

The graph explicitly blocks the following transfers:

- mouse LLC Hoechst/CD31 gradients do not label RCC vessels as functionally perfused;
- the source Cy5 HER2 calibration does not convert RCC Xenium channels to molecules per cell;
- the SKOV3 trastuzumab reference does not provide an RCC concentration field; and
- literature IgG FRAP values constrain sensitivity only and do not identify RCC diffusivity.

## Output boundary

The graph status is:

```text
EVIDENCE_GRAPH_COMPLETE_ABSOLUTE_RCC_REACHABILITY_NOT_COMPUTED
```

The following remain `NOT_COMPUTED`:

```text
reachable_fraction
penetration_depth
expression_reach_gap
model_pharmacological_concordance
```

Machine-readable outputs are under `results/evidence_synthesis_v0.8/evidence/`.
