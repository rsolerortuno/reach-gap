# Critical evaluation — v0.8.0

## Overall assessment: 8.8/10

Version 0.8 is a meaningful analytical increment. It does not add a new biological assay; instead, it
turns the heterogeneous evidence accumulated through v0.7.1 into a reviewable dependency system and
extracts one new same-tissue result that is supportable without pretending that absolute antibody
reachability is known.

### Biologics discovery hiring manager — 8.8/10

The project now shows not only spatial analysis and mechanistic modelling, but also experimental-design
reasoning. The most useful output is the measurement-priority audit: it identifies the administered-drug
field, functional perfusion and surface-antigen calibration as the main missing pieces rather than hiding
them behind a single confidence score.

### Methods reviewer — 8.6/10

Strengths are explicit evidence classes, transfer-blocking edges, a fixed scoring rubric, six-definition
geometry sensitivity and leave-one-component-out analysis. The principal limitation is that the 40.5
readiness score depends on declared weights. It is transparent and useful for audit, but it is not a
validated scientific scale.

### Skeptical statistician — 8.3/10

The target rank probabilities are Monte Carlo frequencies over an explicitly defined design space, not
sampling-based inferential probabilities. The analysis correctly avoids cell-level p-values and does not
fit weights to the observed ranking. VISTA's stability is strong within the proxy, but there is one RCC
section and four targets, so external generalisation is untested.

### Software engineer — 9.3/10

The new functionality is typed, deterministic, CLI-accessible and independently validated. Compact SVG
figures avoid adding a plotting runtime dependency. Hash manifests and tamper tests make the result
package auditable.

### Reproducibility auditor — 9.2/10

Every v0.8 compact output is generated from committed v0.7.1 artifacts, recorded in a SHA-256 manifest
and checked by 70 package invariants. The analysis can be regenerated without the multi-gigabyte raw
bundles.

## What prevents a 9.5/10

- Only one RCC section supports the target ranking.
- Relative target positivity thresholds are within-section and not clinical cut-offs.
- Structural vessel definitions remain unvalidated for functional perfusion in RCC.
- The evidence-readiness weights are expert conventions rather than empirically calibrated utilities.
- No administered-antibody field or same-tissue pharmacological endpoint exists.

A materially stronger next release would require a matched dataset rather than additional software-only
complexity.
