# Critical evaluation — v0.6.0

The project rule remains: a release without concordance against valid measured therapeutic-antibody
distribution cannot score 8/10 or higher.

| Perspective | Score | Main reason it is not higher |
|---|---:|---|
| Biologics discovery hiring manager | 7.8/10 | The software now generalises across formats and exposes unstable proxy rankings, but it still cannot support a dose, format or target decision. |
| Methods reviewer | 7.6/10 | External adapters and ordinal/image benchmarks are useful, but none validates the reaction–diffusion prediction against drug concentration or engagement. |
| Skeptical statistician | 7.8/10 | Exact feature-selection-aware permutation testing and sample-level sensitivity are strong; sample sizes remain small and the SHG subset is only a pilot. |
| Software engineer | 8.2/10 | 68 tests, content-signature validation, strict typing gates for new modules and a non-regression baseline are credible. Legacy typing debt remains substantial. |
| Reproducibility auditor | 8.3/10 | Inputs, checksums, mismatch rates and invalid downloads are exposed. The populated Drive folder was duplicated by the Colab path and must be documented. |

**Weighted overall assessment: 7.8/10.**

## Strongest improvement

The strongest v0.6 result is not a favourable target ranking. It is that the tool detects when such a
ranking is not stable: CosMx target-distance ranks fail to reproduce across samples, and vessel-proxy
definition changes median distances by roughly ninefold. This behaviour is more credible than selecting
the cleanest map.

## Strongest dismissal risk

A reviewer can still dismiss the central pharmacological claim because no valid raw administered-antibody
distribution, perfusion marker, receptor-copy calibration or tissue-specific diffusion measurement is
present. The Bordeau download initially appeared successful but failed content validation. Therefore
pharmacological concordance remains `NOT_COMPUTED`.

## What would move the score above 8

A prespecified out-of-sample comparison between predicted concentration/bound fraction and valid measured
antibody distribution in the same tissue, including dose, time, target, vascular/perfusion channel and
physical scale. Calibrated antigen capacity and multiple independent tumours would then be required to
approach 9/10.
