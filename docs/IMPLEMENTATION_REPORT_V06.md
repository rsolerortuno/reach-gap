# Implementation report — v0.6.0

## Added

- `reach_gap.cosmx`: four-sample CosMx flat-file ingestion, vectorised polygon reconstruction, explicit
  FOV-local distances, three RNA-endothelial proxy definitions, threshold-sensitivity and target-rank
  stability audits.
- `reach_gap.her2_ihc`: AHIHCI parsing, mask-semantics preservation, relative DAB/brown features,
  exact max-statistic permutation testing and denominator-valid sensitivity.
- `reach_gap.shg_collagen`: relative SHG-positive fraction and orientation-coherence extraction with a
  hard abstention from diffusivity.
- `reach_gap.artifact_validation`: PDF/ZIP/GZIP/TIFF/HTML/XML/JSON content-signature detection.
- CLI commands: `prepare-cosmx-external-validation`, `benchmark-her2-ihc`,
  `benchmark-shg-collagen`, `validate-artifact`.
- CI: strict Mypy/Pyright gates for v0.6 modules plus a non-increasing legacy Mypy baseline.

## Real executions

- CosMx: 4 samples, 199,517 metadata cells, 12.7 seconds.
- HER2 IHC: 11 cases, 7.0 seconds.
- SHG pilot: 12 images, approximately 2 seconds.
- Artifact audit: two invalid Bordeau fallback files detected and excluded.

## Validation

- 68 tests pass.
- Coverage 86.83% against an 85% threshold.
- Ruff lint and format checks pass.
- Strict Mypy and Pyright pass for all four v0.6 modules.
- Full-package strict Mypy remains at the disclosed 224-error Python 3.11 baseline.

## Data-location correction

The Colab notebook wrote to `/content/drive/MyDrive/reach-gap-next-stage`, creating a second folder at the
Drive root. The populated data are in that duplicate, not in the initially created nested folder. The
release records this provenance issue rather than treating the empty folder as a failed download.
