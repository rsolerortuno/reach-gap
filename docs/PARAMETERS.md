# Parameter registry

Ranges are broad scenario envelopes, not universal constants or pooled biological estimates. Simulation
defaults lie inside these ranges. Real analyses must replace assay-specific calibration inputs with
verified evidence or abstain.

| Parameter | Range | Unit | Role | Evidence and qualification |
|---|---:|---|---|---|
| IgG free-tissue diffusivity `D0` | 5–30 | µm²/s | Transport | Davies et al., *British Journal of Cancer* 2002, DOI: 10.1038/sj.bjc.6600270; Thurber et al., *Advanced Drug Delivery Reviews* 2008, DOI: 10.1016/j.addr.2008.04.012. Deliberately wider than a single measurement. |
| Fab `D0` | 12–60 | µm²/s | Size-series scenario | Fragment-size dependence follows antibody transport theory. This is a scenario interval, not a pooled estimate. |
| scFv `D0` | 20–100 | µm²/s | Size-series scenario | Same qualification as Fab. |
| Dissociation constant `K_D` | 0.1–10 | nM | Occupancy | Rudnick et al., *Cancer Research* 2011, PMID: 21406401, examined affinity and internalisation effects. The range is a scenario envelope. |
| Calibrated binding capacity field `Bmax(x)` | supplied by input; simulation approximately 0–300 | nM-equivalent tissue volume | Binding-site barrier | Ackerman et al., *Molecular Cancer Therapeutics* 2008, PMID: 18645032, and Thurber et al. support antigen-level effects. Spatial RNA or intensity is not a valid value without calibration. |
| Antigen calibration factor `alpha_B` | 0.25–4 | dimensionless | Propagated calibration uncertainty | Multiplies the supplied `Bmax(x)` field. The broad range represents residual assay-to-capacity uncertainty; it is not fitted to outcome. |
| Internalisation `k_int` | 1e-5–5e-4 | s⁻¹ | Antigen-mediated consumption | Ackerman et al. and Thurber et al. support turnover/internalisation as penetration determinants. Characteristic times are roughly 0.4–19 h. |
| Local free-antibody loss `k_clear` | 1e-8–2e-5 | s⁻¹ | Non-target local loss | Scenario range following the competition between transport and clearance in Thurber et al. It is not systemic clearance. |
| ECM attenuation `beta_E` | 0.2–2.0 | dimensionless | Maps bounded ECM score to diffusivity | Davies et al. and Netti et al., *Cancer Research* 2000, 60:2497–2503, support matrix-dependent hindrance. The mapping coefficient remains uncalibrated. |
| CAF attenuation `beta_F` | 0–1.5 | dimensionless | Maps bounded CAF score to diffusivity | Lu et al., *Nature Biotechnology* 2026, DOI: 10.1038/s41587-026-03152-x, reports FAP-positive CAF/ECM niches associated with reduced delivery. This is not a measured diffusion coefficient. |
| Vessel concentration `c_v` | user-stated; simulation default 20 | nM | Boundary condition | Scenario input, not a systemic PK prediction or a dose-to-tumour conversion. |
| Engagement threshold `tau` | user-stated; simulation default 0.5 | bound fraction | Reachability decision | Operational threshold; always emitted and sensitivity-tested. |

## Excluded defaults

No transcript-to-surface-protein conversion and no universal dose-to-interstitial-concentration
conversion are supplied. A real run missing antigen calibration must abstain. A stated vessel boundary
without perfusion evidence remains conditional, not verified delivery.

## Version 0.7 measured transport sensitivity prior

The committed `results/external_validation/igg_transport_prior/igg_transport_prior.json` narrows the
IgG diffusion scenario to a central 8.7–19.7 µm²/s range and broad 5.4–31.2 µm²/s envelope derived from
four Netti tumour FRAP measurements. Davies' tumour/free ratio of 0.30–0.50 is retained as an independent
relative check. These values refine sensitivity analysis; they do not override the requirement for a
matched tissue measurement.
