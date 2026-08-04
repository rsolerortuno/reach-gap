from __future__ import annotations

from pathlib import Path

from reach_gap.transport_priors import build_igg_transport_prior, cm2_s_1e7_to_um2_s


def test_transport_unit_conversion() -> None:
    assert cm2_s_1e7_to_um2_s(1.0) == 10.0


def test_transport_prior_is_external_and_bounded(tmp_path: Path) -> None:
    result = build_igg_transport_prior(tmp_path)
    diffusion = result["absolute_diffusion_um2_s"]
    assert diffusion["broad_low"] == 5.4
    assert diffusion["broad_high"] == 31.2
    assert result["shg_mapping"]["status"] == "NOT_COMPUTED"
    assert (tmp_path / "igg_transport_prior.json").exists()
