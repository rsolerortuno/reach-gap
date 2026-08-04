"""Literature-derived IgG transport priors with explicit transfer limitations."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class TransportObservation:
    """One reported tumour or in-vitro transport observation."""

    source: str
    system: str
    metric: str
    central_value: float
    lower_value: float | None
    upper_value: float | None
    unit: str
    note: str


OBSERVATIONS: tuple[TransportObservation, ...] = (
    TransportObservation(
        "Netti_2000",
        "MCaIV tumour",
        "IgG diffusion coefficient",
        1.97,
        1.24,
        3.12,
        "1e-7 cm2/s",
        "FRAP; reported interval reproduced from the source table",
    ),
    TransportObservation(
        "Netti_2000",
        "LS174T tumour",
        "IgG diffusion coefficient",
        1.89,
        1.29,
        2.77,
        "1e-7 cm2/s",
        "FRAP; reported interval reproduced from the source table",
    ),
    TransportObservation(
        "Netti_2000",
        "U87 tumour",
        "IgG diffusion coefficient",
        0.87,
        0.73,
        1.03,
        "1e-7 cm2/s",
        "FRAP; reported interval reproduced from the source table",
    ),
    TransportObservation(
        "Netti_2000",
        "HSTS26T tumour",
        "IgG diffusion coefficient",
        0.96,
        0.54,
        1.71,
        "1e-7 cm2/s",
        "FRAP; reported interval reproduced from the source table",
    ),
    TransportObservation(
        "Davies_2002",
        "rhabdomyosarcoma tumour clones",
        "D_tissue / D_free",
        0.40,
        0.30,
        0.50,
        "ratio",
        "Range across tumour clones; collagen and sulphated GAG associated with lower diffusion",
    ),
    TransportObservation(
        "Ramanujan_2002",
        "high-hyaluronan collagen gel",
        "IgG D_gel / D_free",
        0.56,
        0.45,
        0.67,
        "ratio",
        "Reported mean ± spread represented as a conservative interval",
    ),
)


def cm2_s_1e7_to_um2_s(value: float) -> float:
    """Convert a value expressed in 1e-7 cm²/s into µm²/s."""

    return round(value * 10.0, 10)


def build_igg_transport_prior(output_dir: Path) -> dict[str, Any]:
    """Write a conservative external-transfer prior for an IgG-sized molecule."""

    output_dir.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(asdict(observation) for observation in OBSERVATIONS)
    table.to_csv(output_dir / "transport_literature_observations.csv", index=False)

    netti = [observation for observation in OBSERVATIONS if observation.source == "Netti_2000"]
    central_um2_s = [cm2_s_1e7_to_um2_s(observation.central_value) for observation in netti]
    broad_low = min(
        cm2_s_1e7_to_um2_s(observation.lower_value)
        for observation in netti
        if observation.lower_value is not None
    )
    broad_high = max(
        cm2_s_1e7_to_um2_s(observation.upper_value)
        for observation in netti
        if observation.upper_value is not None
    )
    result: dict[str, Any] = {
        "status": "EXTERNAL_IGG_TRANSPORT_PRIOR_NOT_RCC_MEASUREMENT",
        "absolute_diffusion_um2_s": {
            "central_low": min(central_um2_s),
            "central_high": max(central_um2_s),
            "broad_low": broad_low,
            "broad_high": broad_high,
            "basis": "Four tumour FRAP measurements reported by Netti et al.",
        },
        "relative_diffusion_ratio": {
            "low": 0.30,
            "high": 0.50,
            "basis": "Tumour-to-free diffusion range reported across Davies et al. tumour clones",
        },
        "shg_mapping": {
            "status": "NOT_COMPUTED",
            "reasons": [
                "SHG intensity has no universal conversion to collagen concentration",
                "The FRAP studies and SHG benchmark do not measure the same tissue sections",
                "RCC-specific porosity, glycosaminoglycans, pressure and binding were not measured",
            ],
        },
        "recommended_solver_use": {
            "distribution": "loguniform",
            "low_um2_s": broad_low,
            "high_um2_s": broad_high,
            "claim_scope": "sensitivity prior only",
        },
    }
    (output_dir / "igg_transport_prior.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result
