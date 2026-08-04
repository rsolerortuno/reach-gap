"""Generate the committed deterministic reference fields, index and claims."""

from __future__ import annotations

from pathlib import Path

from reach_gap.features import extract_features, save_features
from reach_gap.geometry import GeometryConfig, simulate_geometry
from reach_gap.indexing import compute_index, write_index_outputs
from reach_gap.schemas import ModelParameters
from reach_gap.solver import save_solution, solve_transport_robust

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "simulated"


def main() -> None:
    """Generate one calibrated simulation and propagate its uncertainty."""

    geometry = simulate_geometry(
        GeometryConfig(
            size=42,
            cell_count=900,
            vessel_count=5,
            stroma_level=0.65,
            antigen_level=0.80,
            seed=17,
        )
    )
    features = extract_features(geometry, antigen_calibration_nM_per_signal=300.0)
    RESULTS.mkdir(parents=True, exist_ok=True)
    save_features(features, RESULTS / "features.npz")

    parameters = ModelParameters()
    solution = solve_transport_robust(features, parameters)
    if not solution.converged:
        raise RuntimeError("Committed reference solve did not converge")
    save_solution(solution, RESULTS / "solution.npz")

    output, claims = compute_index(features, parameters, draws=48, seed=17)
    write_index_outputs(output, claims, RESULTS)


if __name__ == "__main__":
    main()
