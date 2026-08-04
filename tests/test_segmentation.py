from __future__ import annotations

import numpy as np
import pandas as pd

from reach_gap.segmentation import _polygon_metrics, segmentation_robustness_summary


def test_polygon_metrics_square() -> None:
    observed = _polygon_metrics("a", np.array([0.0, 1.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0, 1.0]))
    assert observed["boundary_area"] == 1.0
    assert observed["boundary_perimeter"] == 4.0
    assert np.isclose(observed["boundary_compactness"], np.pi / 4.0)


def test_segmentation_summary_exact_area() -> None:
    cells = pd.DataFrame({"cell_id": ["a"], "cell_area": [1.0], "nucleus_area": [1.0]})
    metrics = pd.DataFrame(
        {
            "cell_id": ["a"],
            "boundary_vertices": [4],
            "boundary_area": [1.0],
            "boundary_perimeter": [4.0],
            "boundary_compactness": [np.pi / 4.0],
        }
    )
    _, summary = segmentation_robustness_summary(cells, metrics, metrics)
    assert summary["cell_median_relative_area_error"] == 0.0
    assert summary["nucleus_within_1pct"] == 1.0
