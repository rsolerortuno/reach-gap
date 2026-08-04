from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from reach_gap.shg_collagen import benchmark_shg_collagen, shg_features


def test_shg_features_detect_oriented_fibres() -> None:
    oriented = np.zeros((64, 64), dtype=np.uint8)
    oriented[:, ::8] = 255
    isotropic = np.zeros((64, 64), dtype=np.uint8)
    isotropic[::8, ::8] = 255
    oriented_features = shg_features(oriented)
    isotropic_features = shg_features(isotropic)
    assert (
        oriented_features["orientation_coherence_median"]
        > isotropic_features["orientation_coherence_median"]
    )


def test_shg_benchmark_abstains_from_diffusivity(tmp_path: Path) -> None:
    image = np.zeros((32, 32), dtype=np.uint8)
    image[:, ::4] = 100
    path = tmp_path / "sample.tif"
    tifffile.imwrite(path, image)
    result = benchmark_shg_collagen({"tumour": [path]}, tmp_path / "out")
    assert result["transport_coefficient"]["status"] == "NOT_COMPUTED"
    assert (tmp_path / "out" / "shg_collagen_features.csv").exists()
