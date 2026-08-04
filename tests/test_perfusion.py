from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from reach_gap.perfusion import benchmark_perfusion_tiffs, perfusion_profile_from_rgb


def test_perfusion_profile_detects_near_vessel_enrichment() -> None:
    rgb = np.zeros((80, 80, 3), dtype=np.uint8)
    rgb[:, 10:12, 1] = 255
    for column in range(80):
        distance = abs(column - 10)
        rgb[:, column, 2] = max(0, 220 - 8 * distance)
    profile, summary, vessel, _, _ = perfusion_profile_from_rgb(
        rgb, pixel_size_um=1.0, minimum_component_pixels=1
    )
    assert vessel.any()
    assert summary["distance_hoechst_spearman"] < 0
    assert summary["near_to_far_mean_ratio"] > 1
    assert len(profile) == 5


def test_perfusion_benchmark_writes_outputs(tmp_path: Path) -> None:
    rgb = np.zeros((48, 48, 3), dtype=np.uint8)
    rgb[:, 8:10, 1] = 255
    rgb[..., 2] = np.maximum(0, 180 - 5 * np.abs(np.arange(48) - 8))[None, :]
    path = tmp_path / "image.tif"
    tifffile.imwrite(path, rgb, imagej=True, resolution=(10.0, 10.0), metadata={"unit": "micron"})
    result = benchmark_perfusion_tiffs({"upper": path}, tmp_path / "out", downsample=1)
    assert result["all_images_negative_distance_correlation"]
    assert (tmp_path / "out" / "perfusion_distance_profiles.csv").exists()
