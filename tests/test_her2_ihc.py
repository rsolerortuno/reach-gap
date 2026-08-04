from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PIL import Image

from reach_gap.her2_ihc import benchmark_her2_ihc, image_features, load_her2_cases


def test_image_features_detect_brown_signal() -> None:
    rgb = np.full((10, 10, 3), 255, dtype=np.uint8)
    rgb[:, :5] = np.array([140, 90, 45], dtype=np.uint8)
    mask = np.ones((10, 10), dtype=bool)
    observed = image_features(rgb, mask)
    assert observed["brown_q90"] > 0.0
    assert observed["mask_fraction"] == 1.0


def _write_case(root: Path, case: str, score: int, brown: int, mask_type: str) -> None:
    image_dir = root / "images" / case
    mask_dir = root / "masks" / case
    image_dir.mkdir(parents=True, exist_ok=True)
    mask_dir.mkdir(parents=True, exist_ok=True)
    image = np.full((20, 20, 3), 255, dtype=np.uint8)
    image[:, :] = np.array([180 - brown, 150 - brown // 2, 100 - brown], dtype=np.uint8)
    mask = np.full((20, 20), 255, dtype=np.uint8)
    image_path = image_dir / f"ID{case}_HER2_1.jpg"
    mask_path = mask_dir / f"ID{case}_HER2_1_mask.png"
    Image.fromarray(image).save(image_path)
    Image.fromarray(mask).save(mask_path)
    payload = {
        "case_id": case,
        "diagnostic_assessment": {
            "HER2/neu": {
                "positive_cell_percentage": 50 if score else 0,
                "her2_score": score,
                "result": "positive" if score >= 2 else "negative",
            }
        },
        "images": [
            {
                "image_id": f"ID{case}_HER2",
                "stain_type": "IHC",
                "marker_name": "HER2/neu",
                "file_path": str(image_path.relative_to(root)),
                "mask_path": str(mask_path.relative_to(root)),
                "mask_type": mask_type,
                "magnification": "x10",
            }
        ],
    }
    (root / f"ID{case}.json").write_text(json.dumps(payload), encoding="utf-8")


def test_her2_mask_semantics_and_abstention(tmp_path: Path) -> None:
    _write_case(tmp_path, "0001", 0, 0, "tumor_cells")
    _write_case(tmp_path, "0002", 3, 60, "only_positive_tumor_cells")
    table = load_her2_cases(tmp_path, downsample=1)
    assert table["mask_is_denominator_valid"].tolist() == [True, False]
    result = benchmark_her2_ihc(tmp_path, tmp_path / "out", downsample=1)
    assert result["status"] == "IHC_SCORE_BENCHMARK_ONLY_NOT_MOLECULE_CALIBRATION"
    assert result["absolute_antigen_density"]["status"] == "NOT_COMPUTED"


def test_max_statistic_permutation_pvalue_is_exact_and_selection_aware() -> None:
    import pandas as pd

    from reach_gap.her2_ihc import max_statistic_permutation_pvalue

    table = pd.DataFrame(
        {
            "her2_score": [0, 0, 1, 2],
            "feature_a": [0.0, 0.1, 0.8, 1.0],
            "feature_b": [0.1, 0.0, 0.7, 0.9],
        }
    )
    pvalue, evaluated = max_statistic_permutation_pvalue(table, ["feature_a", "feature_b"])
    assert evaluated == 12
    assert 0.0 < pvalue <= 1.0
