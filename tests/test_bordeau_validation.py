from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
from PIL import Image

from reach_gap.bordeau_validation import extract_docx_media, split_three_by_three_figure


def test_split_three_by_three_figure() -> None:
    canvas = np.full((94, 124, 3), 255, dtype=np.uint8)
    for row in range(3):
        for column in range(3):
            y0 = 2 + row * 31
            x0 = 2 + column * 41
            canvas[y0 : y0 + 28, x0 : x0 + 38] = np.array([10, 20, 30], dtype=np.uint8)
    panels = split_three_by_three_figure(canvas)
    assert len(panels) == 3
    assert all(len(row) == 3 for row in panels)
    assert panels[0][0].shape[:2] == (28, 38)


def test_extract_docx_media(tmp_path: Path) -> None:
    image_path = tmp_path / "image.png"
    Image.fromarray(np.zeros((10, 10, 3), dtype=np.uint8)).save(image_path)
    docx = tmp_path / "supplement.docx"
    with zipfile.ZipFile(docx, "w") as archive:
        archive.write(image_path, "word/media/image1.png")
    extracted = extract_docx_media(docx, tmp_path / "media")
    assert len(extracted) == 1
    assert extracted[0].exists()


def _synthetic_figure(combo: bool) -> np.ndarray:
    canvas = np.full((94, 124, 3), 255, dtype=np.uint8)
    for row in range(3):
        for column in range(3):
            y0 = 2 + row * 31
            x0 = 2 + column * 41
            panel = canvas[y0 : y0 + 28, x0 : x0 + 38]
            panel[:] = 0
            if row == 0:
                panel[:, 4:6, 0] = 220
                green_column = 25 if combo else 12
                panel[8:20, green_column : green_column + 4, 1] = 220
            elif row == 1:
                panel[8:20, 10:18] = 120
            else:
                panel[:] = 80
                panel[8:20, 10:18] = 0
    return canvas


def test_full_bordeau_benchmark(tmp_path: Path) -> None:
    from reach_gap.bordeau_validation import benchmark_bordeau_supplement

    alone = tmp_path / "image2.png"
    combo = tmp_path / "image3.png"
    Image.fromarray(_synthetic_figure(False)).save(alone)
    Image.fromarray(_synthetic_figure(True)).save(combo)
    docx = tmp_path / "supplement.docx"
    with zipfile.ZipFile(docx, "w") as archive:
        archive.write(alone, "word/media/image2.jpeg")
        archive.write(combo, "word/media/image3.png")
    result = benchmark_bordeau_supplement(docx, tmp_path / "out")
    assert result["status"].startswith("PUBLISHED_ADMINISTERED")
    assert result["model_concordance"]["status"] == "NOT_COMPUTED"
    assert result["representative_figure_analysis"][
        "direction_matches_published_threshold_positive_penetration"
    ]
