from __future__ import annotations

from pathlib import Path

import pytest

from reach_gap.visuals import write_horizontal_bar_svg


def test_write_horizontal_bar_svg(tmp_path: Path) -> None:
    path = write_horizontal_bar_svg(
        ["one", "two"],
        [0.25, 0.75],
        tmp_path / "chart.svg",
        title="Test chart",
        x_label="Fraction",
    )
    text = path.read_text(encoding="utf-8")
    assert "Test chart" in text
    assert "75.0%" in text
    assert text.endswith("</svg>\n")


@pytest.mark.parametrize(
    ("labels", "values", "maximum"),
    [([], [], 1.0), (["one"], [], 1.0), (["one"], [1.1], 1.0), (["one"], [0.1], 0.0)],
)
def test_write_horizontal_bar_svg_rejects_invalid_input(
    tmp_path: Path,
    labels: list[str],
    values: list[float],
    maximum: float,
) -> None:
    with pytest.raises(ValueError):
        write_horizontal_bar_svg(
            labels,
            values,
            tmp_path / "invalid.svg",
            title="Invalid",
            x_label="x",
            maximum=maximum,
        )
