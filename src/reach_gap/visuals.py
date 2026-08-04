"""Dependency-light SVG figures for compact release artifacts."""

from __future__ import annotations

from collections.abc import Sequence
from html import escape
from pathlib import Path


def write_horizontal_bar_svg(
    labels: Sequence[str],
    values: Sequence[float],
    path: Path,
    *,
    title: str,
    x_label: str,
    maximum: float = 1.0,
    value_format: str = ".1%",
) -> Path:
    """Write an accessible monochrome horizontal bar chart as SVG."""

    if len(labels) != len(values):
        raise ValueError("labels and values must have equal length")
    if not labels:
        raise ValueError("at least one bar is required")
    if maximum <= 0.0:
        raise ValueError("maximum must be positive")
    if any(value < 0.0 or value > maximum for value in values):
        raise ValueError("values must lie between zero and maximum")

    width = 960
    label_width = 270
    plot_width = 560
    row_height = 44
    top = 90
    bottom = 80
    height = top + row_height * len(labels) + bottom
    bar_height = 24

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        f'<title id="title">{escape(title)}</title>',
        f'<desc id="desc">{escape(x_label)} shown as horizontal bars.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:.0f}" y="38" text-anchor="middle" '
        'font-family="sans-serif" font-size="24">'
        f"{escape(title)}</text>",
    ]
    x0 = label_width
    for tick in range(5):
        fraction = tick / 4
        x = x0 + plot_width * fraction
        value = maximum * fraction
        lines.extend(
            [
                f'<line x1="{x:.1f}" y1="65" x2="{x:.1f}" y2="{height - bottom + 8}" '
                'stroke="#d0d0d0" stroke-width="1"/>',
                f'<text x="{x:.1f}" y="{height - bottom + 30}" text-anchor="middle" '
                'font-family="sans-serif" font-size="13">'
                f"{format(value, value_format)}</text>",
            ]
        )

    for index, (label, value) in enumerate(zip(labels, values, strict=True)):
        y = top + index * row_height
        bar_width = plot_width * value / maximum
        lines.extend(
            [
                f'<text x="{x0 - 12}" y="{y + 17}" text-anchor="end" '
                'font-family="sans-serif" font-size="15">'
                f"{escape(label)}</text>",
                f'<rect x="{x0}" y="{y}" width="{bar_width:.2f}" height="{bar_height}" '
                'fill="#404040"/>',
                f'<text x="{x0 + bar_width + 8:.2f}" y="{y + 17}" '
                'font-family="sans-serif" font-size="14">'
                f"{escape(format(value, value_format))}</text>",
            ]
        )

    lines.extend(
        [
            f'<text x="{x0 + plot_width / 2:.1f}" y="{height - 18}" text-anchor="middle" '
            'font-family="sans-serif" font-size="15">'
            f"{escape(x_label)}</text>",
            "</svg>",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
