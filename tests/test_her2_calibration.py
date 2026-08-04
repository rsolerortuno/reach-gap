from __future__ import annotations

import zipfile
from pathlib import Path

import pandas as pd

from reach_gap.her2_calibration import (
    fit_source_protocol_calibration,
    load_quantitative_her2_workbook,
)


def _write_minimal_xlsx(path: Path) -> None:
    workbook = (
        '<?xml version="1.0"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets>'
        "</workbook>"
    )
    relationships = (
        '<?xml version="1.0"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Target="worksheets/sheet1.xml" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"/>'
        "</Relationships>"
    )
    headers = [
        "Replicate 1 Tumor ID",
        "Replicate 1 Raw Cy5 MFI",
        "Replicate 1 Result (HER2 receptors/cell)",
        "Replicate 2 Tumor ID",
        "Replicate 2 Raw Cy5 MFI",
        "Replicate 2 Result (HER2 receptors/cell)",
        "CAP/CLIA HER2 Score",
    ]
    rows = [headers]
    for index, raw in enumerate((5.0, 10.0, 20.0), start=1):
        receptors = 10_000.0 * raw**0.8
        rows.append(
            [f"T{index}", raw, receptors, f"T{index}", raw * 1.1, receptors * 1.05, index - 1]
        )

    def cell(column: int, row: int, value: object) -> str:
        reference = f"{chr(ord('A') + column)}{row}"
        if isinstance(value, str):
            return f'<c r="{reference}" t="inlineStr"><is><t>{value}</t></is></c>'
        return f'<c r="{reference}"><v>{value}</v></c>'

    sheet_rows = []
    for row_index, values in enumerate(rows, start=1):
        sheet_rows.append(
            f'<row r="{row_index}">'
            + "".join(cell(i, row_index, value) for i, value in enumerate(values))
            + "</row>"
        )
    worksheet = (
        '<?xml version="1.0"?><worksheet '
        'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        "<sheetData>" + "".join(sheet_rows) + "</sheetData></worksheet>"
    )
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", relationships)
        archive.writestr("xl/worksheets/sheet1.xml", worksheet)


def test_quantitative_her2_xlsx_and_calibration(tmp_path: Path) -> None:
    path = tmp_path / "calibration.xlsx"
    _write_minimal_xlsx(path)
    table = load_quantitative_her2_workbook(path)
    assert len(table) == 3
    model, pairs = fit_source_protocol_calibration(table)
    assert len(pairs) == 6
    assert model.log_log_r2 > 0.95
    assert model.predict(model.raw_mfi_min) is not None
    assert model.predict(model.raw_mfi_max * 2) is None


def test_calibration_requires_pairs() -> None:
    table = pd.DataFrame(
        {
            "Replicate 1 Raw Cy5 MFI": [1.0],
            "Replicate 1 Result (HER2 receptors/cell)": [10_000.0],
        }
    )
    try:
        fit_source_protocol_calibration(table)
    except ValueError as exc:
        assert "four" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("Expected ValueError")


def test_full_her2_calibration_benchmark(tmp_path: Path) -> None:
    from reach_gap.her2_calibration import benchmark_her2_receptor_calibration

    path = tmp_path / "calibration.xlsx"
    _write_minimal_xlsx(path)
    result = benchmark_her2_receptor_calibration(path, tmp_path / "out")
    assert result["status"].startswith("SOURCE_PROTOCOL")
    assert result["uncensored_replicate_pairs"] == 6
    assert result["xenium_transfer"]["status"] == "NOT_COMPUTED"
    assert (tmp_path / "out" / "her2_receptor_calibration.json").exists()
