"""Source-protocol HER2 receptor calibration from a quantitative-IHC workbook."""

from __future__ import annotations

import json
import math
import re
import time
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, cast
from xml.etree import ElementTree

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

_XLSX_NAMESPACE = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


@dataclass(frozen=True)
class Her2CalibrationModel:
    """Empirical relation valid only for the source Cy5 quantitative-IHC protocol."""

    raw_mfi_min: float
    raw_mfi_max: float
    receptor_lloq: float
    receptor_uloq: float
    log10_slope: float
    log10_intercept: float
    log_log_r2: float

    def predict(self, raw_cy5_mfi: float) -> float | None:
        """Predict receptor copies only inside the observed source-protocol MFI range."""

        if raw_cy5_mfi < self.raw_mfi_min or raw_cy5_mfi > self.raw_mfi_max:
            return None
        prediction = 10.0 ** (self.log10_slope * math.log10(raw_cy5_mfi) + self.log10_intercept)
        if prediction < self.receptor_lloq or prediction > self.receptor_uloq:
            return None
        return float(prediction)


def _column_index(reference: str) -> int:
    letters = "".join(character for character in reference if character.isalpha())
    value = 0
    for character in letters:
        value = value * 26 + (ord(character.upper()) - ord("A") + 1)
    return value - 1


def _shared_strings(archive: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
    strings: list[str] = []
    for item in root.findall("a:si", _XLSX_NAMESPACE):
        strings.append("".join(node.text or "" for node in item.findall(".//a:t", _XLSX_NAMESPACE)))
    return strings


def read_first_xlsx_sheet(path: Path) -> list[list[Any]]:
    """Read the first worksheet using only the XLSX ZIP/XML standard library."""

    with zipfile.ZipFile(path) as archive:
        workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
        relationships = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        relation_map = {
            relation.attrib["Id"]: relation.attrib["Target"] for relation in relationships
        }
        first_sheet = workbook.find("a:sheets/a:sheet", _XLSX_NAMESPACE)
        if first_sheet is None:
            raise ValueError("Workbook contains no worksheets")
        relationship_id = first_sheet.attrib[
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
        ]
        target = relation_map[relationship_id]
        worksheet_name = target.lstrip("/")
        if not worksheet_name.startswith("xl/"):
            worksheet_name = f"xl/{worksheet_name}"
        worksheet = ElementTree.fromstring(archive.read(worksheet_name))
        strings = _shared_strings(archive)

    rows: list[list[Any]] = []
    for row in worksheet.findall(".//a:sheetData/a:row", _XLSX_NAMESPACE):
        cells: dict[int, Any] = {}
        for cell in row.findall("a:c", _XLSX_NAMESPACE):
            reference = cell.attrib.get("r", "A1")
            cell_type = cell.attrib.get("t")
            value_node = cell.find("a:v", _XLSX_NAMESPACE)
            inline_node = cell.find("a:is/a:t", _XLSX_NAMESPACE)
            value: Any = None
            if inline_node is not None:
                value = inline_node.text or ""
            elif value_node is not None:
                raw = value_node.text or ""
                if cell_type == "s":
                    value = strings[int(raw)]
                else:
                    try:
                        value = float(raw)
                    except ValueError:
                        value = raw
            cells[_column_index(reference)] = value
        if cells:
            width = max(cells) + 1
            rows.append([cells.get(index) for index in range(width)])
    return rows


def _normalise_header(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def _numeric(value: Any) -> float | None:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        number = float(value)
        return number if math.isfinite(number) else None
    text = str(value or "").replace(",", "").strip()
    if not text or text.casefold() in {"nc", "nan", "none"} or text.startswith(("<", ">")):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def load_quantitative_her2_workbook(path: Path) -> pd.DataFrame:
    """Load replicate-level Cy5 and receptor-copy observations from the workbook."""

    rows = read_first_xlsx_sheet(path)
    header_index = next(
        (
            index
            for index, row in enumerate(rows)
            if any(_normalise_header(value) == "Replicate 1 Tumor ID" for value in row)
        ),
        None,
    )
    if header_index is None:
        raise ValueError("Could not locate the quantitative HER2 table header")
    headers = [_normalise_header(value) for value in rows[header_index]]
    records: list[dict[str, Any]] = []
    for row in rows[header_index + 1 :]:
        padded = row + [None] * max(0, len(headers) - len(row))
        record = {header: padded[index] for index, header in enumerate(headers) if header}
        if not record.get("Replicate 1 Tumor ID"):
            continue
        records.append(record)
    return pd.DataFrame(records)


def _resolve_column(table: pd.DataFrame, *tokens: str) -> str:
    normalised_tokens = [re.sub(r"[^a-z0-9]+", "", token.casefold()) for token in tokens]
    for column in table.columns:
        normalised = re.sub(r"[^a-z0-9]+", "", str(column).casefold())
        if all(token in normalised for token in normalised_tokens):
            return str(column)
    raise KeyError(f"Could not resolve a column containing {tokens}")


def fit_source_protocol_calibration(
    table: pd.DataFrame,
) -> tuple[Her2CalibrationModel, pd.DataFrame]:
    """Fit a log-log empirical relation on uncensored replicate observations."""

    rows: list[dict[str, Any]] = []
    try:
        score_column: str | None = _resolve_column(table, "CAP", "CLIA", "HER2", "Score")
    except KeyError:
        score_column = None
    for replicate in (1, 2):
        try:
            raw_column = _resolve_column(table, f"Replicate {replicate}", "Raw Cy5 MFI")
            result_column = _resolve_column(
                table, f"Replicate {replicate}", "Result", "HER2 receptors"
            )
            id_column = _resolve_column(table, f"Replicate {replicate}", "Tumor ID")
        except KeyError:
            continue
        for _, record in table.iterrows():
            raw = _numeric(record.get(raw_column))
            receptors = _numeric(record.get(result_column))
            if raw is None or receptors is None or raw <= 0 or receptors <= 0:
                continue
            rows.append(
                {
                    "replicate": replicate,
                    "tumor_id": str(record.get(id_column)),
                    "raw_cy5_mfi": raw,
                    "her2_receptors_per_cell": receptors,
                    "cap_clia_her2_score": _numeric(record.get(score_column))
                    if score_column is not None
                    else None,
                }
            )
    pairs = pd.DataFrame(rows)
    if len(pairs) < 4:
        raise ValueError("At least four uncensored replicate pairs are required")
    raw_values = pairs["raw_cy5_mfi"].to_numpy(dtype=np.float64)
    receptor_values = pairs["her2_receptors_per_cell"].to_numpy(dtype=np.float64)
    slope, intercept = np.polyfit(np.log10(raw_values), np.log10(receptor_values), 1)
    fitted = 10.0 ** (slope * np.log10(raw_values) + intercept)
    residual = receptor_values - fitted
    total = receptor_values - np.mean(receptor_values)
    r2 = 1.0 - float(np.sum(residual**2) / np.sum(total**2))
    model = Her2CalibrationModel(
        raw_mfi_min=float(np.min(raw_values)),
        raw_mfi_max=float(np.max(raw_values)),
        receptor_lloq=10_375.0,
        receptor_uloq=178_649.0,
        log10_slope=float(slope),
        log10_intercept=float(intercept),
        log_log_r2=r2,
    )
    return model, pairs


def _plot_calibration(
    pairs: pd.DataFrame, model: Her2CalibrationModel, output_dir: Path
) -> list[str]:
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover
        return []
    output_dir.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(7, 5))
    plt.scatter(pairs["raw_cy5_mfi"], pairs["her2_receptors_per_cell"])
    grid = np.geomspace(model.raw_mfi_min, model.raw_mfi_max, 200)
    prediction = 10.0 ** (model.log10_slope * np.log10(grid) + model.log10_intercept)
    plt.plot(grid, prediction)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Raw Cy5 MFI (source protocol)")
    plt.ylabel("HER2 receptors per cell")
    plt.title("Source-protocol quantitative HER2 calibration")
    plt.tight_layout()
    path = output_dir / "her2_source_protocol_calibration.png"
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return [str(path)]


def benchmark_her2_receptor_calibration(workbook_path: Path, output_dir: Path) -> dict[str, Any]:
    """Audit and fit the published quantitative-IHC receptor calibration table."""

    started = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)
    table = load_quantitative_her2_workbook(workbook_path)
    model, pairs = fit_source_protocol_calibration(table)
    correlation = cast(
        tuple[float, float],
        spearmanr(pairs["raw_cy5_mfi"], pairs["her2_receptors_per_cell"]),
    )
    correlation_statistic = correlation[0]
    correlation_pvalue = correlation[1]
    pairs.to_csv(output_dir / "her2_uncensored_replicate_pairs.csv", index=False)

    score_rows: list[dict[str, Any]] = []
    for score, group in pairs.dropna(subset=["cap_clia_her2_score"]).groupby(
        "cap_clia_her2_score", sort=True
    ):
        numeric_score = _numeric(score)
        if numeric_score is None:
            continue
        score_rows.append(
            {
                "cap_clia_her2_score": int(numeric_score),
                "observations": len(group),
                "median_receptors_per_cell": float(group["her2_receptors_per_cell"].median()),
                "minimum_receptors_per_cell": float(group["her2_receptors_per_cell"].min()),
                "maximum_receptors_per_cell": float(group["her2_receptors_per_cell"].max()),
            }
        )
    pd.DataFrame(score_rows).to_csv(output_dir / "her2_score_group_summary.csv", index=False)
    figures = _plot_calibration(pairs, model, output_dir / "figures")

    result: dict[str, Any] = {
        "status": "SOURCE_PROTOCOL_HER2_RECEPTOR_CALIBRATION_NOT_XENIUM_TRANSFER",
        "tumours_in_workbook": len(table),
        "uncensored_replicate_pairs": len(pairs),
        "raw_cy5_receptor_spearman": correlation_statistic,
        "raw_cy5_receptor_spearman_pvalue": correlation_pvalue,
        "model": asdict(model),
        "xenium_transfer": {
            "status": "NOT_COMPUTED",
            "reasons": [
                "The calibration is specific to the source Cy5 quantitative-IHC workflow",
                "No shared calibrator was imaged in both the source assay and Xenium protein assay",
                "CAP/CLIA score categories are not receptor-copy measurements for the RCC section",
            ],
        },
        "runtime_seconds": time.time() - started,
        "figures": figures,
    }
    (output_dir / "her2_receptor_calibration.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result
