"""H&E and pathology-annotation preparation for geometry-only reach-gap runs.

The outputs from this module are deliberately not model-ready target accessibility
features. H&E morphology can support tissue geometry and quality-control priors, but
it cannot establish target positivity, surface antigen density, functional perfusion,
or therapeutic-antibody concentration. The generated claims therefore abstain from
all target-specific reachability metrics.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeAlias, cast

import numpy as np
from numpy.typing import NDArray
from scipy import ndimage

FloatArray: TypeAlias = NDArray[np.float64]
BoolArray: TypeAlias = NDArray[np.bool_]
UInt8Array: TypeAlias = NDArray[np.uint8]


@dataclass(frozen=True)
class OMEImageInfo:
    """Minimal auditable OME-TIFF geometry metadata."""

    full_height_px: int
    full_width_px: int
    channels: int
    levels: int
    physical_size_x_um: float
    physical_size_y_um: float
    selected_level: int
    selected_height_px: int
    selected_width_px: int

    @property
    def selected_pixel_size_x_um(self) -> float:
        """Physical pixel width at the selected pyramid level."""

        return self.physical_size_x_um * self.full_width_px / self.selected_width_px

    @property
    def selected_pixel_size_y_um(self) -> float:
        """Physical pixel height at the selected pyramid level."""

        return self.physical_size_y_um * self.full_height_px / self.selected_height_px


@dataclass(frozen=True)
class LumenCandidate:
    """One H&E-derived low-confidence intratissue lumen candidate."""

    label: int
    centroid_x_px: float
    centroid_y_px: float
    area_px: int
    equivalent_diameter_um: float
    circularity: float
    interior_brightness: float
    ring_hematoxylin: float
    confidence: float


def sha256_file(path: Path, *, chunk_size: int = 16 * 1024 * 1024) -> str:
    """Compute a streaming SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path, *, chunk_size: int = 16 * 1024 * 1024) -> str:
    """Compute a streaming MD5 digest for provider transfer verification."""

    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _find_part(parts_dir: Path, expected_name: str) -> Path:
    exact = parts_dir / expected_name
    if exact.exists():
        return exact
    matches = sorted(parts_dir.glob(f"{expected_name}.*"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one split part for {expected_name}; found {len(matches)}: {matches}"
        )
    return matches[0]


def reassemble_split_file(
    manifest_path: Path,
    parts_dir: Path,
    output_path: Path,
    *,
    verify_part_sha256: bool = True,
    verify_source_md5: bool = True,
) -> dict[str, Any]:
    """Reassemble a bytewise split file with per-part and source verification."""

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    parts = sorted(manifest["parts"], key=lambda item: int(item["index"]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    source_md5 = hashlib.md5()
    source_sha256 = hashlib.sha256()
    part_records: list[dict[str, Any]] = []
    with output_path.open("wb") as target:
        for expected_index, part in enumerate(parts, start=1):
            if int(part["index"]) != expected_index:
                raise ValueError("Split manifest indices are not contiguous and one-based")
            path = _find_part(parts_dir, str(part["name"]))
            observed_size = path.stat().st_size
            if observed_size != int(part["size"]):
                raise ValueError(
                    f"Split part size mismatch for {path.name}: {observed_size} != {part['size']}"
                )
            part_sha256 = hashlib.sha256()
            with path.open("rb") as source:
                while chunk := source.read(8 * 1024 * 1024):
                    part_sha256.update(chunk)
                    source_md5.update(chunk)
                    source_sha256.update(chunk)
                    target.write(chunk)
            observed_part_sha256 = part_sha256.hexdigest()
            if verify_part_sha256 and observed_part_sha256 != str(part["sha256"]):
                raise ValueError(f"Split part SHA-256 mismatch for {path.name}")
            part_records.append(
                {
                    "index": expected_index,
                    "path": str(path),
                    "bytes": observed_size,
                    "sha256": observed_part_sha256,
                }
            )
    expected_size = int(manifest["source"]["size"])
    observed_source_md5 = source_md5.hexdigest()
    expected_source_md5 = str(manifest["source"].get("md5Checksum", ""))
    if output_path.stat().st_size != expected_size:
        raise ValueError(
            f"Reassembled size mismatch: {output_path.stat().st_size} != {expected_size}"
        )
    if verify_source_md5 and expected_source_md5 and observed_source_md5 != expected_source_md5:
        raise ValueError("Reassembled source MD5 does not match the provider checksum")
    return {
        "status": "VERIFIED",
        "output": str(output_path),
        "bytes": output_path.stat().st_size,
        "md5": observed_source_md5,
        "sha256": source_sha256.hexdigest(),
        "parts": part_records,
    }


def _physical_sizes_from_ome(ome_xml: str | None) -> tuple[float, float]:
    if not ome_xml:
        raise ValueError("OME-TIFF has no OME-XML metadata")
    from xml.etree import ElementTree

    root = ElementTree.fromstring(ome_xml)
    namespace = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    pixels = root.find(".//ome:Pixels", namespace)
    if pixels is None:
        raise ValueError("OME-XML has no Pixels element")
    x = pixels.attrib.get("PhysicalSizeX")
    y = pixels.attrib.get("PhysicalSizeY")
    x_unit = pixels.attrib.get("PhysicalSizeXUnit")
    y_unit = pixels.attrib.get("PhysicalSizeYUnit")
    if x is None or y is None or x_unit not in {"µm", "um"} or y_unit not in {"µm", "um"}:
        raise ValueError("OME-XML lacks physical pixel sizes in micrometres")
    return float(x), float(y)


def read_ome_pyramid_level(path: Path, *, level: int = 4) -> tuple[UInt8Array, OMEImageInfo]:
    """Read one bounded OME-TIFF pyramid level rather than the full-resolution image."""

    from tifffile import TiffFile

    with TiffFile(path) as tif:
        if not tif.is_ome or not tif.series:
            raise ValueError("Input must be an OME-TIFF with at least one image series")
        series = tif.series[0]
        if level < 0 or level >= len(series.levels):
            raise ValueError(f"Pyramid level {level} outside 0..{len(series.levels) - 1}")
        full = series.levels[0]
        selected = series.levels[level]
        array = np.asarray(selected.asarray())
        if array.ndim != 3 or array.shape[-1] != 3 or array.dtype != np.uint8:
            raise ValueError(f"Expected uint8 RGB image; observed {array.shape} {array.dtype}")
        physical_x, physical_y = _physical_sizes_from_ome(tif.ome_metadata)
        info = OMEImageInfo(
            full_height_px=int(full.shape[0]),
            full_width_px=int(full.shape[1]),
            channels=int(full.shape[2]),
            levels=len(series.levels),
            physical_size_x_um=physical_x,
            physical_size_y_um=physical_y,
            selected_level=level,
            selected_height_px=int(selected.shape[0]),
            selected_width_px=int(selected.shape[1]),
        )
    return np.asarray(array, dtype=np.uint8), info


def load_geojson(path: Path) -> dict[str, Any]:
    """Read a GeoJSON FeatureCollection with basic structural validation."""

    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("type") != "FeatureCollection" or not isinstance(payload.get("features"), list):
        raise ValueError("Pathology annotation must be a GeoJSON FeatureCollection")
    return cast(dict[str, Any], payload)


def _iter_polygon_rings(geometry: Mapping[str, Any]) -> Iterable[Sequence[Sequence[float]]]:
    geometry_type = str(geometry.get("type"))
    coordinates = geometry.get("coordinates")
    if not isinstance(coordinates, Sequence):
        raise ValueError("Annotation geometry coordinates must be a sequence")
    if geometry_type == "Polygon":
        for ring in coordinates:
            if not isinstance(ring, Sequence):
                raise ValueError("Polygon ring must be a sequence")
            yield cast(Sequence[Sequence[float]], ring)
    elif geometry_type == "MultiPolygon":
        for polygon in coordinates:
            if not isinstance(polygon, Sequence):
                raise ValueError("MultiPolygon member must be a sequence")
            for ring in polygon:
                if not isinstance(ring, Sequence):
                    raise ValueError("Polygon ring must be a sequence")
                yield cast(Sequence[Sequence[float]], ring)
    else:
        raise ValueError(f"Unsupported annotation geometry: {geometry_type}")


def rasterize_annotations(
    geojson: Mapping[str, Any],
    *,
    output_shape: tuple[int, int],
    full_shape: tuple[int, int],
) -> tuple[NDArray[np.int16], dict[int, str]]:
    """Rasterize annotation polygons into a compact categorical mask."""

    from PIL import Image, ImageDraw

    height, width = output_shape
    full_height, full_width = full_shape
    if height <= 0 or width <= 0 or full_height <= 0 or full_width <= 0:
        raise ValueError("Image dimensions must be positive")
    scale_x = width / full_width
    scale_y = height / full_height
    names = sorted(
        {
            str(feature.get("properties", {}).get("name", "UNNAMED"))
            for feature in geojson["features"]
        }
    )
    label_by_name = {name: index + 1 for index, name in enumerate(names)}
    name_by_label = {value: key for key, value in label_by_name.items()}
    canvas = Image.new("I", (width, height), 0)
    draw = ImageDraw.Draw(canvas)
    for feature in geojson["features"]:
        name = str(feature.get("properties", {}).get("name", "UNNAMED"))
        label = label_by_name[name]
        geometry = feature["geometry"]
        rings = list(_iter_polygon_rings(geometry))
        if not rings:
            continue
        outer = [(float(x) * scale_x, float(y) * scale_y) for x, y, *_ in rings[0]]
        draw.polygon(outer, fill=label)
        for hole in rings[1:]:
            inner = [(float(x) * scale_x, float(y) * scale_y) for x, y, *_ in hole]
            draw.polygon(inner, fill=0)
    return np.asarray(canvas, dtype=np.int16), name_by_label


def tissue_mask(rgb: UInt8Array) -> BoolArray:
    """Estimate tissue occupancy using optical density and conservative morphology."""

    rgb_float = np.asarray(rgb, dtype=np.float64)
    height, width = rgb.shape[:2]
    border_width = max(4, min(height, width) // 50)
    border = np.concatenate(
        [
            rgb_float[:border_width].reshape(-1, 3),
            rgb_float[-border_width:].reshape(-1, 3),
            rgb_float[:, :border_width].reshape(-1, 3),
            rgb_float[:, -border_width:].reshape(-1, 3),
        ]
    )
    background_rgb = np.median(border, axis=0)
    border_distance = np.linalg.norm(border - background_rgb, axis=1)
    distance_threshold = max(12.0, float(np.quantile(border_distance, 0.95)) + 4.0)
    colour_distance = np.linalg.norm(rgb_float - background_rgb, axis=2)
    raw = colour_distance > distance_threshold
    cleaned = ndimage.binary_closing(raw, iterations=2)
    cleaned = ndimage.binary_opening(cleaned, iterations=1)
    labels_raw, count = cast(tuple[NDArray[Any], int], ndimage.label(cleaned))
    labels = np.asarray(labels_raw, dtype=np.int64)
    if count == 0:
        return np.zeros(raw.shape, dtype=np.bool_)
    sizes = np.bincount(labels.ravel())
    keep = sizes >= max(64, int(0.00005 * raw.size))
    keep[0] = False
    return np.asarray(keep[labels], dtype=np.bool_)


def he_stain_channels(rgb: UInt8Array) -> tuple[FloatArray, FloatArray]:
    """Return normalized hematoxylin and eosin optical-density channels."""

    from skimage.color import rgb2hed

    hed = rgb2hed(np.asarray(rgb, dtype=np.float64) / 255.0)
    hematoxylin = np.maximum(np.asarray(hed[..., 0], dtype=np.float64), 0.0)
    eosin = np.maximum(np.asarray(hed[..., 1], dtype=np.float64), 0.0)
    for channel in (hematoxylin, eosin):
        finite = channel[np.isfinite(channel)]
        high = float(np.quantile(finite, 0.995)) if finite.size else 1.0
        if high > 0:
            channel /= high
        np.clip(channel, 0.0, 1.0, out=channel)
    return hematoxylin, eosin


def _component_perimeter(component: BoolArray) -> int:
    eroded = ndimage.binary_erosion(component)
    return int(np.count_nonzero(component & ~eroded))


def detect_lumen_candidates(
    rgb: UInt8Array,
    tissue: BoolArray,
    hematoxylin: FloatArray,
    *,
    pixel_size_x_um: float,
    pixel_size_y_um: float,
    exclusion_mask: BoolArray | None = None,
) -> tuple[BoolArray, list[LumenCandidate]]:
    """Detect unvalidated bright intratissue spaces that could include vascular lumina.

    This is intentionally a candidate generator, not a vessel classifier. Tears,
    tubules, ducts, adipocyte spaces and processing artefacts can satisfy these rules.
    """

    if rgb.shape[:2] != tissue.shape or tissue.shape != hematoxylin.shape:
        raise ValueError("RGB, tissue and stain arrays must have matching spatial shapes")
    brightness = np.mean(np.asarray(rgb, dtype=np.float64), axis=2) / 255.0
    filled_tissue = ndimage.binary_fill_holes(tissue)
    holes = filled_tissue & ~tissue & (brightness > 0.82)
    if exclusion_mask is not None:
        if exclusion_mask.shape != holes.shape:
            raise ValueError("exclusion_mask has the wrong shape")
        holes &= ~exclusion_mask
    labels_raw, _count = cast(tuple[NDArray[Any], int], ndimage.label(holes))
    labels = np.asarray(labels_raw, dtype=np.int64)
    accepted = np.zeros_like(holes, dtype=np.bool_)
    candidates: list[LumenCandidate] = []
    pixel_area_um2 = pixel_size_x_um * pixel_size_y_um
    height, width = holes.shape
    objects = ndimage.find_objects(labels)
    for label, object_slice in enumerate(objects, start=1):
        if object_slice is None:
            continue
        y_slice, x_slice = object_slice
        y0 = max(0, int(y_slice.start) - 2)
        y1 = min(height, int(y_slice.stop) + 2)
        x0 = max(0, int(x_slice.start) - 2)
        x1 = min(width, int(x_slice.stop) + 2)
        local_labels = labels[y0:y1, x0:x1]
        component = local_labels == label
        area_px = int(np.count_nonzero(component))
        if area_px < 4:
            continue
        area_um2 = area_px * pixel_area_um2
        diameter_um = 2.0 * math.sqrt(area_um2 / math.pi)
        if diameter_um < 12.0 or diameter_um > 350.0:
            continue
        perimeter_px = _component_perimeter(component)
        if perimeter_px == 0:
            continue
        circularity = float(4.0 * math.pi * area_px / (perimeter_px**2))
        if circularity < 0.05:
            continue
        local_tissue = tissue[y0:y1, x0:x1]
        ring = ndimage.binary_dilation(component, iterations=2) & ~component & local_tissue
        if np.count_nonzero(ring) < 4:
            continue
        local_h = hematoxylin[y0:y1, x0:x1]
        local_brightness = brightness[y0:y1, x0:x1]
        ring_h = float(np.mean(local_h[ring]))
        interior_brightness = float(np.mean(local_brightness[component]))
        confidence = float(
            np.clip(
                0.40 * min(circularity / 0.65, 1.0)
                + 0.35 * min(ring_h / 0.30, 1.0)
                + 0.25 * min(max(interior_brightness - 0.80, 0.0) / 0.20, 1.0),
                0.0,
                1.0,
            )
        )
        if confidence < 0.25:
            continue
        y_values, x_values = np.nonzero(component)
        candidates.append(
            LumenCandidate(
                label=label,
                centroid_x_px=float(np.mean(x_values) + x0),
                centroid_y_px=float(np.mean(y_values) + y0),
                area_px=area_px,
                equivalent_diameter_um=diameter_um,
                circularity=circularity,
                interior_brightness=interior_brightness,
                ring_hematoxylin=ring_h,
                confidence=confidence,
            )
        )
        accepted[y0:y1, x0:x1] |= component
    return accepted, candidates


def annotation_area_summary(
    geojson: Mapping[str, Any], *, pixel_size_x_um: float, pixel_size_y_um: float
) -> dict[str, Any]:
    """Calculate exact polygon areas in image-pixel and physical units."""

    from shapely.geometry import shape

    grouped_px: dict[str, float] = defaultdict(float)
    feature_rows: list[dict[str, Any]] = []
    for feature in geojson["features"]:
        name = str(feature.get("properties", {}).get("name", "UNNAMED"))
        geometry = shape(feature["geometry"])
        if not geometry.is_valid:
            raise ValueError(f"Invalid pathology geometry for {name}")
        area_px = float(geometry.area)
        area_mm2 = area_px * pixel_size_x_um * pixel_size_y_um / 1_000_000.0
        grouped_px[name] += area_px
        feature_rows.append(
            {
                "id": str(feature.get("id", "")),
                "name": name,
                "geometry_type": geometry.geom_type,
                "area_px2": area_px,
                "area_mm2": area_mm2,
                "bounds_px": [float(value) for value in geometry.bounds],
            }
        )
    grouped = {
        name: {
            "area_px2": area_px,
            "area_mm2": area_px * pixel_size_x_um * pixel_size_y_um / 1_000_000.0,
        }
        for name, area_px in sorted(grouped_px.items())
    }
    return {"features": feature_rows, "by_name": grouped}


def _array_summary(values: FloatArray, mask: BoolArray) -> dict[str, float | None]:
    selected = values[mask & np.isfinite(values)]
    if selected.size == 0:
        return {"mean": None, "median": None, "q10": None, "q90": None}
    return {
        "mean": float(np.mean(selected)),
        "median": float(np.median(selected)),
        "q10": float(np.quantile(selected, 0.10)),
        "q90": float(np.quantile(selected, 0.90)),
    }


def region_morphology_summary(
    annotation_mask: NDArray[np.int16],
    name_by_label: Mapping[int, str],
    tissue: BoolArray,
    hematoxylin: FloatArray,
    eosin: FloatArray,
    lumen_mask: BoolArray,
    *,
    pixel_area_um2: float,
) -> dict[str, Any]:
    """Summarize stain and lumen-candidate morphology by annotated region."""

    summaries: dict[str, Any] = {}
    labels = [0, *sorted(name_by_label)]
    for label in labels:
        region = annotation_mask == label
        region_name = "UNANNOTATED" if label == 0 else name_by_label[label]
        tissue_region = region & tissue
        region_px = int(np.count_nonzero(region))
        tissue_px = int(np.count_nonzero(tissue_region))
        lumen_px = int(np.count_nonzero(region & lumen_mask))
        summaries[region_name] = {
            "region_pixels_at_analysis_level": region_px,
            "tissue_pixels_at_analysis_level": tissue_px,
            "tissue_fraction": tissue_px / region_px if region_px else None,
            "tissue_area_mm2_estimate": tissue_px * pixel_area_um2 / 1_000_000.0,
            "hematoxylin": _array_summary(hematoxylin, tissue_region),
            "eosin": _array_summary(eosin, tissue_region),
            "lumen_candidate_pixels": lumen_px,
            "lumen_candidate_fraction_of_tissue": lumen_px / tissue_px if tissue_px else None,
        }
    return summaries


def _normalize_to_uint8(values: FloatArray, mask: BoolArray | None = None) -> UInt8Array:
    finite_mask = np.isfinite(values)
    if mask is not None:
        finite_mask &= mask
    selected = values[finite_mask]
    if selected.size == 0:
        return np.zeros(values.shape, dtype=np.uint8)
    low = float(np.quantile(selected, 0.01))
    high = float(np.quantile(selected, 0.99))
    if high <= low:
        return np.zeros(values.shape, dtype=np.uint8)
    scaled = np.clip((values - low) / (high - low), 0.0, 1.0)
    return np.asarray(np.round(255.0 * scaled), dtype=np.uint8)


def _preview(image: Any, *, maximum_dimension: int = 1600, nearest: bool = False) -> Any:
    """Return a compact image copy for committed QC artefacts."""

    from PIL import Image

    preview = image.copy()
    resampling = Image.Resampling.NEAREST if nearest else Image.Resampling.LANCZOS
    preview.thumbnail((maximum_dimension, maximum_dimension), resampling)
    return preview


def write_qc_images(
    output_dir: Path,
    rgb: UInt8Array,
    annotation_mask: NDArray[np.int16],
    name_by_label: Mapping[int, str],
    tissue: BoolArray,
    hematoxylin: FloatArray,
    eosin: FloatArray,
    lumen_mask: BoolArray,
) -> dict[str, str]:
    """Write compact QC images at the bounded analysis pyramid level."""

    from PIL import Image, ImageDraw

    output_dir.mkdir(parents=True, exist_ok=True)
    original_path = output_dir / "he_analysis_level.jpg"
    _preview(Image.fromarray(rgb)).save(original_path, quality=90)
    tissue_path = output_dir / "tissue_mask.png"
    _preview(Image.fromarray(np.asarray(tissue, dtype=np.uint8) * 255), nearest=True).save(
        tissue_path
    )
    h_path = output_dir / "hematoxylin_proxy.png"
    e_path = output_dir / "eosin_proxy.png"
    _preview(Image.fromarray(_normalize_to_uint8(hematoxylin, tissue))).save(h_path)
    _preview(Image.fromarray(_normalize_to_uint8(eosin, tissue))).save(e_path)

    palette = [
        (230, 25, 75),
        (60, 180, 75),
        (255, 225, 25),
        (0, 130, 200),
        (245, 130, 48),
        (145, 30, 180),
        (70, 240, 240),
        (240, 50, 230),
    ]
    overlay = Image.fromarray(rgb).convert("RGBA")
    overlay_array = np.asarray(overlay).copy()
    for index, label in enumerate(sorted(name_by_label)):
        selected = annotation_mask == label
        color = palette[index % len(palette)]
        overlay_array[selected, :3] = np.asarray(color, dtype=np.uint8)
        overlay_array[selected, 3] = 95
    base = np.asarray(Image.fromarray(rgb).convert("RGBA"), dtype=np.float64)
    alpha = overlay_array[..., 3:4].astype(np.float64) / 255.0
    blended = base.copy()
    blended[..., :3] = alpha * overlay_array[..., :3] + (1.0 - alpha) * base[..., :3]
    blended[..., 3] = 255
    annotation_path = output_dir / "pathology_overlay.png"
    _preview(Image.fromarray(np.asarray(blended, dtype=np.uint8))).save(annotation_path)

    lumen_overlay = Image.fromarray(rgb).convert("RGB")
    draw = ImageDraw.Draw(lumen_overlay)
    boundary = lumen_mask & ~ndimage.binary_erosion(lumen_mask)
    y_values, x_values = np.nonzero(boundary)
    for x, y in zip(x_values.tolist(), y_values.tolist(), strict=True):
        draw.point((x, y), fill=(0, 255, 0))
    lumen_path = output_dir / "unvalidated_lumen_candidates.png"
    _preview(lumen_overlay).save(lumen_path)
    return {
        "he_analysis_level": str(original_path),
        "tissue_mask": str(tissue_path),
        "hematoxylin_proxy": str(h_path),
        "eosin_proxy": str(e_path),
        "pathology_overlay": str(annotation_path),
        "unvalidated_lumen_candidates": str(lumen_path),
    }


def _candidate_records(
    candidates: Sequence[LumenCandidate], info: OMEImageInfo
) -> list[dict[str, Any]]:
    return [
        {
            "candidate_id": f"lumen_{index:06d}",
            "classification": "UNVALIDATED_LUMEN_CANDIDATE",
            "x_analysis_px": candidate.centroid_x_px,
            "y_analysis_px": candidate.centroid_y_px,
            "x_full_resolution_px": candidate.centroid_x_px
            * info.full_width_px
            / info.selected_width_px,
            "y_full_resolution_px": candidate.centroid_y_px
            * info.full_height_px
            / info.selected_height_px,
            "equivalent_diameter_um": candidate.equivalent_diameter_um,
            "circularity": candidate.circularity,
            "interior_brightness": candidate.interior_brightness,
            "ring_hematoxylin": candidate.ring_hematoxylin,
            "candidate_confidence": candidate.confidence,
        }
        for index, candidate in enumerate(candidates, start=1)
    ]


def prepare_he_pathology_rcc(
    *,
    he_path: Path,
    annotation_path: Path,
    alignment_path: Path | None,
    output_dir: Path,
    analysis_level: int = 4,
) -> dict[str, Any]:
    """Prepare real RCC H&E geometry while enforcing target-level abstention."""

    output_dir.mkdir(parents=True, exist_ok=True)
    rgb, info = read_ome_pyramid_level(he_path, level=analysis_level)
    geojson = load_geojson(annotation_path)
    annotation_mask, name_by_label = rasterize_annotations(
        geojson,
        output_shape=rgb.shape[:2],
        full_shape=(info.full_height_px, info.full_width_px),
    )
    tissue = tissue_mask(rgb)
    hematoxylin, eosin = he_stain_channels(rgb)
    excluded_names = {"Adipose tissue", "Necrosis", "Hemorrhage"}
    exclusion_mask = np.zeros(tissue.shape, dtype=np.bool_)
    for label, name in name_by_label.items():
        if name in excluded_names:
            exclusion_mask |= annotation_mask == label
    lumen_mask, candidates = detect_lumen_candidates(
        rgb,
        tissue,
        hematoxylin,
        pixel_size_x_um=info.selected_pixel_size_x_um,
        pixel_size_y_um=info.selected_pixel_size_y_um,
        exclusion_mask=exclusion_mask,
    )

    annotation_summary = annotation_area_summary(
        geojson,
        pixel_size_x_um=info.physical_size_x_um,
        pixel_size_y_um=info.physical_size_y_um,
    )
    pixel_area_um2 = info.selected_pixel_size_x_um * info.selected_pixel_size_y_um
    morphology = region_morphology_summary(
        annotation_mask,
        name_by_label,
        tissue,
        hematoxylin,
        eosin,
        lumen_mask,
        pixel_area_um2=pixel_area_um2,
    )
    qc_images = write_qc_images(
        output_dir / "qc",
        rgb,
        annotation_mask,
        name_by_label,
        tissue,
        hematoxylin,
        eosin,
        lumen_mask,
    )
    candidate_rows = _candidate_records(candidates, info)
    candidate_path = output_dir / "unvalidated_lumen_candidates.json"
    candidate_path.write_text(json.dumps(candidate_rows, indent=2), encoding="utf-8")

    alignment: list[list[float]] | None = None
    if alignment_path is not None:
        matrix = np.loadtxt(alignment_path, delimiter=",")
        if matrix.shape != (3, 3) or not np.isfinite(matrix).all():
            raise ValueError("Alignment file must contain a finite 3x3 affine matrix")
        alignment = np.asarray(matrix, dtype=np.float64).tolist()

    image_width_mm = info.full_width_px * info.physical_size_x_um / 1000.0
    image_height_mm = info.full_height_px * info.physical_size_y_um / 1000.0
    real_result = {
        "status": "REAL_HE_PATHOLOGY_GEOMETRY_PREPARED_TARGET_INDEX_NOT_COMPUTED",
        "dataset": "10x Xenium FFPE human renal cell carcinoma gene and protein dataset",
        "evidence_level": "REAL_HISTOLOGY_AND_PATHOLOGY_ANNOTATIONS_ONLY",
        "image": {
            "full_shape_yxs": [info.full_height_px, info.full_width_px, info.channels],
            "pyramid_levels": info.levels,
            "analysis_level": info.selected_level,
            "analysis_shape_yxs": [
                info.selected_height_px,
                info.selected_width_px,
                info.channels,
            ],
            "physical_size_um_per_pixel": [
                info.physical_size_x_um,
                info.physical_size_y_um,
            ],
            "analysis_size_um_per_pixel": [
                info.selected_pixel_size_x_um,
                info.selected_pixel_size_y_um,
            ],
            "physical_extent_mm": [image_width_mm, image_height_mm],
        },
        "annotations": annotation_summary,
        "morphology_by_region": morphology,
        "tissue": {
            "analysis_level_pixels": int(np.count_nonzero(tissue)),
            "estimated_area_mm2": float(np.count_nonzero(tissue) * pixel_area_um2 / 1_000_000.0),
        },
        "lumen_candidates": {
            "classification": "UNVALIDATED_LUMEN_CANDIDATE",
            "count": len(candidates),
            "area_fraction_of_tissue": float(
                float(np.count_nonzero(lumen_mask)) / float(max(int(np.count_nonzero(tissue)), 1))
            ),
            "warning": (
                "H&E bright intratissue spaces are not validated vessels and may include "
                "tubules, ducts, tears or other artefacts. They are excluded "
                "from mechanistic indexing."
            ),
            "table": str(candidate_path),
        },
        "alignment_matrix_supplied": alignment,
        "qc_images": qc_images,
        "absolute_index": {
            "status": "NOT_COMPUTED",
            "abstention_code": "INSUFFICIENT_EVIDENCE",
            "reasons": [
                "TARGET_EXPRESSION_NOT_AVAILABLE",
                "SURFACE_ANTIGEN_CALIBRATION_NOT_AVAILABLE",
                "FUNCTIONALLY_PERFUSED_VESSELS_NOT_IDENTIFIED",
                "CELL_COORDINATES_NOT_AVAILABLE",
            ],
        },
    }
    result_path = output_dir / "he_pathology_result.json"
    result_path.write_text(json.dumps(real_result, indent=2), encoding="utf-8")

    claims = {
        "permitted": [
            "The complete supplied H&E OME-TIFF was reconstructed and transfer-verified.",
            "Pathologist-provided tumour, necrosis, immune-infiltration, "
            "haemorrhage, adipose and blood-vessel annotations were mapped "
            "in H&E pixel space.",
            "Tissue occupancy and relative H&E stain morphology were computed "
            "at a bounded pyramid level.",
            "Bright intratissue spaces were generated only as unvalidated lumen candidates for QC.",
        ],
        "conditional": [
            "H&E morphology may define geometry priors after independent validation "
            "against molecular or endothelial markers.",
            "The supplied affine matrix can align H&E and Xenium spaces once "
            "cell coordinates are available.",
        ],
        "unsupported": [
            "Target-positive cell fraction.",
            "Reachable fraction.",
            "Expression-reach gap.",
            "Dominant antibody-transport barrier.",
            "Penetration depth.",
            "Functional vascular perfusion.",
            "Clinical efficacy or programme outcome.",
        ],
        "abstention": real_result["absolute_index"],
    }
    claims_path = output_dir / "claims.json"
    claims_path.write_text(json.dumps(claims, indent=2), encoding="utf-8")

    source_paths = [he_path, annotation_path]
    if alignment_path is not None:
        source_paths.append(alignment_path)
    provenance = {
        "sources": [
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "md5": md5_file(path),
            }
            for path in source_paths
        ],
        "parameters": {
            "analysis_level": analysis_level,
            "tissue_background_colour_distance_min": 12.0,
            "lumen_min_equivalent_diameter_um": 12.0,
            "lumen_max_equivalent_diameter_um": 350.0,
            "lumen_min_confidence": 0.25,
        },
        "outputs": {
            "result": str(result_path),
            "claims": str(claims_path),
            "candidate_table": str(candidate_path),
            "qc_images": qc_images,
        },
    }
    manifest_path = output_dir / "processing_manifest.json"
    manifest_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    real_result["result_path"] = str(result_path)
    real_result["claims_path"] = str(claims_path)
    real_result["manifest_path"] = str(manifest_path)
    return real_result
