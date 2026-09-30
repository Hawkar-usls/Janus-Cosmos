#!/usr/bin/env python3
"""LUCI x Palomar JPFM star-morphology transfer and sky-overlap pilot.

This runner deliberately keeps two questions separate:
1) Can the frozen JPFM-2F-B star-morphology measurement contract be transported
   to real LUCI/LUCIFER near-IR images?
2) Do the exact 64 frozen Palomar pilot sources have any public LUCI archival
   observations close enough on the sky to permit a direct decades-later IR
   counterpart test?

A morphology-transfer PASS is not a Palomar-candidate confirmation. An archive
centre overlap is only an acquisition opportunity until the coordinate is
proven to land inside the actual LUCI FITS WCS footprint.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import statistics
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy import ndimage
from scipy.spatial import ConvexHull

from janus_cosmos.luci import read_luci_fits_image

REGISTRY_COMMIT = "7890cd5c8f4650c02dd439dbf96f09bc45638654"
POSS_COMMIT = "4005e200541b321ead3d6608f0162a14430ef1c2"
REGISTRY_MANIFEST = (
    "https://raw.githubusercontent.com/Hawkar-usls/janus-meta-registry/"
    f"{REGISTRY_COMMIT}/data/JANUS-PALOMAR-JPFM-2F-A-BLIND-STRUCTURAL-CLUSTER-MANIFEST-RUN-001.csv.gz"
)
POSS_BASE = f"https://raw.githubusercontent.com/jannefi/poss1-plate-slice/{POSS_COMMIT}/results/s0-642-20260814"
S0_URL = f"{POSS_BASE}/stage_S0.csv.gz"
STRUCTURAL_GZ_SHA256 = "166f5e6621ed2b065b7981b3c8208670f3c989b1394bd559c9005ab1fa6d07d9"
STRUCTURAL_CSV_SHA256 = "34b0ccde7c3683d07626774e52dac0a197451f729242204e59aae81397bdbc2e"
S0_GZ_SHA256 = "f19cf987756c62a68f55a472992d860e73ae63b3a4664189092b0e1fda77f7bb"
S0_CSV_SHA256 = "2ff92f2210acb387ef9ef4b88d561595d3883e9aab27065042627272b96590f0"
TAP_SYNC = "https://archive.lbto.org/tap/sync"

SAMPLE_N = 64
PER_CLUSTER = 4
MEASURE_RADIUS = 10
SOURCE_EDGE = 14
MAX_SOURCES_PER_FRAME = 64
ARCHIVE_CENTER_SEARCH_RADIUS_DEG = 0.20


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch_bytes(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "JANUS-COSMOS-LUCI-PALOMAR/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def require_hash(label: str, data: bytes, expected: str) -> None:
    got = sha256_bytes(data)
    if got != expected:
        raise RuntimeError(f"{label} sha256 mismatch: got={got} expected={expected}")


def tap_query(adql: str, timeout: int = 120) -> Table:
    body = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "QUERY": adql, "FORMAT": "votable"}).encode()
    req = urllib.request.Request(
        TAP_SYNC,
        data=body,
        headers={"User-Agent": "JANUS-COSMOS-LUCI-PALOMAR/1.0", "Content-Type": "application/x-www-form-urlencoded"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return Table.read(io.BytesIO(r.read()), format="votable")


def _s(value) -> str:
    try:
        if getattr(value, "mask", False):
            return ""
    except Exception:
        pass
    return "" if value is None else str(value).strip()


def _f(value, default=float("nan")) -> float:
    try:
        x = float(value)
        return x if math.isfinite(x) else default
    except Exception:
        return default


def load_frozen_palomar_sample() -> tuple[list[dict], dict]:
    mgz = fetch_bytes(REGISTRY_MANIFEST)
    require_hash("Palomar structural manifest gzip", mgz, STRUCTURAL_GZ_SHA256)
    mcsv = gzip.decompress(mgz)
    require_hash("Palomar structural manifest csv", mcsv, STRUCTURAL_CSV_SHA256)
    rows = list(csv.DictReader(io.StringIO(mcsv.decode("utf-8"))))
    if len(rows) != 122820:
        raise RuntimeError(f"structural manifest row invariant failed: {len(rows)}")

    groups: dict[int, list[dict]] = defaultdict(list)
    for r in rows:
        groups[int(float(r["structural_cluster"]))].append(r)

    chosen: list[dict] = []
    for cluster in range(16):
        g = groups[cluster]
        scores = [float(r["anomaly_score"]) for r in g]
        med = statistics.median(scores)
        typical = sorted(g, key=lambda r: (abs(float(r["anomaly_score"]) - med), r["src_id"]))[:2]
        used = {r["src_id"] for r in typical}
        unusual = sorted(
            (r for r in g if r["src_id"] not in used),
            key=lambda r: (-float(r["anomaly_score"]), r["src_id"]),
        )[:2]
        for r in typical:
            chosen.append({"src_id": r["src_id"], "structural_cluster": cluster, "sample_role": "typical"})
        for r in unusual:
            chosen.append({"src_id": r["src_id"], "structural_cluster": cluster, "sample_role": "unusual"})

    if len(chosen) != SAMPLE_N or len({r["src_id"] for r in chosen}) != SAMPLE_N:
        raise RuntimeError("Palomar deterministic 64-source selection invariant failed")

    s0gz = fetch_bytes(S0_URL)
    require_hash("POSS S0 gzip", s0gz, S0_GZ_SHA256)
    s0csv = gzip.decompress(s0gz)
    require_hash("POSS S0 csv", s0csv, S0_CSV_SHA256)
    wanted = {r["src_id"] for r in chosen}
    positions = {}
    for r in csv.DictReader(io.StringIO(s0csv.decode("utf-8"))):
        sid = r.get("src_id", "")
        if sid in wanted:
            positions[sid] = (float(r["ra"]), float(r["dec"]))
    if len(positions) != SAMPLE_N:
        missing = sorted(wanted - positions.keys())[:5]
        raise RuntimeError(f"Palomar sample position join incomplete: {len(positions)}/64 missing={missing}")
    for r in chosen:
        r["ra_deg"], r["dec_deg"] = positions[r["src_id"]]

    return chosen, {
        "registry_commit": REGISTRY_COMMIT,
        "poss_commit": POSS_COMMIT,
        "structural_manifest_gzip_sha256": sha256_bytes(mgz),
        "structural_manifest_csv_sha256": sha256_bytes(mcsv),
        "stage_S0_gzip_sha256": sha256_bytes(s0gz),
        "stage_S0_csv_sha256": sha256_bytes(s0csv),
        "sample_n": SAMPLE_N,
        "selection": "JPFM-2F-B frozen deterministic 2 typical + 2 unusual per structural cluster",
    }


def angular_sep_deg(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    a1, d1, a2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    s = math.sin((d2 - d1) / 2) ** 2 + math.cos(d1) * math.cos(d2) * math.sin((a2 - a1) / 2) ** 2
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(max(0.0, s)))))


def query_luci_near_palomar(sample: list[dict]) -> tuple[list[dict], dict]:
    clauses = [
        f"1=CONTAINS(lbt.lbt.s_point,CIRCLE('ICRS',{r['ra_deg']:.10f},{r['dec_deg']:.10f},{ARCHIVE_CENTER_SEARCH_RADIUS_DEG:.6f}))"
        for r in sample
    ]
    adql = (
        "SELECT TOP 5000 lbt.luci.instrument,lbt.luci.telescope,lbt.luci.object,lbt.luci.filters,"
        "lbt.luci.gratname,lbt.luci.imagetyp,lbt.luci.file_name,lbt.luci.date_obs,"
        "lbt.luci.crval1,lbt.luci.crval2,lbt.luci.naxis1,lbt.luci.naxis2,lbt.luci.pixscale,"
        "lbt.lbt.file_url,lbt.lbt.policy "
        "FROM lbt.luci JOIN lbt.lbt ON lbt.luci.file_name=lbt.lbt.file_name "
        "WHERE lbt.lbt.policy='FREE' AND lbt.luci.imagetyp='SCIENCE' "
        "AND (lbt.luci.gratname='Mirror' OR lbt.luci.gratname='mirror') AND (" + " OR ".join(clauses) + ")"
    )
    rows = tap_query(adql)
    candidates: list[dict] = []
    for row in rows:
        cra, cdec = _f(row["crval1"]), _f(row["crval2"])
        nx, ny, scale = _f(row["naxis1"]), _f(row["naxis2"]), _f(row["pixscale"])
        if not all(math.isfinite(x) for x in (cra, cdec, nx, ny, scale)) or scale <= 0:
            continue
        # LUCI pixscale is archived in arcsec/pixel. Use half-diagonal as a conservative
        # rectangular-footprint radius, then require actual FITS WCS before a direct counterpart claim.
        half_diag_deg = 0.5 * math.hypot(nx * scale, ny * scale) / 3600.0
        if not (0 < half_diag_deg < 0.5):
            continue
        for src in sample:
            sep = angular_sep_deg(src["ra_deg"], src["dec_deg"], cra, cdec)
            if sep <= half_diag_deg:
                candidates.append({
                    "src_id": src["src_id"],
                    "structural_cluster": src["structural_cluster"],
                    "sample_role": src["sample_role"],
                    "ra_deg": src["ra_deg"],
                    "dec_deg": src["dec_deg"],
                    "archive_file_name": _s(row["file_name"]),
                    "file_url": _s(row["file_url"]),
                    "instrument": _s(row["instrument"]),
                    "target": _s(row["object"]),
                    "filters": _s(row["filters"]),
                    "date_obs": _s(row["date_obs"]),
                    "centre_sep_deg": sep,
                    "conservative_half_diagonal_deg": half_diag_deg,
                    "status": "CENTER_GEOMETRY_OVERLAP__FITS_WCS_REQUIRED",
                })
    return candidates, {
        "tap_service": TAP_SYNC,
        "search_radius_deg": ARCHIVE_CENTER_SEARCH_RADIUS_DEG,
        "tap_rows_returned": len(rows),
        "conservative_geometry_overlaps": len(candidates),
        "unique_palomar_sources_with_overlap": len({x["src_id"] for x in candidates}),
    }


def robust_sigma(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return 1.0
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    s = 1.4826 * mad
    if not math.isfinite(s) or s <= 0:
        s = float(np.std(x))
    return s if math.isfinite(s) and s > 0 else 1.0


def radial_profile(signal: np.ndarray, cy: float, cx: float, radius: int = MEASURE_RADIUS) -> np.ndarray:
    yy, xx = np.indices(signal.shape, dtype=float)
    rr = np.hypot(xx - cx, yy - cy)
    edges = np.arange(0.0, radius + 0.5001, 0.5)
    prof = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        vals = signal[(rr >= lo) & (rr < hi)]
        prof.append(float(np.mean(vals)) if vals.size else 0.0)
    a = np.asarray(prof, dtype=float)
    mx = float(np.max(a)) if len(a) else 0.0
    return a / mx if mx > 0 else a


def measure_source(signal: np.ndarray, y: int, x: int, sigma_bg: float) -> dict | None:
    y0, y1 = y - MEASURE_RADIUS, y + MEASURE_RADIUS + 1
    x0, x1 = x - MEASURE_RADIUS, x + MEASURE_RADIUS + 1
    if y0 < 0 or x0 < 0 or y1 > signal.shape[0] or x1 > signal.shape[1]:
        return None
    sub = np.asarray(signal[y0:y1, x0:x1], dtype=float)
    sy, sx = np.indices(sub.shape, dtype=float)
    py = px = MEASURE_RADIUS
    peak = float(sub[py, px])
    if not math.isfinite(peak) or peak <= 0:
        return None
    rr0 = np.hypot(sx - px, sy - py)
    threshold = max(2.0 * sigma_bg, 0.03 * peak)
    weights = np.where((sub > threshold) & (rr0 <= MEASURE_RADIUS), sub, 0.0)
    sw = float(weights.sum())
    if sw <= 0:
        return None
    cx = float((weights * sx).sum() / sw); cy = float((weights * sy).sum() / sw)
    dx, dy = sx - cx, sy - cy
    cov = np.array([
        [float((weights * dx * dx).sum() / sw), float((weights * dx * dy).sum() / sw)],
        [float((weights * dx * dy).sum() / sw), float((weights * dy * dy).sum() / sw)],
    ])
    vals, vecs = np.linalg.eigh(cov)
    vals = np.maximum(vals, 1e-8)
    fmin, fmax = 2.354820045 * math.sqrt(float(vals[0])), 2.354820045 * math.sqrt(float(vals[1]))
    fwhm = math.sqrt(fmin * fmax)
    elong = fmax / fmin if fmin > 0 else float("nan")
    major_vec = vecs[:, 1]
    orient = math.degrees(math.atan2(float(major_vec[1]), float(major_vec[0])))
    rr = np.hypot(sx - cx, sy - cy)
    f2 = float(np.clip(sub[rr <= 2.0], 0, None).sum())
    f6 = float(np.clip(sub[rr <= 6.0], 0, None).sum())
    concentration = f2 / f6 if f6 > 0 else float("nan")

    mask = (sub >= max(3.0 * sigma_bg, 0.30 * peak)) & (rr <= MEASURE_RADIUS)
    lab, _ = ndimage.label(mask)
    iy = int(np.clip(round(cy), 0, sub.shape[0] - 1)); ix = int(np.clip(round(cx), 0, sub.shape[1] - 1))
    target_label = int(lab[iy, ix])
    comp = (lab == target_label) if target_label > 0 else np.zeros_like(mask)
    area = int(comp.sum())
    perimeter_mask = comp & ~ndimage.binary_erosion(comp)
    perimeter = int(perimeter_mask.sum())
    circularity = 4.0 * math.pi * area / (perimeter * perimeter) if perimeter > 0 else float("nan")
    boundary = np.argwhere(perimeter_mask)
    circle_dev = convex_defect = float("nan")
    if len(boundary) >= 3:
        by, bx = boundary[:, 0].astype(float), boundary[:, 1].astype(float)
        br = np.hypot(bx - cx, by - cy)
        mr = float(np.mean(br))
        circle_dev = float(np.std(br) / mr) if mr > 0 else float("nan")
        try:
            hull_area = float(ConvexHull(np.column_stack([bx, by])).volume)
            convex_defect = max(0.0, min(1.0, (hull_area - area) / hull_area)) if hull_area > 0 else float("nan")
        except Exception:
            pass
    return {
        "x": int(x), "y": int(y), "peak_snr": peak / sigma_bg,
        "moment_fwhm_px": fwhm, "moment_fwhm_major_px": fmax, "moment_fwhm_minor_px": fmin,
        "elongation": elong, "orientation_deg": orient,
        "core_concentration_r2_over_r6": concentration, "threshold_area_px": area,
        "circularity": circularity, "convex_defect_fraction": convex_defect,
        "circle_deviation": circle_dev, "profile": radial_profile(np.clip(sub, 0, None), cy, cx),
    }


def detect_and_measure(image: np.ndarray) -> tuple[list[dict], dict]:
    arr = np.asarray(image, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size < 1024:
        return [], {"reason": "insufficient_finite_pixels"}
    bg = float(np.median(finite)); sigma = robust_sigma(finite)
    signal = np.where(np.isfinite(arr), arr - bg, 0.0)
    # LUCI detector images are direct-positive. If the archive product is inverted,
    # pick the polarity with the stronger high tail but record the choice.
    pos_tail = float(np.percentile(signal, 99.9)); neg_tail = float(-np.percentile(signal, 0.1))
    polarity = "direct_positive"
    if neg_tail > 1.5 * max(pos_tail, 1e-9):
        signal = -signal
        polarity = "inverted_by_tail_test"
    signal = np.clip(signal, 0, None)
    sm = ndimage.gaussian_filter(signal, 1.0)
    maxima = sm == ndimage.maximum_filter(sm, size=7, mode="nearest")
    maxima &= signal >= 7.0 * sigma
    maxima[:SOURCE_EDGE, :] = False; maxima[-SOURCE_EDGE:, :] = False
    maxima[:, :SOURCE_EDGE] = False; maxima[:, -SOURCE_EDGE:] = False
    ys, xs = np.nonzero(maxima)
    order = sorted(zip(ys, xs), key=lambda p: (-float(signal[p[0], p[1]]), int(p[0]), int(p[1])))
    measured = []
    for y, x in order:
        if any((x - q["x"]) ** 2 + (y - q["y"]) ** 2 < 12 ** 2 for q in measured):
            continue
        m = measure_source(signal, int(y), int(x), sigma)
        if m is not None and math.isfinite(m["moment_fwhm_px"]) and math.isfinite(m["elongation"]):
            measured.append(m)
        if len(measured) >= MAX_SOURCES_PER_FRAME:
            break

    for m in measured:
        refs = [
            q for q in measured if q is not m and 0.5 * m["peak_snr"] <= q["peak_snr"] <= 2.0 * m["peak_snr"]
        ][:12]
        m["reference_count"] = len(refs)
        if len(refs) >= 3:
            ref_f = float(np.median([q["moment_fwhm_px"] for q in refs]))
            m["reference_median_fwhm_px"] = ref_f
            m["fwhm_ratio"] = m["moment_fwhm_px"] / ref_f if ref_f > 0 else float("nan")
            n = min([len(m["profile"])] + [len(q["profile"]) for q in refs])
            refp = np.mean(np.vstack([q["profile"][:n] for q in refs]), axis=0)
            diff = (m["profile"][:n] - refp) * np.where(refp > 0.1, refp, 0.0)
            diff[:2] = 0.0
            m["profile_diff"] = float(math.sqrt(float(np.mean(diff * diff))))
        else:
            m["reference_median_fwhm_px"] = m["fwhm_ratio"] = m["profile_diff"] = float("nan")
        m["profile"] = [float(x) for x in m["profile"]]

    finite_fraction = (
        sum(math.isfinite(m["moment_fwhm_px"]) and math.isfinite(m["elongation"]) for m in measured) / len(measured)
        if measured else 0.0
    )
    ref_fraction = sum(m["reference_count"] >= 3 for m in measured) / len(measured) if measured else 0.0
    gate = len(measured) >= 12 and finite_fraction >= 0.95 and ref_fraction >= 0.50
    return measured, {
        "background_median": bg, "robust_sigma": sigma, "polarity": polarity,
        "source_count": len(measured), "finite_fwhm_elongation_fraction": finite_fraction,
        "local_reference_ge3_fraction": ref_fraction,
        "frame_transfer_gate_pass": bool(gate),
    }


def download_to_cache(url: str, cache: Path) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    name = url.rsplit("/", 1)[-1] or hashlib.sha256(url.encode()).hexdigest() + ".fits"
    path = cache / name
    if not path.exists():
        path.write_bytes(fetch_bytes(url, timeout=180))
    return path


def load_manifests(paths: list[Path]) -> list[dict]:
    targets = []
    seen = set()
    for path in paths:
        doc = json.loads(path.read_text(encoding="utf-8"))
        for t in doc.get("targets", []):
            key = t.get("target")
            if key in seen:
                continue
            seen.add(key); targets.append(t)
    return targets


def run_morphology_transfer(manifest_paths: list[Path], cache_dir: Path) -> dict:
    frames = []
    source_rows = []
    for target in load_manifests(manifest_paths):
        for item in target.get("filters", []):
            url = str(item.get("url", ""))
            if not url:
                continue
            path = download_to_cache(url, cache_dir)
            image, provenance = read_luci_fits_image(path, require_imaging=True, expected_instrument=item.get("instrument"))
            sources, qa = detect_and_measure(image)
            frame_id = f"{target['target']}::{item.get('filter','UNKNOWN')}::{item.get('archive_file_name',path.name)}"
            frames.append({
                "frame_id": frame_id, "target": target["target"], "filter": item.get("filter"),
                "archive_file_name": item.get("archive_file_name", path.name), "sha256": sha256_bytes(path.read_bytes()),
                "provenance": provenance, **qa,
                "metric_medians": {
                    k: float(np.nanmedian([m[k] for m in sources])) if sources else float("nan")
                    for k in ["moment_fwhm_px", "elongation", "core_concentration_r2_over_r6", "circularity", "fwhm_ratio", "profile_diff"]
                },
            })
            for m in sources:
                source_rows.append({"frame_id": frame_id, **m})
    pass_count = sum(f["frame_transfer_gate_pass"] for f in frames)
    return {
        "frames": frames,
        "sources": source_rows,
        "frame_count": len(frames),
        "frame_gate_pass_count": pass_count,
        "frame_gate_pass_fraction": pass_count / len(frames) if frames else 0.0,
        "outcome": "TRANSFER_FEASIBILITY_PASS" if frames and pass_count / len(frames) >= 0.75 else "FAIL_CLOSED_TRANSFER_FEASIBILITY",
        "palomar_metric_contract": [
            "moment_fwhm_px", "moment_fwhm_major_px", "moment_fwhm_minor_px", "elongation", "orientation_deg",
            "core_concentration_r2_over_r6", "threshold_area_px", "circularity", "convex_defect_fraction",
            "circle_deviation", "reference_count", "reference_median_fwhm_px", "fwhm_ratio", "profile_diff",
        ],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", action="append", required=True, help="LUCI archive manifest; may be repeated")
    ap.add_argument("--output", default="results/luci_palomar_transfer/receipt.json")
    ap.add_argument("--cache-dir", default=".cache/luci_palomar_transfer")
    args = ap.parse_args()

    sample, sample_binding = load_frozen_palomar_sample()
    try:
        overlaps, overlap_meta = query_luci_near_palomar(sample)
        overlap_status = "OVERLAP_PREFLIGHT_COMPLETE"
    except Exception as exc:
        overlaps, overlap_meta = [], {"error": f"{type(exc).__name__}: {exc}"}
        overlap_status = "FAIL_CLOSED_ARCHIVE_OVERLAP_QUERY"

    morphology = run_morphology_transfer([Path(x) for x in args.manifest], Path(args.cache_dir))
    receipt = {
        "schema": "janus.cosmos.luci_palomar_transfer.v1",
        "experiment_id": "LUCI-JPFM-2F-TRANSFER-001",
        "status": "EXECUTED",
        "instrument_boundary": "LUCI/LUCIFER_ONLY",
        "palomar_parent": "JPFM-2F-B star-morphology contract + deterministic 64-source sample",
        "frozen_palomar_binding": sample_binding,
        "morphology_transfer": morphology,
        "direct_sky_overlap_preflight": {
            "status": overlap_status,
            **overlap_meta,
            "overlaps": overlaps,
            "interpretation": (
                "A conservative centre-geometry overlap is only an acquisition opportunity. "
                "A direct Palomar-to-LUCI counterpart claim additionally requires the frozen Palomar coordinate "
                "to land inside the downloaded LUCI FITS WCS footprint and a preregistered local-source test."
            ),
        },
        "claim_ceiling": (
            "METHOD_TRANSFER_AND_ARCHIVE_OVERLAP_PREFLIGHT_ONLY__NO_PALOMAR_ANOMALY_CONFIRMATION__"
            "NO_UAP_ORIGIN_IDENTIFICATION__NO_CAUSALITY"
        ),
    }
    out = Path(args.output); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "outcome": morphology["outcome"],
        "frames": morphology["frame_count"],
        "frame_pass_fraction": morphology["frame_gate_pass_fraction"],
        "palomar_sources_with_archive_overlap": overlap_meta.get("unique_palomar_sources_with_overlap", 0),
        "overlap_status": overlap_status,
        "receipt": str(out),
    }, indent=2))
    return 0 if morphology["outcome"] == "TRANSFER_FEASIBILITY_PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
