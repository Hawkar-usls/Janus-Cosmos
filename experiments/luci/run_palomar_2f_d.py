#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy.spatial import cKDTree

from janus_cosmos.luci import read_luci_fits_image
from janus_cosmos.luci_psf import detect_psf_sources
from janus_cosmos.luci_psf_r1 import measure_psf_at, psf_relative_injection_recovery_gate
from janus_cosmos.pipeline import EventWriter, download_source, sha256_file
from experiments.luci.run_palomar_2f_c import (
    POSS_COMMIT,
    S0_CSV_SHA256,
    S0_GZ_SHA256,
    S0_URL,
    TAP_SYNC,
    _f,
    _s,
    angular_sep_deg,
    exact_wcs_pixel,
    tap_query,
)

EXPECTED_S0_N = 122820
ARCHIVE_BATCH = 5000
MAX_ARCHIVE_BATCHES = 50
MAX_FITS_DOWNLOAD_FILES = 250
MAX_CONTROLS = 20
MIN_CONTROLS = 8
CONTROL_RADIUS_PX = 300.0
CONTROL_SNR_RATIO_MIN = 0.5
CONTROL_SNR_RATIO_MAX = 2.0
PARENT_R1_RESULT = "data/luci_palomar/LUCI-PALOMAR-JPFM-2F-C-R1-RUN-001.json"
PARENT_R1_ARTIFACT_SHA256 = "2d7b8d4e233908b94a89c9b7dfb4e8c91c3e146a92e1fb3982434f3b213e4f71"


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _fetch(url: str, timeout: int = 180) -> bytes:
    req = urllib.request.Request(
        url, headers={"User-Agent": "JANUS-COSMOS-LUCI-JPFM-2F-D/1.0"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _require_parent_r1(repo_root: Path) -> dict:
    p = repo_root / PARENT_R1_RESULT
    d = json.loads(p.read_text(encoding="utf-8"))
    if d.get("scientific_status") != "PSF_GATE_VALIDATED__NO_EXACT_LUCI_OVERLAP_IN_FROZEN_640":
        raise RuntimeError("parent R1 scientific status mismatch")
    v = d.get("psf_injection_recovery_validation", {})
    if not v.get("all_frames_pass") or int(v.get("passed_frame_count", 0)) != 6:
        raise RuntimeError("parent R1 PSF validation not admitted")
    if d.get("artifact", {}).get("sha256") != PARENT_R1_ARTIFACT_SHA256:
        raise RuntimeError("parent R1 artifact SHA mismatch")
    return d


def build_full_frozen_s0_coordinates(out_csv: Path) -> dict:
    gz = _fetch(S0_URL)
    if _sha(gz) != S0_GZ_SHA256:
        raise RuntimeError(f"S0 gzip hash mismatch: {_sha(gz)} != {S0_GZ_SHA256}")
    raw = gzip.decompress(gz)
    if _sha(raw) != S0_CSV_SHA256:
        raise RuntimeError(f"S0 csv hash mismatch: {_sha(raw)} != {S0_CSV_SHA256}")

    rows = []
    seen = set()
    for r in csv.DictReader(io.StringIO(raw.decode("utf-8"))):
        sid = str(r.get("src_id", "")).strip()
        if not sid or sid in seen:
            raise RuntimeError(f"missing/duplicate src_id in S0: {sid!r}")
        ra, dec = float(r["ra"]), float(r["dec"])
        if not (math.isfinite(ra) and math.isfinite(dec) and 0.0 <= ra < 360.0 and -90.0 <= dec <= 90.0):
            raise RuntimeError(f"invalid coordinate for {sid}: {ra},{dec}")
        rows.append({"src_id": sid, "ra_deg": ra, "dec_deg": dec})
        seen.add(sid)

    if len(rows) != EXPECTED_S0_N:
        raise RuntimeError(f"full S0 row count changed: {len(rows)} != {EXPECTED_S0_N}")
    rows.sort(key=lambda x: x["src_id"])

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["src_id", "ra_deg", "dec_deg"])
        w.writeheader()
        w.writerows(rows)
    return {
        "rows": rows,
        "path": str(out_csv),
        "sha256": sha256_file(out_csv),
        "sample_n": len(rows),
        "poss_commit": POSS_COMMIT,
        "stage_S0_gzip_sha256": _sha(gz),
        "stage_S0_csv_sha256": _sha(raw),
    }


def _adql_quote(s: str) -> str:
    return "'" + s.replace("'", "''") + "'"


def query_full_luci_inventory() -> list[dict]:
    cols = (
        "lbt.luci.instrument,lbt.luci.telescope,lbt.luci.object,lbt.luci.filters,"
        "lbt.luci.gratname,lbt.luci.imagetyp,lbt.luci.file_name,lbt.luci.date_obs,"
        "lbt.luci.crval1,lbt.luci.crval2,lbt.luci.crpix1,lbt.luci.crpix2,"
        "lbt.luci.cd1_1,lbt.luci.cd1_2,lbt.luci.cd2_1,lbt.luci.cd2_2,"
        "lbt.luci.ctype1,lbt.luci.ctype2,lbt.luci.naxis1,lbt.luci.naxis2,"
        "lbt.luci.pixscale,lbt.lbt.file_url,lbt.lbt.policy"
    )
    base = (
        " FROM lbt.luci JOIN lbt.lbt ON lbt.luci.file_name=lbt.lbt.file_name "
        "WHERE lbt.lbt.policy='FREE' AND lbt.luci.imagetyp='SCIENCE' "
        "AND (lbt.luci.gratname='Mirror' OR lbt.luci.gratname='mirror') "
    )
    by_name: dict[str, dict] = {}
    last = ""
    for batch_idx in range(MAX_ARCHIVE_BATCHES):
        more = "" if not last else f"AND lbt.luci.file_name > {_adql_quote(last)} "
        adql = (
            f"SELECT TOP {ARCHIVE_BATCH} {cols}{base}{more}"
            "ORDER BY lbt.luci.file_name"
        )
        tab = tap_query(adql, timeout=240)
        if not len(tab):
            break
        names = []
        for row in tab:
            fname = _s(row["file_name"])
            if not fname:
                continue
            names.append(fname)
            rec = {
                "file_name": fname,
                "file_url": _s(row["file_url"]),
                "instrument": _s(row["instrument"]),
                "telescope": _s(row["telescope"]),
                "target": _s(row["object"]),
                "filters": _s(row["filters"]),
                "date_obs": _s(row["date_obs"]),
                "crval1": _f(row["crval1"]),
                "crval2": _f(row["crval2"]),
                "crpix1": _f(row["crpix1"]),
                "crpix2": _f(row["crpix2"]),
                "cd1_1": _f(row["cd1_1"]),
                "cd1_2": _f(row["cd1_2"]),
                "cd2_1": _f(row["cd2_1"]),
                "cd2_2": _f(row["cd2_2"]),
                "ctype1": _s(row["ctype1"]),
                "ctype2": _s(row["ctype2"]),
                "naxis1": _f(row["naxis1"]),
                "naxis2": _f(row["naxis2"]),
                "pixscale": _f(row["pixscale"]),
            }
            by_name[fname] = rec
        if not names:
            raise RuntimeError("LUCI inventory pagination returned no usable filenames")
        new_last = max(names)
        if new_last <= last:
            raise RuntimeError("LUCI inventory pagination did not advance")
        last = new_last
        if len(tab) < ARCHIVE_BATCH:
            break
    else:
        raise RuntimeError("LUCI archive inventory exceeded pagination safety bound")
    return sorted(by_name.values(), key=lambda x: x["file_name"])


INVENTORY_FIELDS = [
    "file_name", "file_url", "instrument", "telescope", "target", "filters", "date_obs",
    "crval1", "crval2", "crpix1", "crpix2",
    "cd1_1", "cd1_2", "cd2_1", "cd2_2",
    "ctype1", "ctype2", "naxis1", "naxis2", "pixscale",
]


def write_inventory(rows: list[dict], path: Path) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=INVENTORY_FIELDS)
        w.writeheader()
        w.writerows(rows)
    return {"path": str(path), "sha256": sha256_file(path), "row_count": len(rows)}


def _unit_vectors(rows: list[dict], rakey: str, deckey: str) -> np.ndarray:
    ra = np.deg2rad(np.asarray([float(x[rakey]) for x in rows]))
    dec = np.deg2rad(np.asarray([float(x[deckey]) for x in rows]))
    c = np.cos(dec)
    return np.column_stack((c*np.cos(ra), c*np.sin(ra), np.sin(dec)))


def _frame_half_diag_deg(frame: dict) -> float | None:
    base = [frame["crval1"], frame["crval2"], frame["naxis1"], frame["naxis2"]]
    if not all(math.isfinite(float(v)) for v in base):
        return None
    nx, ny = float(frame["naxis1"]), float(frame["naxis2"])
    if nx <= 0 or ny <= 0:
        return None
    cd = [frame["cd1_1"], frame["cd1_2"], frame["cd2_1"], frame["cd2_2"]]
    if all(math.isfinite(float(v)) for v in cd):
        sx = 3600.0 * math.hypot(float(frame["cd1_1"]), float(frame["cd2_1"]))
        sy = 3600.0 * math.hypot(float(frame["cd1_2"]), float(frame["cd2_2"]))
    else:
        scale = float(frame["pixscale"])
        if not math.isfinite(scale) or scale <= 0:
            return None
        sx = sy = scale
    if sx <= 0 or sy <= 0:
        return None
    half_diag = 0.5 * math.hypot(nx*sx, ny*sy) / 3600.0
    if not (0.0 < half_diag < 0.5):
        return None
    return half_diag


def archive_first_spatial_crossmatch(sources: list[dict], inventory: list[dict]) -> tuple[list[dict], dict]:
    xyz = _unit_vectors(sources, "ra_deg", "dec_deg")
    tree = cKDTree(xyz)
    pairs = []
    valid_frames = 0
    max_half_diag = 0.0
    for frame in inventory:
        half_diag = _frame_half_diag_deg(frame)
        if half_diag is None:
            continue
        cra, cdec = float(frame["crval1"]), float(frame["crval2"])
        if not (0.0 <= cra < 360.0 and -90.0 <= cdec <= 90.0):
            continue
        valid_frames += 1
        max_half_diag = max(max_half_diag, half_diag)
        r = math.radians(half_diag)
        center = _unit_vectors([{"ra": cra, "dec": cdec}], "ra", "dec")[0]
        chord = 2.0 * math.sin(r/2.0)
        for idx in tree.query_ball_point(center, chord):
            src = sources[int(idx)]
            sep = angular_sep_deg(src["ra_deg"], src["dec_deg"], cra, cdec)
            if sep <= half_diag:
                pairs.append({
                    "src_id": src["src_id"],
                    "ra_deg": src["ra_deg"],
                    "dec_deg": src["dec_deg"],
                    "file_name": frame["file_name"],
                    "file_url": frame["file_url"],
                    "instrument": frame["instrument"],
                    "target": frame["target"],
                    "filters": frame["filters"],
                    "date_obs": frame["date_obs"],
                    "center_sep_deg": sep,
                    "half_diagonal_deg": half_diag,
                })
    pairs.sort(key=lambda x: (x["src_id"], x["file_name"]))
    return pairs, {
        "inventory_rows": len(inventory),
        "valid_geometry_frames": valid_frames,
        "max_half_diagonal_deg": max_half_diag,
        "coarse_pair_count": len(pairs),
        "coarse_unique_sources": len({x["src_id"] for x in pairs}),
        "coarse_unique_files": len({x["file_name"] for x in pairs}),
    }


PAIR_FIELDS = [
    "src_id", "ra_deg", "dec_deg", "file_name", "file_url", "instrument",
    "target", "filters", "date_obs", "center_sep_deg", "half_diagonal_deg",
]


def write_pairs(rows: list[dict], path: Path) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=PAIR_FIELDS)
        w.writeheader()
        w.writerows(rows)
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "pair_count": len(rows),
        "unique_sources": len({x["src_id"] for x in rows}),
        "unique_files": len({x["file_name"] for x in rows}),
    }


def metadata_wcs_contains(frame: dict, ra: float, dec: float) -> tuple[bool, dict]:
    required = [
        "crval1", "crval2", "crpix1", "crpix2", "cd1_1", "cd1_2", "cd2_1", "cd2_2",
        "naxis1", "naxis2",
    ]
    if not all(math.isfinite(float(frame[k])) for k in required):
        return False, {"reason": "INCOMPLETE_ARCHIVE_WCS"}
    if not frame.get("ctype1") or not frame.get("ctype2"):
        return False, {"reason": "MISSING_CTYPE"}
    hdr = fits.Header()
    hdr["NAXIS"] = 2
    hdr["NAXIS1"] = int(round(float(frame["naxis1"])))
    hdr["NAXIS2"] = int(round(float(frame["naxis2"])))
    hdr["CTYPE1"] = frame["ctype1"]
    hdr["CTYPE2"] = frame["ctype2"]
    hdr["CRVAL1"] = float(frame["crval1"])
    hdr["CRVAL2"] = float(frame["crval2"])
    hdr["CRPIX1"] = float(frame["crpix1"])
    hdr["CRPIX2"] = float(frame["crpix2"])
    hdr["CD1_1"] = float(frame["cd1_1"])
    hdr["CD1_2"] = float(frame["cd1_2"])
    hdr["CD2_1"] = float(frame["cd2_1"])
    hdr["CD2_2"] = float(frame["cd2_2"])
    try:
        w = WCS(hdr).celestial
        x, y = w.world_to_pixel_values(float(ra), float(dec))
        x, y = float(x), float(y)
    except Exception as e:
        return False, {"reason": "ARCHIVE_WCS_ERROR", "error_type": type(e).__name__}
    nx, ny = int(hdr["NAXIS1"]), int(hdr["NAXIS2"])
    inside = math.isfinite(x) and math.isfinite(y) and (-0.5 <= x < nx-0.5) and (-0.5 <= y < ny-0.5)
    return bool(inside), {"x": x, "y": y, "naxis1": nx, "naxis2": ny}


def metadata_wcs_preflight(coarse: list[dict], inventory: list[dict]) -> list[dict]:
    by_file = {x["file_name"]: x for x in inventory}
    out = []
    for p in coarse:
        inside, wm = metadata_wcs_contains(by_file[p["file_name"]], p["ra_deg"], p["dec_deg"])
        if inside:
            out.append({**p, "metadata_wcs_x": wm["x"], "metadata_wcs_y": wm["y"]})
    out.sort(key=lambda x: (x["src_id"], x["file_name"]))
    return out


META_PAIR_FIELDS = PAIR_FIELDS + ["metadata_wcs_x", "metadata_wcs_y"]


def write_meta_pairs(rows: list[dict], path: Path) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=META_PAIR_FIELDS)
        w.writeheader()
        w.writerows(rows)
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "pair_count": len(rows),
        "unique_sources": len({x["src_id"] for x in rows}),
        "unique_files": len({x["file_name"] for x in rows}),
    }


def _robust_z(value: float, vals: list[float]) -> float | None:
    if not vals:
        return None
    a = np.asarray(vals, dtype=float)
    med = float(np.median(a))
    mad = float(np.median(np.abs(a-med)))
    scale = 1.4826*mad
    if not math.isfinite(scale) or scale <= 0:
        return None
    return float((value-med)/scale)


def counterpart_with_matched_controls(image: np.ndarray, x: float, y: float) -> dict:
    target = measure_psf_at(image, y, x, search_radius_px=2)
    if target is None:
        return {
            "status": "NO_R1_PSF_SOURCE_WITHIN_2P5PX",
            "counterpart_present": False,
            "measurement_method": "R1 local coordinate classifier",
        }

    all_src = detect_psf_sources(image, max_sources=4096)
    candidates = []
    for q in all_src:
        d = math.hypot(q.x-x, q.y-y)
        if d < 10.0 or d > CONTROL_RADIUS_PX:
            continue
        ratio = q.peak_snr / target.peak_snr if target.peak_snr > 0 else float("inf")
        if CONTROL_SNR_RATIO_MIN <= ratio <= CONTROL_SNR_RATIO_MAX:
            candidates.append((abs(math.log(ratio)), d, q))
    candidates.sort(key=lambda z: (z[0], z[1], z[2].x, z[2].y))
    controls = [z[2] for z in candidates[:MAX_CONTROLS]]

    result = {
        "status": "COUNTERPART_DETECTED",
        "counterpart_present": True,
        "measurement_method": "R1 local coordinate classifier",
        "source": target.to_dict(),
        "matched_control_contract": {
            "radius_px": CONTROL_RADIUS_PX,
            "peak_snr_ratio_min": CONTROL_SNR_RATIO_MIN,
            "peak_snr_ratio_max": CONTROL_SNR_RATIO_MAX,
            "min_controls": MIN_CONTROLS,
            "max_controls": MAX_CONTROLS,
        },
        "matched_control_count": len(controls),
    }
    if len(controls) < MIN_CONTROLS:
        result["morphology_status"] = "INSUFFICIENT_MATCHED_LOCAL_CONTROLS"
        return result

    fw = [q.fwhm_geom_px for q in controls]
    el = [q.elongation for q in controls]
    result["morphology_status"] = "MATCHED_LOCAL_CONTROL_COMPARISON_AVAILABLE"
    result["matched_controls"] = {
        "median_fwhm_geom_px": float(np.median(fw)),
        "median_elongation": float(np.median(el)),
        "target_robust_z_fwhm": _robust_z(target.fwhm_geom_px, fw),
        "target_robust_z_elongation": _robust_z(target.elongation, el),
        "control_sources": [q.to_dict() for q in controls],
    }
    return result


def run_fits_exact_chain(
    frozen_meta_pairs: list[dict],
    *,
    cache_dir: Path,
    events: EventWriter,
) -> tuple[list[dict], bool]:
    files = sorted({x["file_name"] for x in frozen_meta_pairs})
    if len(files) > MAX_FITS_DOWNLOAD_FILES:
        return [], True

    by_file: dict[str, list[dict]] = defaultdict(list)
    for p in frozen_meta_pairs:
        by_file[p["file_name"]].append(p)

    out = []
    for file_idx, fname in enumerate(files):
        group = by_file[fname]
        path, source_meta = download_source(
            group[0]["file_url"], cache_dir, events,
            target="PALOMAR_2F_D_EXACT_OVERLAP",
            filter_name=group[0].get("filters", "UNKNOWN"),
        )
        image = None
        image_meta = None
        gate = None
        for p in group:
            inside, wm = exact_wcs_pixel(path, p["ra_deg"], p["dec_deg"])
            rec = {
                **p,
                "file_sha256": source_meta["sha256"],
                "exact_fits_wcs_inside": bool(inside),
                "fits_wcs": wm,
            }
            if inside:
                if image is None:
                    image, image_meta = read_luci_fits_image(
                        path, require_imaging=True, expected_instrument=p.get("instrument")
                    )
                    gate = psf_relative_injection_recovery_gate(
                        image, seed=20262000 + file_idx
                    )
                rec["image_meta"] = image_meta
                rec["overlap_frame_injection_gate"] = gate
                if gate.get("passed"):
                    rec["counterpart_test"] = counterpart_with_matched_controls(
                        image, wm["x"], wm["y"]
                    )
                else:
                    rec["counterpart_test"] = {
                        "status": "BLOCKED_BY_OVERLAP_FRAME_R1_INJECTION_GATE"
                    }
            out.append(rec)
    out.sort(key=lambda x: (x["src_id"], x["file_name"]))
    return out, False


def main() -> int:
    ap = argparse.ArgumentParser(
        description="LUCI-PALOMAR-JPFM-2F-D: exhaustive archive-first LUCI coverage of frozen POSS-I S0"
    )
    ap.add_argument("--output-dir", default="results/luci_palomar_2f_d")
    ap.add_argument("--cache-dir", default=".cache/luci_palomar_2f_d")
    args = ap.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    parent = _require_parent_r1(repo_root)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache_dir)
    events = EventWriter(out / "events.jsonl")

    full = build_full_frozen_s0_coordinates(out / "frozen_poss_s0_122820.csv")
    events.emit("full_s0_frozen", sample_n=full["sample_n"], sha256=full["sha256"])

    inventory_rows = query_full_luci_inventory()
    inventory = write_inventory(inventory_rows, out / "frozen_luci_public_imaging_inventory.csv")
    events.emit("luci_inventory_frozen", row_count=inventory["row_count"], sha256=inventory["sha256"])

    coarse_rows, coarse_meta = archive_first_spatial_crossmatch(full["rows"], inventory_rows)
    coarse = write_pairs(coarse_rows, out / "frozen_coarse_overlap_pairs.csv")
    events.emit("coarse_overlap_set_frozen", pair_count=coarse["pair_count"], sha256=coarse["sha256"])

    meta_rows = metadata_wcs_preflight(coarse_rows, inventory_rows)
    meta = write_meta_pairs(meta_rows, out / "frozen_metadata_wcs_overlap_pairs.csv")
    events.emit("metadata_wcs_overlap_set_frozen", pair_count=meta["pair_count"], sha256=meta["sha256"])

    exact_results, cap_exceeded = run_fits_exact_chain(
        meta_rows, cache_dir=cache / "exact_fits", events=events
    )
    exact_inside = [x for x in exact_results if x.get("exact_fits_wcs_inside")]
    gate_pass = [x for x in exact_inside if x.get("overlap_frame_injection_gate", {}).get("passed")]
    gate_fail = [x for x in exact_inside if not x.get("overlap_frame_injection_gate", {}).get("passed")]
    counterparts = [x for x in gate_pass if x.get("counterpart_test", {}).get("counterpart_present")]
    no_counterparts = [x for x in gate_pass if not x.get("counterpart_test", {}).get("counterpart_present")]
    morph_available = [
        x for x in counterparts
        if x.get("counterpart_test", {}).get("morphology_status") == "MATCHED_LOCAL_CONTROL_COMPARISON_AVAILABLE"
    ]

    if cap_exceeded:
        status = "BLOCKED"
        scientific_status = "EXHAUSTIVE_COVERAGE_FROZEN__FITS_DOWNLOAD_CAP_EXCEEDED__NO_PIXEL_INSPECTION"
    elif not exact_inside:
        status = "PASS"
        scientific_status = "EXHAUSTIVE_FULL_S0_COVERAGE__NO_EXACT_PUBLIC_LUCI_FITS_OVERLAP"
    elif gate_fail:
        status = "BLOCKED"
        scientific_status = "EXACT_LUCI_OVERLAP_FOUND__ONE_OR_MORE_OVERLAP_FRAMES_FAILED_R1_GATE"
    else:
        status = "PASS"
        scientific_status = "EXACT_LUCI_OVERLAP_FOUND__R1_VALIDATED_COUNTERPART_CHAIN_EXECUTED"

    receipt = {
        "schema": "janus.cosmos.luci_palomar.jpfm_2f_d.receipt.v1",
        "experiment_id": "LUCI-PALOMAR-JPFM-2F-D-EXHAUSTIVE-COVERAGE",
        "status": status,
        "scientific_status": scientific_status,
        "parent_r1": {
            "path": PARENT_R1_RESULT,
            "artifact_sha256": PARENT_R1_ARTIFACT_SHA256,
            "validation_pass_fraction": parent["psf_injection_recovery_validation"]["pass_fraction"],
        },
        "frozen_full_poss_s0": {k: v for k, v in full.items() if k != "rows"},
        "frozen_luci_inventory": inventory,
        "archive_first_crossmatch": {
            **coarse_meta,
            "frozen_coarse_pairs": coarse,
            "frozen_metadata_wcs_pairs": meta,
            "method": "3D unit-vector cKDTree against LUCI frame centers followed by frame-specific half-diagonal filter; archive inventory frozen before crossmatch",
        },
        "exact_fits_chain": {
            "download_file_cap": MAX_FITS_DOWNLOAD_FILES,
            "cap_exceeded": cap_exceeded,
            "exact_fits_wcs_pair_count": len(exact_inside),
            "exact_fits_wcs_unique_sources": len({x["src_id"] for x in exact_inside}),
            "exact_fits_wcs_unique_files": len({x["file_name"] for x in exact_inside}),
            "r1_gate_pass_pair_count": len(gate_pass),
            "r1_gate_fail_pair_count": len(gate_fail),
            "counterpart_present_count": len(counterparts),
            "counterpart_absent_count": len(no_counterparts),
            "matched_morphology_available_count": len(morph_available),
            "results": exact_results,
        },
        "freeze_order": [
            "verify fixed POSS-I S0 gzip/csv hashes",
            "write and hash all 122820 Palomar coordinates",
            "query LUCI archive independently of Palomar coordinates",
            "write and hash LUCI public FREE SCIENCE Mirror inventory",
            "compute coarse spatial pairs and hash them",
            "compute archive-metadata WCS preflight pairs and hash them",
            "only then download candidate FITS",
            "exact FITS WCS containment",
            "per-overlap-frame R1 injection-recovery",
            "IR source/no-source using local R1 coordinate classifier",
            "PSF morphology vs preregistered same-frame SNR-matched local controls",
        ],
        "chain": "frozen full Palomar S0 -> frozen archive-first LUCI inventory -> frozen coarse pairs -> frozen metadata-WCS pairs -> exact FITS WCS -> overlap-frame R1 injection/recovery -> IR source/no-source -> PSF morphology -> matched local controls",
        "claim_ceiling": "INDEPENDENT_NEAR_IR_COUNTERPART_TEST_ONLY__NO_ANOMALY_OR_UAP_ORIGIN_CLAIM__NO_CAUSALITY",
    }
    rp = out / "receipt.json"
    rp.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": status,
        "scientific_status": scientific_status,
        "full_s0_n": full["sample_n"],
        "full_s0_sha256": full["sha256"],
        "luci_inventory_n": inventory["row_count"],
        "luci_inventory_sha256": inventory["sha256"],
        "coarse_pairs": coarse["pair_count"],
        "metadata_wcs_pairs": meta["pair_count"],
        "exact_fits_wcs_pairs": len(exact_inside),
        "r1_gate_fail_pairs": len(gate_fail),
        "counterpart_present": len(counterparts),
        "counterpart_absent": len(no_counterparts),
        "matched_morphology_available": len(morph_available),
        "receipt": str(rp),
    }, indent=2))
    return 0 if status == "PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
