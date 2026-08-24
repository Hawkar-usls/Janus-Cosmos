#!/usr/bin/env python3
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

from janus_cosmos.luci import read_luci_fits_image
from janus_cosmos.luci_psf import detect_psf_sources, injection_recovery_gate
from janus_cosmos.pipeline import EventWriter, _source_for_item, download_source, load_manifest, sha256_file

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
CORPUS_PER_CLUSTER = 40
CORPUS_N = 16 * CORPUS_PER_CLUSTER
TYPICAL_PER_CLUSTER = 20
UNUSUAL_PER_CLUSTER = 20
SEARCH_RADIUS_DEG = 0.20
QUERY_CHUNK = 32
MAX_EXACT_CANDIDATE_FILES = 50


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _fetch(url: str, timeout: int = 180) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "JANUS-COSMOS-LUCI-JPFM-2F-C/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _require(label: str, b: bytes, expected: str) -> None:
    got = _sha(b)
    if got != expected:
        raise RuntimeError(f"{label} hash mismatch: got={got} expected={expected}")


def _s(v) -> str:
    try:
        if getattr(v, "mask", False):
            return ""
    except Exception:
        pass
    return "" if v is None else str(v).strip()


def _f(v, default=float("nan")) -> float:
    try:
        x = float(v)
        return x if math.isfinite(x) else default
    except Exception:
        return default


def tap_query(adql: str, timeout: int = 180) -> Table:
    body = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "QUERY": adql, "FORMAT": "votable"}).encode()
    req = urllib.request.Request(
        TAP_SYNC, data=body,
        headers={"User-Agent": "JANUS-COSMOS-LUCI-JPFM-2F-C/1.0", "Content-Type": "application/x-www-form-urlencoded"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return Table.read(io.BytesIO(r.read()), format="votable")


def build_frozen_expanded_corpus(out_csv: Path) -> dict:
    mgz = _fetch(REGISTRY_MANIFEST); _require("structural gzip", mgz, STRUCTURAL_GZ_SHA256)
    mcsv = gzip.decompress(mgz); _require("structural csv", mcsv, STRUCTURAL_CSV_SHA256)
    rows = list(csv.DictReader(io.StringIO(mcsv.decode("utf-8"))))
    groups: dict[int, list[dict]] = defaultdict(list)
    for r in rows:
        groups[int(float(r["structural_cluster"]))].append(r)
    chosen = []
    for cluster in range(16):
        g = groups[cluster]
        med = statistics.median(float(r["anomaly_score"]) for r in g)
        typical = sorted(g, key=lambda r: (abs(float(r["anomaly_score"]) - med), r["src_id"]))[:TYPICAL_PER_CLUSTER]
        used = {r["src_id"] for r in typical}
        unusual = sorted((r for r in g if r["src_id"] not in used), key=lambda r: (-float(r["anomaly_score"]), r["src_id"]))[:UNUSUAL_PER_CLUSTER]
        chosen += [{"src_id": r["src_id"], "structural_cluster": cluster, "sample_role": "typical"} for r in typical]
        chosen += [{"src_id": r["src_id"], "structural_cluster": cluster, "sample_role": "unusual"} for r in unusual]
    if len(chosen) != CORPUS_N or len({x["src_id"] for x in chosen}) != CORPUS_N:
        raise RuntimeError("expanded Palomar corpus invariant failed")

    s0gz = _fetch(S0_URL); _require("S0 gzip", s0gz, S0_GZ_SHA256)
    s0csv = gzip.decompress(s0gz); _require("S0 csv", s0csv, S0_CSV_SHA256)
    wanted = {x["src_id"] for x in chosen}; pos = {}
    for r in csv.DictReader(io.StringIO(s0csv.decode("utf-8"))):
        if r.get("src_id") in wanted:
            pos[r["src_id"]] = (float(r["ra"]), float(r["dec"]))
    if len(pos) != CORPUS_N:
        raise RuntimeError(f"expanded position join incomplete: {len(pos)}/{CORPUS_N}")
    for x in chosen:
        x["ra_deg"], x["dec_deg"] = pos[x["src_id"]]
    chosen.sort(key=lambda x: (x["structural_cluster"], x["sample_role"], x["src_id"]))

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["src_id", "structural_cluster", "sample_role", "ra_deg", "dec_deg"])
        w.writeheader(); w.writerows(chosen)
    return {
        "rows": chosen,
        "path": str(out_csv),
        "sha256": sha256_file(out_csv),
        "sample_n": len(chosen),
        "selection": "20 median-near typical + 20 highest structural-anomaly unusual per each of 16 frozen clusters; LUCI outcomes unavailable at selection time",
        "registry_commit": REGISTRY_COMMIT,
        "poss_commit": POSS_COMMIT,
        "structural_manifest_gzip_sha256": _sha(mgz),
        "structural_manifest_csv_sha256": _sha(mcsv),
        "stage_S0_gzip_sha256": _sha(s0gz),
        "stage_S0_csv_sha256": _sha(s0csv),
    }


def angular_sep_deg(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    a1, d1, a2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    z = math.sin((d2-d1)/2)**2 + math.cos(d1)*math.cos(d2)*math.sin((a2-a1)/2)**2
    return math.degrees(2*math.asin(min(1.0, math.sqrt(max(0.0, z)))))


def query_conservative_overlaps(corpus: list[dict]) -> tuple[list[dict], dict]:
    overlaps = []
    tap_rows = 0
    for start in range(0, len(corpus), QUERY_CHUNK):
        chunk = corpus[start:start+QUERY_CHUNK]
        clauses = [f"1=CONTAINS(lbt.lbt.s_point,CIRCLE('ICRS',{x['ra_deg']:.10f},{x['dec_deg']:.10f},{SEARCH_RADIUS_DEG:.6f}))" for x in chunk]
        adql = (
            "SELECT TOP 5000 lbt.luci.instrument,lbt.luci.telescope,lbt.luci.object,lbt.luci.filters,"
            "lbt.luci.gratname,lbt.luci.imagetyp,lbt.luci.file_name,lbt.luci.date_obs,lbt.luci.crval1,lbt.luci.crval2,"
            "lbt.luci.naxis1,lbt.luci.naxis2,lbt.luci.pixscale,lbt.lbt.file_url,lbt.lbt.policy "
            "FROM lbt.luci JOIN lbt.lbt ON lbt.luci.file_name=lbt.lbt.file_name "
            "WHERE lbt.lbt.policy='FREE' AND lbt.luci.imagetyp='SCIENCE' "
            "AND (lbt.luci.gratname='Mirror' OR lbt.luci.gratname='mirror') AND (" + " OR ".join(clauses) + ")"
        )
        rows = tap_query(adql); tap_rows += len(rows)
        for row in rows:
            cra, cdec = _f(row["crval1"]), _f(row["crval2"])
            nx, ny, scale = _f(row["naxis1"]), _f(row["naxis2"]), _f(row["pixscale"])
            if not all(math.isfinite(v) for v in (cra,cdec,nx,ny,scale)) or scale <= 0:
                continue
            half_diag = 0.5 * math.hypot(nx*scale, ny*scale)/3600.0
            if not (0 < half_diag < 0.5):
                continue
            for src in chunk:
                sep = angular_sep_deg(src["ra_deg"], src["dec_deg"], cra, cdec)
                if sep <= half_diag:
                    overlaps.append({
                        **src, "file_name": _s(row["file_name"]), "file_url": _s(row["file_url"]),
                        "instrument": _s(row["instrument"]), "target": _s(row["object"]), "filters": _s(row["filters"]),
                        "date_obs": _s(row["date_obs"]), "center_sep_deg": sep, "half_diagonal_deg": half_diag,
                    })
    dedup = {}
    for x in overlaps:
        dedup[(x["src_id"], x["file_name"])] = x
    out = sorted(dedup.values(), key=lambda x: (x["src_id"], x["file_name"]))
    return out, {"tap_rows_returned_across_chunks": tap_rows, "conservative_overlap_pairs": len(out), "unique_sources": len({x['src_id'] for x in out})}


def exact_wcs_pixel(path: Path, ra: float, dec: float) -> tuple[bool, dict]:
    with fits.open(path, memmap=False, lazy_load_hdus=True) as hdul:
        for idx, hdu in enumerate(hdul):
            hdr = hdu.header
            nx, ny = int(hdr.get("NAXIS1", 0) or 0), int(hdr.get("NAXIS2", 0) or 0)
            if nx <= 0 or ny <= 0:
                continue
            try:
                w = WCS(hdr).celestial
                if not w.has_celestial:
                    continue
                x, y = w.world_to_pixel_values(float(ra), float(dec))
                x, y = float(x), float(y)
                if math.isfinite(x) and math.isfinite(y):
                    inside = (-0.5 <= x < nx-0.5) and (-0.5 <= y < ny-0.5)
                    return bool(inside), {"hdu": idx, "x": x, "y": y, "naxis1": nx, "naxis2": ny}
            except Exception:
                continue
    return False, {"reason": "NO_USABLE_CELESTIAL_WCS"}


def local_counterpart(image: np.ndarray, x: float, y: float) -> dict:
    sources = detect_psf_sources(image)
    if not sources:
        return {"status": "NO_PSF_SOURCES_IN_FRAME", "counterpart_present": False, "source_count": 0}
    ds = sorted(((q.x-x)**2 + (q.y-y)**2, q) for q in sources)
    d2, q = ds[0]
    if d2 > 3.0**2:
        return {"status": "NO_IR_SOURCE_WITHIN_3PX", "counterpart_present": False, "source_count": len(sources), "nearest_distance_px": math.sqrt(d2)}
    controls = [s for _, s in ds[1:] if (s.x-x)**2 + (s.y-y)**2 <= 300.0**2][:20]
    out = {"status": "COUNTERPART_DETECTED", "counterpart_present": True, "source_count": len(sources), "distance_px": math.sqrt(d2), "source": q.to_dict(), "local_control_count": len(controls)}
    if len(controls) < 5:
        out["morphology_status"] = "INSUFFICIENT_LOCAL_CONTROLS"
        return out
    def rz(value: float, vals: list[float]) -> float | None:
        med = float(np.median(vals)); mad = float(np.median(np.abs(np.asarray(vals)-med))); scale = 1.4826*mad
        return None if scale <= 0 or not math.isfinite(scale) else float((value-med)/scale)
    fw = [s.fwhm_geom_px for s in controls]; el = [s.elongation for s in controls]
    out["morphology_status"] = "LOCAL_CONTROL_COMPARISON_AVAILABLE"
    out["local_controls"] = {"median_fwhm_geom_px": float(np.median(fw)), "median_elongation": float(np.median(el)), "target_robust_z_fwhm": rz(q.fwhm_geom_px, fw), "target_robust_z_elongation": rz(q.elongation, el)}
    return out


def validate_manifest_frames(paths: list[Path], cache: Path, events: EventWriter) -> dict:
    frames = []
    for mp in paths:
        manifest = load_manifest(mp)
        for target in manifest["targets"]:
            for item in target["filters"]:
                filt = str(item.get("filter", "UNKNOWN")); src = _source_for_item(item)
                path, meta = download_source(src, cache, events, target=target["target"], filter_name=filt)
                image, imeta = read_luci_fits_image(path, require_imaging=True, expected_instrument=item.get("instrument"))
                gate = injection_recovery_gate(image, seed=20260815 + len(frames))
                frames.append({"target": target["target"], "filter": filt, "file_sha256": meta["sha256"], "image_meta": imeta, "injection_gate": gate})
    pass_n = sum(x["injection_gate"]["passed"] for x in frames)
    return {"frame_count": len(frames), "passed_frame_count": pass_n, "pass_fraction": pass_n/len(frames) if frames else 0.0, "all_frames_pass": bool(frames and pass_n == len(frames)), "frames": frames}


def main() -> int:
    ap = argparse.ArgumentParser(description="LUCI-PALOMAR-JPFM-2F-C: PSF-aware injection/recovery + expanded frozen Palomar WCS test")
    ap.add_argument("--validation-manifest", action="append", required=True)
    ap.add_argument("--output-dir", default="results/luci_palomar_2f_c")
    ap.add_argument("--cache-dir", default=".cache/luci_palomar_2f_c")
    args = ap.parse_args()
    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache_dir); events = EventWriter(out / "events.jsonl")

    corpus_meta = build_frozen_expanded_corpus(out / "frozen_palomar_640.csv")
    # Temporal firewall: corpus and SHA are materialized before the first LUCI-overlap query.
    events.emit("palomar_corpus_frozen", sample_n=corpus_meta["sample_n"], sha256=corpus_meta["sha256"])
    validation = validate_manifest_frames([Path(x) for x in args.validation_manifest], cache / "validation", events)
    conservative, overlap_meta = query_conservative_overlaps(corpus_meta["rows"])

    exact_results = []
    unique_files = sorted({x["file_name"] for x in conservative})
    cap_exceeded = len(unique_files) > MAX_EXACT_CANDIDATE_FILES
    if not cap_exceeded:
        by_file: dict[str, list[dict]] = defaultdict(list)
        for x in conservative: by_file[x["file_name"]].append(x)
        for fname in unique_files:
            group = by_file[fname]; url = group[0]["file_url"]
            path, meta = download_source(url, cache / "overlaps", events, target="PALOMAR_OVERLAP", filter_name=group[0].get("filters", "UNKNOWN"))
            for cand in group:
                inside, wmeta = exact_wcs_pixel(path, cand["ra_deg"], cand["dec_deg"])
                rec = {**cand, "exact_wcs_inside": inside, "wcs": wmeta, "file_sha256": meta["sha256"]}
                if inside:
                    image, imeta = read_luci_fits_image(path, require_imaging=True, expected_instrument=cand.get("instrument"))
                    gate = injection_recovery_gate(image, seed=20261815 + len(exact_results))
                    rec["overlap_frame_injection_gate"] = gate
                    rec["image_meta"] = imeta
                    if gate["passed"]:
                        rec["counterpart_test"] = local_counterpart(image, wmeta["x"], wmeta["y"])
                    else:
                        rec["counterpart_test"] = {"status": "BLOCKED_BY_INJECTION_RECOVERY_GATE"}
                exact_results.append(rec)

    exact_inside = [x for x in exact_results if x.get("exact_wcs_inside")]
    admitted_counterparts = [x for x in exact_inside if x.get("overlap_frame_injection_gate", {}).get("passed") and x.get("counterpart_test", {}).get("counterpart_present")]
    status = "PASS"
    scientific = "NO_EXACT_LUCI_OVERLAP_IN_EXPANDED_FROZEN_CORPUS"
    if cap_exceeded:
        status = "BLOCKED"; scientific = "OVERLAP_CANDIDATE_CAP_EXCEEDED__NO_INSPECTION"
    elif not validation["all_frames_pass"]:
        status = "BLOCKED"; scientific = "PSF_INJECTION_RECOVERY_VALIDATION_FAILED"
    elif exact_inside:
        scientific = "EXACT_OVERLAP_AVAILABLE__COUNTERPART_CHAIN_EXECUTED"

    receipt = {
        "schema": "janus.cosmos.luci_palomar.jpfm_2f_c.receipt.v1",
        "experiment_id": "LUCI-PALOMAR-JPFM-2F-C",
        "status": status,
        "scientific_status": scientific,
        "frozen_palomar_corpus": {k:v for k,v in corpus_meta.items() if k != "rows"},
        "validation_injection_recovery": validation,
        "archive_overlap_preflight": {**overlap_meta, "search_radius_deg": SEARCH_RADIUS_DEG, "query_chunk_size": QUERY_CHUNK, "exact_candidate_file_cap": MAX_EXACT_CANDIDATE_FILES, "candidate_cap_exceeded": cap_exceeded},
        "exact_wcs_overlap_pair_count": len(exact_inside),
        "exact_wcs_unique_palomar_sources": len({x["src_id"] for x in exact_inside}),
        "admitted_ir_counterpart_count": len(admitted_counterparts),
        "exact_results": exact_results,
        "chain": "Palomar coordinate -> exact LUCI FITS WCS containment -> frame injection-recovery gate -> IR source/no-source -> PSF morphology -> matched local controls",
        "firewall": "Expanded Palomar corpus selection is frozen from Palomar-only structural data before LUCI archive outcomes are queried.",
        "claim_ceiling": "INDEPENDENT_NEAR_IR_COUNTERPART_TEST_ONLY__NO_ANOMALY_OR_UAP_ORIGIN_CLAIM__NO_CAUSALITY",
    }
    rp = out / "receipt.json"; rp.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "scientific_status": scientific, "validation_pass_fraction": validation["pass_fraction"], "palomar_corpus_n": CORPUS_N, "conservative_pairs": len(conservative), "exact_wcs_pairs": len(exact_inside), "counterparts": len(admitted_counterparts), "receipt": str(rp)}, indent=2))
    return 0 if status == "PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
