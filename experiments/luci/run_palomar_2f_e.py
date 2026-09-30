#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

from janus_cosmos.luci import read_luci_fits_image
from janus_cosmos.luci_psf import _inject_gaussian, robust_background
from janus_cosmos.luci_psf_r1 import measure_psf_at, psf_relative_injection_recovery_gate
from janus_cosmos.pipeline import EventWriter, download_source, sha256_file
from experiments.luci.run_palomar_2f_c import exact_wcs_pixel
from experiments.luci.run_palomar_2f_d import counterpart_with_matched_controls

PARENT_RUN_ID = 31888129868
PARENT_ARTIFACT_ID = 9247843315
PARENT_ARTIFACT_ZIP_SHA256 = "1240a600dfb0189243dfb3188ab53dcc8ad6f7b270236c56756b49c2e4fc6184"
PARENT_META_WCS_SHA256 = "aa57e64deca11bb2afd09364fa8b72837e92fc3c2187357c5041b273023754c6"
PARENT_META_WCS_PAIRS = 918
PARENT_META_WCS_SOURCES = 64
PARENT_META_WCS_FILES = 798

EXPECTED_SKY_PAIRS = 443
EXPECTED_SKY_SOURCES = 42
EXPECTED_SKY_FILES = 403
HEADER_SHARD_SIZE = 225  # deterministic header-only replay; each shard remains <= old 250-file cap
MAX_REPRESENTATIVES_PER_SOURCE = 3
LOCAL_SNR_GRID = (8.0, 12.0)

PAIR_FIELDS = [
    "src_id", "ra_deg", "dec_deg", "file_name", "file_url", "instrument",
    "target", "filters", "date_obs", "center_sep_deg", "half_diagonal_deg",
    "metadata_wcs_x", "metadata_wcs_y",
]

EXACT_FIELDS = PAIR_FIELDS + [
    "file_sha256", "exact_hdu", "exact_x", "exact_y", "exact_naxis1", "exact_naxis2",
]
REP_FIELDS = EXACT_FIELDS + ["selection_pool", "temporal_role"]


def require_sha(path: Path, expected: str, label: str) -> None:
    got = sha256_file(path)
    if got != expected:
        raise RuntimeError(f"{label} SHA mismatch: {got} != {expected}")


def _read_pairs(path: Path) -> list[dict]:
    rows = list(csv.DictReader(path.open("r", encoding="utf-8", newline="")))
    if len(rows) != PARENT_META_WCS_PAIRS:
        raise RuntimeError(f"parent metadata-WCS row count changed: {len(rows)}")
    for r in rows:
        r["ra_deg"] = float(r["ra_deg"])
        r["dec_deg"] = float(r["dec_deg"])
        r["center_sep_deg"] = float(r["center_sep_deg"])
        r["half_diagonal_deg"] = float(r["half_diagonal_deg"])
        r["metadata_wcs_x"] = float(r["metadata_wcs_x"])
        r["metadata_wcs_y"] = float(r["metadata_wcs_y"])
    if len({r["src_id"] for r in rows}) != PARENT_META_WCS_SOURCES:
        raise RuntimeError("parent metadata-WCS unique source count changed")
    if len({r["file_name"] for r in rows}) != PARENT_META_WCS_FILES:
        raise RuntimeError("parent metadata-WCS unique file count changed")
    return rows


def exposure_admissible(row: dict) -> bool:
    s = str(row.get("filters", "") or "").lower()
    return "blind" not in s and "pv lens" not in s


def standard_broadband(filters: str) -> bool:
    return bool(re.fullmatch(r"clear\s+(?:z|j|h|k|ks)", str(filters or "").strip(), flags=re.I))


def freeze_rows(rows: list[dict], path: Path, fields: list[str]) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "pair_count": len(rows),
        "unique_sources": len({r["src_id"] for r in rows}),
        "unique_files": len({r["file_name"] for r in rows}),
    }


def deterministic_header_shards(rows: list[dict]) -> list[list[str]]:
    files = sorted({r["file_name"] for r in rows})
    return [files[i:i+HEADER_SHARD_SIZE] for i in range(0, len(files), HEADER_SHARD_SIZE)]


def exact_header_replay(
    rows: list[dict],
    *,
    cache_dir: Path,
    events: EventWriter,
) -> tuple[list[dict], dict[str, Path], list[dict]]:
    by_file: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_file[r["file_name"]].append(r)

    shards = deterministic_header_shards(rows)
    if any(len(s) > 250 for s in shards):
        raise RuntimeError("header shard exceeds inherited 250-file cap")

    cached: dict[str, Path] = {}
    exact: list[dict] = []
    file_receipts: list[dict] = []

    # Critical firewall: this stage calls exact_wcs_pixel only. No image data array is read.
    for shard_idx, names in enumerate(shards):
        events.emit("header_shard_started", shard_index=shard_idx, file_count=len(names))
        for fname in names:
            group = by_file[fname]
            first = group[0]
            path, meta = download_source(
                first["file_url"], cache_dir, events,
                target="PALOMAR_2F_E_HEADER_ONLY",
                filter_name=str(first.get("filters", "UNKNOWN")),
            )
            cached[fname] = path
            file_receipts.append({
                "file_name": fname,
                "file_url": first["file_url"],
                "sha256": meta["sha256"],
                "bytes": meta["bytes"],
                "header_shard": shard_idx,
            })
            for r in group:
                inside, wm = exact_wcs_pixel(path, r["ra_deg"], r["dec_deg"])
                if inside:
                    exact.append({
                        **r,
                        "file_sha256": meta["sha256"],
                        "exact_hdu": int(wm["hdu"]),
                        "exact_x": float(wm["x"]),
                        "exact_y": float(wm["y"]),
                        "exact_naxis1": int(wm["naxis1"]),
                        "exact_naxis2": int(wm["naxis2"]),
                    })
        events.emit("header_shard_completed", shard_index=shard_idx, file_count=len(names))

    exact.sort(key=lambda r: (r["src_id"], r["date_obs"], r["file_name"]))
    file_receipts.sort(key=lambda r: r["file_name"])
    return exact, cached, file_receipts


def temporal_representatives(exact: list[dict]) -> list[dict]:
    by_source: dict[str, list[dict]] = defaultdict(list)
    for r in exact:
        by_source[r["src_id"]].append(r)

    selected: list[dict] = []
    for sid in sorted(by_source):
        all_rows = by_source[sid]
        broad = [r for r in all_rows if standard_broadband(r.get("filters", ""))]
        pool = broad if broad else all_rows
        pool_name = "STANDARD_BROADBAND" if broad else "ANY_ADMISSIBLE_SKY_FILTER"
        # one row per file; exact pair keys are src_id,file_name but be explicit
        unique = {r["file_name"]: r for r in pool}
        seq = sorted(unique.values(), key=lambda r: ((r.get("date_obs") or "9999"), r["file_name"]))
        if not seq:
            continue
        idxs = sorted({0, len(seq)//2, len(seq)-1})
        roles_by_idx = {0: "earliest", len(seq)//2: "median", len(seq)-1: "latest"}
        for i in idxs:
            selected.append({**seq[i], "selection_pool": pool_name, "temporal_role": roles_by_idx[i]})
    selected.sort(key=lambda r: (r["src_id"], r["temporal_role"], r["file_name"]))
    if len(selected) > MAX_REPRESENTATIVES_PER_SOURCE * len(by_source):
        raise RuntimeError("temporal representative cardinality invariant failed")
    return selected


def local_coordinate_injection(image: np.ndarray, x: float, y: float, fwhm_px: float) -> dict:
    a = np.asarray(image, dtype=float)
    h, w = a.shape
    if not (12.0 <= x < w-12.0 and 12.0 <= y < h-12.0):
        return {"passed": False, "reason": "COORDINATE_TOO_CLOSE_TO_EDGE"}
    med, sigma = robust_background(a)
    base_fwhm = max(2.0, min(8.0, float(fwhm_px)))
    trials = []
    for snr in LOCAL_SNR_GRID:
        z = np.array(a, copy=True)
        _inject_gaussian(z, y, x, base_fwhm, snr*sigma)
        q = measure_psf_at(z, y, x)
        trials.append({"snr": snr, "recovered": q is not None})
    return {
        "passed": all(t["recovered"] for t in trials),
        "background_median": med,
        "background_sigma": sigma,
        "injection_fwhm_px": base_fwhm,
        "trials": trials,
    }


def pixel_replay(
    selected: list[dict],
    cached: dict[str, Path],
) -> list[dict]:
    by_file: dict[str, list[dict]] = defaultdict(list)
    for r in selected:
        by_file[r["file_name"]].append(r)

    results: list[dict] = []
    for file_idx, fname in enumerate(sorted(by_file)):
        group = by_file[fname]
        path = cached[fname]
        image, image_meta = read_luci_fits_image(
            path, require_imaging=True, expected_instrument=group[0].get("instrument")
        )
        gate = psf_relative_injection_recovery_gate(image, seed=20262100 + file_idx)
        for r in group:
            rec = {
                **r,
                "image_meta": image_meta,
                "overlap_frame_r1_gate": gate,
            }
            if not gate.get("passed"):
                rec["counterpart_test"] = {"status": "BLOCKED_BY_OVERLAP_FRAME_R1_GATE"}
                results.append(rec)
                continue

            target = measure_psf_at(image, r["exact_y"], r["exact_x"])
            if target is not None:
                rec["counterpart_test"] = counterpart_with_matched_controls(
                    image, r["exact_x"], r["exact_y"]
                )
                rec["local_coordinate_injection"] = {"not_required_for_presence": True}
            else:
                base_fwhm = gate.get("injection_base_fwhm_px", 2.0)
                local = local_coordinate_injection(
                    image, r["exact_x"], r["exact_y"], base_fwhm
                )
                rec["local_coordinate_injection"] = local
                if local.get("passed"):
                    rec["counterpart_test"] = {
                        "status": "NO_R1_PSF_SOURCE_WITH_LOCAL_SENSITIVITY_PASS",
                        "counterpart_present": False,
                    }
                else:
                    rec["counterpart_test"] = {
                        "status": "NO_SOURCE_RESULT_BLOCKED_BY_LOCAL_SENSITIVITY",
                        "counterpart_present": None,
                    }
            results.append(rec)
    results.sort(key=lambda r: (r["src_id"], r["temporal_role"], r["file_name"]))
    return results


def main() -> int:
    ap = argparse.ArgumentParser(description="2F-E staged exact-FITS-WCS replay and representative near-IR counterpart test")
    ap.add_argument("--parent-artifact-dir", required=True)
    ap.add_argument("--output-dir", default="results/luci_palomar_2f_e")
    ap.add_argument("--cache-dir", default=".cache/luci_palomar_2f_e")
    args = ap.parse_args()

    parent_dir = Path(args.parent_artifact_dir)
    parent_pairs = parent_dir / "frozen_metadata_wcs_overlap_pairs.csv"
    parent_receipt = parent_dir / "receipt.json"
    if not parent_pairs.exists() or not parent_receipt.exists():
        raise RuntimeError("parent 2F-D artifact contents missing")
    require_sha(parent_pairs, PARENT_META_WCS_SHA256, "2F-D metadata-WCS pair set")
    pd = json.loads(parent_receipt.read_text(encoding="utf-8"))
    if pd.get("scientific_status") != "EXHAUSTIVE_COVERAGE_FROZEN__FITS_DOWNLOAD_CAP_EXCEEDED__NO_PIXEL_INSPECTION":
        raise RuntimeError("parent 2F-D scientific status mismatch")

    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache_dir)
    events = EventWriter(out / "events.jsonl")

    parent_rows = _read_pairs(parent_pairs)
    sky = sorted((r for r in parent_rows if exposure_admissible(r)), key=lambda r: (r["src_id"], r["file_name"]))
    if len(sky) != EXPECTED_SKY_PAIRS or len({r['src_id'] for r in sky}) != EXPECTED_SKY_SOURCES or len({r['file_name'] for r in sky}) != EXPECTED_SKY_FILES:
        raise RuntimeError("exposure-admissibility audit cardinality changed")
    sky_meta = freeze_rows(sky, out / "frozen_exposure_admissible_pairs.csv", PAIR_FIELDS)
    events.emit("exposure_admissible_set_frozen", **sky_meta)

    shards = deterministic_header_shards(sky)
    shard_manifest = {
        "schema": "janus.cosmos.luci_palomar.header_shards.v1",
        "parent_metadata_wcs_sha256": PARENT_META_WCS_SHA256,
        "header_only": True,
        "pixel_access_before_exact_freeze": False,
        "shards": [{"index": i, "file_count": len(s), "files": s} for i,s in enumerate(shards)],
    }
    shard_path = out / "frozen_header_shards.json"
    shard_path.write_text(json.dumps(shard_manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    shard_sha = sha256_file(shard_path)
    events.emit("header_shards_frozen", shard_count=len(shards), sha256=shard_sha)

    exact, cached, file_receipts = exact_header_replay(sky, cache_dir=cache / "headers", events=events)
    exact_meta = freeze_rows(exact, out / "frozen_exact_fits_wcs_pairs.csv", EXACT_FIELDS)
    events.emit("exact_fits_wcs_set_frozen", **exact_meta)

    selected = temporal_representatives(exact)
    selected_meta = freeze_rows(selected, out / "frozen_temporal_representatives.csv", REP_FIELDS)
    events.emit("temporal_representatives_frozen", **selected_meta)

    # Pixel access starts only here, after exact and representative set SHA values exist.
    results = pixel_replay(selected, cached) if selected else []

    gate_fail = [r for r in results if not r.get("overlap_frame_r1_gate", {}).get("passed")]
    present = [r for r in results if r.get("counterpart_test", {}).get("counterpart_present") is True]
    absent = [r for r in results if r.get("counterpart_test", {}).get("counterpart_present") is False]
    local_block = [r for r in results if r.get("counterpart_test", {}).get("counterpart_present") is None]
    morph = [r for r in present if r.get("counterpart_test", {}).get("morphology_status") == "MATCHED_LOCAL_CONTROL_COMPARISON_AVAILABLE"]

    if gate_fail:
        status = "BLOCKED"
        scientific_status = "EXACT_FITS_WCS_FROZEN__REPRESENTATIVE_PIXEL_REPLAY_PARTIALLY_BLOCKED_BY_R1_GATE"
    elif not exact:
        status = "PASS"
        scientific_status = "STAGED_HEADER_REPLAY_COMPLETE__NO_EXACT_FITS_WCS_OVERLAP_IN_ADMISSIBLE_SET"
    else:
        status = "PASS"
        scientific_status = "EXACT_FITS_WCS_OVERLAPS_CONFIRMED__REPRESENTATIVE_R1_COUNTERPART_CHAIN_EXECUTED"

    receipt = {
        "schema": "janus.cosmos.luci_palomar.jpfm_2f_e.receipt.v1",
        "experiment_id": "LUCI-PALOMAR-JPFM-2F-E-STAGED-EXACT-WCS",
        "status": status,
        "scientific_status": scientific_status,
        "parent": {
            "workflow_run_id": PARENT_RUN_ID,
            "artifact_id": PARENT_ARTIFACT_ID,
            "artifact_zip_sha256": PARENT_ARTIFACT_ZIP_SHA256,
            "metadata_wcs_sha256": PARENT_META_WCS_SHA256,
        },
        "exposure_admissibility": {
            **sky_meta,
            "rule": "exclude any archive filters string containing blind or PV lens, case-insensitive",
        },
        "header_replay": {
            "shard_size": HEADER_SHARD_SIZE,
            "shard_count": len(shards),
            "shard_manifest_sha256": shard_sha,
            "pixel_access_before_exact_set_freeze": False,
            "downloaded_unique_files": len(file_receipts),
            "exact_fits_wcs": exact_meta,
            "file_receipts": file_receipts,
        },
        "representative_selection": {
            **selected_meta,
            "rule": "per exact-overlap source: prefer standard broad-band archive configurations matching clear z/J/H/K/Ks when available; otherwise use any exposure-admissible exact overlap; select deterministic earliest, median and latest unique files",
            "selection_before_pixel_access": True,
        },
        "pixel_replay": {
            "representative_pair_count": len(results),
            "unique_sources": len({r['src_id'] for r in results}),
            "unique_files": len({r['file_name'] for r in results}),
            "r1_gate_fail_count": len(gate_fail),
            "counterpart_present_count": len(present),
            "counterpart_absent_with_local_sensitivity_pass_count": len(absent),
            "no_source_local_sensitivity_block_count": len(local_block),
            "matched_morphology_available_count": len(morph),
            "results": results,
        },
        "claim_ceiling": "REPRESENTATIVE_EXACT_OVERLAP_NEAR_IR_COUNTERPART_TEST_ONLY__NOT_FULL_ARCHIVE_PHOTOMETRIC_COMPLETENESS__NO_ANOMALY_OR_UAP_ORIGIN_CLAIM__NO_CAUSALITY",
    }
    rp = out / "receipt.json"
    rp.write_text(json.dumps(receipt, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({
        "status": status,
        "scientific_status": scientific_status,
        "exposure_admissible_pairs": sky_meta["pair_count"],
        "exact_fits_wcs_pairs": exact_meta["pair_count"],
        "exact_fits_wcs_sources": exact_meta["unique_sources"],
        "exact_fits_wcs_files": exact_meta["unique_files"],
        "representative_pairs": selected_meta["pair_count"],
        "representative_sources": selected_meta["unique_sources"],
        "r1_gate_fail": len(gate_fail),
        "counterpart_present": len(present),
        "counterpart_absent_sensitive": len(absent),
        "local_sensitivity_block": len(local_block),
        "morphology_available": len(morph),
        "receipt": str(rp),
    }, indent=2))
    return 0 if status == "PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
