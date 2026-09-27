#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime
from pathlib import Path

UPSTREAM_COMMIT = "4005e200541b321ead3d6608f0162a14430ef1c2"
RELEASE = "results/s0-642-20260814"
RAW_BASE = f"https://raw.githubusercontent.com/jannefi/poss1-plate-slice/{UPSTREAM_COMMIT}"
STAGE_URL = f"{RAW_BASE}/{RELEASE}/stage_S0.csv.gz"
TILE_URL = f"{RAW_BASE}/{RELEASE}/tile_manifest.csv.gz"
REPAIRED_URL = f"{RAW_BASE}/{RELEASE}/repaired_astrometry_tiles.csv"
PLATE_MANIFEST_URL = f"{RAW_BASE}/data/plate_manifest.csv"
IRSA_BASE = "https://irsa.ipac.caltech.edu/data/DSS/images/dss1red"
STUDY_START = date(1949, 11, 19)
STUDY_END = date(1957, 4, 28)
EXPECTED = {
    "stage_rows": 122820,
    "tile_rows": 31458,
    "plate_count": 642,
    "plates_in_window": 640,
    "candidate_rows_in_window": 122511,
    "unique_observation_days": 312,
    "stage_sha256_uncompressed": "2ff92f2210acb387ef9ef4b88d561595d3883e9aab27065042627272b96590f0",
    "tile_sha256_uncompressed": "5dcb90dc5d98550e5a60246aced2b097922a267c69e81f27d45d16a288142a99",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def get(url: str, *, max_bytes: int | None = None, range_header: str | None = None) -> bytes:
    headers = {"User-Agent": "JANUS-COSMOS-PALOMAR-DAY-MATERIALIZER/1.0"}
    if range_header:
        headers["Range"] = range_header
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read() if max_bytes is None else r.read(max_bytes)


def parse_csv_bytes(raw: bytes) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))


def gunzip_checked(blob: bytes, expected: str, label: str) -> bytes:
    raw = gzip.decompress(blob)
    got = sha256(raw)
    if got != expected:
        raise RuntimeError(f"{label} uncompressed SHA mismatch: {got} != {expected}")
    return raw


def parse_fits_header(data: bytes) -> dict[str, str]:
    out: dict[str, str] = {}
    for pos in range(0, len(data) - 79, 80):
        card = data[pos:pos + 80].decode("ascii", errors="ignore")
        key = card[:8].strip()
        if key == "END":
            break
        if not key or card[8:10] != "= ":
            continue
        value = card[10:80].split("/", 1)[0].strip()
        if value.startswith("'") and "'" in value[1:]:
            value = value[1:value.find("'", 1)]
        out[key] = value.strip()
    return out


def parse_obs_date(value: str) -> date:
    v = (value or "").strip().strip("'")
    for fmt in ("%Y-%m-%d", "%d/%m/%y", "%m/%d/%y", "%d/%m/%Y", "%m/%d/%Y"):
        try:
            d = datetime.strptime(v[:10], fmt).date()
            if d.year >= 2000 and "/" in v and len(v.split("/")[-1][:2]) == 2:
                d = d.replace(year=d.year - 100)
            return d
        except ValueError:
            pass
    raise ValueError(f"cannot parse DATE-OBS={value!r}")


def as_float(v: str | None) -> float | None:
    try:
        return float(str(v).strip().replace("D", "E"))
    except Exception:
        return None


def fetch_plate_meta(plate_id: str) -> dict[str, object]:
    url = f"{IRSA_BASE}/dss1red_{plate_id}.fits"
    raw = get(url, max_bytes=131072, range_header="bytes=0-131071")
    hdr = parse_fits_header(raw)
    if not hdr:
        raise RuntimeError(f"no FITS header for {plate_id}")
    region = (hdr.get("REGION") or plate_id).strip()
    if region and region != plate_id:
        raise RuntimeError(f"plate identity mismatch {plate_id} vs REGION={region}")
    obs = parse_obs_date(hdr.get("DATE-OBS", ""))
    exposure = as_float(hdr.get("EXPOSURE"))
    if exposure is None or exposure <= 0:
        raise RuntimeError(f"missing/nonpositive EXPOSURE for {plate_id}")
    return {
        "plate_id": plate_id,
        "date": obs.isoformat(),
        "exposure_minutes": exposure,
        "date_obs_raw": hdr.get("DATE-OBS", ""),
        "emulsion": hdr.get("EMULSION", ""),
        "filter": hdr.get("FILTER", ""),
        "telescop": hdr.get("TELESCOP", ""),
        "source_url": url,
    }


def fetch_all_plate_meta(plates: list[str]) -> list[dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=12, thread_name_prefix="irsa") as pool:
        futures = {pool.submit(fetch_plate_meta, p): p for p in plates}
        for i, f in enumerate(as_completed(futures), 1):
            p = futures[f]
            try:
                out[p] = f.result()
            except Exception as e:
                errors.append(f"{p}:{type(e).__name__}:{e}")
            if i % 50 == 0 or i == len(plates):
                print(f"IRSA_HEADERS {i}/{len(plates)} errors={len(errors)}", flush=True)
    if errors:
        raise RuntimeError("plate metadata fail-closed: " + " | ".join(sorted(errors)[:20]))
    return [out[p] for p in sorted(plates)]


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    stage_gz = get(STAGE_URL)
    tile_gz = get(TILE_URL)
    repaired_raw = get(REPAIRED_URL)
    plate_manifest_raw = get(PLATE_MANIFEST_URL)

    stage_raw = gunzip_checked(stage_gz, EXPECTED["stage_sha256_uncompressed"], "stage_S0")
    tile_raw = gunzip_checked(tile_gz, EXPECTED["tile_sha256_uncompressed"], "tile_manifest")
    stage_rows = parse_csv_bytes(stage_raw)
    tile_rows = parse_csv_bytes(tile_raw)
    repaired_rows = parse_csv_bytes(repaired_raw)
    plate_manifest = parse_csv_bytes(plate_manifest_raw)

    assert len(stage_rows) == EXPECTED["stage_rows"], len(stage_rows)
    assert len(tile_rows) == EXPECTED["tile_rows"], len(tile_rows)
    assert len(plate_manifest) == EXPECTED["plate_count"], len(plate_manifest)

    plate_center: dict[str, tuple[float, float]] = {}
    for r in plate_manifest:
        p = (r.get("plate_id") or "").strip()
        ra, dec = as_float(r.get("ra_deg")), as_float(r.get("dec_deg"))
        if not p or ra is None or dec is None:
            raise RuntimeError(f"bad plate manifest row: {r}")
        plate_center[p] = (ra % 360.0, dec)

    tile_to_plate: dict[str, str] = {}
    tiles_per_plate: Counter[str] = Counter()
    for r in tile_rows:
        tile = (r.get("tile_id") or "").strip()
        plate = (r.get("plate_id") or "").strip()
        if not tile or not plate:
            raise RuntimeError("blank tile_id/plate_id")
        prior = tile_to_plate.setdefault(tile, plate)
        if prior != plate:
            raise RuntimeError(f"ambiguous tile mapping {tile}: {prior} vs {plate}")
        tiles_per_plate[plate] += 1

    candidates_per_plate: Counter[str] = Counter()
    for r in stage_rows:
        tile = (r.get("tile_id") or "").strip()
        plate = tile_to_plate.get(tile)
        if not plate:
            raise RuntimeError(f"stage row references unknown tile {tile}")
        candidates_per_plate[plate] += 1

    repaired_tiles_per_plate: Counter[str] = Counter()
    repaired_before_per_plate: Counter[str] = Counter()
    repaired_after_per_plate: Counter[str] = Counter()
    for r in repaired_rows:
        if (r.get("action") or "").strip().lower() != "repaired":
            continue
        tile = (r.get("tile_id") or "").strip()
        plate = (r.get("plate") or tile_to_plate.get(tile) or "").strip()
        if not plate:
            continue
        repaired_tiles_per_plate[plate] += 1
        repaired_before_per_plate[plate] += int(float(r.get("rows_before") or 0))
        repaired_after_per_plate[plate] += int(float(r.get("rows_after") or 0))

    plates = sorted(plate_center)
    meta = fetch_all_plate_meta(plates)
    meta_by_plate = {str(r["plate_id"]): r for r in meta}

    in_window = [p for p in plates if STUDY_START <= date.fromisoformat(str(meta_by_plate[p]["date"])) <= STUDY_END]
    if len(in_window) != EXPECTED["plates_in_window"]:
        raise RuntimeError(f"plates in window {len(in_window)} != {EXPECTED['plates_in_window']}")
    if sum(candidates_per_plate[p] for p in in_window) != EXPECTED["candidate_rows_in_window"]:
        raise RuntimeError("candidate rows in window mismatch")

    grouped: defaultdict[str, list[str]] = defaultdict(list)
    for p in in_window:
        grouped[str(meta_by_plate[p]["date"])].append(p)
    if len(grouped) != EXPECTED["unique_observation_days"]:
        raise RuntimeError(f"unique observation days {len(grouped)} != {EXPECTED['unique_observation_days']}")

    rows = []
    for ds in sorted(grouped):
        ps = sorted(grouped[ds])
        candidate_count = sum(candidates_per_plate[p] for p in ps)
        tile_count = sum(tiles_per_plate[p] for p in ps)
        plate_exposure = sum(float(meta_by_plate[p]["exposure_minutes"]) for p in ps)
        exposure_tile_min = sum(float(meta_by_plate[p]["exposure_minutes"]) * tiles_per_plate[p] for p in ps)
        repaired_count = sum(repaired_tiles_per_plate[p] for p in ps)
        repaired_before = sum(repaired_before_per_plate[p] for p in ps)
        repaired_after = sum(repaired_after_per_plate[p] for p in ps)
        centers = [{"plate_id": p, "ra_deg": round(plate_center[p][0], 6), "dec_deg": round(plate_center[p][1], 6)} for p in ps]
        rows.append({
            "date": ds,
            "observed_flag": True,
            "candidate_count": candidate_count,
            "plate_count": len(ps),
            "tile_count": tile_count,
            "exposure": {
                "plate_exposure_minutes_sum": round(plate_exposure, 6),
                "processed_tile_minutes": round(exposure_tile_min, 6),
                "contract": "processed tile-time denominator; not exact non-overlapping geometric sky area-time"
            },
            "sky_position": {
                "plate_centers_ra_dec": centers,
                "ra_deg_min": round(min(plate_center[p][0] for p in ps), 6),
                "ra_deg_max": round(max(plate_center[p][0] for p in ps), 6),
                "dec_deg_min": round(min(plate_center[p][1] for p in ps), 6),
                "dec_deg_max": round(max(plate_center[p][1] for p in ps), 6)
            },
            "detector_quality": {
                "known_repaired_astrometry_tiles": repaired_count,
                "rows_before_repair_on_known_repaired_tiles": repaired_before,
                "rows_after_repair_on_known_repaired_tiles": repaired_after,
                "plate_ids_with_known_repairs": [p for p in ps if repaired_tiles_per_plate[p]],
            },
            "artifact_state": "KNOWN_REPAIRED_ASTROMETRY_TILE_PRESENT" if repaired_count else "NO_BOUND_REPAIR_LOG_TILE_ON_THIS_DAY__NOT_ARTIFACT_FREE_CLAIM",
            "plate_ids": ps,
        })

    total_candidates = sum(r["candidate_count"] for r in rows)
    total_tiles = sum(r["tile_count"] for r in rows)
    assert total_candidates == EXPECTED["candidate_rows_in_window"]
    assert len(rows) == EXPECTED["unique_observation_days"]

    out = {
        "schema": "janus.cosmos.space.palomar_observation_day_ledger.v0.1",
        "artifact_id": "JANUS-COSMOS-PALOMAR-OBSERVATION-DAY-LEDGER-v0.1",
        "status": "MATERIALIZED_SOURCE_BOUND_312_DAY_LEDGER",
        "study_window": [STUDY_START.isoformat(), STUDY_END.isoformat()],
        "representation_contract": {
            "row_unit": "unique UTC observation date",
            "zero_law": "absence of a row outside this ledger must never be interpreted as an observed zero",
            "observed_zero_law": "an observed zero is legal only if an observation-day row exists with observed_flag=true and candidate_count=0",
            "sky_contract": "RA/Dec are angular plate-center coordinates on the celestial sphere, not physical distances",
            "artifact_contract": "repair-log absence is not proof of artifact absence"
        },
        "source_freeze": {
            "upstream_repository": "jannefi/poss1-plate-slice",
            "upstream_commit": UPSTREAM_COMMIT,
            "release": RELEASE,
            "stage_url": STAGE_URL,
            "tile_manifest_url": TILE_URL,
            "repaired_astrometry_url": REPAIRED_URL,
            "plate_manifest_url": PLATE_MANIFEST_URL,
            "irsa_plate_base": IRSA_BASE,
            "stage_sha256_gzip": sha256(stage_gz),
            "stage_sha256_uncompressed": sha256(stage_raw),
            "tile_manifest_sha256_gzip": sha256(tile_gz),
            "tile_manifest_sha256_uncompressed": sha256(tile_raw),
            "repaired_astrometry_sha256": sha256(repaired_raw),
            "plate_manifest_sha256": sha256(plate_manifest_raw),
            "irsa_headers_parsed": len(meta)
        },
        "counts": {
            "rows": len(rows),
            "plates_in_window": len(in_window),
            "candidate_rows_in_window": total_candidates,
            "processed_tiles_in_window": total_tiles,
            "all_observed_days_candidate_positive": all(r["candidate_count"] > 0 for r in rows)
        },
        "rows": rows,
        "canonical_seal": "OBSERVATION OPPORTUNITY IS MATERIALIZED BEFORE ANY EARTH↔SPACE JOIN. NO OBSERVATION IS NEVER SILENTLY TURNED INTO ZERO."
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("PALOMAR_DAY_LEDGER", json.dumps(out["counts"], sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
