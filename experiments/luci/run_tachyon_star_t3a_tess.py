#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
from astropy.io import fits

SECTOR_MANIFEST_SHA256 = "9c2fade3139a194fb14b1971817520f6531fff94aebe0243a10103bcf386a849"
EXPECTED_SOURCES = 42
CUTOUT_X = 11
CUTOUT_Y = 11
ALLOWED_VALUE_COLUMNS = ("TIME", "QUALITY", "CADENCENO")
FORBIDDEN_VALUE_COLUMNS = ("FLUX", "FLUX_ERR", "FLUX_BKG", "FLUX_BKG_ERR")
GAP_LOW = 0.5
GAP_HIGH = 1.5
CLAIM = "TESS_METADATA_ONLY_CADENCE_PREPOINTING__NO_TARGET_FLUX_INSPECTION__NO_TRANSIENT_DETECTION__NO_FTL__NO_RETROCAUSALITY__NO_TACHYON_IDENTITY"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _fetch(url: str, *, timeout: int = 600, retries: int = 2) -> bytes:
    err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "Janus-Cosmos-TachyonStar-T3A/1.0", "Accept": "application/zip,application/json,*/*"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:
            err = f"{type(exc).__name__}: {exc}"
            if attempt < retries:
                time.sleep(3 * attempt)
    raise RuntimeError(f"GET failed after {retries} attempts: {url}: {err}")


def tesscut_url(row: dict) -> str:
    return "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?" + urllib.parse.urlencode({
        "ra": row["ra_deg"],
        "dec": row["dec_deg"],
        "x": CUTOUT_X,
        "y": CUTOUT_Y,
        "units": "px",
        "sector": int(row["sector"]),
    })


def _extract_single_fits(payload: bytes, target_dir: Path, safe: str) -> tuple[Path, str]:
    target_dir.mkdir(parents=True, exist_ok=True)
    bio = io.BytesIO(payload)
    if not zipfile.is_zipfile(bio):
        sample = payload[:300].decode("utf-8", errors="replace")
        raise RuntimeError(f"TESSCut response is not ZIP: {sample}")
    bio.seek(0)
    with zipfile.ZipFile(bio) as z:
        names = [n for n in z.namelist() if n.lower().endswith(".fits") and not n.endswith("/")]
        if len(names) != 1:
            raise RuntimeError(f"expected exactly one FITS for fixed sector, got {len(names)}: {names}")
        name = names[0]
        data = z.read(name)
    out = target_dir / f"{safe}.fits"
    out.write_bytes(data)
    return out, name


def read_time_quality_only(path: Path) -> dict:
    """Read only TIME/QUALITY/CADENCENO values. Science array columns are forbidden in T3A."""
    with fits.open(path, mode="readonly", memmap=True, lazy_load_hdus=True) as hdul:
        chosen = None
        for i, hdu in enumerate(hdul):
            cols = list(getattr(getattr(hdu, "columns", None), "names", []) or [])
            upper = {str(c).upper() for c in cols}
            if "TIME" in upper and "QUALITY" in upper:
                chosen = (i, hdu, cols)
                break
        if chosen is None:
            raise RuntimeError("no binary table containing TIME and QUALITY")
        hdu_index, hdu, cols = chosen
        upper_map = {str(c).upper(): str(c) for c in cols}
        # Firewall: only these columns are dereferenced. FLUX-like columns are never read here.
        data = hdu.data
        time_values = np.asarray(data[upper_map["TIME"]], dtype=float).copy()
        quality_values = np.asarray(data[upper_map["QUALITY"]], dtype=np.int64).copy()
        if "CADENCENO" in upper_map:
            cadence_values = np.asarray(data[upper_map["CADENCENO"]], dtype=np.int64).copy()
        else:
            cadence_values = np.arange(len(time_values), dtype=np.int64)
        return {
            "hdu": int(hdu_index),
            "columns": cols,
            "time": time_values,
            "quality": quality_values,
            "cadenceno": cadence_values,
        }


def eligible_triples(time_values: np.ndarray, quality: np.ndarray, expected_cadence_s: float) -> list[dict]:
    t = np.asarray(time_values, dtype=float)
    q = np.asarray(quality, dtype=np.int64)
    if len(t) != len(q) or len(t) < 3:
        return []
    expected_days = float(expected_cadence_s) / 86400.0
    out = []
    for b in range(1, len(t) - 1):
        a, c = b - 1, b + 1
        if not (np.isfinite(t[a]) and np.isfinite(t[b]) and np.isfinite(t[c])):
            continue
        if not (int(q[a]) == 0 and int(q[b]) == 0 and int(q[c]) == 0):
            continue
        pre = float(t[b] - t[a])
        post = float(t[c] - t[b])
        if not (GAP_LOW * expected_days <= pre <= GAP_HIGH * expected_days):
            continue
        if not (GAP_LOW * expected_days <= post <= GAP_HIGH * expected_days):
            continue
        out.append({
            "a_row": a,
            "b_row": b,
            "c_row": c,
            "a_time": float(t[a]),
            "b_time": float(t[b]),
            "c_time": float(t[c]),
            "delta_pre_s": pre * 86400.0,
            "delta_post_s": post * 86400.0,
        })
    return out


def choose_prepointed_triple(time_values: np.ndarray, quality: np.ndarray, cadenceno: np.ndarray, expected_cadence_s: float) -> dict | None:
    triples = eligible_triples(time_values, quality, expected_cadence_s)
    finite = np.asarray(time_values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not triples or finite.size == 0:
        return None
    mid = float(np.median(finite))
    best = min(triples, key=lambda r: (abs(float(r["b_time"]) - mid), int(r["b_row"])))
    best = dict(best)
    cv = np.asarray(cadenceno, dtype=np.int64)
    best.update({
        "a_cadenceno": int(cv[best["a_row"]]),
        "b_cadenceno": int(cv[best["b_row"]]),
        "c_cadenceno": int(cv[best["c_row"]]),
        "sector_time_median": mid,
        "eligible_triple_count": len(triples),
    })
    return best


def _read_rows(path: Path) -> list[dict]:
    return list(csv.DictReader(path.open("r", encoding="utf-8", newline="")))


def _write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector-manifest", default="data/tachyon_star/JANUS-TACHYON-STAR-T3A-TESS-SECTOR-MANIFEST.csv")
    ap.add_argument("--output-dir", default="results/tachyon_star_t3a_tess")
    ap.add_argument("--cache-dir", default=".cache/tachyon_star_t3a_tess")
    args = ap.parse_args()

    sector_manifest = Path(args.sector_manifest)
    out = Path(args.output_dir)
    cache = Path(args.cache_dir)
    out.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)

    if sha256_file(sector_manifest) != SECTOR_MANIFEST_SHA256:
        raise RuntimeError(f"frozen TESS sector manifest SHA mismatch: {sha256_file(sector_manifest)}")
    rows = _read_rows(sector_manifest)
    if len(rows) != EXPECTED_SOURCES or len({r["src_id"] for r in rows}) != EXPECTED_SOURCES:
        raise RuntimeError("TESS sector manifest cardinality changed")

    selected = []
    blocked = []
    raw_receipts = []
    for i, row in enumerate(rows):
        sid = row["src_id"]
        safe = f"{i:02d}_{hashlib.sha256(sid.encode()).hexdigest()[:12]}"
        url = tesscut_url(row)
        try:
            payload = _fetch(url)
            zip_path = cache / f"{safe}.zip"
            zip_path.write_bytes(payload)
            fits_path, source_name = _extract_single_fits(payload, cache, safe)
            meta = read_time_quality_only(fits_path)
            forbidden_present = sorted(set(FORBIDDEN_VALUE_COLUMNS) & {str(c).upper() for c in meta["columns"]})
            triple = choose_prepointed_triple(
                meta["time"], meta["quality"], meta["cadenceno"], float(row["ffi_cadence_s"])
            )
            raw_receipts.append({
                "src_id": sid,
                "sector": int(row["sector"]),
                "tesscut_url": url,
                "zip_sha256": sha256_file(zip_path),
                "fits_sha256": sha256_file(fits_path),
                "source_fits_name": source_name,
                "table_hdu": meta["hdu"],
                "table_columns": meta["columns"],
                "forbidden_science_columns_present_but_not_dereferenced": forbidden_present,
                "row_count": int(len(meta["time"])),
                "flux_values_inspected": False,
            })
            if triple is None:
                blocked.append({"src_id": sid, "sector": int(row["sector"]), "reason": "NO_METADATA_ONLY_QUALITY_ABC_TRIPLE"})
            else:
                selected.append({
                    **row,
                    "tesscut_zip_sha256": sha256_file(zip_path),
                    "tesscut_fits_sha256": sha256_file(fits_path),
                    "source_fits_name": source_name,
                    "table_hdu": meta["hdu"],
                    **triple,
                    "flux_values_inspected": False,
                })
        except Exception as exc:
            blocked.append({
                "src_id": sid,
                "sector": int(row["sector"]),
                "reason": "TESSCUT_OR_METADATA_REPLAY_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
            })
        time.sleep(0.25)

    fields = [
        "src_id", "ra_deg", "dec_deg", "sector", "sectorName", "camera", "ccd", "ffi_cadence_s", "selection_rule",
        "tesscut_zip_sha256", "tesscut_fits_sha256", "source_fits_name", "table_hdu",
        "a_row", "b_row", "c_row", "a_cadenceno", "b_cadenceno", "c_cadenceno",
        "a_time", "b_time", "c_time", "delta_pre_s", "delta_post_s", "sector_time_median", "eligible_triple_count",
        "flux_values_inspected",
    ]
    selected.sort(key=lambda r: r["src_id"])
    _write_csv(out / "tess_prepointed_cadence_manifest.csv", selected, fields)
    (out / "raw_metadata_receipts.json").write_text(json.dumps(raw_receipts, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (out / "blocked.json").write_text(json.dumps(blocked, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    status = "PASS_PREPOINTED_CADENCE_MANIFEST_COMPLETE" if len(selected) == EXPECTED_SOURCES and not blocked else "BLOCKED_PREPOINTED_CADENCE_MANIFEST_INCOMPLETE"
    rec = {
        "schema": "janus.cosmos.tachyon_star.t3a.tess_cadence_prepoint.receipt.v1",
        "experiment_id": "JANUS-TACHYON-STAR-T3A-TESS-CADENCE-NATIVE-PREPOINTED-DISCOVERY",
        "status": status,
        "sector_manifest_sha256": SECTOR_MANIFEST_SHA256,
        "sources_expected": EXPECTED_SOURCES,
        "sources_prepointed": len(selected),
        "sources_blocked": len(blocked),
        "cadence_200s_prepointed": sum(int(r["ffi_cadence_s"]) == 200 for r in selected),
        "cadence_600s_prepointed": sum(int(r["ffi_cadence_s"]) == 600 for r in selected),
        "target_pixel_bytes_downloaded": True,
        "flux_values_inspected": False,
        "selection_uses_only": list(ALLOWED_VALUE_COLUMNS),
        "forbidden_science_value_columns": list(FORBIDDEN_VALUE_COLUMNS),
        "prepointed_manifest_sha256": sha256_file(out / "tess_prepointed_cadence_manifest.csv"),
        "raw_metadata_receipts_sha256": sha256_file(out / "raw_metadata_receipts.json"),
        "blocked": blocked,
        "claim_ceiling": CLAIM,
    }
    (out / "receipt.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(rec, indent=2, sort_keys=True))
    return 0 if status == "PASS_PREPOINTED_CADENCE_MANIFEST_COMPLETE" else 3


if __name__ == "__main__":
    raise SystemExit(main())
