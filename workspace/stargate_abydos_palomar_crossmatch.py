from __future__ import annotations

import gzip
import hashlib
import io
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

TARGET_RA = 223.415064157
TARGET_DEC = 33.979315670
RADII_ARCSEC = [2, 5, 10, 30, 60, 300]
SOURCE_URL = "https://raw.githubusercontent.com/jannefi/poss1-plate-slice/4005e200541b321ead3d6608f0162a14430ef1c2/results/s0-642-20260814/stage_S0.csv.gz"
EXPECTED_ROWS = 122820
EXPECTED_SHA = "2ff92f2210acb387ef9ef4b88d561595d3883e9aab27065042627272b96590f0"
OUT = Path("data/stargate/STARGATE-ABYDOS-PALOMAR-SPATIAL-CROSSMATCH-v1-LATEST-RECEIPT.json")

EVENTS = {
    "PALOMAR_NINE_1950_FIELD_CENTER": {"ra_deg": 212.9291, "dec_deg": 26.8311, "date": "1950-04-12"},
    "PALOMAR_TRIPLE_1952_CENTER": {"ra_deg": 319.5433333333333, "dec_deg": 50.37872222222222, "date": "1952-07-19"},
}


def angular_sep_deg(ra1, dec1, ra2, dec2):
    r1 = np.deg2rad(ra1); d1 = np.deg2rad(dec1)
    r2 = np.deg2rad(ra2); d2 = np.deg2rad(dec2)
    x = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.rad2deg(np.arccos(np.clip(x, -1.0, 1.0)))


def scalar_sep_deg(a_ra, a_dec, b_ra, b_dec):
    return float(angular_sep_deg(a_ra, a_dec, b_ra, b_dec))


def main():
    r = requests.get(SOURCE_URL, timeout=90)
    r.raise_for_status()
    raw = gzip.decompress(r.content)
    digest = hashlib.sha256(raw).hexdigest()
    df = pd.read_csv(io.BytesIO(raw))
    if digest != EXPECTED_SHA:
        raise RuntimeError(f"uncompressed sha mismatch: {digest}")
    if len(df) != EXPECTED_ROWS:
        raise RuntimeError(f"row count mismatch: {len(df)}")
    if not {"ra", "dec"}.issubset(df.columns):
        raise RuntimeError(f"missing ra/dec columns: {list(df.columns)}")

    seps_deg = angular_sep_deg(TARGET_RA, TARGET_DEC, df["ra"].to_numpy(float), df["dec"].to_numpy(float))
    seps_arcsec = seps_deg * 3600.0
    order = np.argsort(seps_arcsec)
    nearest_rows = []
    for idx in order[:20]:
        row = df.iloc[int(idx)]
        nearest_rows.append({
            "separation_arcsec": float(seps_arcsec[int(idx)]),
            "ra_deg": float(row["ra"]),
            "dec_deg": float(row["dec"]),
            "src_id": None if pd.isna(row.get("src_id")) else str(row.get("src_id")),
            "tile_id": None if pd.isna(row.get("tile_id")) else str(row.get("tile_id")),
            "object_id": None if pd.isna(row.get("object_id")) else str(row.get("object_id")),
        })

    counts = {str(rad): int(np.sum(seps_arcsec <= rad)) for rad in RADII_ARCSEC}
    nearest = nearest_rows[0]
    event_controls = {}
    for name, e in EVENTS.items():
        s = scalar_sep_deg(TARGET_RA, TARGET_DEC, e["ra_deg"], e["dec_deg"])
        event_controls[name] = {**e, "separation_deg": s, "separation_arcmin": s * 60.0}

    payload = {
        "schema": "janus.cosmos.stargate_abydos.palomar_spatial_crossmatch.receipt.v1",
        "experiment_id": "STARGATE-ABYDOS-PALOMAR-SPATIAL-CROSSMATCH-v1",
        "status": "COMPLETE",
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_target": {"ra_deg": TARGET_RA, "dec_deg": TARGET_DEC, "frame": "ICRS"},
        "source": {
            "url": SOURCE_URL,
            "rows": int(len(df)),
            "uncompressed_sha256": digest,
            "catalog_semantics": "independent public POSS-I candidate reconstruction; not identical to request-only 107875-row VASCO catalog"
        },
        "radius_counts_arcsec": counts,
        "nearest_open_cohort_candidate": nearest,
        "nearest_20_open_cohort_candidates": nearest_rows,
        "published_event_controls": event_controls,
        "summary": {
            "exact_2arcsec_match": counts["2"] > 0,
            "within_5arcsec": counts["5"] > 0,
            "within_30arcsec": counts["30"] > 0,
            "nearest_separation_arcsec": nearest["separation_arcsec"],
            "same_field_as_1950_nine": event_controls["PALOMAR_NINE_1950_FIELD_CENTER"]["separation_deg"] <= (10.0/60.0),
            "same_field_as_1952_triple": event_controls["PALOMAR_TRIPLE_1952_CENTER"]["separation_deg"] <= (3.0/60.0)
        },
        "firewall": {
            "nearest_candidate_is_verified_physical_transient": False,
            "nuclear_temporal_association_is_spatial_identity": False,
            "nuclear_association_proves_light_echo": False,
            "solar_reflection_hypothesis_equals_nuclear_flash_reflection": False,
            "claim_ceiling": "SPATIAL_CROSSMATCH_ONLY"
        }
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()
