#!/usr/bin/env python3
"""Live, deterministic MAST discovery for the expanded Janus Cosmos corpus.

Discovery only: no OCR, face, semantic, cipher, or candidate scoring.
Selects one reproducible HST IMAGE FITS product per target/filter when available.
Transient MAST failures are recorded per target and retried instead of aborting
an entire corpus discovery run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from pathlib import Path

from astroquery.mast import Observations

FILTERS = ("F435W", "F555W", "F814W")
DEFAULT_TARGETS = [
    "NGC1365", "NGC1425", "NGC1637", "NGC2841", "NGC3031", "NGC3627", "NGC4321",
    "M51", "M81", "M82", "NGC253", "NGC4258", "NGC4565", "NGC5194", "NGC5195",
    "NGC2403", "NGC6946", "NGC5457", "NGC3351", "NGC3621", "NGC6744", "NGC7793",
]
DISCOVERY_RETRIES = 3
DISCOVERY_BACKOFF_SECONDS = 3.0
MAST_TIMEOUT_SECONDS = 120

try:
    Observations.TIMEOUT = MAST_TIMEOUT_SECONDS
except Exception:
    pass


def safe(v):
    return "" if v is None else str(v)


def mast_query_name(target: str) -> str:
    """Use catalogue-friendly spacing for NGC designations."""
    if re.fullmatch(r"NGC\d+", target, flags=re.I):
        return re.sub(r"^NGC", "NGC ", target, flags=re.I)
    return target


def product_score(row):
    filename = safe(row["productFilename"]).lower()
    desc = safe(row["description"]).lower()
    group = safe(row["productGroupDescription"]).lower()
    ptype = safe(row["productType"]).lower()
    sci = 1 if ptype == "science" else 0
    minrec = 1 if "minimum recommended" in group else 0
    mosaic = 1 if any(k in filename for k in ("mosaic", "_drz", "_drc", "_sci")) else 0
    fits = 1 if filename.endswith((".fits", ".fits.gz")) else 0
    try:
        size = int(row["size"])
    except Exception:
        size = 0
    return (sci, minrec, mosaic, fits, size, filename, desc)


def discover(target: str):
    query_name = mast_query_name(target)
    last_error = None
    for attempt in range(1, DISCOVERY_RETRIES + 1):
        try:
            obs = Observations.query_object(query_name, radius="0.05 deg", obs_collection="HST", dataproduct_type="IMAGE")
            if len(obs) == 0:
                return [], {"target": target, "query_name": query_name, "status": "NO_OBSERVATIONS", "attempts": attempt}
            products = Observations.get_unique_product_list(obs)
            rows = []
            for filt in FILTERS:
                candidates = []
                for row in products:
                    filters = safe(row["filters"])
                    filename = safe(row["productFilename"])
                    if filt not in filters.split(";") and filt not in filters.split(",") and filt not in filters:
                        continue
                    if not re.search(r"\.fits(?:\.gz)?$", filename, re.I):
                        continue
                    if safe(row["dataproduct_type"]).lower() not in ("image", ""):
                        continue
                    candidates.append(row)
                if not candidates:
                    rows.append({"target": target, "filter": filt, "discovery_status": "NOT_FOUND"})
                    continue
                row = max(candidates, key=product_score)
                rows.append({
                    "target": target,
                    "query_name": query_name,
                    "filter": filt,
                    "band": {"F435W": "B", "F555W": "V", "F814W": "I"}[filt],
                    "discovery_status": "FOUND",
                    "dataURI": safe(row["dataURI"]),
                    "obs_collection": safe(row["obs_collection"]),
                    "obs_id": safe(row["obs_id"]),
                    "obsid": safe(row["obsid"]),
                    "productFilename": safe(row["productFilename"]),
                    "productType": safe(row["productType"]),
                    "productGroupDescription": safe(row["productGroupDescription"]),
                    "size": int(row["size"]) if safe(row["size"]).isdigit() else None,
                })
            return rows, {"target": target, "query_name": query_name, "status": "OK", "attempts": attempt}
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            if attempt < DISCOVERY_RETRIES:
                time.sleep(DISCOVERY_BACKOFF_SECONDS * attempt)
    return [], {
        "target": target,
        "query_name": query_name,
        "status": "MAST_QUERY_ERROR",
        "attempts": DISCOVERY_RETRIES,
        "error": last_error,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="*", default=DEFAULT_TARGETS)
    ap.add_argument("--out", default="data/hst_expanded_live_manifest.json")
    args = ap.parse_args()

    products = []
    target_status = []
    for target in sorted(dict.fromkeys(args.targets)):
        rows, status = discover(target)
        products.extend(rows)
        target_status.append(status)
        print(json.dumps(status, sort_keys=True), flush=True)

    targets = []
    for target in sorted(dict.fromkeys(args.targets)):
        found = [p for p in products if p["target"] == target and p["discovery_status"] == "FOUND"]
        if len(found) >= 2:
            targets.append({
                "target": target,
                "class": "galaxy",
                "filters": sorted([
                    {"filter": p["filter"], "band": p["band"], "dataURI": p["dataURI"], "obs_id": p["obs_id"], "productFilename": p["productFilename"]}
                    for p in found
                ], key=lambda x: x["filter"]),
            })

    query_errors = [x for x in target_status if x["status"] == "MAST_QUERY_ERROR"]
    receipt = {
        "schema": "janus.cosmos.hst.live_mast_manifest.v0.4",
        "status": "LIVE_MAST_DISCOVERY",
        "source": "MAST / STScI",
        "selection": "Deterministic per-target/filter HST IMAGE FITS selection from live MAST products; targets with fewer than two requested filters are excluded from scoring.",
        "requested_targets": len(args.targets),
        "selected_targets": len(targets),
        "selected_products": sum(len(t["filters"]) for t in targets),
        "filters": list(FILTERS),
        "discovery_retries": DISCOVERY_RETRIES,
        "mast_timeout_seconds": MAST_TIMEOUT_SECONDS,
        "query_error_count": len(query_errors),
        "target_status": target_status,
        "blind_restrictions": {"ocr": False, "face_search": False, "semantic_analysis": False, "cipher_search": False, "post_hoc_tuning": False},
        "targets": targets,
        "discovery_log": products,
    }
    raw = json.dumps(receipt, sort_keys=True, separators=(",", ":")).encode()
    receipt["manifest_sha256"] = hashlib.sha256(raw).hexdigest()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"selected_targets": len(targets), "selected_products": receipt["selected_products"], "query_error_count": len(query_errors), "manifest_sha256": receipt["manifest_sha256"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
