#!/usr/bin/env python3
"""Reproduce a narrow MAST archival inventory, not an IFU nonexistence proof."""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path
from urllib.parse import urlencode


RA = 53.13127
DEC = -27.78908
RADIUS_DEG = 0.0041667  # ~15 arcsec, not the source physical extent
REQUEST = {
    "service": "Mast.Caom.Cone",
    "params": {"ra": RA, "dec": DEC, "radius": RADIUS_DEG},
    "format": "json", "pagesize": 5000,
}
QUERY_URL = "https://mast.stsci.edu/api/v0/invoke?" + urlencode(
    {"request": json.dumps(REQUEST, separators=(",", ":"))}
)


def audit(raw):
    data = json.loads(raw)
    assert data["status"] == "COMPLETE"
    rows = data["data"]
    assert data["paging"]["rowsFiltered"] == len(rows)
    assert data["paging"]["pagesFiltered"] == 1
    names = Counter(str(x.get("instrument_name")) for x in rows)
    spectral = [x for x in rows if "NIRSPEC" in str(x.get("instrument_name", "")).upper()]
    ifu = [x for x in rows if "IFU" in str(x.get("instrument_name", "")).upper()]
    jwst = [x for x in rows if x.get("obs_collection") == "JWST"]
    return {
        "schema": "janus-mocss-mast-15arcsec-observation-inventory-v1",
        "archive_query": "Mast.Caom.Cone",
        "exact_request": REQUEST,
        "replay_query_url": QUERY_URL,
        "snapshot_response_sha256": hashlib.sha256(raw).hexdigest(),
        "status": "MAST_COMPLETE_CATALOG_INVENTORY_NOT_GEOMETRIC_IFU_SKY_COVERAGE_PROOF",
        "source_ra_dec_deg": [RA, DEC],
        "cone_radius_arcsec": RADIUS_DEG * 3600,
        "total_catalog_records": len(rows),
        "jwst_collection_records": len(jwst),
        "nirspec_labelled_records": len(spectral),
        "ifu_labelled_records": len(ifu),
        "labelled_instrument_counts": dict(sorted(names.items())),
        "ifu_inventory_mode": "substring IFU in instrument_name; not a test of actual aperture footprint or hidden IFU in generic mode labels",
        "line_16004_recovery": None,
        "additional_IFU_physical_spectra": None,
        "directional_gas_trail": None,
        "scientific_gate": "NO_IFU_LABELLED_PRODUCT_IN_THIS_INDEXED_CONE__NOT_NO_IFU_DATA_ANYWHERE__NOT_AN_ASTROPHYSICAL_FLUX_LIMIT",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", required=True, type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    result = audit(args.snapshot.read_bytes())
    s = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(s)
    else:
        print(s, end="")
    print("MAST_NIRSPEC_IFU_CONE_INVENTORY_REPRODUCED_NO_LINE_ID")


if __name__ == "__main__":
    main()
