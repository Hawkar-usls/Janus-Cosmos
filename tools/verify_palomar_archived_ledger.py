#!/usr/bin/env python3
"""Verify already archived Palomar 312-day ledger offline, without pretending to
refetch mutable plate metadata from IRSA or the upstream main branch.

This gate does NOT prove that remote source files are currently unchanged.
Explicit workflow_dispatch keeps the existing network materializer available.
"""
from __future__ import annotations

import argparse
from datetime import date
import json
from pathlib import Path

EXPECTED_REPOSITORY = "jannefi/poss1-plate-slice"
EXPECTED_COMMIT = "4005e200541b321ead3d6608f0162a14430ef1c2"
EXPECTED_STAGE = "2ff92f2210acb387ef9ef4b88d561595d3883e9aab27065042627272b96590f0"
EXPECTED_TILES = "5dcb90dc5d98550e5a60246aced2b097922a267c69e81f27d45d16a288142a99"
EXPECTED_COUNTS = {"rows": 312, "plates_in_window": 640, "candidate_rows_in_window": 122511, "processed_tiles_in_window": 31360}
SCHEMA = "janus.cosmos.space.palomar_observation_day_ledger.v0.1"
STATUS = "MATERIALIZED_SOURCE_BOUND_312_DAY_LEDGER"
WINDOW = ["1949-11-19", "1957-04-28"]


class ArchiveHold(ValueError):
    pass


def require(value, code):
    if not value:
        raise ArchiveHold(code)


def verify(doc):
    require(isinstance(doc, dict) and doc.get("schema") == SCHEMA
            and doc.get("status") == STATUS, "PALOMAR_ARCHIVE_SCHEMA")
    require(doc.get("study_window") == WINDOW, "PALOMAR_ARCHIVE_WINDOW")
    frozen = doc.get("source_freeze")
    require(isinstance(frozen, dict) and
            frozen.get("upstream_repository") == EXPECTED_REPOSITORY and
            frozen.get("upstream_commit") == EXPECTED_COMMIT and
            frozen.get("stage_sha256_uncompressed") == EXPECTED_STAGE and
            frozen.get("tile_manifest_sha256_uncompressed") == EXPECTED_TILES and
            frozen.get("irsa_headers_parsed") == 642, "PALOMAR_ARCHIVE_SOURCE_BINDING")
    # The upstream plate manifest was intentionally referenced by main in v1.1;
    # its historical digest is retained, but it is NOT treated as live-frozen.
    manifest = frozen.get("plate_manifest_sha256")
    require(isinstance(manifest, str) and len(manifest) == 64 and
            all(c in "0123456789abcdef" for c in manifest),
            "PALOMAR_ARCHIVE_MANIFEST_HASH")
    rows = doc.get("rows")
    counts = doc.get("counts")
    require(isinstance(rows, list) and isinstance(counts, dict) and
            all(type(counts.get(k)) is int and counts[k] == v
                for k, v in EXPECTED_COUNTS.items()) and
            len(rows) == EXPECTED_COUNTS["rows"],
            "PALOMAR_ARCHIVE_COUNTS")
    previous = ""
    seen = set()
    total_candidates = total_tiles = total_plates = 0
    for row in rows:
        require(isinstance(row, dict) and row.get("observed_flag") is True,
                "PALOMAR_ARCHIVE_OBSERVATION")
        ds = row.get("date")
        require(isinstance(ds, str) and WINDOW[0] <= ds <= WINDOW[1]
                and ds > previous, "PALOMAR_ARCHIVE_DATE_ORDER")
        try:
            require(date.fromisoformat(ds).isoformat() == ds, "PALOMAR_ARCHIVE_DATE")
        except ValueError as exc:
            raise ArchiveHold("PALOMAR_ARCHIVE_DATE") from exc
        previous = ds
        plates = row.get("plate_ids")
        require(isinstance(plates, list) and len(plates) > 0 and
                len(plates) == row.get("plate_count") and
                all(isinstance(p, str) and p and p not in seen for p in plates),
                "PALOMAR_ARCHIVE_PLATES")
        seen.update(plates)
        positions = (row.get("sky_position") or {}).get("plate_centers_ra_dec")
        require(isinstance(positions, list) and
                {p.get("plate_id") for p in positions if isinstance(p, dict)} == set(plates)
                and len(positions) == len(plates), "PALOMAR_ARCHIVE_SKY_BINDING")
        for position in positions:
            ra, dec = position.get("ra_deg"), position.get("dec_deg")
            require(type(ra) in (int, float) and type(dec) in (int, float)
                    and 0 <= ra < 360 and -90 <= dec <= 90,
                    "PALOMAR_ARCHIVE_SKY_COORDINATES")
        candidates, tiles = row.get("candidate_count"), row.get("tile_count")
        require(type(candidates) is int and candidates > 0 and
                type(tiles) is int and tiles > 0, "PALOMAR_ARCHIVE_ROW_COUNTS")
        total_candidates += candidates
        total_tiles += tiles
        total_plates += len(plates)
    require(total_candidates == EXPECTED_COUNTS["candidate_rows_in_window"]
            and total_tiles == EXPECTED_COUNTS["processed_tiles_in_window"]
            and total_plates == EXPECTED_COUNTS["plates_in_window"]
            and counts.get("all_observed_days_candidate_positive") is True,
            "PALOMAR_ARCHIVE_AGGREGATES")
    require(doc.get("canonical_seal", "").startswith("OBSERVATION OPPORTUNITY IS MATERIALIZED"),
            "PALOMAR_ARCHIVE_SEAL")
    return {"status": "ARCHIVE_VERIFIED_OFFLINE_NOT_REMOTE_REFRESH",
            "rows": len(rows), "plates": total_plates,
            "upstream_current_state_verified": False, "authority_delta": 0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ledger", type=Path, required=True)
    args = ap.parse_args()
    try:
        raw = args.ledger.read_bytes()
        require(0 < len(raw) <= 1_000_000, "PALOMAR_ARCHIVE_SIZE")
        doc = json.loads(raw, parse_constant=lambda _: (_ for _ in ()).throw(ArchiveHold("PALOMAR_NONFINITE_JSON")))
        report = verify(doc)
    except (ArchiveHold, OSError, ValueError, TypeError, KeyError) as exc:
        print("PALOMAR_ARCHIVE=HOLD_" + type(exc).__name__)
        return 2
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
