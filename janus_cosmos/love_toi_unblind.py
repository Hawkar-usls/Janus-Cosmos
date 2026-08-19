from __future__ import annotations

import argparse
import csv
import io
import json
import math
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

TAP_SYNC = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"


def _query_toi(tic_ids: list[str], retries: int = 3) -> tuple[list[dict], str | None]:
    ids = ",".join(str(int(x)) for x in tic_ids)
    query = (
        "select toi,tid,tfopwg_disp,pl_orbper,pl_trandurh,pl_trandep,pl_rade,"
        "st_tmag,st_dist,st_teff,rowupdate,ra,dec "
        f"from toi where tid in ({ids})"
    )
    url = TAP_SYNC + "?" + urllib.parse.urlencode({"query": query, "format": "csv"})
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Janus-Cosmos-LOVE-TOI-UNBLIND/1.0"})
            with urllib.request.urlopen(req, timeout=120) as response:
                text = response.read().decode("utf-8")
            rows = []
            for r in csv.DictReader(io.StringIO(text)):
                rows.append({
                    "toi": r.get("toi") or None,
                    "tic_id": str(int(float(r["tid"]))) if r.get("tid") else None,
                    "disposition": r.get("tfopwg_disp") or None,
                    "period_days": float(r["pl_orbper"]) if r.get("pl_orbper") else None,
                    "duration_hours": float(r["pl_trandurh"]) if r.get("pl_trandurh") else None,
                    "depth_ppm": float(r["pl_trandep"]) if r.get("pl_trandep") else None,
                    "radius_earth": float(r["pl_rade"]) if r.get("pl_rade") else None,
                    "tmag": float(r["st_tmag"]) if r.get("st_tmag") else None,
                    "distance_pc": float(r["st_dist"]) if r.get("st_dist") else None,
                    "teff_k": float(r["st_teff"]) if r.get("st_teff") else None,
                    "rowupdate": r.get("rowupdate") or None,
                    "ra_deg": float(r["ra"]) if r.get("ra") else None,
                    "dec_deg": float(r["dec"]) if r.get("dec") else None,
                })
            return rows, None
        except Exception as exc:
            last_error = f"attempt {attempt}/{retries}: {type(exc).__name__}: {exc}"
            if attempt < retries:
                time.sleep(5 * attempt)
    return [], last_error


def _relative_period_delta(a: float | None, b: float | None) -> float | None:
    if a is None or b is None or a <= 0 or b <= 0:
        return None
    return abs(a - b) / ((a + b) / 2.0)


def run(raw_receipt_path: Path, prereg_path: Path, output_path: Path) -> dict:
    raw = json.loads(raw_receipt_path.read_text(encoding="utf-8"))
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    frozen = prereg["frozen_sources"]
    labels = {str(x["tic_id"]): f"TARGET_A_R{x['rank']}" for x in frozen}
    toi_rows, error = _query_toi(list(labels))
    if error is not None:
        result = {
            "schema": "janus.cosmos.love.toi_unblind.result.v1",
            "experiment_id": "LOVE-G1-007-TOI-UNBLIND-v1",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "status": "BLOCKED_CATALOG_QUERY_ERROR",
            "error": error,
            "love_gate": {"candidate_activated": False, "reserved_codename": "Love"},
            "claim_ceiling": "CATALOG_UNBLIND_FAILED__NO_NEW_OBJECT_CLAIM__LOVE_SEALED",
        }
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return result

    by_tic: dict[str, list[dict]] = {}
    for row in toi_rows:
        row["blind_label"] = labels.get(str(row["tic_id"]))
        by_tic.setdefault(str(row["tic_id"]), []).append(row)

    hit = raw["raw_hit"]
    hit_tic = str(hit["tic_id"])
    blind_period = float(hit["combined_bls"]["best_period_days"])
    matches = []
    for row in by_tic.get(hit_tic, []):
        delta = _relative_period_delta(blind_period, row.get("period_days"))
        enriched = dict(row)
        enriched["blind_bls_period_days"] = blind_period
        enriched["relative_period_delta"] = delta
        enriched["period_matches_within_1pct"] = bool(delta is not None and delta <= 0.01)
        matches.append(enriched)

    same_period = [x for x in matches if x["period_matches_within_1pct"]]
    rediscovered = bool(same_period)
    classification = (
        "KNOWN_TESS_PROJECT_SIGNAL_REDISCOVERED"
        if rediscovered
        else "NO_MATCHING_TOI_PERIOD__REMAINS_UNCONFIRMED_AND_REQUIRES_MORE_CATALOG_VETTING"
    )

    result = {
        "schema": "janus.cosmos.love.toi_unblind.result.v1",
        "experiment_id": "LOVE-G1-007-TOI-UNBLIND-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "source_raw_receipt": str(raw_receipt_path),
        "authoritative_catalog": "NASA Exoplanet Archive TESS Project Candidates / TOI table",
        "catalog_rows_for_frozen_tics": toi_rows,
        "blind_hit": {
            "blind_label": hit["blind_label"],
            "tic_id": hit_tic,
            "gaia_source_id": hit["gaia_source_id"],
            "blind_bls_period_days": blind_period,
            "catalog_matches_same_tic": matches,
            "matching_period_rows": same_period,
            "classification": classification,
        },
        "science_interpretation": {
            "new_planet_candidate_from_this_run": False if rediscovered else None,
            "positive_control_value": (
                "The blind JANUS BLS gate independently recovered a period already present in the authoritative TOI catalog; "
                "this validates sensitivity to a real cataloged TESS transit-like signal but is not a new discovery."
                if rediscovered else
                "No same-period TOI row was found; this does not establish novelty and requires broader catalog and false-positive checks."
            ),
        },
        "love_gate": {
            "reserved_codename": "Love",
            "candidate_activated": False,
            "reason": (
                "The blind hit is a rediscovery of an already cataloged TESS Object of Interest, so it is excluded by the novelty firewall."
                if rediscovered else
                "Novelty is not yet established; independent replication and false-positive rejection are still required."
            ),
            "next_gate": (
                "RETURN_TO_BLIND_SEARCH__DO_NOT_TUNE_ON_REDISCOVERED_TOI"
                if rediscovered else
                "BROADER_CATALOG_AND_FALSE_POSITIVE_VETTING"
            ),
        },
        "claim_ceiling": "TOI_UNBLIND_CLASSIFICATION_ONLY__NO_NEW_EXOPLANET_DISCOVERY__LOVE_SEALED",
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description="Unblind a blind TESS BLS hit against the authoritative NASA TOI table.")
    ap.add_argument("--raw-receipt", default="data/love/LOVE-G1-007-TESS-BLS-v1-RUN-001-RAW-RECEIPT.json")
    ap.add_argument("--prereg", default="data/love/LAVE_TO_LOVE_TESS_BLS_PREREG.json")
    ap.add_argument("--output", default="results/love_toi_unblind/love-toi-unblind-result.json")
    args = ap.parse_args()
    result = run(Path(args.raw_receipt), Path(args.prereg), Path(args.output))
    print(json.dumps({
        "status": result["status"],
        "classification": result.get("blind_hit", {}).get("classification"),
        "love_candidate_activated": result["love_gate"]["candidate_activated"],
        "next_gate": result["love_gate"].get("next_gate"),
    }, indent=2))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
