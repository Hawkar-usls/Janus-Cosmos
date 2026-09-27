#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import build_observer_index as base

SITE_REGISTRY = "SITE_COORDINATE_REGISTRY_v0.2.json"
DAY_LEDGER = "PALOMAR_OBSERVATION_DAY_LEDGER_v0.1.json"


def walk_bridge_rows(x, source_path: str, aliases, day_by_date, out):
    if isinstance(x, dict):
        if x.get("bridge_eligible") is True and isinstance(x.get("date"), str) and isinstance(x.get("site"), str):
            site = x["site"].strip()
            reg = aliases.get(site.casefold())
            probe = {"date": x.get("date")}
            state = base.match_state(probe, day_by_date)
            row = {
                "row_id": x.get("row_id"),
                "site": site,
                "date": x.get("date"),
                "source_path": source_path,
                "source_grade": x.get("source_grade"),
                "classification_state": x.get("classification_state"),
                "palomar_match_state": state,
            }
            if reg:
                row["lat"] = reg.get("lat")
                row["lng"] = reg.get("lon")
            if state == "EXACT_DAY_JOIN_AVAILABLE":
                d = day_by_date.get(x.get("date")) or {}
                row["palomar_day"] = {
                    "candidate_count": d.get("candidate_count"),
                    "plate_count": d.get("plate_count"),
                    "tile_count": d.get("tile_count"),
                    "artifact_state": d.get("artifact_state"),
                    "exposure": d.get("exposure"),
                }
            out.append(row)
        for v in x.values():
            walk_bridge_rows(v, source_path, aliases, day_by_date, out)
    elif isinstance(x, list):
        for v in x:
            walk_bridge_rows(v, source_path, aliases, day_by_date, out)


def primary_bridge_cohort(earth_root: Path, aliases, day_by_date):
    rows = []
    for path in sorted(earth_root.rglob("*.json")):
        if base.skip_path(path) or path.name == SITE_REGISTRY:
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        rel = str(path.relative_to(earth_root.parent)).replace("\\", "/")
        walk_bridge_rows(data, rel, aliases, day_by_date, rows)
    unique = []
    seen = set()
    for r in rows:
        key = (r.get("row_id"), r.get("site"), r.get("date"), r.get("source_path"))
        if key in seen:
            continue
        seen.add(key)
        unique.append(r)
    counts = {}
    for r in unique:
        s = r["palomar_match_state"]
        counts[s] = counts.get(s, 0) + 1
    return unique, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--meta-root", type=Path)
    ap.add_argument("--plate-manifest", type=Path)
    ap.add_argument("--day-ledger", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    repo = args.repo_root.resolve()
    aliases, bootstrap_sites = base.parse_registry(repo / "domains/earth" / SITE_REGISTRY)
    day = base.load_day_ledger(args.day_ledger)

    earth, sky, unmapped, scan_receipts = [], [], [], []
    suppressed_earth, suppressed_sky = [], []
    roots = [(repo / "domains", "JANUS_COSMOS")]
    if args.meta_root and args.meta_root.exists():
        roots.append((args.meta_root.resolve(), "JANUS_META"))

    for root, family in roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob("*.json")):
            if base.skip_path(path):
                continue
            e, s, u, receipt = base.scan_json_file(path, root, family, aliases)
            rel = str(path.relative_to(root)).replace("\\", "/")

            suppress_all_render = False
            suppress_reason = None
            if family == "JANUS_COSMOS" and rel.startswith("bridge/"):
                suppress_all_render = True
                suppress_reason = "BRIDGE_METADATA_DESCRIBES_RELATIONS__NOT_NEW_EARTH_OR_SKY_ENTITIES"
            elif family == "JANUS_COSMOS" and path.name == SITE_REGISTRY:
                suppress_all_render = True
                suppress_reason = "SITE_REGISTRY_IS_COORDINATE_AUTHORITY__NOT_AN_EVENT_LAYER"

            if suppress_all_render:
                if e:
                    suppressed_earth.append({"source_family": family, "path": rel, "suppressed_coordinate_rows": len(e), "reason": suppress_reason})
                if s:
                    suppressed_sky.append({"source_family": family, "path": rel, "suppressed_coordinate_rows": len(s), "reason": suppress_reason})
            else:
                earth.extend(e)
                unmapped.extend(u)
                if family == "JANUS_COSMOS" and path.name == DAY_LEDGER:
                    if s:
                        suppressed_sky.append({
                            "source_family": family,
                            "path": rel,
                            "suppressed_coordinate_rows": len(s),
                            "reason": "OBSERVATION_DAY_METADATA_DUPLICATES_CANONICAL_POSS1_PLATE_CENTER_LAYER",
                        })
                else:
                    sky.extend(s)

            if receipt.get("state") != "OK":
                scan_receipts.append({"source_family": family, "path": rel, "receipt": receipt})

    earth = base.dedupe(earth, "earth")
    sky = base.dedupe(sky, "sky")
    plates = base.parse_plates(args.plate_manifest)
    enriched_plate_count = base.enrich_plates_from_days(plates, day)
    for p in earth:
        p["palomar_match_state"] = base.match_state(p, day)

    presentation_states = {}
    for p in earth:
        s = p["palomar_match_state"]
        presentation_states[s] = presentation_states.get(s, 0) + 1

    cohort, cohort_states = primary_bridge_cohort(repo / "domains/earth", aliases, day)

    out = {
        "schema": "janus.cosmos.observer.autosync_index.v0.5",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "authority_boundary": "GENERATED_PRESENTATION_INDEX__SOURCE_JSON_REMAINS_CANONICAL",
        "autosync": {
            "janus_cosmos_scan": True,
            "janus_meta_scan": bool(args.meta_root and args.meta_root.exists()),
            "protected_exclusions": ["cousteau", "kusto"],
            "plate_manifest_rows": len(plates),
            "palomar_day_rows": len(day),
        },
        "representation_quotient": {
            "laws": [
                "METADATA_ABOUT_AN_ENTITY != NEW_ENTITY",
                "DIFFERENT_REPRESENTATIONS_OF_SAME_EVENT -> ONE_RENDERED_ENTITY",
                "QUOTIENT_BEFORE_RENDER__QUOTIENT_BEFORE_COUNT",
            ],
            "suppressed_generic_earth_layers": suppressed_earth,
            "suppressed_generic_sky_layers": suppressed_sky,
            "suppressed_earth_coordinate_rows": sum(x["suppressed_coordinate_rows"] for x in suppressed_earth),
            "suppressed_sky_coordinate_rows": sum(x["suppressed_coordinate_rows"] for x in suppressed_sky),
            "canonical_palomar_render_layer": "POSS1_PUBLIC_MANIFEST plate centers",
            "palomar_day_ledger_role": "TEMPORAL_AND_DETECTOR_METADATA_ENRICHMENT__NOT_SECOND_SKY_POINT_LAYER",
            "bridge_domain_role": "RELATION_AND_TEST_METADATA__NOT_SECOND_EARTH_EVENT_LAYER",
            "site_registry_role": "COORDINATE_ALIAS_AUTHORITY__NOT_EVENT_LAYER",
            "plate_centers_enriched_with_observation_day": enriched_plate_count,
        },
        "counts": {
            "earth_points": len(earth),
            "earth_unmapped_source_rows": len(unmapped),
            "janus_sky_targets": len(sky),
            "poss1_plate_centers": len(plates),
            "all_space_points": len(sky) + len(plates),
            "primary_bridge_test_rows": len(cohort),
        },
        "earth_points": earth,
        "earth_unmapped": unmapped[:2000],
        "janus_sky_targets": sky,
        "poss1_plate_centers": plates,
        "bridge": {
            "palomar_window": ["1949-11-19", "1957-04-28"],
            "palomar_unique_observation_days_expected": 312,
            "palomar_day_rows_materialized": len(day),
            "presentation_match_state_counts": presentation_states,
            "primary_test_rows": cohort,
            "primary_test_state_counts": cohort_states,
            "exact_day_join_gate": "OPEN" if day and len(day) >= 312 else "BLOCKED_312_DAY_TABLE",
            "laws": [
                "NO_OBSERVATION != ZERO",
                "TEMPORAL_NONOVERLAP != PHYSICAL_NONCORRELATION",
                "YEAR_PRECISION != DAY_PRECISION",
                "REPORT_COUNT != PHYSICAL_RATE",
                "PRESENTATION_ROWS != INFERENTIAL_COHORT",
            ],
        },
        "bootstrap_site_registry_count": len(bootstrap_sites),
        "scan_receipts": scan_receipts,
        "canonical_seal": "THE OBSERVER MAY DISCOVER MANY REPRESENTATIONS, BUT IT RENDERS AND COUNTS DECISION-RELEVANT ENTITIES IN THEIR CONTRACTED QUOTIENT. BRIDGE METADATA DOES NOT BECOME A SECOND EARTH EVENT."
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("OBSERVER_INDEX", json.dumps(out["counts"], sort_keys=True))
    print("PRESENTATION_BRIDGE_STATES", json.dumps(presentation_states, sort_keys=True))
    print("PRIMARY_BRIDGE_STATES", json.dumps(cohort_states, sort_keys=True))
    print("REPRESENTATION_QUOTIENT", json.dumps(out["representation_quotient"], sort_keys=True))


if __name__ == "__main__":
    main()
