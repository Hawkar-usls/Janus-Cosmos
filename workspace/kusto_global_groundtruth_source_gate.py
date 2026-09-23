#!/usr/bin/env python3
from __future__ import annotations

import json
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "workspace" / "kusto_global_groundtruth_out"
OUT.mkdir(parents=True, exist_ok=True)

PREREG = ROOT / "data/cousteau/JANUS-KUSTO-GLOBAL-GROUNDTRUTH-BENCHMARK-PREREG-2026-09-23-v1.0.json"
REGISTRY = ROOT / "data/cousteau/JANUS-KUSTO-GLOBAL-GROUNDTRUTH-SOURCE-REGISTRY-2026-09-23-v1.0.json"
LAYER = "https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
UA = {"User-Agent": "JANUS-KUSTO-groundtruth-source-gate/1.0"}

prereg = json.loads(PREREG.read_text())
registry = json.loads(REGISTRY.read_text())

def get_json(url, params=None, tries=6):
    last = None
    for i in range(tries):
        try:
            r = requests.get(url, params=params, headers=UA, timeout=90)
            if r.status_code == 429:
                time.sleep(min(20, 2 ** i))
                continue
            r.raise_for_status()
            d = r.json()
            if isinstance(d, dict) and d.get("error"):
                raise RuntimeError(d["error"])
            return d
        except Exception as e:
            last = e
            if i + 1 == tries:
                raise
            time.sleep(min(20, 2 ** i))
    raise last

def source_probe(url):
    try:
        r = requests.get(url, headers=UA, timeout=45, allow_redirects=True, stream=True)
        return {
            "url": url,
            "status_code": r.status_code,
            "resolved_url": r.url,
            "reachable": bool(200 <= r.status_code < 400),
            "content_type": r.headers.get("content-type"),
        }
    except Exception as e:
        return {"url": url, "reachable": False, "error": type(e).__name__ + ": " + str(e)}

def survey_query(sid):
    return get_json(LAYER, {
        "where": f"SURVEY_ID='{sid}'",
        "outFields": "SURVEY_ID,PLATFORM,INSTRUMENT,SURVEY_YEAR,START_TIME,END_TIME,DOWNLOAD_URL,SOURCE,NCEI_ID",
        "returnGeometry": "false",
        "f": "json",
    }).get("features", [])

def point_query(lon, lat):
    return get_json(LAYER, {
        "geometry": f"{lon},{lat}",
        "geometryType": "esriGeometryPoint",
        "inSR": "4326",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "SURVEY_ID,PLATFORM,INSTRUMENT,SURVEY_YEAR,START_TIME,END_TIME,DOWNLOAD_URL,SOURCE,NCEI_ID",
        "returnGeometry": "false",
        "f": "json",
    }).get("features", [])

# Truth is used here only for source/coverage binding. This stage NEVER reads bathymetric depth values
# and is not a detector/scoring stage.
shark = next(s for s in registry["sources"] if s["id"] == "NOAA_EX1903L2_DIVE07_SHARK_ROCK")
p = shark["truth_coordinate"]
point_hits = point_query(p["lon"], p["lat"])

survey_ids_to_bind = ["EX1106", "EX1903L2", "EX1903L1"]
survey_bindings = {}
for sid in survey_ids_to_bind:
    feats = survey_query(sid)
    survey_bindings[sid] = [
        dict(f.get("attributes") or {})
        for f in feats
    ]

# Only probe a compact, frozen authority subset so failures do not become scientific negatives.
probe_urls = []
for s in registry["sources"]:
    u = s.get("source_url")
    if u:
        probe_urls.append(u)
    probe_urls.extend(s.get("source_urls") or [])
probe_urls = sorted(set(probe_urls))
probes = [source_probe(u) for u in probe_urls]

point_surveys = sorted({
    str((f.get("attributes") or {}).get("SURVEY_ID"))
    for f in point_hits
    if (f.get("attributes") or {}).get("SURVEY_ID")
})

out = {
    "artifact_id": "JANUS-KUSTO-GLOBAL-GROUNDTRUTH-SOURCE-GATE-RUN-2026-09-23-v1.0",
    "prereg": prereg["artifact_id"],
    "registry": registry["artifact_id"],
    "stage": "G0_SOURCE_AND_PROVENANCE_INVENTORY",
    "depth_values_read": False,
    "detector_run": False,
    "truth_scoring_run": False,
    "ncei_shark_rock_exact_point": {
        "coordinate": p,
        "covering_survey_ids": point_surveys,
        "covering_feature_count": len(point_hits),
        "features": [dict(f.get("attributes") or {}) for f in point_hits],
        "interpretation": "Exact NCEI multibeam footprint intersection only; footprint intersection is not ping/beam insonification proof."
    },
    "ncei_named_survey_bindings": survey_bindings,
    "authority_url_probes": probes,
    "source_probe_reachable_count": sum(1 for x in probes if x.get("reachable")),
    "source_probe_total": len(probes),
    "hard_rules_preserved": [
        "NO_DEPTH_VALUES_READ",
        "NO_DETECTOR_THRESHOLD_CHANGE",
        "FOOTPRINT_INTERSECTION_NE_PING_BEAM_COVERAGE",
        "HTTP_FAILURE_NE_DATA_ABSENCE",
        "TRUTH_USED_FOR_BINDING_ONLY_NOT_DISCOVERY"
    ],
    "next_stage": "G1_FREEZE_TARGET_HIDDEN_INPUT_WINDOWS_AND_NATIVE_RESOLUTION_SCALE_LADDERS",
    "claim_ceiling": "SOURCE_AND_COVERAGE_BINDING_ONLY"
}

path = OUT / "JANUS-KUSTO-GLOBAL-GROUNDTRUTH-SOURCE-GATE-RUN-2026-09-23-v1.0.json"
path.write_text(json.dumps(out, indent=2))
print(json.dumps({
    "artifact_id": out["artifact_id"],
    "shark_rock_covering_surveys": point_surveys,
    "shark_rock_covering_feature_count": len(point_hits),
    "named_survey_binding_counts": {k: len(v) for k, v in survey_bindings.items()},
    "source_probe_reachable": [out["source_probe_reachable_count"], out["source_probe_total"]],
    "depth_values_read": False,
}, indent=2))
