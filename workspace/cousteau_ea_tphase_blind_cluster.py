#!/usr/bin/env python3
"""JANUS Echo Cousteau — blind clustering of the Equatorial Atlantic T-phase catalog.

Scientific contract:
1) acquire the authoritative MGDS catalog;
2) freeze clustering parameters and cluster products WITHOUT the LOVE–EDEM anchor;
3) hash/freeze the blind phase;
4) only then reveal the frozen anchor and score distances.

The raw MGDS catalog is intentionally not committed. The output records provenance,
source SHA-256, parsing metadata, blind clusters, and post-reveal scores.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests
from lxml import etree
from sklearn.cluster import DBSCAN

EARTH_RADIUS_KM = 6371.0088
MGDS_FILESERVER = "https://www.marine-geo.org/services/FileServer"
MGDS_DATASET_PAGE = "https://www.marine-geo.org/tools/files/30497"
MGDS_DOI_PAGE = "https://www.marine-geo.org/doi/10.26022/IEDA/330497"
SOURCE_DOI = "10.26022/IEDA/330497"
EXPECTED_EVENT_COUNT = 6843

# PRE-REGISTERED BEFORE ANCHOR REVEAL. Do not tune after seeing anchor distances.
DBSCAN_GRID = [
    {"eps_km": 25.0, "min_samples": 5},
    {"eps_km": 25.0, "min_samples": 10},
    {"eps_km": 25.0, "min_samples": 20},
    {"eps_km": 50.0, "min_samples": 5},
    {"eps_km": 50.0, "min_samples": 10},
    {"eps_km": 50.0, "min_samples": 20},
    {"eps_km": 100.0, "min_samples": 5},
    {"eps_km": 100.0, "min_samples": 10},
    {"eps_km": 100.0, "min_samples": 20},
]
# Diagnostic look-elsewhere domain; intentionally frozen before anchor reveal.
NULL_DOMAIN = {"lat_min": -15.0, "lat_max": 20.0, "lon_min": -50.0, "lon_max": 10.0}
NULL_SAMPLES = 20000
NULL_SEED = 119520


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def canonical_hash(obj: Any) -> str:
    b = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return sha256_bytes(b)


def haversine_km(lat1, lon1, lat2, lon2):
    lat1 = np.radians(lat1); lon1 = np.radians(lon1)
    lat2 = np.radians(lat2); lon2 = np.radians(lon2)
    dlat = lat2 - lat1; dlon = lon2 - lon1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2.0)**2
    return EARTH_RADIUS_KM * 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))


def spherical_center(lat: np.ndarray, lon: np.ndarray) -> tuple[float, float]:
    latr = np.radians(lat); lonr = np.radians(lon)
    x = np.cos(latr) * np.cos(lonr)
    y = np.cos(latr) * np.sin(lonr)
    z = np.sin(latr)
    x0, y0, z0 = float(np.mean(x)), float(np.mean(y)), float(np.mean(z))
    lon0 = math.degrees(math.atan2(y0, x0))
    hyp = math.hypot(x0, y0)
    lat0 = math.degrees(math.atan2(z0, hyp))
    return lat0, lon0


def text_of_element(el) -> str:
    pieces = [str(v) for v in el.attrib.values() if v is not None]
    for child in el.iter():
        if child.text and child.text.strip(): pieces.append(child.text.strip())
        pieces.extend(str(v) for v in child.attrib.values() if v is not None)
    return " ".join(pieces)


def discover_mgds_file(session: requests.Session) -> dict[str, Any]:
    attempts = []
    query_types = ["Earthquake:Catalog:Microseismicity", "Seismic:Passive"]
    for dtype in query_types:
        params = {
            "minlatitude": -20,
            "maxlatitude": 20,
            "minlongitude": -60,
            "maxlongitude": 15,
            "starttime": "2011-01-01",
            "endtime": "2016-01-01",
            "format": "full_info",
            "data_type": dtype,
        }
        r = session.get(MGDS_FILESERVER, params=params, timeout=90)
        attempts.append({"url": r.url, "status": r.status_code, "bytes": len(r.content), "data_type": dtype})
        if r.status_code != 200 or not r.content.strip():
            continue
        try:
            root = etree.fromstring(r.content)
        except Exception:
            continue
        candidates = []
        for el in root.iter():
            if etree.QName(el).localname.lower() != "file":
                continue
            blob = text_of_element(el)
            low = blob.lower()
            attrs = dict(el.attrib)
            score = 0
            for needle, weight in [
                ("ea_hydroacoustics", 20), ("330497", 20), ("t-phase", 12),
                ("hydroacoustic", 8), ("microseismic", 5), ("catalog", 3), (".txt", 2)
            ]:
                if needle in low: score += weight
            candidates.append({"score": score, "attrs": attrs, "blob": blob[:4000]})
        candidates.sort(key=lambda x: x["score"], reverse=True)
        if candidates and candidates[0]["score"] > 0:
            c = candidates[0]
            url = c["attrs"].get("download") or c["attrs"].get("url")
            uid = (c["attrs"].get("data_uid") or c["attrs"].get("dataUID") or
                   c["attrs"].get("uid") or c["attrs"].get("id"))
            if not url and uid:
                url = f"https://www.marine-geo.org/services/FileDownloadServer?data_uid={uid}"
            if url:
                if url.startswith("http://"):
                    url = "https://" + url[len("http://"):]
                return {"download_url": url, "candidate": c, "attempts": attempts}

    # Fallback: inspect landing pages for direct FileDownloadServer URLs or data UIDs.
    for page in [MGDS_DATASET_PAGE, MGDS_DOI_PAGE]:
        r = session.get(page, timeout=90)
        attempts.append({"url": r.url, "status": r.status_code, "bytes": len(r.content), "fallback_page": True})
        text = r.text
        urls = re.findall(r'https?://[^"\'<> ]*FileDownloadServer[^"\'<> ]*', text, flags=re.I)
        if urls:
            url = urls[0].replace("&amp;", "&")
            if url.startswith("http://"): url = "https://" + url[7:]
            return {"download_url": url, "candidate": {"score": 1, "attrs": {}, "blob": "landing-page URL"}, "attempts": attempts}
        m = re.search(r'(?:data[_-]?uid|file[_-]?uid)["\' :=]+([A-Za-z0-9._-]+)', text, flags=re.I)
        if m:
            uid = m.group(1)
            return {"download_url": f"https://www.marine-geo.org/services/FileDownloadServer?data_uid={uid}",
                    "candidate": {"score": 1, "attrs": {"uid": uid}, "blob": "landing-page UID"}, "attempts": attempts}

    raise RuntimeError("MGDS catalog file could not be discovered from FileServer or landing pages")


def download_catalog(session: requests.Session, url: str) -> tuple[bytes, dict[str, Any]]:
    r = session.get(url, timeout=120, allow_redirects=True)
    r.raise_for_status()
    b = r.content
    if len(b) < 1000:
        raise RuntimeError(f"download too small ({len(b)} bytes); final_url={r.url}; prefix={b[:300]!r}")
    # Reject obvious HTML error/login page.
    prefix = b[:500].lower()
    if b"<html" in prefix or b"<!doctype html" in prefix:
        raise RuntimeError(f"download returned HTML instead of ASCII catalog; final_url={r.url}")
    return b, {"requested_url": url, "final_url": r.url, "bytes": len(b), "content_type": r.headers.get("content-type")}


def normalize_col(c: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(c).strip().lower())


def parse_catalog(raw: bytes) -> tuple[pd.DataFrame, dict[str, Any]]:
    text = raw.decode("utf-8", errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise RuntimeError("catalog is empty after decode")

    # Try common delimiters and header conventions, including a leading '#'.
    variants = []
    clean = "\n".join((ln[1:].lstrip() if i == 0 and ln.lstrip().startswith("#") else ln) for i, ln in enumerate(lines))
    for sep in [r"\s+", "\t", ",", r"\s*,\s*"]:
        try:
            df = pd.read_csv(io.StringIO(clean), sep=sep, engine="python")
            if len(df) > 0 and len(df.columns) >= 2:
                variants.append((df, sep))
        except Exception:
            pass
    if not variants:
        raise RuntimeError("could not parse catalog with whitespace/tab/comma strategies")

    lat_alias = {"lat", "latitude", "eventlat", "eventlatitude", "originlat", "originlatitude"}
    lon_alias = {"lon", "long", "longitude", "eventlon", "eventlong", "eventlongitude", "originlon", "originlongitude"}

    best = None
    for df, sep in variants:
        norm = {normalize_col(c): c for c in df.columns}
        lat_col = next((norm[a] for a in lat_alias if a in norm), None)
        lon_col = next((norm[a] for a in lon_alias if a in norm), None)
        if lat_col is not None and lon_col is not None:
            score = -abs(len(df) - EXPECTED_EVENT_COUNT)
            if best is None or score > best[0]: best = (score, df, sep, lat_col, lon_col)
    if best is None:
        # Last resort: infer two numeric columns whose ranges look like lat/lon.
        for df, sep in variants:
            numeric = {}
            for c in df.columns:
                s = pd.to_numeric(df[c], errors="coerce")
                if s.notna().mean() > 0.9: numeric[c] = s
            for lc, ls in numeric.items():
                if not ((ls >= -90) & (ls <= 90)).mean() > 0.98: continue
                for oc, os_ in numeric.items():
                    if oc == lc: continue
                    if not ((os_ >= -180) & (os_ <= 180)).mean() > 0.98: continue
                    # Equatorial Atlantic catalog should have lat/lon spread and negative longitudes.
                    if ls.std() > 0.2 and os_.std() > 0.2:
                        score = -abs(len(df) - EXPECTED_EVENT_COUNT)
                        if best is None or score > best[0]: best = (score, df, sep, lc, oc)
    if best is None:
        raise RuntimeError(f"could not identify latitude/longitude columns; columns={[str(c) for c in variants[0][0].columns]}")

    _, df, sep, lat_col, lon_col = best
    lat = pd.to_numeric(df[lat_col], errors="coerce")
    lon = pd.to_numeric(df[lon_col], errors="coerce")
    mask = lat.between(-90, 90) & lon.between(-180, 180)
    out = pd.DataFrame({"lat": lat[mask].astype(float), "lon": lon[mask].astype(float)}).reset_index(drop=True)
    meta = {
        "delimiter_strategy": sep,
        "raw_row_count": int(len(df)),
        "valid_coordinate_count": int(len(out)),
        "lat_column": str(lat_col),
        "lon_column": str(lon_col),
        "columns": [str(c) for c in df.columns],
        "first_header_line": lines[0][:500],
    }
    return out, meta


def blind_cluster(coords: pd.DataFrame) -> dict[str, Any]:
    """ANCHOR-FREE function by contract."""
    X_deg = coords[["lat", "lon"]].to_numpy(dtype=float)
    X_rad = np.radians(X_deg)
    configs = []
    for cfg in DBSCAN_GRID:
        model = DBSCAN(
            eps=cfg["eps_km"] / EARTH_RADIUS_KM,
            min_samples=cfg["min_samples"],
            metric="haversine",
            algorithm="ball_tree",
            n_jobs=-1,
        )
        labels = model.fit_predict(X_rad)
        cluster_ids = sorted(int(x) for x in set(labels) if x >= 0)
        clusters = []
        for cid in cluster_ids:
            idx = np.flatnonzero(labels == cid)
            pts = X_deg[idx]
            clat, clon = spherical_center(pts[:, 0], pts[:, 1])
            d = haversine_km(pts[:, 0], pts[:, 1], clat, clon)
            clusters.append({
                "cluster_id": cid,
                "n": int(len(idx)),
                "center_lat": round(float(clat), 7),
                "center_lon": round(float(clon), 7),
                "radius_p50_km": round(float(np.quantile(d, 0.50)), 3),
                "radius_p95_km": round(float(np.quantile(d, 0.95)), 3),
                "radius_max_km": round(float(np.max(d)), 3),
            })
        clusters.sort(key=lambda c: (-c["n"], c["cluster_id"]))
        configs.append({
            "eps_km": cfg["eps_km"],
            "min_samples": cfg["min_samples"],
            "cluster_count": len(clusters),
            "clustered_events": int(np.sum(labels >= 0)),
            "noise_events": int(np.sum(labels < 0)),
            "noise_fraction": round(float(np.mean(labels < 0)), 6),
            "clusters": clusters,
        })
    blind = {
        "anchor_visible": False,
        "metric": "DBSCAN_HAVERSINE_GREAT_CIRCLE",
        "earth_radius_km": EARTH_RADIUS_KM,
        "parameter_grid": DBSCAN_GRID,
        "event_count": int(len(coords)),
        "catalog_bbox": {
            "lat_min": float(coords.lat.min()), "lat_max": float(coords.lat.max()),
            "lon_min": float(coords.lon.min()), "lon_max": float(coords.lon.max()),
        },
        "configs": configs,
    }
    blind["freeze_sha256"] = canonical_hash(blind)
    return blind


def sample_area_uniform_latlon(rng: np.random.Generator, n: int, dom: dict[str, float]):
    lon = rng.uniform(dom["lon_min"], dom["lon_max"], n)
    s0, s1 = math.sin(math.radians(dom["lat_min"])), math.sin(math.radians(dom["lat_max"]))
    lat = np.degrees(np.arcsin(rng.uniform(s0, s1, n)))
    return lat, lon


def reveal_and_score(coords: pd.DataFrame, blind: dict[str, Any], anchor_lat: float, anchor_lon: float) -> dict[str, Any]:
    # This is the FIRST stage that is allowed to access anchor coordinates.
    direct = haversine_km(coords.lat.to_numpy(), coords.lon.to_numpy(), anchor_lat, anchor_lon)
    nearest_idx = int(np.argmin(direct))
    reveal = {
        "anchor_lat": anchor_lat,
        "anchor_lon": anchor_lon,
        "blind_freeze_sha256_verified": canonical_hash({k: v for k, v in blind.items() if k != "freeze_sha256"}) == blind["freeze_sha256"],
        "nearest_event": {
            "distance_km": round(float(direct[nearest_idx]), 3),
            "lat": round(float(coords.lat.iloc[nearest_idx]), 7),
            "lon": round(float(coords.lon.iloc[nearest_idx]), 7),
            "catalog_row_zero_based": nearest_idx,
        },
        "configs": [],
    }
    rng = np.random.default_rng(NULL_SEED)
    null_lat, null_lon = sample_area_uniform_latlon(rng, NULL_SAMPLES, NULL_DOMAIN)

    for cfg in blind["configs"]:
        centers = cfg["clusters"]
        if not centers:
            reveal["configs"].append({"eps_km": cfg["eps_km"], "min_samples": cfg["min_samples"], "nearest_cluster": None})
            continue
        clat = np.array([c["center_lat"] for c in centers], dtype=float)
        clon = np.array([c["center_lon"] for c in centers], dtype=float)
        d = haversine_km(clat, clon, anchor_lat, anchor_lon)
        j = int(np.argmin(d))
        nearest = centers[j]
        anchor_d = float(d[j])

        # Diagnostic rectangular look-elsewhere null only; not an ocean-masked formal p-value.
        # Chunk to avoid a large N_null x N_cluster matrix.
        null_min = np.full(NULL_SAMPLES, np.inf)
        for k in range(len(centers)):
            dk = haversine_km(null_lat, null_lon, clat[k], clon[k])
            null_min = np.minimum(null_min, dk)
        percentile = float(np.mean(null_min <= anchor_d))
        reveal["configs"].append({
            "eps_km": cfg["eps_km"],
            "min_samples": cfg["min_samples"],
            "nearest_cluster": {
                **nearest,
                "anchor_to_center_km": round(anchor_d, 3),
                "anchor_inside_cluster_p95_radius": bool(anchor_d <= nearest["radius_p95_km"]),
                "anchor_inside_cluster_max_radius": bool(anchor_d <= nearest["radius_max_km"]),
            },
            "diagnostic_rectangular_null": {
                "domain": NULL_DOMAIN,
                "samples": NULL_SAMPLES,
                "seed": NULL_SEED,
                "fraction_random_points_with_equal_or_smaller_nearest_cluster_distance": round(percentile, 6),
                "formal_p_value": False,
                "warning": "RECTANGULAR_DOMAIN_NOT_OCEAN_MASKED__DIAGNOSTIC_ONLY",
            },
        })
    return reveal


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--status-output", required=False)
    args = ap.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    status_path = Path(args.status_output) if args.status_output else out.with_name(out.stem + "-STATUS.json")

    status = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-001-2026-08-21-v1.0",
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "RUNNING",
    }
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")

    try:
        s = requests.Session()
        s.headers.update({"User-Agent": "Janus-Echo-Cousteau/1.0 scientific reproducibility audit"})
        discovery = discover_mgds_file(s)
        raw, dlmeta = download_catalog(s, discovery["download_url"])
        coords, parse_meta = parse_catalog(raw)

        source_hash = sha256_bytes(raw)
        blind = blind_cluster(coords)  # no anchor passed or available to this function

        # Gate: blind object is now frozen. Only below this line is the project anchor revealed.
        FROZEN_ANCHOR_LAT = -3.865418
        FROZEN_ANCHOR_LON = 3.854924
        reveal = reveal_and_score(coords, blind, FROZEN_ANCHOR_LAT, FROZEN_ANCHOR_LON)

        expected_match = len(coords) == EXPECTED_EVENT_COUNT
        nearest_event_km = reveal["nearest_event"]["distance_km"]
        nearest_cluster_distances = [
            c["nearest_cluster"]["anchor_to_center_km"]
            for c in reveal["configs"] if c.get("nearest_cluster")
        ]
        any_p95 = any(
            c.get("nearest_cluster") and c["nearest_cluster"]["anchor_inside_cluster_p95_radius"]
            for c in reveal["configs"]
        )
        verdict = (
            "ANCHOR_OVERLAPS_BLIND_CLUSTER_P95__REQUIRES_LOOK_ELSEWHERE_AND_TECTONIC_CONTROL"
            if any_p95 else
            "NO_BLIND_CLUSTER_P95_OVERLAP_WITH_FROZEN_ANCHOR"
        )

        result = {
            "artifact_id": "JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-001-2026-08-21-v1.0",
            "research_branch": "Janus-Echo-Кусто",
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "purpose": "Blindly cluster the public 2011-2015 Equatorial Atlantic T-phase event catalog before revealing the LOVE-EDEM frozen anchor.",
            "source": {
                "doi": SOURCE_DOI,
                "dataset": "EA_Hydroacoustics",
                "dataset_page": MGDS_DATASET_PAGE,
                "expected_event_count": EXPECTED_EVENT_COUNT,
                "parsed_valid_coordinate_count": int(len(coords)),
                "expected_count_exact_match": expected_match,
                "raw_catalog_sha256": source_hash,
                "raw_catalog_committed": False,
                "license": "CC BY-NC-SA 3.0",
                "download": dlmeta,
                "discovery": discovery,
                "parse": parse_meta,
            },
            "preregistration": {
                "anchor_hidden_during_clustering": True,
                "clustering_parameters_frozen_before_anchor_reveal": True,
                "dbscan_grid": DBSCAN_GRID,
                "metric": "haversine",
                "null_domain": NULL_DOMAIN,
                "null_samples": NULL_SAMPLES,
                "null_seed": NULL_SEED,
            },
            "blind_phase": blind,
            "post_reveal": reveal,
            "summary": {
                "nearest_catalog_event_to_anchor_km": nearest_event_km,
                "nearest_blind_cluster_center_across_grid_km": round(min(nearest_cluster_distances), 3) if nearest_cluster_distances else None,
                "anchor_inside_any_blind_cluster_p95_radius": any_p95,
                "verdict": verdict,
                "semantic_status": "UNCONFIRMED",
            },
            "hard_rules": [
                "BLIND_CLUSTER_BEFORE_ANCHOR_REVEAL",
                "DO_NOT_RETUNE_EPS_OR_MIN_SAMPLES_AFTER_REVEAL",
                "MID_ATLANTIC_RIDGE_SEISMICITY_IS_A_MANDATORY_TECTONIC_CONTROL",
                "RECTANGULAR_LOOK_ELSEWHERE_NULL_IS_DIAGNOSTIC_NOT_FORMAL",
                "DISTANCE_TO_CLUSTER_IS_NOT_CAUSATION",
                "NO_RECENTERING",
                "NO_UNDERWATER_PYRAMID_DETECTED_YET",
            ],
            "status": "BLIND_CLUSTER_RUN_COMPLETE",
        }
        out.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
        status.update({
            "status": "SUCCESS",
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "result_path": str(out),
            "source_sha256": source_hash,
            "parsed_event_count": int(len(coords)),
            "verdict": verdict,
        })
        status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        print(json.dumps(result["summary"], indent=2))
        return 0
    except Exception as e:
        status.update({
            "status": "BLOCKED_DATA_ACQUISITION_OR_PARSE",
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "error_type": type(e).__name__,
            "error": str(e),
        })
        status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        print(json.dumps(status, indent=2), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
