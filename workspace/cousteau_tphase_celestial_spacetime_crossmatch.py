#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from cousteau_ea_tphase_blind_cluster_v7 import acquire_exact_file, parse_exact
from cousteau_ea_tphase_blind_cluster import sha256_bytes

EARTH_RADIUS_KM = 6371.0088
LOVE_RA = 204.30267916666668
LOVE_DEC = -36.78240527777778
EDEM_RA = 139.22409686590188
EDEM_DEC = 30.26038779947318
BISECTOR_RA = 170.38390326101873
BISECTOR_DEC = -3.8654180644718967
RELAY_ALT = 44.70660585548666
N_PERM = 2000
NULL_SEED = 170383903
TOP_N = 25


def wrap180(x):
    return (np.asarray(x, dtype=float) + 180.0) % 360.0 - 180.0


def parse_time_code(code: str) -> datetime:
    s = str(code).strip()
    if len(s) != 14 or not s.isdigit():
        raise ValueError(f"invalid source_time_code {s!r}")
    year = int(s[:4])
    doy = int(s[4:7])
    hour = int(s[7:9])
    minute = int(s[9:11])
    sec_tenths = int(s[11:14])
    second = sec_tenths / 10.0
    dt = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy - 1, hours=hour, minutes=minute, seconds=second)
    return dt


def jd_from_datetime(dt: datetime) -> float:
    return dt.timestamp() / 86400.0 + 2440587.5


def gmst_deg(dt: datetime) -> float:
    jd = jd_from_datetime(dt)
    T = (jd - 2451545.0) / 36525.0
    g = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T*T - (T*T*T) / 38710000.0
    return g % 360.0


def haversine_angle_deg(lat1, lon1, lat2, lon2):
    lat1 = np.radians(np.asarray(lat1, dtype=float)); lon1 = np.radians(np.asarray(lon1, dtype=float))
    lat2 = np.radians(np.asarray(lat2, dtype=float)); lon2 = np.radians(np.asarray(lon2, dtype=float))
    dlat = lat2-lat1; dlon = lon2-lon1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2.0)**2
    return np.degrees(2.0*np.arcsin(np.minimum(1.0, np.sqrt(a))))


def altitude_deg(lat_deg, lon_deg, gmst, ra_deg, dec_deg):
    lat = np.radians(np.asarray(lat_deg, dtype=float))
    dec = math.radians(dec_deg)
    lst = (np.asarray(gmst, dtype=float) + np.asarray(lon_deg, dtype=float)) % 360.0
    H = np.radians(wrap180(lst - ra_deg))
    sin_alt = np.sin(lat)*math.sin(dec) + np.cos(lat)*math.cos(dec)*np.cos(H)
    return np.degrees(np.arcsin(np.clip(sin_alt, -1.0, 1.0)))


def row_payload(full, i, gmst, sub_lon, sep, love_alt, edem_alt):
    r = full.iloc[int(i)]
    return {
        "rank": None,
        "catalog_row_zero_based": int(i),
        "source_line": int(r.source_line),
        "source_time_code": str(r.source_time_code),
        "event_time_utc": parse_time_code(str(r.source_time_code)).isoformat().replace('+00:00','Z'),
        "event_lat": float(r.lat),
        "event_lon": float(r.lon),
        "n_hydrophones": int(r.n_hydrophones),
        "source_magnitude_db": float(r.source_magnitude_db),
        "reported_lat_error_deg": float(r.lat_error_deg),
        "reported_lon_error_deg": float(r.lon_error_deg),
        "gmst_deg": round(float(gmst[i]), 9),
        "bisector_subpoint_lat_deg": BISECTOR_DEC,
        "bisector_subpoint_lon_deg_east": round(float(sub_lon[i]), 9),
        "event_to_bisector_subpoint_angular_deg": round(float(sep[i]), 9),
        "event_to_bisector_subpoint_km": round(float(math.radians(float(sep[i]))*EARTH_RADIUS_KM), 3),
        "bisector_altitude_at_event_deg": round(90.0-float(sep[i]), 9),
        "love_alt_deg": round(float(love_alt[i]), 9),
        "edem_alt_deg": round(float(edem_alt[i]), 9),
        "love_edem_altitude_difference_abs_deg": round(abs(float(love_alt[i]-edem_alt[i])), 9),
        "mean_love_edem_alt_deg": round(float((love_alt[i]+edem_alt[i])/2.0), 9),
        "mean_altitude_residual_from_relay_44_7066_deg": round(abs(float((love_alt[i]+edem_alt[i])/2.0)-RELAY_ALT), 9),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True)
    ap.add_argument('--status-output', required=True)
    args = ap.parse_args()
    out = Path(args.output); status_out = Path(args.status_output)
    out.parent.mkdir(parents=True, exist_ok=True)

    started = datetime.now(timezone.utc)
    archive, gz, raw, member, trace = acquire_exact_file()
    full, parse_meta = parse_exact(raw)
    n = len(full)
    if n != 5943:
        raise RuntimeError(f"authoritative row count drift: expected current 5943, got {n}")

    times = [parse_time_code(str(x)) for x in full.source_time_code]
    gmst = np.array([gmst_deg(t) for t in times], dtype=float)
    sub_lon = wrap180(BISECTOR_RA - gmst)
    event_lat = full.lat.to_numpy(float); event_lon = full.lon.to_numpy(float)
    sep = haversine_angle_deg(event_lat, event_lon, BISECTOR_DEC, sub_lon)
    love_alt = altitude_deg(event_lat, event_lon, gmst, LOVE_RA, LOVE_DEC)
    edem_alt = altitude_deg(event_lat, event_lon, gmst, EDEM_RA, EDEM_DEC)

    order = np.argsort(sep)
    top = []
    for rank, i in enumerate(order[:TOP_N], start=1):
        p = row_payload(full, int(i), gmst, sub_lon, sep, love_alt, edem_alt)
        p['rank'] = rank
        top.append(p)

    observed_min = float(sep[order[0]])
    rng = np.random.default_rng(NULL_SEED)
    null_mins = np.empty(N_PERM, dtype=float)
    # Freeze both marginals; permute only timestamp-derived longitude assignments among fixed event positions.
    for k in range(N_PERM):
        perm = rng.permutation(n)
        d = haversine_angle_deg(event_lat, event_lon, BISECTOR_DEC, sub_lon[perm])
        null_mins[k] = float(np.min(d))
    global_p = (1.0 + float(np.sum(null_mins <= observed_min))) / (N_PERM + 1.0)

    # Secondary rank: actual local LOVE/EDEM relay symmetry. This is descriptive and not a new discovery statistic.
    symmetry_score = np.abs(love_alt-edem_alt) + np.abs((love_alt+edem_alt)/2.0 - RELAY_ALT)
    sym_order = np.argsort(symmetry_score)
    top_sym = []
    for rank, i in enumerate(sym_order[:TOP_N], start=1):
        p = row_payload(full, int(i), gmst, sub_lon, sep, love_alt, edem_alt)
        p['rank'] = rank
        p['relay_symmetry_score_deg'] = round(float(symmetry_score[i]), 9)
        top_sym.append(p)

    result = {
        "artifact_id": "JANUS-ECHO-COUSTEAU-TPHASE-CELESTIAL-SPACETIME-CROSSMATCH-RUN-001-2026-08-21-v1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "BLIND_ALL_EVENT_DATE_LOCATION_TO_FROZEN_LOVE_EDEM_CELESTIAL_GEOMETRY_CROSSMATCH",
        "source": {
            "doi": "10.26022/IEDA/330497",
            "dataset_uid": 30497,
            "file_uid": "2504732",
            "member": member,
            "authoritative_rows": n,
            "catalog_ascii_sha256": sha256_bytes(raw),
            "archive_sha256": sha256_bytes(archive),
            "parse": parse_meta,
            "acquisition_trace": trace,
        },
        "frozen_celestial_model": {
            "LOVE": {"ra_deg_icrs_j2000": LOVE_RA, "dec_deg_icrs_j2000": LOVE_DEC},
            "EDEM": {"ra_deg_icrs_j2000": EDEM_RA, "dec_deg_icrs_j2000": EDEM_DEC, "identity_confirmed": False},
            "BISECTOR": {"ra_deg_icrs_j2000": BISECTOR_RA, "dec_deg_icrs_j2000": BISECTOR_DEC},
            "relay_equal_altitude_target_deg": RELAY_ALT,
            "earth_subpoint_rule": "lat=bisector_dec; lon_east=wrap180(bisector_RA-GMST(event_UTC))",
        },
        "time_parser": {
            "format": "YYYYDDDHHMMSSS",
            "final_three_digits": "seconds_in_tenths",
            "gmst_formula": "Meeus-style polynomial using UTC-derived JD; sufficient for km-scale screening, not precision astrometry",
        },
        "primary_metric": "great_circle_angular_distance_between_event_location_and_time_specific_bisector_zenith_subpoint",
        "top_matches_by_bisector_subpoint_distance": top,
        "top_matches_by_love_edem_relay_symmetry_score": top_sym,
        "look_elsewhere": {
            "null": "permute_timestamp_derived_bisector_subpoint_longitudes_among_the_same_5943_event_locations",
            "permutations": N_PERM,
            "seed": NULL_SEED,
            "observed_min_angular_deg": round(observed_min, 9),
            "observed_min_km": round(math.radians(observed_min)*EARTH_RADIUS_KM, 3),
            "null_min_quantiles_deg": {
                "q01": round(float(np.quantile(null_mins,0.01)),9),
                "q05": round(float(np.quantile(null_mins,0.05)),9),
                "q50": round(float(np.quantile(null_mins,0.50)),9),
                "q95": round(float(np.quantile(null_mins,0.95)),9),
            },
            "empirical_global_p": round(global_p, 9),
            "formal_claim": "DIAGNOSTIC_GLOBAL_TIME_LOCATION_PAIRING_NULL__NOT_CAUSAL_PROOF",
        },
        "scientific_interpretation": {
            "target_evidence_can_increase_from_best_match_alone": False,
            "tectonic_mar_control_remains_mandatory": True,
            "historical_anomaly_overlay_must_be_separate": True,
            "claim_ceiling": "CROSSMATCH_RANKING_AND_LOOK_ELSEWHERE_DIAGNOSTIC_ONLY",
        },
        "hard_rules": [
            "ALL_5943_EVENTS_SCORED_BEFORE_HISTORICAL_CASE_SELECTION",
            "NO_POSTHOC_RETARGET_TO_BEST_LONGITUDE",
            "BEST_OF_5943_REQUIRES_LOOK_ELSEWHERE_CONTROL",
            "EDEM_IDENTITY_UNCONFIRMED",
            "DISTANCE_OR_ALIGNMENT_IS_NOT_CAUSATION",
            "MID_ATLANTIC_RIDGE_TECTONIC_CONTROL_PRECEDES_TARGET_INTERPRETATION",
            "NO_UNDERWATER_PYRAMID_DETECTED_YET",
        ],
        "status": "RUN_COMPLETE",
    }
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    status = {
        "artifact_id": result['artifact_id'],
        "status": "SUCCESS",
        "authoritative_rows": n,
        "best_event_time_utc": top[0]['event_time_utc'],
        "best_event_lat": top[0]['event_lat'],
        "best_event_lon": top[0]['event_lon'],
        "best_event_distance_km": top[0]['event_to_bisector_subpoint_km'],
        "best_event_angular_deg": top[0]['event_to_bisector_subpoint_angular_deg'],
        "empirical_global_p": result['look_elsewhere']['empirical_global_p'],
        "target_evidence": "NOT_PROMOTED_BY_CROSSMATCH_ALONE",
    }
    status_out.write_text(json.dumps(status, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
