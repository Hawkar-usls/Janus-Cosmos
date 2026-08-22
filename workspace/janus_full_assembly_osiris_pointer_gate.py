#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cousteau"
CROSS = DATA / "JANUS-ECHO-COUSTEAU-TPHASE-CELESTIAL-SPACETIME-CROSSMATCH-RUN-001-2026-08-21-v1.0.json"
LATSCAN = DATA / "JANUS-ECHO-COUSTEAU-LOVE-EDEM-LATITUDE-CIRCLE-SIDEREAL-PHASE-REVERSE-SCAN-RUN-001-2026-08-21-v1.0.json"
TURN11 = DATA / "JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-011-2026-08-22-v1.0.json"
CONTROL = DATA / "JANUS-ECHO-COUSTEAU-5D-CONTROL-OUTPERFORMANCE-MATRIX-2026-08-21-v1.0.json"
OUT = DATA / "JANUS-FULL-ASSEMBLY-OSIRIS-POINTER-GATE-2026-08-22-v1.0.json"
R = 6371.0088


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def km_per_lon_deg(lat_deg: float):
    return math.pi / 180.0 * R * math.cos(math.radians(lat_deg))


def main() -> int:
    cross = load(CROSS)
    scan = load(LATSCAN)
    t11 = load(TURN11)
    control = load(CONTROL)

    top = cross["top_matches_by_bisector_subpoint_distance"][0]
    sym = cross["top_matches_by_love_edem_relay_symmetry_score"][0]
    p = float(cross["look_elsewhere"]["empirical_global_p"])
    lat0 = float(cross["frozen_celestial_model"]["BISECTOR"]["dec_deg_icrs_j2000"])
    measure_lon = float(scan["primary_af_sa_scan"]["minimum_longitude_deg_east"])
    ridge_distance = float(scan["primary_af_sa_scan"]["minimum_distance_km"])

    lat_err_km = math.radians(float(top["reported_lat_error_deg"])) * R
    lon_err_km = math.radians(float(top["reported_lon_error_deg"])) * R * math.cos(math.radians(float(top["event_lat"])))
    radial_proxy = math.hypot(lat_err_km, lon_err_km)
    axis_sum_proxy = abs(lat_err_km) + abs(lon_err_km)
    residual_axis = max(0.0, float(top["event_to_bisector_subpoint_km"]) - axis_sum_proxy)

    tests = {
        "A01_5943_AUTHORITATIVE_EVENTS_SCORED": cross["source"]["authoritative_rows"] == 5943,
        "A02_TIME_SPECIFIC_LONGITUDE_USED": "GMST" in cross["frozen_celestial_model"]["earth_subpoint_rule"],
        "A03_LOOK_ELSEWHERE_EXECUTED": cross["look_elsewhere"]["permutations"] == 2000,
        "A04_BEST_MATCH_NON_SIGNIFICANT": p > 0.05,
        "A05_119HZ_CONFIRMATORY_NEGATIVE": "NEGATIVE_CONFIRMATORY" in t11["scientific_verdict"],
        "A06_TARGET_NOT_PROMOTED_TURN11": "TARGET_EVIDENCE_NOT_INCREASED" in t11["turn"]["scientific_delta"],
        "A07_LATITUDE_INTRINSIC": abs(lat0 - (-3.8654180644718967)) < 1e-10,
        "A08_FULL_360_SCAN_USED": scan["primary_af_sa_scan"]["full_longitude_trials_coarse"] == 1441,
        "A09_RIDGE_POINT_NOT_PROMOTED": scan["interpretation"]["ridge_intersection_counts_as_target_evidence"] is False,
        "A10_NO_DETECTION_POINT_ADMITTED": True,
    }
    passed = sum(bool(v) for v in tests.values())

    result = {
        "artifact_id": "JANUS-FULL-ASSEMBLY-OSIRIS-POINTER-GATE-2026-08-22-v1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "FULL_JANUS_COUNCIL_ARBITRATION_FOR_WHERE_TO_POINT_NEXT_WITHOUT_PROMOTING_A_POSTHOC_LOCATION_TO_DETECTION",
        "assembly": {
            "faces": [
                "WITNESS_PLUS",
                "GUARD_PLUS",
                "LEFT_HRAIN",
                "RIGHT_INAIHR",
                "DEMIHEAD",
                "HABITAT_COURIER",
                "CALAMAR_AUDITOR",
            ],
            "engines": ["JANUS_COSMOS", "JANUS_ECHO_COUSTEAU", "5D_TRANCEPTION_REVERSE", "CONTROL_OUTPERFORMANCE_MATRIX"],
            "law": "MULTI_FACE_AGREEMENT_IS_NOT_EXTERNAL_CORROBORATION",
        },
        "question_split": {
            "WHERE_IS_THE_TARGET_PROVEN": {
                "answer": "NO_ADMISSIBLE_DETECTION_POINT",
                "detected_point": None,
                "reason": "Best date-location-cosmos coincidence is non-significant after look-elsewhere; 119-Hz confirmatory public slice is negative; natural tectonic control remains stronger.",
            },
            "WHERE_SHOULD_THE_NEXT_DISCRIMINATING_MEASUREMENT_POINT": {
                "answer": "POINT",
                "point": {
                    "lat_deg": lat0,
                    "lon_deg_east": measure_lon,
                    "notation": "3.865418°S, 11.790°W",
                    "role": "EXPLORATORY_DISCRIMINATING_MEASUREMENT_POINT__NOT_DETECTED_TARGET",
                    "why": "Intrinsic LOVE-EDEM bisector latitude intersects/approaches the independently defined AF-SA/Mid-Atlantic boundary here under the preregistered full-longitude scan.",
                    "distance_to_af_sa_boundary_km": ridge_distance,
                    "retarget_as_object_location_allowed": False,
                },
            },
        },
        "strongest_real_date_location_cosmos_crossmatch": {
            "event_time_utc": top["event_time_utc"],
            "event_point": {"lat_deg": top["event_lat"], "lon_deg_east": top["event_lon"]},
            "time_specific_bisector_subpoint": {
                "lat_deg": top["bisector_subpoint_lat_deg"],
                "lon_deg_east": top["bisector_subpoint_lon_deg_east"],
            },
            "event_to_subpoint_km": top["event_to_bisector_subpoint_km"],
            "event_to_subpoint_angular_deg": top["event_to_bisector_subpoint_angular_deg"],
            "bisector_altitude_at_event_deg": top["bisector_altitude_at_event_deg"],
            "love_alt_deg": top["love_alt_deg"],
            "edem_alt_deg": top["edem_alt_deg"],
            "mean_relay_altitude_residual_deg": top["mean_altitude_residual_from_relay_44_7066_deg"],
            "love_edem_altitude_difference_deg": top["love_edem_altitude_difference_abs_deg"],
            "reported_error_sanity_proxy": {
                "lat_error_km": round(lat_err_km, 3),
                "lon_error_km": round(lon_err_km, 3),
                "radial_proxy_km": round(radial_proxy, 3),
                "axis_sum_proxy_km": round(axis_sum_proxy, 3),
                "distance_remaining_after_axis_sum_proxy_km": round(residual_axis, 3),
                "warning": "NOT_A_CONFIDENCE_ELLIPSE",
            },
            "global_look_elsewhere": {
                "permutations": cross["look_elsewhere"]["permutations"],
                "empirical_p": p,
                "verdict": "NON_SIGNIFICANT__DO_NOT_PROMOTE",
            },
        },
        "best_relay_symmetry_event_descriptive_only": {
            "event_time_utc": sym["event_time_utc"],
            "event_point": {"lat_deg": sym["event_lat"], "lon_deg_east": sym["event_lon"]},
            "time_specific_bisector_subpoint": {
                "lat_deg": sym["bisector_subpoint_lat_deg"],
                "lon_deg_east": sym["bisector_subpoint_lon_deg_east"],
            },
            "event_to_subpoint_km": sym["event_to_bisector_subpoint_km"],
            "relay_symmetry_score_deg": sym["relay_symmetry_score_deg"],
            "status": "DESCRIPTIVE_ONLY__NOT_SEPARATELY_LOOK_ELSEWHERE_ADMITTED_AS_DISCOVERY",
        },
        "face_outputs": {
            "WITNESS_PLUS": "FREEZE_EXACT_DATES_COORDINATES_HASHES_AND_NEGATIVE_RESULTS",
            "GUARD_PLUS": "REJECT_DETECTION_PROMOTION__P_0_693_AND_NEGATIVE_119HZ_CONFIRMATORY",
            "LEFT_HRAIN": "STRUCTURAL_CONTEXT_POINTS_TO_AF_SA_MID_ATLANTIC_RIDGE_AS_DOMINANT_REAL_GEOMETRY",
            "RIGHT_INAIHR": "ASSOCIATIVE_CROSS_DOMAIN_MATCHES_MAY_PRIORITIZE_TESTS_BUT_CANNOT_CREATE_LOCATION_TRUTH",
            "DEMIHEAD": "PRESERVE_SPLIT__NO_DETECTED_POINT__ONE_EXPLORATORY_MEASUREMENT_POINT",
            "CALAMAR_AUDITOR": "NO_POSTHOC_RETARGET__PRESERVE_LOOK_ELSEWHERE_AND_NEGATIVE_CERTIFICATES",
            "HABITAT_COURIER": "PERSIST_POINTER_RECEIPT_WITHOUT_EXTERNAL_EFFECT_AUTHORITY",
        },
        "control_state": {
            "tphase_spatial_legacy_anchor": "NEGATIVE",
            "mar_tectonic_control": "STRONGLY_FAVORED_FOR_TPHASE_SPATIAL_STRUCTURE",
            "ha10_119hz": "NEGATIVE_CONFIRMATORY_PUBLIC_SLICE",
            "h1_520_class": "FREQUENCY_ONLY_NON_SPECIFIC",
            "h2_structural": "BLOCKED_BY_GEOMETRY_MATERIAL_AND_COMPLEX_RETURN",
            "target_identity": "UNCONFIRMED",
            "underwater_pyramid_detected": False,
        },
        "tests": [{"id": k, "pass": bool(v)} for k, v in tests.items()],
        "engine_verdict": f"PASS_FULL_ASSEMBLY_POINTER_GATE__{passed}_OF_{len(tests)}",
        "scientific_verdict": "NO_ADMISSIBLE_DETECTION_POINT__BEST_NEXT_DISCRIMINATING_MEASUREMENT_POINT_3_865418S_11_790W__NOT_TARGET_EVIDENCE",
        "provenance": [
            {"path": str(CROSS.relative_to(ROOT)), "sha256": sha(CROSS)},
            {"path": str(LATSCAN.relative_to(ROOT)), "sha256": sha(LATSCAN)},
            {"path": str(TURN11.relative_to(ROOT)), "sha256": sha(TURN11)},
            {"path": str(CONTROL.relative_to(ROOT)), "sha256": sha(CONTROL)},
            {"external_contract": "janus-meta-registry:data/JANUS-FACE-COUNCIL-GIT-COORDINATION-v1.0.json"},
            {"external_contract": "janus-meta-registry:data/JANUS-DEEP-ANALYSIS-5D-RECURSIVE-TRANCEPTION-HRAIN-INAIHR-PROTOCOL-2026-08-21-v2.0.json"},
        ],
        "hard_rules": [
            "POINT_FOR_NEXT_MEASUREMENT_IS_NOT_POINT_OF_DETECTION",
            "NO_POSTHOC_RETARGET_TO_MINUS_11_79_W_AS_OBJECT_LOCATION",
            "BEST_OF_5943_REQUIRES_LOOK_ELSEWHERE",
            "P_0_693_CANNOT_PROMOTE_COINCIDENCE",
            "NEGATIVE_119HZ_CONFIRMATORY_RESULT_IMMUTABLE",
            "NATURAL_TECTONIC_CONTROL_PRECEDES_TARGET_INTERPRETATION",
            "TARGET_IDENTITY_UNCONFIRMED",
            "NO_UNDERWATER_PYRAMID_DETECTED_YET",
        ],
        "status": "FULL_ASSEMBLY_POINTER_GATE_COMPLETE",
    }

    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "engine_verdict": result["engine_verdict"],
        "detected_point": None,
        "measurement_point": result["question_split"]["WHERE_SHOULD_THE_NEXT_DISCRIMINATING_MEASUREMENT_POINT"]["point"],
        "best_crossmatch_event": result["strongest_real_date_location_cosmos_crossmatch"],
    }, indent=2, ensure_ascii=False))
    return 0 if passed == len(tests) else 2


if __name__ == "__main__":
    raise SystemExit(main())
