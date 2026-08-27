#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
BASE_RECEIPT = OUT / "vega_occulted_planet_receipt.json"
EVIDENCE = ROOT / "evidence_snapshot_v1_1.json"
RV_AUDIT = OUT / "vega_tres_rv_audit.json"
RECEIPT = OUT / "vega_occulted_planet_receipt_v1_1.json"
TOPA = OUT / "topa_queue_v1_1.json"
SPIDER = OUT / "spider_queue_v1_1.json"


def canonical_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def aligned_true_mass(min_mass_earth: float, inclination_deg: float) -> float:
    return min_mass_earth / math.sin(math.radians(inclination_deg))


def build() -> dict[str, Any]:
    base = load_json(BASE_RECEIPT)
    evidence = load_json(EVIDENCE)
    rv = load_json(RV_AUDIT) if RV_AUDIT.exists() else None

    min_mass = 20.0
    aligned_mass_earth = aligned_true_mass(min_mass, 6.5)
    aligned_mass_jup = aligned_mass_earth / 317.83

    rv_live = {
        "available": rv is not None,
        "status": rv["interpretation"]["status"] if rv else "NOT_FETCHED",
        "row_count": rv.get("row_count") if rv else None,
        "delta_bic_fixed_2p43": rv.get("diagnostics", {}).get("delta_bic_rotation_to_rotation_plus_2p43") if rv else None,
        "candidate_amplitude_full_m_s": rv.get("diagnostics", {}).get("candidate_amplitude_full_m_s") if rv else None,
        "candidate_amplitude_early_m_s": rv.get("diagnostics", {}).get("candidate_amplitude_early_m_s") if rv else None,
        "candidate_amplitude_late_m_s": rv.get("diagnostics", {}).get("candidate_amplitude_late_m_s") if rv else None,
        "claim_boundary": "Fixed-period re-fit is a diagnostic reproduction only; it is not an independent confirmation and does not supersede Hurt et al. activity-aware analysis."
    }

    receipt: dict[str, Any] = {
        "schema": "janus.cosmos.vega_occulted_planet.receipt.v1.1",
        "version": "1.1.0",
        "target": "Vega",
        "parent_v1_freeze_sha256": base["freeze_sha256"],
        "evidence_snapshot_sha256": canonical_hash(evidence),
        "parent_grid": base["grid"],
        "published_evidence_partition": {
            "H1A_2P43D_RV_CANDIDATE": {
                "status": "OPEN_UNCONFIRMED",
                "period_days": 2.43,
                "rv_semiamplitude_m_s_published": 6.0,
                "minimum_mass_earth_published_approx": min_mass,
                "true_mass_if_i_6p5deg_earth": aligned_mass_earth,
                "true_mass_if_i_6p5deg_jupiter": aligned_mass_jup,
                "note": "If orbital inclination is near Vega's spin inclination, the same m*sin(i) implies a much larger true mass. Alignment is an assumption, not a measurement of the candidate orbit."
            },
            "H1B_3_TO_5_AU_NEPTUNE_SHEPHERD": {
                "status": "NOT_EXCLUDED_BY_CURRENT_PUBLISHED_SUMMARY",
                "reason": [
                    "MIRI explicitly allows a possible modest/Neptune-size shepherd near the 3-5 au warm-disk inner edge.",
                    "JWST/NIRCam encoded direct-imaging anchors begin near 7.7 au.",
                    "Published RV sensitivity is strongly orientation-dependent for nearly pole-on Vega; low/moderate masses at several au are not hard-excluded by the encoded summary."
                ]
            }
        },
        "live_rv_audit": rv_live,
        "astrometry_gate": {
            "status": "OPEN_NO_HARD_GAIA_CONSTRAINT",
            "reason": "Hurt et al. state that Vega is far beyond the bright limit for standard Gaia processing and do not provide a Vega-specific reliable precision. Missing/uncertain Gaia astrometry is not evidence for a companion."
        },
        "topa_decision": {
            "rank_1": "Acquire/extend activity-aware radial velocities resolving 0.676 d stellar rotation and 2.43 d candidate simultaneously.",
            "rank_2": "Build a 3-5 au disk-dynamical injection/recovery grid for Neptune-to-sub-Saturn masses.",
            "rank_3": "Seek a provenance-backed Vega-specific astrometric acceleration/proper-motion constraint; otherwise leave astrometry open.",
            "do_not_do": "Do not merge H1A and H1B into one planet hypothesis and do not treat a fixed-period re-fit as confirmation."
        },
        "verdict": "TWO_OPEN_INNER_PLANET_SUBHYPOTHESES",
        "claim_ceiling": "H1A_UNCONFIRMED_AND_H1B_NOT_EXCLUDED_ONLY"
    }
    receipt["freeze_sha256"] = canonical_hash(receipt)

    topa = {
        "schema": "janus.cosmos.topa_queue.vega.v1.1",
        "receipt_hash": receipt["freeze_sha256"],
        "ranked_tests": [
            {
                "rank": 1,
                "id": "RV_ACTIVITY_AWARE_2P43_RETEST",
                "discriminates": ["H1A_2P43D_RV_CANDIDATE", "STELLAR_ACTIVITY_ALIAS"],
                "success": "Independent/extended RV series retains coherent 2.43 d signal after activity-aware modeling.",
                "failure": "Signal loses coherence/significance with improved activity model or new epochs."
            },
            {
                "rank": 2,
                "id": "INNER_EDGE_3_TO_5_AU_DYNAMICAL_INJECTION_RECOVERY",
                "discriminates": ["H1B_3_TO_5_AU_NEPTUNE_SHEPHERD", "DRAG_DOMINATED_NO_PLANET_REQUIRED"],
                "success": "A frozen low-mass planet family reproduces the inner edge/gap while preserving the observed smoothness better than no-planet controls.",
                "failure": "No-planet dust dynamics match as well or better, or required masses conflict with RV/imaging limits."
            },
            {
                "rank": 3,
                "id": "VEGA_SPECIFIC_ASTROMETRY_PROVENANCE_GATE",
                "discriminates": ["MASSIVE_INNER_COMPANION", "LOW_MASS_OR_NO_COMPANION"],
                "success": "A documented acceleration/proper-motion anomaly with covariance and bright-star calibration exists.",
                "failure": "No reliable Vega-specific astrometric constraint can be established; gate remains open rather than negative."
            }
        ]
    }
    spider = {
        "schema": "janus.cosmos.spider_queue.vega.v1.1",
        "receipt_hash": receipt["freeze_sha256"],
        "requests": [
            {
                "priority": 1,
                "id": "TRES_VIZIER_MACHINE_READABLE",
                "status": "INGESTED" if rv else "QUEUED",
                "source": "VizieR J/AJ/161/157/table2",
                "expected_rows": 1524,
                "purpose": "Reproducible RV ingest and diagnostic fixed-period audit."
            },
            {
                "priority": 1,
                "id": "NEW_OR_INDEPENDENT_VEGA_RV",
                "status": "QUEUED",
                "purpose": "Resolve stellar rotation/activity from the 2.43 d candidate with new epochs or independent instrumentation."
            },
            {
                "priority": 2,
                "id": "MIRI_NIRCAM_ALMA_INNER_EDGE_PRODUCTS",
                "status": "QUEUED",
                "purpose": "Construct 3-5 au low-mass dynamical injection/recovery against the smooth disk."
            },
            {
                "priority": 3,
                "id": "VEGA_BRIGHT_STAR_ASTROMETRY",
                "status": "QUEUED_WITH_CAVEAT",
                "purpose": "Only accept a hard constraint with explicit Vega-specific bright-star calibration/covariance."
            }
        ],
        "rule": "Missing data are not evidence; every accepted product must carry provenance."
    }

    OUT.mkdir(exist_ok=True)
    RECEIPT.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    TOPA.write_text(json.dumps(topa, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    SPIDER.write_text(json.dumps(spider, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return receipt


def self_test() -> None:
    r1 = build()
    r2 = build()
    assert r1["freeze_sha256"] == r2["freeze_sha256"]
    assert r1["published_evidence_partition"]["H1A_2P43D_RV_CANDIDATE"]["status"] == "OPEN_UNCONFIRMED"
    assert r1["published_evidence_partition"]["H1B_3_TO_5_AU_NEPTUNE_SHEPHERD"]["status"] == "NOT_EXCLUDED_BY_CURRENT_PUBLISHED_SUMMARY"
    assert 0.5 < r1["published_evidence_partition"]["H1A_2P43D_RV_CANDIDATE"]["true_mass_if_i_6p5deg_jupiter"] < 0.7
    print("VEGA EVIDENCE SYNTHESIS v1.1 PASS")
    print("verdict =", r1["verdict"])
    print("freeze_sha256 =", r1["freeze_sha256"])


if __name__ == "__main__":
    self_test()
