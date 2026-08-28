#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
PARENT = OUT / "vega_occulted_planet_receipt_v1_1.json"
QP = OUT / "vega_tres_rv_quasiperiodic_activity_audit.json"
SHEPHERD = OUT / "vega_3to5au_shepherd_analytic_grid.json"
MANIFEST = OUT / "vega_3to5au_injection_recovery_manifest.json"

RECEIPT = OUT / "vega_occulted_planet_receipt_v1_2.json"
TOPA = OUT / "topa_queue_v1_2.json"
SPIDER = OUT / "spider_queue_v1_2.json"


def canonical_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def classify_h1a(qp_status: str) -> str:
    if qp_status == "2P43_SURVIVES_HURT_LIKE_QP_ACTIVITY_WITH_SPECIFICITY":
        return "OPEN_UNCONFIRMED_ACTIVITY_AWARE_SURVIVOR"
    if qp_status == "HURT_LIKE_QP_ACTIVITY_ABSORBS_2P43":
        return "WEAKENED_ACTIVITY_COMPATIBLE"
    if qp_status == "WEAK_2P43_PREFERENCE_AFTER_HURT_LIKE_QP_ACTIVITY":
        return "OPEN_WEAK_ACTIVITY_AWARE_PREFERENCE"
    return "OPEN_UNCONFIRMED_MIXED_ACTIVITY_SPECIFICITY"


def build() -> dict:
    parent = load(PARENT)
    qp = load(QP)
    shepherd = load(SHEPHERD)
    manifest = load(MANIFEST)

    qps = qp["specificity"]
    h1a_status = classify_h1a(qps["status"])
    h1b_status = (
        "ANALYTICALLY_FEASIBLE_REAL_IMAGE_RECOVERY_PENDING"
        if shepherd["summary"]["analytic_feasible_neptune_like_pole_on_rows"] > 0
        else "NO_FROZEN_GRID_EDGE_MATCH_REAL_IMAGE_RECOVERY_STILL_REQUIRED"
    )

    # Optimizer outputs can drift at the 1e-6 level (and the derived coherence
    # proxy by a small fraction of a day) across scipy/cpu runs while leaving
    # the scientific classification unchanged. Raw values remain in the QP
    # report; only scientifically meaningful quantized summaries enter the
    # evidence receipt hash so freeze_sha256 is reproducible.
    qp_summary = {
        "delta_bic_candidate_vs_activity": round(float(qps["delta_bic_candidate_vs_activity"]), 3),
        "delta_bic_control_vs_activity": round(float(qps["delta_bic_control_vs_activity"]), 3),
        "delta_bic_candidate_over_control": round(float(qps["delta_bic_candidate_over_control"]), 3),
        "candidate_amplitude_m_s": round(
            float(qp["models"]["quasiperiodic_activity_plus_2p43"]["global_period_amplitudes_m_s"]["2.43"]),
            2,
        ),
        "activity_coherence_proxy_days": int(
            round(float(qp["models"]["quasiperiodic_activity_only"]["activity_hyperparameters"]["coherence_timescale_proxy_days"]))
        ),
    }

    receipt = {
        "schema": "janus.cosmos.vega_occulted_planet.receipt.v1.2.1",
        "version": "1.2.1",
        "target": "Vega",
        "parent_v1_1_freeze_sha256": parent["freeze_sha256"],
        "freeze_numeric_policy": {
            "purpose": "Prevent insignificant optimizer/platform microdrift from changing the evidence hash.",
            "delta_bic_decimals": 3,
            "rv_amplitude_decimals": 2,
            "coherence_proxy_days_resolution": 1,
            "raw_values_location": "vega_tres_rv_quasiperiodic_activity_audit.json",
        },
        "H1A_2P43D_RV_CANDIDATE": {
            "status": h1a_status,
            "quasiperiodic_specificity_status": qps["status"],
            **qp_summary,
            "model_relation_to_hurt_2021": qp["model"]["relationship_to_hurt_2021"],
            "claim": "The 2.43-day candidate remains unconfirmed regardless of model preference in this audit.",
        },
        "H1B_3_TO_5_AU_NEPTUNE_SHEPHERD": {
            "status": h1b_status,
            "analytic_grid_freeze_sha256": shepherd["freeze_sha256"],
            "injection_recovery_manifest_freeze_sha256": manifest["freeze_sha256"],
            "analytic_feasible_neptune_like_pole_on_rows": shepherd["summary"]["analytic_feasible_neptune_like_pole_on_rows"],
            "rv_k_range_m_s_among_feasible": [
                round(float(shepherd["summary"]["minimum_rv_k_m_s_among_feasible"]), 6),
                round(float(shepherd["summary"]["maximum_rv_k_m_s_among_feasible"]), 6),
            ],
            "claim": "Analytic chaotic-zone/RV feasibility is not image-domain recovery and is not evidence of a planet.",
        },
        "topa_decision": {
            "rank_1": "Acquire independent/new Vega RV epochs and repeat the activity-aware 2.43 d test with an independent instrument or extended baseline.",
            "rank_2": "Obtain provenance-backed MIRI warm-disk products and execute the frozen 3-5 au image-domain injection/recovery manifest.",
            "rank_3": "Only then combine H1A/H1B constraints in a joint architecture model; do not assume they are the same planet.",
        },
        "verdict": "TWO_OPEN_INNER_PLANET_SUBHYPOTHESES_ACTIVITY_AUDITED",
        "claim_ceiling": "H1A_UNCONFIRMED_AND_H1B_ANALYTIC_FEASIBILITY_ONLY",
        "firewall": [
            "No current output is a planet detection.",
            "H1A and H1B remain separate hypotheses.",
            "The Hurt-like GP is a scalable approximation, not an exact reimplementation.",
            "The 3-5 au branch is preregistered for real image-domain injection/recovery but has not yet run on a MIRI image.",
        ],
    }
    receipt["freeze_sha256"] = canonical_hash(receipt)

    topa = {
        "schema": "janus.cosmos.topa_queue.vega.v1.2.1",
        "receipt_hash": receipt["freeze_sha256"],
        "ranked_tests": [
            {
                "rank": 1,
                "id": "INDEPENDENT_RV_2P43_ACTIVITY_AWARE_REPLICATION",
                "status": "WAITING_FOR_NEW_OR_INDEPENDENT_RV",
                "success": "A coherent 2.43 d component survives an activity-aware model in independent/new RV data.",
                "failure": "The 2.43 d preference disappears, loses coherence, or tracks stellar activity diagnostics.",
            },
            {
                "rank": 2,
                "id": "MIRI_3_TO_5_AU_REAL_IMAGE_INJECTION_RECOVERY",
                "status": manifest["status"],
                "manifest_hash": manifest["freeze_sha256"],
                "success": manifest["frozen_test"]["success_rule"],
                "failure": manifest["frozen_test"]["failure_rule"],
            },
        ],
    }

    spider = {
        "schema": "janus.cosmos.spider_queue.vega.v1.2.1",
        "receipt_hash": receipt["freeze_sha256"],
        "requests": [
            {
                "priority": 1,
                "id": "INDEPENDENT_OR_EXTENDED_VEGA_RV",
                "status": "QUEUED",
                "purpose": "Independent confirmation/refutation of the 2.43 d candidate under activity-aware modeling.",
            },
            {
                "priority": 2,
                "id": "JWST_MIRI_WARM_DISK_SCIENCE_PRODUCT",
                "status": "QUEUED",
                "purpose": "Execute the frozen 3-5 au disk injection/recovery manifest.",
                "required_manifest_hash": manifest["freeze_sha256"],
            },
        ],
        "rule": "Missing data remain missing data; no queue state may be interpreted as positive evidence.",
    }

    OUT.mkdir(exist_ok=True)
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    TOPA.write_text(json.dumps(topa, indent=2) + "\n", encoding="utf-8")
    SPIDER.write_text(json.dumps(spider, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    for path in (PARENT, QP, SHEPHERD, MANIFEST):
        if not path.exists():
            raise SystemExit(f"missing prerequisite: {path.name}")
    r1 = build()
    r2 = build()
    assert r1["freeze_sha256"] == r2["freeze_sha256"]
    assert r1["claim_ceiling"] == "H1A_UNCONFIRMED_AND_H1B_ANALYTIC_FEASIBILITY_ONLY"
    assert r1["verdict"] == "TWO_OPEN_INNER_PLANET_SUBHYPOTHESES_ACTIVITY_AUDITED"
    print("VEGA EVIDENCE SYNTHESIS v1.2.1 PASS")
    print("H1A_status =", r1["H1A_2P43D_RV_CANDIDATE"]["status"])
    print("H1B_status =", r1["H1B_3_TO_5_AU_NEPTUNE_SHEPHERD"]["status"])
    print("quantized_QP =", r1["H1A_2P43D_RV_CANDIDATE"])
    print("freeze_sha256 =", r1["freeze_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
