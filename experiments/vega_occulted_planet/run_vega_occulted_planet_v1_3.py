#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
PARENT = OUT / "vega_occulted_planet_receipt_v1_2.json"
MIRI_PROV = OUT / "vega_miri_final_provenance.json"
MIRI_GATE = OUT / "vega_miri_resolvability_gate.json"
MIRI_PROFILES = OUT / "vega_miri_radial_profiles.json"

RECEIPT = OUT / "vega_occulted_planet_receipt_v1_3.json"
TOPA = OUT / "topa_queue_v1_3.json"
SPIDER = OUT / "spider_queue_v1_3.json"


def canonical_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def q(value: Any, digits: int) -> float | None:
    if value is None:
        return None
    return round(float(value), digits)


def build() -> dict:
    parent = load(PARENT)
    prov = load(MIRI_PROV)
    gate = load(MIRI_GATE)
    profiles = load(MIRI_PROFILES)

    products = {
        p["filter"]: {
            "sha256": p["sha256"],
            "bytes": p["bytes"],
            "filename": p["filename"],
        }
        for p in prov["products"]
    }

    receipt = {
        "schema": "janus.cosmos.vega_occulted_planet.receipt.v1.3.1",
        "version": "1.3.1",
        "target": "Vega",
        "parent_v1_2_freeze_sha256": parent["freeze_sha256"],
        "freeze_numeric_policy": {
            "pixel_scale_arcsec_decimals": 6,
            "diffraction_arcsec_decimals": 6,
            "finite_support_fraction_decimals": 3,
            "raw_values_location": "vega_miri_resolvability_gate.json",
        },
        "H1A_2P43D_RV_CANDIDATE": parent["H1A_2P43D_RV_CANDIDATE"],
        "H1B_3_TO_5_AU_NEPTUNE_SHEPHERD": {
            **parent["H1B_3_TO_5_AU_NEPTUNE_SHEPHERD"],
            "real_miri_products_ingested": True,
            "miri_products": products,
            "miri_resolvability_gate": gate["gate"],
            "miri_resolvability_status": gate["interpretation"]["status"],
            "status": "ANALYTICALLY_FEASIBLE_DIRECT_3_TO_5_AU_IMAGE_RECOVERY_RESOLUTION_LIMITED",
            "corrected_next_test": "Forward-model the observable radial surface-brightness profile and SED under no-planet PR-drag and shepherd families; do not claim localized 3-5 AU image recovery when the edge is at/below MIRI diffraction scales.",
            "claim": "The authors' final MIRI products are now ingested, but their angular resolution does not turn the 3-5 AU shepherd hypothesis into a directly imaged planet test.",
        },
        "real_data_gate": {
            "miri_provenance_status": prov["status"],
            "target_edge_au": [q(x, 3) for x in gate["target_warm_disk_inner_edge_au"]],
            "target_edge_arcsec": [q(x, 6) for x in gate["target_warm_disk_inner_edge_arcsec"]],
            "filters": [
                {
                    "filter": f["filter"],
                    "pixel_scale_arcsec": q(f["pixel_scale_arcsec"], 6),
                    "diffraction_fwhm_proxy_arcsec": q(
                        f["diffraction"]["fwhm_proxy_1p03_lambda_over_D_arcsec"], 6
                    ),
                    "diffraction_rayleigh_arcsec": q(
                        f["diffraction"]["rayleigh_1p22_lambda_over_D_arcsec"], 6
                    ),
                    "finite_support_fraction_3_to_5_au": q(
                        f["finite_support_fraction_3_to_5_au"], 3
                    ),
                    "resolvability": f["resolvability"],
                }
                for f in gate["filters"]
            ],
            "radial_profiles_schema": profiles["schema"],
        },
        "topa_decision": {
            "rank_1": "Independent/new RV replication remains the decisive H1A gate because the Hurt-like GP absorbs the 2.43 d preference in the current TRES series.",
            "rank_2": "For H1B, pivot from localized 3-5 AU image injection to a joint radial-profile + SED forward model using the ingested final MIRI reductions, with no-planet PR-drag controls.",
            "rank_3": "Use any future shorter-wavelength/higher-angular-resolution or interferometric constraints to attack the 3-5 AU edge directly.",
        },
        "verdict": "H1A_WEAKENED_H1B_REAL_MIRI_INGESTED_BUT_DIRECT_EDGE_RECOVERY_RESOLUTION_LIMITED",
        "claim_ceiling": "NO_PLANET_DETECTION_H1A_ACTIVITY_COMPATIBLE_H1B_PROFILE_SED_TEST_PENDING",
        "firewall": [
            "The 2.43 d RV candidate is not confirmed and is weakened by the quasi-periodic activity model.",
            "The 3-5 AU shepherd remains an analytic/dynamical hypothesis, not a directly imaged object.",
            "Real MIRI FITS ingestion improves provenance but does not overcome diffraction or coronagraphic limits.",
            "Radial profile bins smaller than the PSF are descriptive oversampling and must not be counted as independent resolution elements.",
        ],
    }
    receipt["freeze_sha256"] = canonical_hash(receipt)

    topa = {
        "schema": "janus.cosmos.topa_queue.vega.v1.3.1",
        "receipt_hash": receipt["freeze_sha256"],
        "ranked_tests": [
            {
                "rank": 1,
                "id": "INDEPENDENT_RV_ACTIVITY_AWARE_REPLICATION",
                "status": "WAITING_FOR_NEW_OR_INDEPENDENT_RV",
            },
            {
                "rank": 2,
                "id": "MIRI_RADIAL_PROFILE_SED_FORWARD_MODEL",
                "status": "READY_WITH_REAL_MIRI_IMAGES_SED_PHOTOMETRY_MODEL_PENDING",
                "controls": [
                    "NO_PLANET_PR_DRAG_ONLY",
                    "LOW_MASS_SHEPHERD_FAMILIES",
                    "PROFILE_PERTURBATION_PLACEBOS",
                ],
                "success": "A frozen shepherd family improves held-out radial-profile/SED likelihood over no-planet and placebo controls without violating RV/direct-imaging constraints.",
                "failure": "No-planet PR-drag or placebo models match as well or better, or required masses conflict with independent constraints.",
            },
        ],
    }
    spider = {
        "schema": "janus.cosmos.spider_queue.vega.v1.3.1",
        "receipt_hash": receipt["freeze_sha256"],
        "requests": [
            {
                "priority": 1,
                "id": "INDEPENDENT_OR_EXTENDED_VEGA_RV",
                "status": "QUEUED",
            },
            {
                "priority": 2,
                "id": "VEGA_INNER_DISK_SED_PHOTOMETRY_MACHINE_READABLE",
                "status": "QUEUED",
                "purpose": "Combine unresolved 3-5 AU inner-edge information with the real MIRI radial profiles instead of pretending direct spatial recovery.",
            },
        ],
        "rule": "Resolution limits are part of the evidence model; they may not be bypassed by over-resolving radial bins.",
    }

    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    TOPA.write_text(json.dumps(topa, indent=2) + "\n", encoding="utf-8")
    SPIDER.write_text(json.dumps(spider, indent=2) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    for path in (PARENT, MIRI_PROV, MIRI_GATE, MIRI_PROFILES):
        if not path.exists():
            raise SystemExit(f"missing prerequisite: {path.name}")
    r1 = build()
    r2 = build()
    assert r1["freeze_sha256"] == r2["freeze_sha256"]
    assert r1["claim_ceiling"] == "NO_PLANET_DETECTION_H1A_ACTIVITY_COMPATIBLE_H1B_PROFILE_SED_TEST_PENDING"
    print("VEGA EVIDENCE SYNTHESIS v1.3.1 PASS")
    print("verdict =", r1["verdict"])
    print("H1A_status =", r1["H1A_2P43D_RV_CANDIDATE"]["status"])
    print("H1B_status =", r1["H1B_3_TO_5_AU_NEPTUNE_SHEPHERD"]["status"])
    print("MIRI_gate =", r1["H1B_3_TO_5_AU_NEPTUNE_SHEPHERD"]["miri_resolvability_gate"])
    print("freeze_sha256 =", r1["freeze_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
