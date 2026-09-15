#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

PARENT_TRIAL_RESULTS_SHA256 = "25763a6ad7db994717e0e698cc2fa1612672baf63e3f8d7b0aa61be6ace45162"
TARGET_IDS = ("TS-T3B-TESS-15", "TS-T3B-TESS-39")
Z_MIN = 8.0
DELTA_MIN = 4.0
CLAIM = "TESS_RESIDUAL_SENSITIVITY_CHARACTERIZATION_ONLY__NO_PARENT_NEGATIVE_UPGRADE__NO_NEW_PIXEL_READ__NO_FTL__NO_RETROCAUSALITY__NO_TACHYON_IDENTITY"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def required_injected_snr(raw: dict) -> float:
    vals = [
        0.0,
        Z_MIN - float(raw["contrast_z"]),
        DELTA_MIN - float(raw["b_minus_a_sigma"]),
        DELTA_MIN - float(raw["b_minus_c_sigma"]),
    ]
    if not all(math.isfinite(v) for v in vals):
        raise ValueError("nonfinite raw candidate metrics")
    return float(max(vals))


def predicted_candidate(raw: dict, injected_snr: float) -> bool:
    x = float(injected_snr)
    return bool(
        float(raw["contrast_z"]) + x >= Z_MIN
        and float(raw["b_minus_a_sigma"]) + x >= DELTA_MIN
        and float(raw["b_minus_c_sigma"]) + x >= DELTA_MIN
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trial-results", required=True)
    ap.add_argument("--output-dir", default="results/tachyon_star_t3c_tess")
    args = ap.parse_args()

    src = Path(args.trial_results)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if sha256_file(src) != PARENT_TRIAL_RESULTS_SHA256:
        raise RuntimeError("frozen T3B trial-results SHA mismatch")
    rows = json.loads(src.read_text(encoding="utf-8"))
    by_id = {r["trial_id"]: r for r in rows}
    if set(TARGET_IDS) - set(by_id):
        raise RuntimeError("frozen unresolved T3B target missing")

    receipts = []
    inconsistent = []
    for tid in TARGET_IDS:
        row = by_id[tid]
        if row.get("classification") != "UNRESOLVED_TRIAL":
            raise RuntimeError(f"parent target {tid} is no longer unresolved")
        raw = row.get("raw") or {}
        req = required_injected_snr(raw)
        ceil_req = int(math.ceil(req - 1e-12))
        frozen_inj = {float(q["injected_snr"]): bool(q["candidate_recovered"]) for q in row.get("injections", [])}
        predicted_10 = predicted_candidate(raw, 10.0)
        predicted_12 = predicted_candidate(raw, 12.0)
        consistent = (
            frozen_inj.get(10.0) is False
            and frozen_inj.get(12.0) is True
            and predicted_10 is False
            and predicted_12 is True
            and 10.0 < req <= 12.0
        )
        if not consistent:
            inconsistent.append(tid)
        receipts.append({
            "trial_id": tid,
            "src_id": row.get("src_id"),
            "parent_classification": row.get("classification"),
            "raw_contrast_z": float(raw["contrast_z"]),
            "raw_b_minus_a_sigma": float(raw["b_minus_a_sigma"]),
            "raw_b_minus_c_sigma": float(raw["b_minus_c_sigma"]),
            "required_injected_snr_exact": req,
            "smallest_integer_injected_snr_meeting_candidate_gate": ceil_req,
            "frozen_10sigma_recovered": frozen_inj.get(10.0),
            "frozen_12sigma_recovered": frozen_inj.get(12.0),
            "analytic_10sigma_candidate": predicted_10,
            "analytic_12sigma_candidate": predicted_12,
            "analytic_consistent_with_frozen_bracket": consistent,
            "parent_reclassification": "NONE",
        })

    status = "PASS_RESIDUAL_SENSITIVITY_CHARACTERIZED" if not inconsistent else "BLOCKED_ANALYTIC_INCONSISTENCY"
    rec = {
        "schema": "janus.cosmos.tachyon_star.t3c.tess_residual_sensitivity_bound.receipt.v1",
        "experiment_id": "JANUS-TACHYON-STAR-T3C-TESS-RESIDUAL-SENSITIVITY-BOUND",
        "status": status,
        "parent_trial_results_sha256": PARENT_TRIAL_RESULTS_SHA256,
        "parent_qualified_null_trials_immutable": 40,
        "parent_candidate_trials": 0,
        "targeted_unresolved_trials": 2,
        "new_pixel_read": False,
        "parent_trials_reclassified": 0,
        "targets": receipts,
        "inconsistent_trial_ids": inconsistent,
        "claim_ceiling": CLAIM,
    }
    (out / "receipt.json").write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(rec, indent=2, sort_keys=True))
    return 0 if status == "PASS_RESIDUAL_SENSITIVITY_CHARACTERIZED" else 3


if __name__ == "__main__":
    raise SystemExit(main())
