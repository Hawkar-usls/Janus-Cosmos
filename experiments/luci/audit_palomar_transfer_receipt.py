#!/usr/bin/env python3
"""Corrective fail-closed audit for LUCI x Palomar morphology-transfer receipts.

This audit was added *after* the first real LUCI transfer receipt exposed
obviously degenerate pixel-shape estimates. It is therefore a corrective
implementation audit, not a preregistered astrophysical gate, and it may not
be used to rescue or strengthen the original morphology PASS.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path


def finite(x) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def audit(receipt: dict) -> dict:
    sources = receipt.get("morphology_transfer", {}).get("sources", [])
    frames = receipt.get("morphology_transfer", {}).get("frames", [])
    by_frame = defaultdict(list)
    for row in sources:
        by_frame[str(row.get("frame_id", ""))].append(row)

    frame_audits = []
    for frame in frames:
        frame_id = str(frame.get("frame_id", ""))
        rows = by_frame.get(frame_id, [])
        n = len(rows)
        undersampled = 0
        single_pixel = 0
        circularity_gt1 = 0
        nonfinite_required = 0
        for r in rows:
            minor = r.get("moment_fwhm_minor_px")
            area = r.get("threshold_area_px")
            circ = r.get("circularity")
            required = [r.get("moment_fwhm_px"), minor, r.get("moment_fwhm_major_px"), r.get("elongation")]
            if any(not finite(x) for x in required):
                nonfinite_required += 1
            if finite(minor) and float(minor) < 1.0:
                # A sub-pixel minor-axis FWHM is not spatially resolved by the detector;
                # treating it as measured morphology is an implementation failure.
                undersampled += 1
            if finite(area) and float(area) <= 1.0:
                # One-pixel support cannot constrain a 2-D shape tensor.
                single_pixel += 1
            if finite(circ) and float(circ) > 1.000001:
                # For 4*pi*A/P^2, values >1 violate the continuous isoperimetric bound.
                # The pilot's pixel-perimeter estimator therefore cannot be admitted as
                # a calibrated circularity measurement when this occurs.
                circularity_gt1 += 1

        def frac(k: int) -> float:
            return k / n if n else 0.0

        invalid = bool(
            n == 0
            or frac(undersampled) > 0.10
            or frac(single_pixel) > 0.10
            or circularity_gt1 > 0
            or nonfinite_required > 0
        )
        frame_audits.append({
            "frame_id": frame_id,
            "source_count": n,
            "undersampled_minor_fwhm_lt_1px_count": undersampled,
            "undersampled_minor_fwhm_lt_1px_fraction": frac(undersampled),
            "single_pixel_shape_support_count": single_pixel,
            "single_pixel_shape_support_fraction": frac(single_pixel),
            "circularity_gt_one_count": circularity_gt1,
            "nonfinite_required_shape_count": nonfinite_required,
            "implementation_valid": not invalid,
        })

    invalid_frames = [x for x in frame_audits if not x["implementation_valid"]]
    return {
        "schema": "janus.cosmos.luci_palomar_transfer.corrective_audit.v1",
        "audit_type": "POST_OUTCOME_CORRECTIVE_IMPLEMENTATION_AUDIT",
        "status": "IMPLEMENTATION_INVALID_MORPHOLOGY_MEASUREMENT" if invalid_frames else "IMPLEMENTATION_SANITY_PASS",
        "original_morphology_outcome": receipt.get("morphology_transfer", {}).get("outcome"),
        "original_frame_pass_fraction": receipt.get("morphology_transfer", {}).get("frame_gate_pass_fraction"),
        "frame_count": len(frame_audits),
        "invalid_frame_count": len(invalid_frames),
        "frame_audits": frame_audits,
        "direct_sky_overlap_preflight_preserved": receipt.get("direct_sky_overlap_preflight", {}),
        "decision": (
            "The original morphology-transfer PASS is void for scientific interpretation because the source-shape "
            "implementation produced unresolved/single-pixel and/or mathematically invalid shape estimates. "
            "The independent archive sky-overlap preflight is not changed by this morphology audit."
        ),
        "claim_ceiling": "CORRECTIVE_IMPLEMENTATION_AUDIT_ONLY__NO_LUCI_PALOMAR_MORPHOLOGY_CLAIM__NO_ANOMALY_CONFIRMATION",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("receipt")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    receipt = json.loads(Path(args.receipt).read_text(encoding="utf-8"))
    result = audit(receipt)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "invalid_frame_count": result["invalid_frame_count"],
        "frame_count": result["frame_count"],
        "palomar_sources_with_archive_overlap": result["direct_sky_overlap_preflight_preserved"].get("unique_palomar_sources_with_overlap", 0),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
