#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
MANIFEST_PATH = ROOT / "manifest.json"
OUT = ROOT / "output"


def logspace(lo: float, hi: float, steps: int) -> list[float]:
    if lo <= 0 or hi <= 0 or steps < 2:
        raise ValueError("logspace requires positive bounds and >=2 steps")
    a = math.log10(lo)
    b = math.log10(hi)
    return [10 ** (a + (b - a) * i / (steps - 1)) for i in range(steps)]


def interp_logx(x: float, anchors: list[dict[str, float]]) -> float | None:
    anchors = sorted(anchors, key=lambda r: r["separation_au"])
    if x < anchors[0]["separation_au"]:
        return None
    if x >= anchors[-1]["separation_au"]:
        return anchors[-1]["max_jupiter_masses"]
    lx = math.log(x)
    for left, right in zip(anchors[:-1], anchors[1:]):
        x0 = left["separation_au"]
        x1 = right["separation_au"]
        if x0 <= x <= x1:
            t = (lx - math.log(x0)) / (math.log(x1) - math.log(x0))
            return left["max_jupiter_masses"] + t * (
                right["max_jupiter_masses"] - left["max_jupiter_masses"]
            )
    raise AssertionError("interpolation fell through")


def classify_cell(a_au: float, mass_earth: float, m: dict[str, Any]) -> dict[str, Any]:
    rules = m["frozen_rules"]
    pubs = m["published_constraints"]
    mjup_earth = float(rules["jupiter_mass_in_earth_masses"])
    saturn_earth = float(rules["saturn_mass_in_earth_masses"])

    hard_reasons: list[str] = []
    soft_flags: list[str] = []

    # Published qualitative MIRI constraint encoded as a hard exclusion for this arm.
    if rules["reject_if_mass_ge_saturn_and_a_ge_10_au"]:
        if a_au >= 10.0 and mass_earth >= saturn_earth:
            hard_reasons.append("MIRI_SMOOTH_DISK_SATURN_PLUS_OUTSIDE_10_AU")

    # Model-dependent NIRCam direct-imaging upper limits.
    anchors = pubs["nircam"]["contrast_mass_limits"]
    limit_mjup = interp_logx(a_au, anchors)
    if limit_mjup is not None and mass_earth >= limit_mjup * mjup_earth:
        hard_reasons.append("NIRCAM_MODEL_DEPENDENT_DIRECT_IMAGING_LIMIT")

    # This is deliberately a soft flag because the <6 Earth-mass value is tied
    # to a specific shepherding interpretation, not a universal exclusion law.
    if 70.0 <= a_au <= 100.0 and mass_earth > float(
        pubs["miri"]["outer_belt_shepherd_example_max_earth_masses"]
    ):
        soft_flags.append("OUTER_BELT_SHEPHERD_INTERPRETATION_DISFAVORED")

    if 3.0 <= a_au <= 5.0 and 5.0 <= mass_earth <= 30.0:
        soft_flags.append("INNER_EDGE_NEPTUNE_CLASS_PRIORITY_REGION")

    if a_au < anchors[0]["separation_au"]:
        soft_flags.append("INSIDE_ENCODED_NIRCAM_1ARCSEC_ANCHOR")

    return {
        "admissible": not hard_reasons,
        "hard_reasons": hard_reasons,
        "soft_flags": soft_flags,
        "nircam_limit_jupiter_mass": limit_mjup,
    }


def freeze_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def build_topa_queue(receipt: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema": "janus.cosmos.topa_queue.vega.v1",
        "target": "Vega",
        "hypotheses": [
            {
                "id": "H0_DRAG_DOMINATED_NO_PLANET_REQUIRED",
                "priority": 1,
                "falsifiers": [
                    "repeatable non-axisymmetric disk perturbation tied to a Keplerian perturber",
                    "independent astrometric or radial-velocity orbit",
                ],
            },
            {
                "id": "H1_OCCULTED_UNRESOLVED_INNER_PLANET",
                "priority": 1,
                "falsifiers": [
                    "RV/astrometric exclusion covering the remaining inner mass-orbit envelope",
                    "disk-dynamical exclusion of the same envelope",
                ],
            },
            {
                "id": "H2_60_AU_GAP_PLANET",
                "priority": 2,
                "falsifiers": [
                    "multiwavelength gap explained by dust physics without a perturber",
                    "planet-mass constraints incompatible with maintaining disk smoothness",
                ],
            },
            {
                "id": "H3_PSF_REDUCTION_SYSTEMATIC",
                "priority": 1,
                "falsifiers": [
                    "feature replication across independent instruments and reduction pipelines"
                ],
            },
            {
                "id": "H4_ALTERNATIVE_DISK_DYNAMICS",
                "priority": 2,
                "falsifiers": [
                    "dynamical model comparison decisively preferring a planetary perturber"
                ],
            },
        ],
        "receipt_hash": receipt["freeze_sha256"],
        "rule": "TOPA ranks discriminating tests; it does not turn an admissible grid cell into evidence.",
    }


def build_spider_queue(receipt: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema": "janus.cosmos.spider_queue.vega.v1",
        "target": "Vega",
        "requests": [
            {
                "id": "JWST_NIRCAM_F444W",
                "kind": "direct_imaging",
                "request": "public calibrated coronagraphic products, contrast curves, multi-epoch source astrometry",
            },
            {
                "id": "JWST_MIRI_DISK",
                "kind": "disk_morphology",
                "request": "F1550C/F2300C/F2550W calibrated products, masks, PSF/reference products, radial profiles",
            },
            {
                "id": "HST_STIS_VEGA",
                "kind": "independent_coronagraphy",
                "request": "science products, reference PSF metadata, reduction provenance",
            },
            {
                "id": "ALMA_1P34MM",
                "kind": "outer_belt",
                "request": "image/cube products, synthesized beam, radial profile and calibration provenance",
            },
            {
                "id": "VEGA_RV",
                "kind": "radial_velocity",
                "request": "public RV epochs with timestamps, uncertainties, instrument labels and zero-point treatment",
            },
            {
                "id": "VEGA_ASTROMETRY",
                "kind": "astrometry",
                "request": "Hipparcos/Gaia proper-motion anomaly or acceleration constraints with covariance/provenance",
            },
        ],
        "receipt_hash": receipt["freeze_sha256"],
        "rule": "Spider records provenance and negative results; missing data are not evidence.",
    }


def run() -> dict[str, Any]:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    grid = manifest["grid"]
    axes = logspace(
        float(grid["semi_major_axis_au"]["min"]),
        float(grid["semi_major_axis_au"]["max"]),
        int(grid["semi_major_axis_au"]["steps"]),
    )
    masses = logspace(
        float(grid["mass_earth"]["min"]),
        float(grid["mass_earth"]["max"]),
        int(grid["mass_earth"]["steps"]),
    )

    total = 0
    admissible = 0
    inner_priority = 0
    inside_nircam_anchor = 0
    hard_reason_counts: dict[str, int] = {}
    soft_flag_counts: dict[str, int] = {}
    admissible_a: list[float] = []
    admissible_m: list[float] = []

    for a_au in axes:
        for mass_earth in masses:
            total += 1
            c = classify_cell(a_au, mass_earth, manifest)
            for reason in c["hard_reasons"]:
                hard_reason_counts[reason] = hard_reason_counts.get(reason, 0) + 1
            for flag in c["soft_flags"]:
                soft_flag_counts[flag] = soft_flag_counts.get(flag, 0) + 1
            if c["admissible"]:
                admissible += 1
                admissible_a.append(a_au)
                admissible_m.append(mass_earth)
                if "INNER_EDGE_NEPTUNE_CLASS_PRIORITY_REGION" in c["soft_flags"]:
                    inner_priority += 1
                if "INSIDE_ENCODED_NIRCAM_1ARCSEC_ANCHOR" in c["soft_flags"]:
                    inside_nircam_anchor += 1

    receipt: dict[str, Any] = {
        "schema": "janus.cosmos.vega_occulted_planet.receipt.v1",
        "version": manifest["version"],
        "target": manifest["target"],
        "hypothesis": manifest["hypothesis"],
        "manifest_sha256": freeze_hash(manifest),
        "grid": {
            "total_cells": total,
            "admissible_cells": admissible,
            "rejected_cells": total - admissible,
            "admissible_fraction": admissible / total,
            "inner_edge_neptune_priority_cells": inner_priority,
            "admissible_inside_nircam_1arcsec_anchor": inside_nircam_anchor,
            "admissible_envelope": {
                "a_au_min": min(admissible_a),
                "a_au_max": max(admissible_a),
                "mass_earth_min": min(admissible_m),
                "mass_earth_max": max(admissible_m),
            },
        },
        "hard_reason_counts": hard_reason_counts,
        "soft_flag_counts": soft_flag_counts,
        "verdict": "H1_NOT_EXCLUDED_BY_ENCODED_CONSTRAINTS" if admissible else "H1_EXCLUDED_BY_ENCODED_CONSTRAINTS",
        "claim_ceiling": manifest["hypothesis"]["claim_ceiling"],
    }
    receipt["freeze_sha256"] = freeze_hash(receipt)

    OUT.mkdir(exist_ok=True)
    (OUT / "vega_occulted_planet_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (OUT / "topa_queue.json").write_text(
        json.dumps(build_topa_queue(receipt), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (OUT / "spider_queue.json").write_text(
        json.dumps(build_spider_queue(receipt), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return receipt


def self_test() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    # Inner Neptune-class example should survive the encoded direct-imaging rule.
    c1 = classify_cell(4.0, 17.0, manifest)
    assert c1["admissible"] is True
    assert "INNER_EDGE_NEPTUNE_CLASS_PRIORITY_REGION" in c1["soft_flags"]

    # A Jupiter-mass body at 60 au is incompatible with the encoded smooth-disk rule.
    c2 = classify_cell(60.0, 317.83, manifest)
    assert c2["admissible"] is False
    assert "MIRI_SMOOTH_DISK_SATURN_PLUS_OUTSIDE_10_AU" in c2["hard_reasons"]

    # A 4-Mjup body at ~8 au should be rejected by encoded NIRCam limits.
    c3 = classify_cell(8.0, 4.0 * 317.83, manifest)
    assert c3["admissible"] is False
    assert "NIRCAM_MODEL_DEPENDENT_DIRECT_IMAGING_LIMIT" in c3["hard_reasons"]

    r1 = run()
    r2 = run()
    assert r1["freeze_sha256"] == r2["freeze_sha256"]
    assert r1["grid"]["admissible_cells"] > 0
    print("VEGA OCCULTED PLANET SELF-TEST PASS")
    print("verdict =", r1["verdict"])
    print("admissible_fraction =", r1["grid"]["admissible_fraction"])
    print("freeze_sha256 =", r1["freeze_sha256"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    receipt = run()
    print(json.dumps(receipt, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
