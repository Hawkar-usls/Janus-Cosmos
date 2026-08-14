#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

import demiurge_forge_v2 as forge
import download_sky_v2_1 as downloader
import janus_cosmos_core_v2 as core
import janus_cosmos_specificity_v2_1 as specificity
import janus_cosmos_v2_0 as parent_runtime


ROOT = Path(__file__).resolve().parent


def check(condition, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    print("JANUS COSMOS v2.1.0 DETECTOR SPECIFICITY REPAIR SELF-TEST", flush=True)
    protocol = json.loads((ROOT / "SPECIFICITY_PROTOCOL_v2_1.json").read_text(encoding="utf-8"))

    frozen = parent_runtime.verify_forge()
    check(frozen["genome_sha256"] == protocol["parent_detector"]["genome_sha256"], "parent genome binding")
    check(frozen["freeze_sha256"] == protocol["parent_detector"]["freeze_sha256"], "parent freeze binding")
    receipt = specificity.verify_protocol_sources(protocol)
    check(receipt["real_sky_control_count"] == 20 and receipt["hst_control_count"] == 20, "control cohort count")
    print("[PASS] frozen parent detector + v2.1 protocol/source identity")

    negative_cfg = protocol["negative_parent_certificate"]
    negative_path = ROOT / negative_cfg["path"]
    negative = json.loads(negative_path.read_text(encoding="utf-8"))
    evidence_root = negative_path.parent
    for name, expected_sha in negative["artifact_sha256"].items():
        check(core.sha256_file(evidence_root / name) == expected_sha, f"negative evidence artifact drift: {name}")
    check(negative["scientific_gate"]["specificity_gate"] == "FAIL", "negative certificate gate drift")
    check(not negative["scientific_gate"]["orion_candidate_admitted"], "negative Orion result drift")
    check(not negative["scientific_gate"]["ngc1425_candidate_admitted"], "negative NGC1425 result drift")
    print("[PASS] v2.0.2 negative specificity certificate + 10 artifact hashes")

    regenerated = specificity.regenerate_sky_controls(protocol)
    check(regenerated == protocol["real_sky_controls"]["centers"], "blind control regeneration")
    legacy = json.loads((ROOT / "SKY_MANIFEST_v2_0.json").read_text(encoding="utf-8"))["blind_controls"]["centers"]
    old_coords = {(round(row["ra_deg"], 9), round(row["dec_deg"], 9)) for row in legacy}
    new_coords = {(round(row["ra_deg"], 9), round(row["dec_deg"], 9)) for row in regenerated}
    check(old_coords.isdisjoint(new_coords), "observed v2.0.2 control leaked into fresh v2.1 admission cohort")
    for index, left in enumerate(regenerated):
        for right in regenerated[index + 1 :]:
            separation = specificity.angular_separation_deg(
                (left["ra_deg"], left["dec_deg"]), (right["ra_deg"], right["dec_deg"])
            )
            check(separation >= 12.0 - 1e-10, "blind control separation drift")
    print("[PASS] 20 fresh deterministic real-sky controls; observed legacy controls excluded")

    plan = downloader.plan()
    counts = {}
    for item in plan:
        counts[item["kind"]] = counts.get(item["kind"], 0) + 1
        check(all(url.startswith("https://") for url in item["urls"]), "non-HTTPS frozen source")
    check(
        counts
        == {
            "ORION": 4,
            "REAL_SKY_CONTROL": 80,
            "HST_TARGET": 4,
            "HST_REAL_CONTROL": 80,
        },
        f"unexpected downloader plan: {counts}",
    )
    check(len(plan) == 168, "downloader product count")
    check(all("_wf3" in item["dst"].name for item in plan if item["kind"].startswith("HST")), "HST mosaic leaked into v2.1")
    print("[PASS] downloader plan = 168 HTTPS products; HST locked to WF3 science+weight pairs")

    length = protocol["corridor_null"]["length_pixels_normalized"]
    half_width = protocol["corridor_null"]["half_width_pixels_normalized"]
    specs = [specificity.corridor_spec(("self-test", index), length, half_width) for index in range(64)]
    image = np.random.default_rng(44).normal(size=(128, 128)).astype(np.float32)
    image = np.clip((image - image.min()) / (image.max() - image.min()), 0, 1)
    for spec in specs:
        crop = specificity.extract_corridor(image, spec)
        check(crop.shape == (128, 128) and np.all(np.isfinite(crop)), "corridor extraction")
        cx, cy = spec["center_xy"]
        check(0 <= cx < 128 and 0 <= cy < 128 and 0 <= spec["angle_rad"] < math.pi, "corridor bounds")
    check(len({specificity.canonical_sha256(spec) for spec in specs}) == 64, "corridor sampler collision")
    print("[PASS] deterministic random-position/random-orientation corridor sampler")

    genome = forge.default_genome()
    candidate = specificity.extract_corridor(image, specs[0])
    local_cfg = dict(protocol["corridor_null"])
    local_cfg["null_corridors_per_image"] = 31
    local = specificity.corridor_local_rank(image, candidate, genome, "self-test-file", "SELFTEST", local_cfg)
    check(local["null_count"] == 31 and 0 < local["p_empirical"] <= 1 and np.isfinite(local["candidate_score"]), "corridor local rank")
    print("[PASS] within-image corridor empirical rank")

    raw_a = np.random.default_rng(1).normal(size=(192, 176)).astype(np.float32)
    raw_b = np.random.default_rng(2).normal(size=(192, 176)).astype(np.float32)
    weight_a = np.ones_like(raw_a)
    weight_b = np.ones_like(raw_b)
    weight_a[:12] = 0
    weight_b[:16] = 0
    weight_a[:, :9] = 0
    weight_b[:, :11] = 0
    image_a, image_b, mask = specificity.common_valid_support_pair(
        raw_a, raw_b, weight_a, weight_b, genome, protocol["hst_target"]["common_valid_support"]
    )
    check(image_a.shape == image_b.shape == (128, 128), "mask-aware normalization shape")
    check(mask["common_valid_fraction"] >= 0.98 and mask["mask_gate_pass"], "mask-aware support gate")
    check(mask["native_side"] < max(raw_a.shape), "mask crop did not remove invalid footprint")
    print("[PASS] HST common-valid-support crop removes detector footprint")

    controls = [float(value) for value in range(20)]
    admitted = specificity.real_field_rank(20.5, controls)
    blocked_equal = specificity.real_field_rank(19.0, controls)
    check(admitted["outperforms_all_controls"] and admitted["p_empirical"] == 1 / 21, "real rank admission")
    check(not blocked_equal["outperforms_all_controls"] and blocked_equal["control_exceedances"] == 1, "ties/exceedances must block")
    try:
        specificity.real_field_rank(1.0, controls[:19])
    except RuntimeError:
        pass
    else:
        raise RuntimeError("underfilled real-control cohort escaped fail-closed gate")
    print("[PASS] real-field empirical admission is tie-conservative and fail-closed")

    base = np.random.default_rng(7).normal(size=(128, 128))
    transformed = np.exp(base / 3.0)
    correlation = specificity.morphology_correlation(base, transformed, protocol["morphology_agreement"])
    check(correlation > 0.99, "rank-standardized morphology should ignore monotone photometric transforms")
    print("[PASS] PSF/resolution-matched rank morphology agreement")

    source = (ROOT / "janus_cosmos_v2_1.py").read_text(encoding="utf-8")
    check("synthetic_p_alone_can_admit" not in source, "runtime contains a synthetic-p admission shortcut")
    check("real_field_gates" in source and "corridor_local_family_gate" in source, "runtime omits v2.1 admission gates")
    print("[PASS] runtime admission source scan: no synthetic-p-only shortcut")

    print("SELF-TEST PASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
