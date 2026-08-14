#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import tempfile
import time
from pathlib import Path

import numpy as np

import demiurge_forge_v2 as forge
import download_sky_v2_1 as downloader
import janus_cosmos_core_v2 as core
import janus_cosmos_specificity_v2_1 as specificity
import janus_cosmos_v2_0 as parent_runtime
import janus_cosmos_v2_1 as runtime


ROOT = Path(__file__).resolve().parent


def check(condition, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    print("JANUS COSMOS v2.1.1 DETECTOR SPECIFICITY REPAIR SELF-TEST", flush=True)
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
    check(len(negative["artifact_sha256"]) == 10, "negative certificate artifact inventory drift")
    check(
        all(len(value) == 64 and set(value) <= set("0123456789abcdef") for value in negative["artifact_sha256"].values()),
        "negative certificate contains an invalid artifact hash",
    )
    present_artifacts = 0
    for name, expected_sha in negative["artifact_sha256"].items():
        artifact = evidence_root / name
        if artifact.is_file():
            check(core.sha256_file(artifact) == expected_sha, f"negative evidence artifact drift: {name}")
            present_artifacts += 1
    check(negative["scientific_gate"]["specificity_gate"] == "FAIL", "negative certificate gate drift")
    check(not negative["scientific_gate"]["orion_candidate_admitted"], "negative Orion result drift")
    check(not negative["scientific_gate"]["ngc1425_candidate_admitted"], "negative NGC1425 result drift")
    print(
        f"[PASS] v2.0.2 negative specificity certificate + 10 hash bindings "
        f"({present_artifacts} raw artifacts present in this package)"
    )

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
    check(downloader.MAX_WORKERS == 10, "download worker ceiling drift")
    check(all("_wf3" in item["dst"].name for item in plan if item["kind"].startswith("HST")), "HST mosaic leaked into v2.1")
    print("[PASS] downloader plan = 168 HTTPS products; HST locked to WF3 science+weight pairs")

    check(runtime.resolve_workers(10) == 10, "ten-worker runtime profile rejected")
    check(runtime.resolve_workers(1000) == 10 and runtime.resolve_workers(0) == 1, "worker clamp drift")
    _, legacy_settings = runtime.checkpoint_key(
        "f" * 64,
        protocol["parent_detector"]["genome_sha256"],
        "SELFTEST",
        "phase_iaaft",
        protocol["synthetic_null_diagnostics"]["test_nulls_per_model"],
        protocol["synthetic_null_diagnostics"]["calibration_nulls_per_model"],
        list(protocol["synthetic_null_diagnostics"]["seeds"]),
    )
    check(legacy_settings["version"] == "2.1.0", "v2.1.0 checkpoint identity bridge drift")
    with tempfile.TemporaryDirectory() as temporary:
        receipt = Path(temporary) / "atomic.json"
        runtime.atomic_write_json(receipt, {"status": "PASS", "workers": 10})
        check(json.loads(receipt.read_text(encoding="utf-8"))["workers"] == 10, "atomic checkpoint write")
        check(not list(Path(temporary).glob("*.tmp")), "atomic checkpoint temp leak")
    print("[PASS] bounded ten-worker profile + v2.1.0 resume identity + atomic checkpoint write")

    original_emit = runtime.emit
    original_persist = runtime.persist_running_report
    try:
        runtime.emit = lambda *args, **kwargs: None
        runtime.persist_running_report = lambda report: None
        scheduler_report = {"global_status": {}, "test_rows": {}}

        def delayed_row(item: dict) -> dict:
            time.sleep(float(item["delay"]))
            return {"id": item["id"]}

        scheduled = runtime.parallel_ordered_fields(
            [
                {"id": "FIRST", "delay": 0.03},
                {"id": "SECOND", "delay": 0.01},
                {"id": "THIRD", "delay": 0.00},
            ],
            lambda item: item["id"],
            delayed_row,
            3,
            "test_start",
            "test_complete",
            scheduler_report,
            "test_rows",
        )
        check(list(scheduled) == ["FIRST", "SECOND", "THIRD"], "parallel completion order leaked into report")
    finally:
        runtime.emit = original_emit
        runtime.persist_running_report = original_persist
    print("[PASS] parallel scheduler restores frozen input order after out-of-order completion")

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
    check("ProcessPoolExecutor" in source and "as_completed" in source, "runtime omits bounded process scheduler")
    check('get_context("spawn")' in source, "runtime process scheduler is not Windows-spawn safe")
    check("FROZEN_INPUT_ORDER" in source, "parallel completion order can contaminate report order")
    check("legacy_v2_1_0_model_checkpoints_accepted_when_settings_hash_matches" in source, "v2.1.0 resume bridge missing")
    print("[PASS] runtime admission source scan: no synthetic-p-only shortcut")

    print("SELF-TEST PASS", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
