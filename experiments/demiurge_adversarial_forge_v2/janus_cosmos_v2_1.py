#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
import traceback
from pathlib import Path

import numpy as np

import janus_cosmos_core_v2 as core
import janus_cosmos_specificity_v2_1 as specificity
import janus_cosmos_v2_0 as parent_runtime


VERSION = "2.1.0"
ROOT = Path(__file__).resolve().parent
PROTOCOL = json.loads((ROOT / "SPECIFICITY_PROTOCOL_v2_1.json").read_text(encoding="utf-8"))
DATA = ROOT / "external_data"
OUT = ROOT / "results_v2_1"
CHECKPOINTS = OUT / "checkpoints"
EVENTS = OUT / "janus-cosmos-v2.1-events.jsonl"
REPORT = OUT / "janus-cosmos-v2.1-report.json"
SUMMARY = OUT / "SUMMARY_v2.1.txt"
TERMINAL = OUT / "terminal_v2.1.log"
LOG_HANDLE = None


def log(message: str) -> None:
    print(message, flush=True)
    if LOG_HANDLE is not None:
        LOG_HANDLE.write(message + "\n")
        LOG_HANDLE.flush()


def emit(event: str, **fields) -> None:
    row = {"ts_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "event": event, **fields}
    encoded = json.dumps(row, sort_keys=True, ensure_ascii=False)
    log(encoded)
    with EVENTS.open("a", encoding="utf-8") as handle:
        handle.write(encoded + "\n")


def combined_sha256(paths: list[Path]) -> str:
    rows = [{"path": str(path.relative_to(ROOT)), "sha256": core.sha256_file(path)} for path in paths]
    return hashlib.sha256(specificity.canonical_json(rows).encode("utf-8")).hexdigest()


def checkpoint_key(file_sha: str, genome_sha: str, label: str, model: str, nulls: int, cal: int, seeds: list[int]) -> tuple[str, dict]:
    settings = {
        "version": VERSION,
        "file_sha": file_sha,
        "genome_sha": genome_sha,
        "protocol_sha": specificity.canonical_sha256(PROTOCOL),
        "label": label,
        "variant": "WHOLE",
        "model": model,
        "nulls": nulls,
        "cal": cal,
        "seeds": seeds,
    }
    return hashlib.sha256(specificity.canonical_json(settings).encode("utf-8")).hexdigest(), settings


def run_tail_model(
    x: np.ndarray,
    genome: dict,
    file_sha: str,
    label: str,
    model: str,
    nulls: int,
    cal: int,
    seeds: list[int],
) -> dict:
    key, settings = checkpoint_key(file_sha, core.sha256_bytes(core.canonical_json(genome).encode("utf-8")), label, model, nulls, cal, seeds)
    checkpoint = CHECKPOINTS / f"{label}__WHOLE__{model}.json"
    if checkpoint.exists():
        old = json.loads(checkpoint.read_text(encoding="utf-8"))
        if old.get("settings_hash") == key:
            log(f"[CACHE] {label} WHOLE {model}")
            return old["result"]

    def progress(seed: int, index: int, count: int) -> None:
        log(f"      {label} WHOLE {model} seed={seed}: {index}/{count}")

    emit("model_start", label=label, variant="WHOLE", model=model)
    result = specificity.empirical_test_with_tail(
        x,
        genome,
        model,
        nulls,
        cal,
        seeds,
        (label, "WHOLE"),
        progress=progress,
    )
    checkpoint.write_text(
        json.dumps({"settings_hash": key, "settings": settings, "result": result}, indent=2), encoding="utf-8"
    )
    emit(
        "model_complete",
        label=label,
        variant="WHOLE",
        model=model,
        p=result["p_empirical"],
        tail_ratio_q99=result["tail_ratio_q99"],
    )
    return result


def analyse_normalized(
    x: np.ndarray,
    genome: dict,
    file_sha: str,
    label: str,
    nulls: int,
    cal: int,
    seeds: list[int],
) -> dict:
    models = {
        model: run_tail_model(x, genome, file_sha, label, model, nulls, cal, seeds)
        for model in ("phase_iaaft", "block_shuffle")
    }
    return {"models": models, "band_tail_effect": specificity.band_tail_effect(models)}


def analyse_fits(path: Path, label: str, genome: dict, nulls: int, cal: int, seeds: list[int]) -> tuple[dict, np.ndarray, dict, tuple]:
    if not path.is_file():
        raise RuntimeError(f"missing frozen source product: {path.relative_to(ROOT)}; run download_sky_v2_1.py")
    raw, header, meta = core.read_primary_fits(path)
    x = core.normalize(raw, genome)
    file_sha = core.sha256_file(path)
    analysis = {
        "file": str(path.relative_to(ROOT)),
        "sha256": file_sha,
        "meta": meta,
        **analyse_normalized(x, genome, file_sha, label, nulls, cal, seeds),
    }
    return analysis, x, header, raw.shape


def control_file(field: dict, survey: dict) -> Path:
    filename = f"{field['id'].lower()}_{survey['family'].lower()}_{survey['band'].lower()}.fits".replace("2mass", "tmass")
    return DATA / "controls_v2_1" / filename


def corridor_field_score(bands: dict) -> float:
    # The local p-value remains a gate, but the cross-field rank uses a
    # continuous q99 tail ratio so p-floor saturation cannot create ties by
    # construction.
    values = [float(item["corridor"]["local_rank"]["tail_ratio_q99"]) for item in bands.values()]
    return float(min(values))


def local_corridor_family_gate(bands: dict) -> dict:
    required = int(PROTOCOL["corridor_null"]["required_passing_bands_per_family"])
    result = {}
    for family in ("DSS2", "2MASS"):
        passing = sum(
            bool(item["corridor"]["local_rank"]["passes_local_alpha"])
            for item in bands.values()
            if item["family"] == family
        )
        result[family] = {"passing_bands": passing, "required": required, "pass": passing >= required}
    result["pass"] = all(result[family]["pass"] for family in ("DSS2", "2MASS"))
    return result


def analyse_real_sky_field(
    field: dict,
    genome: dict,
    nulls: int,
    cal: int,
    seeds: list[int],
    is_target: bool,
) -> dict:
    surveys = PROTOCOL["orion_target"]["surveys"] if is_target else PROTOCOL["real_sky_controls"]["surveys"]
    bands = {}
    images: dict[tuple[str, str], np.ndarray] = {}
    corridors: dict[tuple[str, str], np.ndarray] = {}
    candidate_cfg = PROTOCOL["corridor_null"]
    heldout_spec = None
    if not is_target:
        heldout_spec = specificity.corridor_spec(
            (candidate_cfg["candidate_seed_domain"], field["id"]),
            float(candidate_cfg["length_pixels_normalized"]),
            float(candidate_cfg["half_width_pixels_normalized"]),
        )
    for survey in surveys:
        family, band = survey["family"], survey["band"]
        if is_target:
            path = DATA / "orion" / survey["filename"]
            label = f"ORION_{family}_{band}"
        else:
            path = control_file(field, survey)
            label = f"{field['id']}_{family}_{band}"
        analysis, x, header, native_shape = analyse_fits(path, label, genome, nulls, cal, seeds)
        if is_target:
            candidate, diagnostics = core.belt_corridor(
                x,
                header,
                native_shape,
                PROTOCOL["orion_target"]["belt_stars_j2000"],
                half_width=float(candidate_cfg["half_width_pixels_normalized"]),
                margin=8,
            )
        else:
            candidate = specificity.extract_corridor(x, heldout_spec)
            diagnostics = {"heldout_spec": heldout_spec, "designation": "FROZEN_RANDOM_CONTROL_CANDIDATE"}
        local_rank = specificity.corridor_local_rank(x, candidate, genome, analysis["sha256"], label, candidate_cfg)
        analysis["corridor"] = {"diagnostics": diagnostics, "local_rank": local_rank}
        bands[label] = {"family": family, "band": band, "analysis": analysis, "corridor": analysis["corridor"]}
        images[(family, band)] = x
        corridors[(family, band)] = candidate
    whole_analyses = [item["analysis"] for item in bands.values()]
    morphology_cfg = PROTOCOL["morphology_agreement"]
    return {
        "center": field,
        "bands": bands,
        "scores": {
            "whole_detector": specificity.field_tail_effect(whole_analyses),
            "corridor_detector": corridor_field_score(bands),
            "whole_cross_survey_morphology": specificity.orion_cross_survey_agreement(images, morphology_cfg),
            "corridor_cross_survey_morphology": specificity.orion_cross_survey_agreement(corridors, morphology_cfg),
        },
        "corridor_local_family_gate": local_corridor_family_gate(bands),
    }


def hst_paths(field_id: str, control: bool) -> dict[str, Path]:
    root = DATA / "hst" / ("controls_v2_1" if control else "")
    if control:
        root = root / field_id
    else:
        root = DATA / "hst" / field_id
    chip = PROTOCOL["hst_real_controls"]["canonical_chip"]
    paths = {}
    for band in PROTOCOL["hst_target"]["bands"]:
        paths[f"{band}_science"] = root / f"h_{field_id}_{band}_{chip}.fits"
        paths[f"{band}_weight"] = root / f"h_{field_id}_{band}_{chip}_wgt.fits"
    return paths


def analyse_hst_field(
    field_id: str,
    genome: dict,
    nulls: int,
    cal: int,
    seeds: list[int],
    control: bool,
) -> dict:
    paths = hst_paths(field_id, control)
    for path in paths.values():
        if not path.is_file():
            raise RuntimeError(f"missing frozen HST source product: {path.relative_to(ROOT)}; run download_sky_v2_1.py")
    f555_raw, _, f555_meta = core.read_primary_fits(paths["f555_science"])
    f555_weight, _, f555_weight_meta = core.read_primary_fits(paths["f555_weight"])
    f814_raw, _, f814_meta = core.read_primary_fits(paths["f814_science"])
    f814_weight, _, f814_weight_meta = core.read_primary_fits(paths["f814_weight"])
    x555, x814, mask_diagnostics = specificity.common_valid_support_pair(
        f555_raw,
        f814_raw,
        f555_weight,
        f814_weight,
        genome,
        PROTOCOL["hst_target"]["common_valid_support"],
    )
    all_paths = list(paths.values())
    field_sha = combined_sha256(all_paths)
    label_prefix = f"HSTCTRL_{field_id}" if control else "NGC1425"
    analyses = {
        "F555": analyse_normalized(x555, genome, field_sha + ":f555", label_prefix + "_F555_WF3_MASKED", nulls, cal, seeds),
        "F814": analyse_normalized(x814, genome, field_sha + ":f814", label_prefix + "_F814_WF3_MASKED", nulls, cal, seeds),
    }
    analyses["F555"]["source"] = {
        "science": str(paths["f555_science"].relative_to(ROOT)),
        "weight": str(paths["f555_weight"].relative_to(ROOT)),
        "science_meta": f555_meta,
        "weight_meta": f555_weight_meta,
    }
    analyses["F814"]["source"] = {
        "science": str(paths["f814_science"].relative_to(ROOT)),
        "weight": str(paths["f814_weight"].relative_to(ROOT)),
        "science_meta": f814_meta,
        "weight_meta": f814_weight_meta,
    }
    morphology = specificity.morphology_correlation(x555, x814, PROTOCOL["morphology_agreement"])
    return {
        "field_id": field_id,
        "combined_source_sha256": field_sha,
        "mask": mask_diagnostics,
        "bands": analyses,
        "scores": {
            "whole_detector": specificity.field_tail_effect(analyses.values()),
            "cross_filter_morphology": float(morphology),
        },
    }


def rank_gate(target_score: float, controls: list[float]) -> dict:
    rank = specificity.real_field_rank(target_score, controls)
    maximum = int(PROTOCOL["real_field_admission"]["maximum_control_exceedances"])
    rank["maximum_control_exceedances"] = maximum
    rank["pass"] = bool(rank["control_exceedances"] <= maximum)
    return rank


def build_summary(report: dict) -> str:
    lines = [
        "JANUS COSMOS v2.1.0 — DETECTOR SPECIFICITY REPAIR",
        "",
        f"status: {report['status']}",
        f"smoke_only: {report['smoke_only']}",
        f"protocol_sha256: {report['protocol']['protocol_sha256']}",
        f"frozen_genome: {report['frozen_detector']['genome_sha256']}",
    ]
    if report.get("global_status"):
        lines.extend(
            [
                "",
                f"REAL_SKY_CONTROL_FIELDS: {report['global_status'].get('real_sky_control_fields')}",
                f"HST_REAL_CONTROL_FIELDS: {report['global_status'].get('hst_real_control_fields')}",
                f"ORION: {report['orion'].get('status') if report.get('orion') else 'NOT_RUN'}",
                f"NGC1425: {report['ngc1425'].get('status') if report.get('ngc1425') else 'NOT_RUN'}",
            ]
        )
    if report["errors"]:
        lines.extend(["", "ERRORS:", *[f"- {item['error']}" for item in report["errors"]]])
    lines.extend(["", "Claim ceiling: " + PROTOCOL["claim_ceiling"]])
    return "\n".join(lines) + "\n"


def main() -> int:
    global LOG_HANDLE
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--nulls", type=int)
    parser.add_argument("--cal-nulls", type=int)
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    CHECKPOINTS.mkdir(exist_ok=True)
    EVENTS.write_text("", encoding="utf-8")
    LOG_HANDLE = TERMINAL.open("w", encoding="utf-8", buffering=1)
    frozen = parent_runtime.verify_forge()
    protocol_receipt = specificity.verify_protocol_sources(PROTOCOL)
    if frozen["genome_sha256"] != PROTOCOL["parent_detector"]["genome_sha256"] or frozen["freeze_sha256"] != PROTOCOL["parent_detector"]["freeze_sha256"]:
        raise RuntimeError("v2.1 protocol is not bound to the packaged frozen detector")
    if args.self_test:
        log("SELF-TEST PASS: v2.0.2 frozen detector + v2.1 protocol verified")
        LOG_HANDLE.close()
        return 0
    cfg = PROTOCOL["synthetic_null_diagnostics"]
    if not args.smoke and (args.nulls is not None or args.cal_nulls is not None):
        raise RuntimeError("full v2.1 run forbids Monte Carlo overrides; use the frozen protocol")
    nulls = int(args.nulls or (24 if args.smoke else cfg["test_nulls_per_model"]))
    cal = int(args.cal_nulls or (12 if args.smoke else cfg["calibration_nulls_per_model"]))
    seeds = list(cfg["seeds"])
    report = {
        "schema": "janus.cosmos.v2.1.report",
        "version": VERSION,
        "status": "RUNNING",
        "smoke_only": bool(args.smoke),
        "protocol": protocol_receipt,
        "frozen_detector": {"genome_sha256": frozen["genome_sha256"], "freeze_sha256": frozen["freeze_sha256"]},
        "negative_parent_certificate": PROTOCOL["negative_parent_certificate"],
        "synthetic_null_diagnostics": {"test_nulls_per_model": nulls, "calibration_nulls": cal, "seeds": seeds},
        "real_sky_controls": {},
        "hst_real_controls": {},
        "orion": {},
        "ngc1425": {},
        "global_status": {},
        "errors": [],
        "claim_ceiling": PROTOCOL["claim_ceiling"],
    }
    genome = frozen["genome"]
    try:
        sky_fields = PROTOCOL["real_sky_controls"]["centers"][: (2 if args.smoke else None)]
        for field in sky_fields:
            emit("real_sky_control_start", field_id=field["id"])
            report["real_sky_controls"][field["id"]] = analyse_real_sky_field(field, genome, nulls, cal, seeds, False)
            emit("real_sky_control_complete", field_id=field["id"])
        hst_fields = PROTOCOL["hst_real_controls"]["field_ids"][: (2 if args.smoke else None)]
        for field_id in hst_fields:
            emit("hst_real_control_start", field_id=field_id)
            report["hst_real_controls"][field_id] = analyse_hst_field(field_id, genome, nulls, cal, seeds, True)
            emit("hst_real_control_complete", field_id=field_id)
        report["global_status"]["real_sky_control_fields"] = len(report["real_sky_controls"])
        report["global_status"]["hst_real_control_fields"] = len(report["hst_real_controls"])
        if args.smoke:
            report["status"] = "PASS"
            report["global_status"]["admission_disabled"] = True
            report["global_status"]["admission_disabled_reason"] = "SMOKE_ONLY_INCOMPLETE_REAL_CONTROL_COHORT"
        else:
            emit("orion_target_start")
            orion = analyse_real_sky_field(PROTOCOL["orion_target"]["center_j2000"], genome, nulls, cal, seeds, True)
            emit("orion_target_complete")
            sky_controls = list(report["real_sky_controls"].values())
            orion["real_field_gates"] = {
                "whole_detector": rank_gate(orion["scores"]["whole_detector"], [item["scores"]["whole_detector"] for item in sky_controls]),
                "corridor_detector": rank_gate(orion["scores"]["corridor_detector"], [item["scores"]["corridor_detector"] for item in sky_controls]),
                "whole_cross_survey_morphology": rank_gate(
                    orion["scores"]["whole_cross_survey_morphology"]["score"],
                    [item["scores"]["whole_cross_survey_morphology"]["score"] for item in sky_controls],
                ),
                "corridor_cross_survey_morphology": rank_gate(
                    orion["scores"]["corridor_cross_survey_morphology"]["score"],
                    [item["scores"]["corridor_cross_survey_morphology"]["score"] for item in sky_controls],
                ),
            }
            orion_pass = orion["corridor_local_family_gate"]["pass"] and all(
                item["pass"] for item in orion["real_field_gates"].values()
            )
            orion["admitted"] = bool(orion_pass)
            orion["status"] = "SKY_FIXED_MORPHOLOGY_CANDIDATE" if orion_pass else "DETECTOR_SPECIFICITY_BLOCKED"
            report["orion"] = orion

            emit("ngc1425_target_start")
            ngc = analyse_hst_field(PROTOCOL["hst_target"]["id"], genome, nulls, cal, seeds, False)
            emit("ngc1425_target_complete")
            hst_controls = list(report["hst_real_controls"].values())
            ngc["real_field_gates"] = {
                "whole_detector": rank_gate(ngc["scores"]["whole_detector"], [item["scores"]["whole_detector"] for item in hst_controls]),
                "cross_filter_morphology": rank_gate(
                    ngc["scores"]["cross_filter_morphology"],
                    [item["scores"]["cross_filter_morphology"] for item in hst_controls],
                ),
            }
            ngc_pass = ngc["mask"]["mask_gate_pass"] and all(item["pass"] for item in ngc["real_field_gates"].values())
            ngc["admitted"] = bool(ngc_pass)
            ngc["status"] = "HST_CROSS_FILTER_MORPHOLOGY_CANDIDATE" if ngc_pass else "DETECTOR_SPECIFICITY_BLOCKED"
            report["ngc1425"] = ngc
            report["global_status"].update(
                {
                    "admission_disabled": False,
                    "orion_candidate_admitted": bool(orion_pass),
                    "ngc1425_candidate_admitted": bool(ngc_pass),
                    "claim_ceiling_enforced": True,
                }
            )
            report["status"] = "PASS"
    except Exception as error:
        report["status"] = "FAIL"
        report["errors"].append({"error": f"{type(error).__name__}: {error}", "traceback": traceback.format_exc(limit=16)})
        emit("fatal_error", error=str(error))
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    summary = build_summary(report)
    SUMMARY.write_text(summary, encoding="utf-8")
    log(summary.rstrip())
    if LOG_HANDLE:
        LOG_HANDLE.close()
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
