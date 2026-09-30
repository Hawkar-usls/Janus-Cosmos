#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from janus_cosmos.luci import read_luci_fits_image
from janus_cosmos.luci_psf_r1 import psf_relative_injection_recovery_gate
from janus_cosmos.pipeline import EventWriter, _source_for_item, download_source, load_manifest
from experiments.luci.run_palomar_2f_c import (
    MAX_EXACT_CANDIDATE_FILES,
    build_frozen_expanded_corpus,
    exact_wcs_pixel,
    local_counterpart,
    query_conservative_overlaps,
)

FROZEN_CORPUS_SHA256 = "b5211c2224a8db031e4ae87b31a155f5ee7be0b96927123b5b0593ba4390f71d"


def validate_frames(manifest_paths: list[Path], cache_dir: Path, events: EventWriter) -> dict:
    frames = []
    for mp in manifest_paths:
        manifest = load_manifest(mp)
        for target in manifest["targets"]:
            for item in target["filters"]:
                filt = str(item.get("filter", "UNKNOWN"))
                path, source_meta = download_source(
                    _source_for_item(item), cache_dir, events,
                    target=target["target"], filter_name=filt,
                )
                image, image_meta = read_luci_fits_image(
                    path, require_imaging=True, expected_instrument=item.get("instrument")
                )
                gate = psf_relative_injection_recovery_gate(
                    image, seed=20260831 + len(frames)
                )
                frames.append({
                    "target": target["target"],
                    "filter": filt,
                    "file_sha256": source_meta["sha256"],
                    "image_meta": image_meta,
                    "injection_gate": gate,
                })
    passed = sum(bool(x["injection_gate"].get("passed")) for x in frames)
    return {
        "frame_count": len(frames),
        "passed_frame_count": passed,
        "pass_fraction": passed / len(frames) if frames else 0.0,
        "all_frames_pass": bool(frames and passed == len(frames)),
        "frames": frames,
    }


def run_exact_chain(
    conservative: list[dict],
    *,
    global_validation_pass: bool,
    cache_dir: Path,
    events: EventWriter,
) -> tuple[list[dict], bool]:
    unique_files = sorted({x["file_name"] for x in conservative})
    if len(unique_files) > MAX_EXACT_CANDIDATE_FILES:
        return [], True

    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in conservative:
        grouped[row["file_name"]].append(row)

    results = []
    for fname in unique_files:
        group = grouped[fname]
        path, source_meta = download_source(
            group[0]["file_url"], cache_dir, events,
            target="PALOMAR_2F_C_R1_OVERLAP",
            filter_name=group[0].get("filters", "UNKNOWN"),
        )
        image_cache = None
        image_meta_cache = None
        gate_cache = None
        for cand in group:
            inside, wmeta = exact_wcs_pixel(path, cand["ra_deg"], cand["dec_deg"])
            record = {
                **cand,
                "exact_wcs_inside": bool(inside),
                "wcs": wmeta,
                "file_sha256": source_meta["sha256"],
            }
            if inside:
                if not global_validation_pass:
                    record["counterpart_test"] = {
                        "status": "BLOCKED_BY_GLOBAL_VALIDATION_GATE"
                    }
                else:
                    if image_cache is None:
                        image_cache, image_meta_cache = read_luci_fits_image(
                            path,
                            require_imaging=True,
                            expected_instrument=cand.get("instrument"),
                        )
                        gate_cache = psf_relative_injection_recovery_gate(
                            image_cache, seed=20261831 + len(results)
                        )
                    record["image_meta"] = image_meta_cache
                    record["overlap_frame_injection_gate"] = gate_cache
                    if gate_cache.get("passed"):
                        record["counterpart_test"] = local_counterpart(
                            image_cache, wmeta["x"], wmeta["y"]
                        )
                    else:
                        record["counterpart_test"] = {
                            "status": "BLOCKED_BY_OVERLAP_FRAME_INJECTION_GATE"
                        }
            results.append(record)
    return results, False


def main() -> int:
    ap = argparse.ArgumentParser(
        description="LUCI-PALOMAR-JPFM-2F-C-R1: local PSF-relative injection/recovery repair"
    )
    ap.add_argument("--validation-manifest", action="append", required=True)
    ap.add_argument("--output-dir", default="results/luci_palomar_2f_c_r1")
    ap.add_argument("--cache-dir", default=".cache/luci_palomar_2f_c_r1")
    args = ap.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    cache = Path(args.cache_dir)
    events = EventWriter(out / "events.jsonl")

    corpus = build_frozen_expanded_corpus(out / "frozen_palomar_640.csv")
    if corpus["sha256"] != FROZEN_CORPUS_SHA256:
        raise RuntimeError(
            f"frozen corpus SHA changed: {corpus['sha256']} != {FROZEN_CORPUS_SHA256}"
        )
    events.emit(
        "palomar_corpus_r1_sha_verified",
        sample_n=corpus["sample_n"],
        sha256=corpus["sha256"],
    )

    validation = validate_frames(
        [Path(x) for x in args.validation_manifest], cache / "validation", events
    )
    conservative, overlap_meta = query_conservative_overlaps(corpus["rows"])
    exact_results, cap_exceeded = run_exact_chain(
        conservative,
        global_validation_pass=validation["all_frames_pass"],
        cache_dir=cache / "overlaps",
        events=events,
    )

    exact_inside = [x for x in exact_results if x.get("exact_wcs_inside")]
    counterpart_tests = [
        x for x in exact_inside
        if x.get("overlap_frame_injection_gate", {}).get("passed")
        and x.get("counterpart_test", {}).get("counterpart_present")
    ]

    status = "PASS"
    if not validation["all_frames_pass"]:
        status = "BLOCKED"
        scientific = "PSF_RELATIVE_INJECTION_RECOVERY_VALIDATION_FAILED"
    elif cap_exceeded:
        status = "BLOCKED"
        scientific = "OVERLAP_CANDIDATE_CAP_EXCEEDED__NO_INSPECTION"
    elif not exact_inside:
        scientific = "PSF_GATE_VALIDATED__NO_EXACT_LUCI_OVERLAP_IN_FROZEN_640"
    else:
        scientific = "EXACT_OVERLAP_AVAILABLE__COUNTERPART_CHAIN_EXECUTED"

    receipt = {
        "schema": "janus.cosmos.luci_palomar.jpfm_2f_c_r1.receipt.v1",
        "experiment_id": "LUCI-PALOMAR-JPFM-2F-C-R1",
        "status": status,
        "scientific_status": scientific,
        "repair_scope": "measurement repair only; frozen corpus and decision thresholds unchanged from R1 preregistration",
        "frozen_palomar_corpus": {
            k: v for k, v in corpus.items() if k != "rows"
        },
        "required_frozen_corpus_sha256": FROZEN_CORPUS_SHA256,
        "validation_injection_recovery": validation,
        "archive_overlap_preflight": {
            **overlap_meta,
            "exact_candidate_file_cap": MAX_EXACT_CANDIDATE_FILES,
            "candidate_cap_exceeded": cap_exceeded,
        },
        "exact_wcs_overlap_pair_count": len(exact_inside),
        "exact_wcs_unique_palomar_sources": len({x["src_id"] for x in exact_inside}),
        "admitted_ir_counterpart_count": len(counterpart_tests),
        "exact_results": exact_results,
        "chain": "Palomar coordinate -> exact LUCI FITS WCS containment -> PSF-relative frame injection-recovery -> IR source/no-source -> PSF morphology -> matched local controls",
        "firewall": "Palomar-640 SHA must exactly match RUN-001 before LUCI outcomes are queried.",
        "claim_ceiling": "INDEPENDENT_NEAR_IR_COUNTERPART_TEST_ONLY__NO_ANOMALY_OR_UAP_ORIGIN_CLAIM__NO_CAUSALITY",
    }
    rp = out / "receipt.json"
    rp.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": status,
        "scientific_status": scientific,
        "corpus_sha256": corpus["sha256"],
        "validation_pass_fraction": validation["pass_fraction"],
        "conservative_pairs": len(conservative),
        "exact_wcs_pairs": len(exact_inside),
        "counterparts": len(counterpart_tests),
        "receipt": str(rp),
    }, indent=2))
    return 0 if status == "PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
