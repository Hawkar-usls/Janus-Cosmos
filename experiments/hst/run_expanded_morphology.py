from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
from pathlib import Path

import numpy as np

import run_blind_corpus_morphology as morph
from run_blind_corpus_variant import download_variant, base

DEFAULT_SEEDS = (20260810, 20260811, 20260812)


def load_manifest(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("status") != "LIVE_MAST_DISCOVERY":
        raise RuntimeError(f"Expected LIVE_MAST_DISCOVERY manifest, got {data.get('status')!r}")
    return data


def robust_from_filter(result: dict, alpha: float) -> bool:
    return (
        result["morphology_preserving_phase"]["p_empirical"] < alpha
        and result["local_block_shuffle"]["p_empirical"] < alpha
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/hst_expanded_live_manifest.json")
    ap.add_argument("--backend", choices=["mast_api"], default="mast_api")
    ap.add_argument("--nulls", type=int, default=2048)
    ap.add_argument("--seeds", default=",".join(map(str, DEFAULT_SEEDS)))
    ap.add_argument("--output-dir", default="results/expanded_morphology")
    args = ap.parse_args()

    seeds = tuple(int(x) for x in args.seeds.split(",") if x.strip())
    if not seeds:
        raise ValueError("At least one seed is required")
    manifest = load_manifest(Path(args.manifest))
    targets = manifest["targets"]
    target_count = len(targets)
    filter_count = sum(len(t["filters"]) for t in targets)
    # Family-wise control across target/filter decisions and two morphology nulls.
    alpha = 0.05 / max(1, target_count * max(1, filter_count) * 2)

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    base.EVENT_LOG = out_dir / "janus-cosmos-expanded-morphology-events.jsonl"
    base.RECEIPT = out_dir / "janus-cosmos-expanded-morphology-receipt.json"
    base.EVENT_LOG.unlink(missing_ok=True)

    morph.NULLS = args.nulls
    all_seed_results = []
    base.emit("expanded_run_started", schema="janus.cosmos.expanded_morphology_event.v0.1", backend=args.backend, nulls=args.nulls, seeds=list(seeds), alpha_familywise=alpha, target_count=target_count, filter_count=filter_count, semantic_analysis=False, ocr=False, face_search=False, cipher_search=False, post_hoc_tuning=False)

    with tempfile.TemporaryDirectory() as td:
        for seed in seeds:
            rng = np.random.default_rng(seed)
            seed_targets = []
            for target in targets:
                t = {"target": target["target"], "filters": {}, "robust_passing_filters": [], "robust_cross_band_candidate": False}
                for item in target["filters"]:
                    path = Path(td) / f"{target['target']}_{item['filter']}_{seed}.fits"
                    metadata = download_variant(item["dataURI"], path, args.backend)
                    image = base.read_image(path)
                    result = morph.analyze_filter(image, rng, target["target"], item["filter"])
                    result["familywise_alpha"] = alpha
                    result["familywise_candidate"] = robust_from_filter(result, alpha)
                    t["filters"][item["filter"]] = result
                    t.setdefault("source_products", []).append({"filter": item["filter"], "band": item["band"], **metadata})
                    if result["familywise_candidate"]:
                        t["robust_passing_filters"].append(item["filter"])
                t["robust_cross_band_candidate"] = len(t["robust_passing_filters"]) >= 2
                seed_targets.append(t)
            all_seed_results.append({"seed": seed, "targets": seed_targets})

    consensus = []
    for target in targets:
        name = target["target"]
        per_seed = [s["targets"][[t["target"] for t in s["targets"]].index(name)] for s in all_seed_results]
        common_filters = []
        for filt in sorted({f["filter"] for f in target["filters"]}):
            passed_all = all(filt in t["robust_passing_filters"] for t in per_seed)
            if passed_all:
                common_filters.append(filt)
        consensus.append({"target": name, "seed_count": len(per_seed), "consensus_passing_filters": common_filters, "consensus_cross_band_candidate": len(common_filters) >= 2})

    receipt = {
        "schema": "janus.cosmos.hst.expanded_morphology_receipt.v0.1",
        "status": "EXPANDED_MORPHOLOGY_BLIND_TEST",
        "source": manifest["source"],
        "manifest_sha256": manifest["manifest_sha256"],
        "backend": args.backend,
        "target_count": target_count,
        "filter_count": filter_count,
        "nulls_per_null_model": args.nulls,
        "seeds": list(seeds),
        "familywise_alpha": alpha,
        "null_models": {
            "legacy_pixel_permutation": "diagnostic baseline only",
            "morphology_preserving_phase": {"low_freq_fraction": morph.LOW_FREQ_FRACTION, "phase_strength": morph.PHASE_STRENGTH},
            "local_block_shuffle": {"block_size": morph.BLOCK},
        },
        "seed_runs": all_seed_results,
        "consensus": consensus,
        "consensus_candidate_count": sum(1 for x in consensus if x["consensus_cross_band_candidate"]),
        "blind_gate": {"semantic_analysis": False, "ocr": False, "face_search": False, "cipher_search": False, "post_hoc_tuning": False, "cross_filter_required": True},
        "claim_ceiling": "Robust image-level geometric candidates only; no astronomical discovery claim without independent replication and scientific review.",
    }
    base.RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    receipt_sha = hashlib.sha256(base.RECEIPT.read_bytes()).hexdigest()
    base.emit("expanded_run_completed", consensus_candidate_count=receipt["consensus_candidate_count"], receipt_sha256=receipt_sha)
    receipt["event_log"] = {"path": str(base.EVENT_LOG), "sha256": hashlib.sha256(base.EVENT_LOG.read_bytes()).hexdigest()}
    receipt["receipt_sha256"] = receipt_sha
    base.RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"target_count": target_count, "filter_count": filter_count, "nulls": args.nulls, "seeds": list(seeds), "familywise_alpha": alpha, "consensus_candidate_count": receipt["consensus_candidate_count"], "receipt_sha256": receipt_sha}, indent=2))


if __name__ == "__main__":
    main()
