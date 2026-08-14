#!/usr/bin/env python3
"""Compatibility wrapper for the canonical Janus Cosmos v1 pipeline."""
from __future__ import annotations

import argparse
from pathlib import Path

from janus_cosmos.pipeline import parse_seeds, run_pipeline


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/hst_expanded_live_manifest.json")
    ap.add_argument("--nulls", type=int, default=4096)
    ap.add_argument("--seeds", default="20260810,20260811,20260812")
    ap.add_argument("--output-dir", default="results/expanded_morphology")
    ap.add_argument("--exploratory", action="store_true")
    ap.add_argument("--allow-underpowered", action="store_true")
    ap.add_argument("--backend", default="mast_api", help="Deprecated compatibility argument; ignored by canonical v1")
    args = ap.parse_args()

    receipt = run_pipeline(
        manifest_path=Path(args.manifest),
        output_dir=Path(args.output_dir),
        cache_dir=Path(".cache/janus_cosmos"),
        test_nulls=args.nulls,
        seeds=parse_seeds(args.seeds),
        exploratory=args.exploratory,
        allow_underpowered=args.allow_underpowered,
    )
    print("JANUS COSMOS CANONICAL EXPANDED STATUS =", receipt["status"])
    print("candidate_count =", receipt["candidate_count"])
    return 0 if receipt["status"] in {"PASS", "SMOKE_ONLY", "PARTIAL_ANALYSIS", "SMOKE_WITH_ERRORS"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
