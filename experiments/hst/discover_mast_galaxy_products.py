#!/usr/bin/env python3
"""Compatibility wrapper for the canonical Janus Cosmos v1 MAST discovery."""
from __future__ import annotations

import argparse
import json

from janus_cosmos.discovery import DEFAULT_TARGETS, DiscoveryConfig, build_manifest, write_manifest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="*", default=list(DEFAULT_TARGETS))
    ap.add_argument("--out", default="data/hst_expanded_live_manifest.json")
    ap.add_argument("--radius", type=float, default=0.20)
    args = ap.parse_args()

    manifest = build_manifest(
        args.targets,
        cfg=DiscoveryConfig(radius_deg=args.radius),
    )
    write_manifest(args.out, manifest)
    for row in manifest["target_status"]:
        print(json.dumps(row, sort_keys=True), flush=True)
    print(json.dumps({
        "selected_targets": manifest["selected_targets"],
        "selected_products": manifest["selected_products"],
        "output": args.out,
    }, indent=2), flush=True)
    return 0 if manifest["selected_targets"] > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
