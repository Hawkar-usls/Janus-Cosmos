#!/usr/bin/env python3
"""Build a deterministic, provenance-first manifest of public HST galaxy products.

This stage discovers products only. It does not perform OCR, face detection,
semantic search, cipher search, or scientific candidate scoring.
"""
import hashlib
import json
import os
import random
from pathlib import Path

SEED = 20260811
FILTERS = {"F435W", "F555W", "F814W"}
TARGET_N = int(os.getenv("JANUS_TARGET_N", "100"))
OUT = Path(os.getenv("JANUS_MANIFEST_OUT", "janus-cosmos-expanded-mast-manifest.json"))

# Seeded placeholder registry for deterministic expansion. The discovery adapter
# can replace this list with live MAST results without changing the manifest schema.
SEED_TARGETS = [
    "NGC1365", "NGC1425", "NGC1637", "NGC2841", "NGC3031", "NGC3627", "NGC4321",
    "M51", "M81", "M82", "NGC253", "NGC4258", "NGC4565", "NGC5194", "NGC5195",
    "NGC2403", "NGC6946", "NGC5457", "NGC3351", "NGC3621", "NGC6744", "NGC7793"
]


def product_key(item):
    return (item.get("target", ""), item.get("filter", ""), item.get("product_uri", ""))


def main():
    rng = random.Random(SEED)
    targets = list(dict.fromkeys(SEED_TARGETS))
    rng.shuffle(targets)
    targets = targets[:TARGET_N]

    products = []
    for target in sorted(targets):
        for filt in sorted(FILTERS):
            products.append({
                "target": target,
                "filter": filt,
                "product_uri": None,
                "source": "MAST",
                "discovery_status": "PENDING_LIVE_MAST_DISCOVERY"
            })

    products = sorted({product_key(p): p for p in products}.values(), key=product_key)
    receipt = {
        "schema": "janus.cosmos.hst.expanded_mast_manifest.v0.1",
        "status": "MANIFEST_SKELETON",
        "seed": SEED,
        "requested_target_count": TARGET_N,
        "target_count": len(targets),
        "product_count": len(products),
        "filters": sorted(FILTERS),
        "products": products,
        "blind_gate": {
            "ocr": False,
            "face_search": False,
            "semantic_analysis": False,
            "cipher_search": False,
            "post_hoc_tuning": False
        },
        "manifest_sha256": None
    }
    canonical = json.dumps(receipt, sort_keys=True, separators=(",", ":")).encode()
    receipt["manifest_sha256"] = hashlib.sha256(canonical).hexdigest()
    OUT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": receipt["status"], "targets": len(targets), "products": len(products)}, indent=2))


if __name__ == "__main__":
    main()
