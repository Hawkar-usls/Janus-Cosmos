#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
DATA = OUT / "miri_final"
PROV = OUT / "vega_miri_final_provenance.json"

SOURCE_REPO = "merope82/Vega"
SOURCE_COMMIT_REF = "main"
BASE = "https://raw.githubusercontent.com/merope82/Vega/main"
FILES = {
    "F1550C": {
        "filename": "Vega_F1550C_final.fits",
        "expected_bytes": 1059840,
        "github_blob_sha": "d23be7510d60511b7f49c26289aff0a41061a547",
        "wavelength_micron": 15.5,
    },
    "F2300C": {
        "filename": "Vega_F2300C_final.fits",
        "expected_bytes": 1059840,
        "github_blob_sha": "023814db2aac32f02cb6c50bb120ce1487bfe8b2",
        "wavelength_micron": 23.0,
    },
    "F2550W": {
        "filename": "Vega_F2550W_final.fits",
        "expected_bytes": 4204800,
        "github_blob_sha": "b45a9df675e059143c6cb3c0529c6754c6c45aef",
        "wavelength_micron": 25.5,
    },
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    products = []
    for filt, meta in FILES.items():
        url = f"{BASE}/{meta['filename']}"
        req = urllib.request.Request(url, headers={"User-Agent": "Janus-Cosmos-Vega-MIRI-Spider/1.0"})
        with urllib.request.urlopen(req, timeout=120) as response:
            payload = response.read()
            content_type = response.headers.get("Content-Type", "")
        if len(payload) != meta["expected_bytes"]:
            raise RuntimeError(
                f"size drift for {meta['filename']}: got {len(payload)}, expected {meta['expected_bytes']}"
            )
        path = DATA / meta["filename"]
        path.write_bytes(payload)
        products.append(
            {
                "filter": filt,
                "filename": meta["filename"],
                "wavelength_micron": meta["wavelength_micron"],
                "source_url": url,
                "source_repo": SOURCE_REPO,
                "source_ref": SOURCE_COMMIT_REF,
                "github_blob_sha": meta["github_blob_sha"],
                "bytes": len(payload),
                "sha256": sha256_bytes(payload),
                "content_type": content_type,
                "status": "PASS",
            }
        )

    report = {
        "schema": "janus.cosmos.vega.miri_final_provenance.v1.4",
        "source_repo": SOURCE_REPO,
        "source_statement": "Paper-author public repository of final reduced JWST/MIRI Vega debris-disk FITS images.",
        "paper": "Su et al. 2024, ApJ 977, 277, DOI 10.3847/1538-4357/ad8cde",
        "products": products,
        "status": "PASS",
        "claim_firewall": "These are reduced image products and their provenance; downloading them does not by itself test or detect a planet.",
    }
    PROV.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("VEGA MIRI FINAL SPIDER FETCH PASS")
    for p in products:
        print(p["filter"], p["filename"], p["bytes"], p["sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
