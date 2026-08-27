#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
DATA = OUT / "data"
TSV = DATA / "vega_tres_rv.tsv"
PROV = OUT / "vega_tres_rv_provenance.json"

ENDPOINTS = [
    "https://vizier.cds.unistra.fr/viz-bin/asu-tsv",
    "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv",
]
PARAMS = [
    ("-source", "J/AJ/161/157/table2"),
    ("-out", "Seq BJD RVel e_RVel"),
    ("-out.max", "2000"),
]
EXPECTED_ROWS = 1524


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def build_url(endpoint: str) -> str:
    return endpoint + "?" + urllib.parse.urlencode(PARAMS)


def parse_rows(text: str) -> list[tuple[int, float, float, float]]:
    rows: list[tuple[int, float, float, float]] = []
    header_seen = False
    for raw in text.splitlines():
        line = raw.strip("\ufeff\r\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if not header_seen:
            normalized = [p.strip() for p in parts]
            if normalized[:4] == ["Seq", "BJD", "RVel", "e_RVel"]:
                header_seen = True
            continue
        if len(parts) < 4:
            continue
        try:
            seq = int(parts[0].strip())
            bjd = float(parts[1].strip())
            rv = float(parts[2].strip())
            erv = float(parts[3].strip())
        except ValueError:
            # Unit/separator rows are intentionally skipped.
            continue
        if erv <= 0:
            continue
        rows.append((seq, bjd, rv, erv))
    return rows


def fetch(timeout: int = 60) -> dict:
    OUT.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    errors: list[dict[str, str]] = []
    for endpoint in ENDPOINTS:
        url = build_url(endpoint)
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "Janus-Cosmos-Vega-Spider/1.1"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as response:
                payload = response.read()
                content_type = response.headers.get("Content-Type", "")
            text = payload.decode("utf-8", errors="replace")
            rows = parse_rows(text)
            if len(rows) != EXPECTED_ROWS:
                raise RuntimeError(
                    f"VizieR row-count drift: got {len(rows)}, expected {EXPECTED_ROWS}"
                )
            TSV.write_bytes(payload)
            provenance = {
                "schema": "janus.cosmos.vega.spider.tres_rv_provenance.v1.1",
                "source": "VizieR J/AJ/161/157/table2",
                "catalog_doi": "10.26093/cds/vizier.51610157",
                "paper_doi": "10.3847/1538-3881/abdec8",
                "url": url,
                "endpoint": endpoint,
                "content_type": content_type,
                "bytes": len(payload),
                "sha256": sha256_bytes(payload),
                "parsed_rows": len(rows),
                "expected_rows": EXPECTED_ROWS,
                "first_bjd": min(r[1] for r in rows),
                "last_bjd": max(r[1] for r in rows),
                "status": "PASS",
            }
            PROV.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
            return provenance
        except Exception as exc:
            errors.append({"endpoint": endpoint, "error": f"{type(exc).__name__}: {exc}"})
            time.sleep(2)
    failure = {
        "schema": "janus.cosmos.vega.spider.tres_rv_provenance.v1.1",
        "status": "FAIL",
        "errors": errors,
    }
    PROV.write_text(json.dumps(failure, indent=2) + "\n", encoding="utf-8")
    raise RuntimeError("all VizieR endpoints failed: " + json.dumps(errors))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    if args.dry_run:
        for endpoint in ENDPOINTS:
            url = build_url(endpoint)
            assert "J%2FAJ%2F161%2F157%2Ftable2" in url
            assert "-out.max=2000" in url
            print(url)
        print("VEGA TRES SPIDER DRY-RUN PASS")
        return 0
    provenance = fetch()
    print("VEGA TRES SPIDER FETCH PASS")
    print("rows =", provenance["parsed_rows"])
    print("sha256 =", provenance["sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
