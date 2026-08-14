#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
import urllib.parse
from pathlib import Path

import requests


ROOT = Path(__file__).resolve().parent
PROTOCOL = json.loads((ROOT / "SPECIFICITY_PROTOCOL_v2_1.json").read_text(encoding="utf-8"))
DATA = ROOT / "external_data"
PROVENANCE = DATA / "download_provenance_v2_1.json"
ERRORS = DATA / "download_errors_v2_1.json"
HIPS_ENDPOINTS = [
    "https://alasky.cds.unistra.fr/hips-image-services/hips2fits",
    "https://alaskybis.cds.unistra.fr/hips-image-services/hips2fits",
]
SESSION = requests.Session()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fits_ok(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 2880:
        return False
    with path.open("rb") as handle:
        head = handle.read(80)
    return head.startswith(b"SIMPLE") or head.startswith(b"XTENSION")


def hips_url(endpoint: str, hips: str, ra: float, dec: float, fov: float, pixels: int) -> str:
    query = {
        "hips": hips,
        "format": "fits",
        "width": str(pixels),
        "height": str(pixels),
        "fov": str(fov),
        "projection": "TAN",
        "coordsys": "icrs",
        "rotation_angle": "0",
        "ra": f"{ra:.12f}",
        "dec": f"{dec:.12f}",
    }
    return endpoint + "?" + urllib.parse.urlencode(query)


def download(urls: list[str], destination: Path, retries: int = 3) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if fits_ok(destination):
        return {"status": "cached", "bytes": destination.stat().st_size, "sha256": sha256_file(destination), "url": "CACHE"}
    last_error: Exception | None = None
    for url in urls:
        for attempt in range(1, retries + 1):
            partial = destination.with_suffix(destination.suffix + ".part")
            try:
                with SESSION.get(
                    url,
                    stream=True,
                    timeout=(20, 300),
                    headers={"User-Agent": "Janus-Cosmos-v2.1/specificity-repair"},
                ) as response:
                    response.raise_for_status()
                    with partial.open("wb") as handle:
                        for chunk in response.iter_content(1024 * 1024):
                            if chunk:
                                handle.write(chunk)
                if not fits_ok(partial):
                    sample = partial.read_bytes()[:200] if partial.exists() else b""
                    raise RuntimeError("server response is not FITS: " + repr(sample))
                os.replace(partial, destination)
                return {
                    "status": "downloaded",
                    "bytes": destination.stat().st_size,
                    "sha256": sha256_file(destination),
                    "url": url,
                }
            except Exception as error:
                last_error = error
                if partial.exists():
                    partial.unlink()
                if attempt < retries:
                    time.sleep(attempt * 2)
    raise RuntimeError(f"download failed: {last_error}")


def _hips_item(kind: str, identifier: str, destination: Path, survey: dict, center: dict, fov: float, pixels: int) -> dict:
    return {
        "kind": kind,
        "id": identifier,
        "dst": destination,
        "urls": [
            hips_url(endpoint, survey["hips"], center["ra_deg"], center["dec_deg"], fov, pixels)
            for endpoint in HIPS_ENDPOINTS
        ],
        "center": center,
        "survey": survey,
    }


def _hst_products(field_id: str, kind: str, destination_root: Path) -> list[dict]:
    cfg = PROTOCOL["hst_real_controls"]
    base = cfg["base_url"]
    chip = cfg["canonical_chip"]
    rows = []
    for band in PROTOCOL["hst_target"]["bands"]:
        for suffix, role in (("", "SCIENCE"), ("_wgt", "WEIGHT")):
            filename = f"h_{field_id}_{band}_{chip}{suffix}.fits"
            rows.append(
                {
                    "kind": kind,
                    "id": f"HST_{field_id}_{band}_{chip}_{role}",
                    "dst": destination_root / filename,
                    "urls": [base + filename],
                    "field_id": field_id,
                    "band": band,
                    "chip": chip,
                    "role": role,
                }
            )
    return rows


def plan() -> list[dict]:
    rows: list[dict] = []
    target = PROTOCOL["orion_target"]
    center = target["center_j2000"]
    for survey in target["surveys"]:
        rows.append(
            _hips_item(
                "ORION",
                f"ORION_{survey['family']}_{survey['band']}",
                DATA / "orion" / survey["filename"],
                survey,
                center,
                float(target["fov_deg"]),
                int(target["pixels"]),
            )
        )
    controls = PROTOCOL["real_sky_controls"]
    for field in controls["centers"]:
        for survey in controls["surveys"]:
            filename = f"{field['id'].lower()}_{survey['family'].lower()}_{survey['band'].lower()}.fits".replace("2mass", "tmass")
            rows.append(
                _hips_item(
                    "REAL_SKY_CONTROL",
                    f"{field['id']}_{survey['family']}_{survey['band']}",
                    DATA / "controls_v2_1" / filename,
                    survey,
                    field,
                    float(controls["fov_deg"]),
                    int(controls["pixels"]),
                )
            )
    target_id = PROTOCOL["hst_target"]["id"]
    rows.extend(_hst_products(target_id, "HST_TARGET", DATA / "hst" / target_id))
    for field_id in PROTOCOL["hst_real_controls"]["field_ids"]:
        rows.extend(_hst_products(field_id, "HST_REAL_CONTROL", DATA / "hst" / "controls_v2_1" / field_id))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-real-sky-controls", action="store_true")
    parser.add_argument("--skip-hst-controls", action="store_true")
    args = parser.parse_args()
    selected = []
    for item in plan():
        if args.skip_real_sky_controls and item["kind"] == "REAL_SKY_CONTROL":
            continue
        if args.skip_hst_controls and item["kind"] == "HST_REAL_CONTROL":
            continue
        selected.append(item)
    records: list[dict] = []
    errors: list[dict] = []
    for item in selected:
        if args.dry_run:
            print(f"[DRY] {item['id']} -> {item['dst'].relative_to(ROOT)}")
            for url in item["urls"]:
                print("      " + url)
            continue
        try:
            info = download(item["urls"], item["dst"])
            record = {
                "kind": item["kind"],
                "id": item["id"],
                "file": str(item["dst"].relative_to(ROOT)),
                **info,
            }
            for key in ("center", "survey", "field_id", "band", "chip", "role"):
                if key in item:
                    record[key] = item[key]
            records.append(record)
            print(f"[OK] {item['id']} {info['status']} {info['bytes']} bytes", flush=True)
        except Exception as error:
            errors.append({"kind": item["kind"], "id": item["id"], "error": f"{type(error).__name__}: {error}"})
            print(f"[ERROR] {item['id']}: {error}", flush=True)
    if args.dry_run:
        print(f"DRY-RUN PASS: {len(selected)} of {len(plan())} frozen source products planned")
        return 0
    DATA.mkdir(parents=True, exist_ok=True)
    PROVENANCE.write_text(
        json.dumps(
            {
                "schema": "janus.cosmos.download_provenance.v2.1",
                "protocol_sha256": hashlib.sha256(
                    json.dumps(PROTOCOL, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
                ).hexdigest(),
                "records": records,
                "errors": errors,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    ERRORS.write_text(json.dumps(errors, indent=2, ensure_ascii=False), encoding="utf-8")
    print("DOWNLOAD PASS" if not errors else f"DOWNLOAD PARTIAL: {len(errors)} error(s)")
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
