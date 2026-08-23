#!/usr/bin/env python3
"""Acquire public FKt230303 bathymetry products on GitHub Actions and run TOPA v1.1.

Raw public data are kept ephemeral; only checksums, provenance and candidate summaries
are committed. A failed acquisition is also a valid frozen result.
"""
from __future__ import annotations

import hashlib
import html
import json
import os
import re
import shutil
import sys
import tarfile
import zipfile
from pathlib import Path
from urllib.parse import urljoin, urlparse

import requests

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / ".topa_fkt230303_tmp"
OUT = ROOT / "data/tesla-sweep/results/TOPA-FKT230303-ACQUISITION-AND-FIELDGRID-RUN-001.json"
SOURCE_PAGES = [
    "https://www.marine-geo.org/tools/files/32240",  # ESRI ASCII products
    "https://www.marine-geo.org/tools/files/32239",  # NetCDF/GRD products
]
MAX_BYTES = 700 * 1024 * 1024
UA = "JANUS-TOPA-open-science-training/1.0 (+https://github.com/Hawkar-usls/Janus-Cosmos)"


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def safe_name_from_response(r: requests.Response, fallback: str) -> str:
    cd = r.headers.get("content-disposition", "")
    m = re.search(r"filename\*?=(?:UTF-8''|\")?([^\";]+)", cd, flags=re.I)
    if m:
        return Path(m.group(1).strip()).name
    name = Path(urlparse(r.url).path).name
    return name if name and "." in name else fallback


def download(session: requests.Session, url: str, dest_dir: Path, fallback: str):
    rec = {"requested_url": url}
    try:
        with session.get(url, timeout=(30, 240), allow_redirects=True, stream=True) as r:
            rec.update({
                "status_code": r.status_code,
                "final_url": r.url,
                "content_type": r.headers.get("content-type"),
                "content_disposition": r.headers.get("content-disposition"),
            })
            r.raise_for_status()
            name = safe_name_from_response(r, fallback)
            p = dest_dir / name
            total = 0
            with p.open("wb") as f:
                for chunk in r.iter_content(1024 * 1024):
                    if not chunk:
                        continue
                    total += len(chunk)
                    if total > MAX_BYTES:
                        raise RuntimeError(f"payload exceeded safety cap {MAX_BYTES} bytes")
                    f.write(chunk)
            rec.update({"bytes": total, "path": str(p), "sha256": sha256_file(p)})
            return p, rec
    except Exception as e:
        rec["error"] = f"{type(e).__name__}: {e}"
        return None, rec


def looks_html(p: Path, content_type: str | None) -> bool:
    if content_type and "html" in content_type.lower():
        return True
    try:
        return b"<html" in p.read_bytes()[:4096].lower() or b"<!doctype html" in p.read_bytes()[:4096].lower()
    except Exception:
        return False


def extract_links(p: Path, base_url: str):
    text = p.read_text(encoding="utf-8", errors="ignore")
    hrefs = re.findall(r"href\s*=\s*[\"']([^\"']+)[\"']", text, flags=re.I)
    links = []
    for href in hrefs:
        href = html.unescape(href)
        u = urljoin(base_url, href)
        low = u.lower()
        if any(k in low for k in ["download", ".zip", ".grd", ".nc", ".asc", ".gz", ".tgz", "file"]):
            if u not in links:
                links.append(u)
    return links[:40]


def unpack(p: Path, dest: Path):
    dest.mkdir(parents=True, exist_ok=True)
    low = p.name.lower()
    if zipfile.is_zipfile(p):
        with zipfile.ZipFile(p) as z:
            z.extractall(dest)
        return "zip"
    if tarfile.is_tarfile(p):
        with tarfile.open(p) as t:
            t.extractall(dest)
        return "tar"
    if low.endswith(".gz") and not low.endswith(".tar.gz"):
        import gzip
        out = dest / p.stem
        with gzip.open(p, "rb") as src, out.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        return "gzip"
    return "none"


def file_priority(p: Path):
    s = p.name.lower()
    score = 0
    if "puy" in s or "folles" in s: score += 100
    if "topo1m" in s or "1m" in s: score += 40
    if "topo" in s or "bath" in s: score += 20
    if "0.5" in s or "50cm" in s: score -= 10  # reserve ultra-high-res positives for later validation
    if p.suffix.lower() in {".asc", ".grd", ".nc", ".netcdf", ".npy"}: score += 10
    return -score, p.name


def main():
    if WORK.exists(): shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    dl = WORK / "downloads"; dl.mkdir()
    ex = WORK / "extracted"; ex.mkdir()

    sess = requests.Session()
    sess.headers.update({"User-Agent": UA, "Accept": "*/*"})
    acquisition = []
    payloads = []

    for i, page in enumerate(SOURCE_PAGES, 1):
        p, rec = download(sess, page, dl, f"mgds_page_{i}.bin")
        acquisition.append(rec)
        if p is None: continue
        if looks_html(p, rec.get("content_type")):
            links = extract_links(p, rec.get("final_url", page))
            rec["candidate_download_links"] = links
            for j, u in enumerate(links[:12], 1):
                q, qrec = download(sess, u, dl, f"mgds_link_{i}_{j}.bin")
                acquisition.append(qrec)
                if q is not None and not looks_html(q, qrec.get("content_type")):
                    payloads.append(q)
        else:
            payloads.append(p)

    unpack_records = []
    for idx, p in enumerate(payloads):
        target = ex / f"payload_{idx:02d}"
        try:
            kind = unpack(p, target)
            unpack_records.append({"path": str(p), "kind": kind})
            if kind == "none":
                target.mkdir(exist_ok=True)
                shutil.copy2(p, target / p.name)
        except Exception as e:
            unpack_records.append({"path": str(p), "error": f"{type(e).__name__}: {e}"})

    candidates = [p for p in ex.rglob("*") if p.is_file() and p.suffix.lower() in {".asc", ".grd", ".nc", ".netcdf", ".npy", ".txt"}]
    candidates.sort(key=file_priority)

    result = {
        "schema": "JANUS_TOPA_FKT230303_ACQUISITION_FIELDGRID_RUN",
        "version": "1.0",
        "run_environment": "github-actions",
        "branch": os.environ.get("GITHUB_REF_NAME"),
        "source_pages": SOURCE_PAGES,
        "acquisition_attempts": acquisition,
        "unpack_records": unpack_records,
        "discovered_grid_candidates": [str(p.relative_to(WORK)) for p in candidates[:100]],
        "raw_files_committed": false,
        "raw_data_policy": "public source bytes remain ephemeral; preserve source URL, response metadata and checksums only",
        "benchmark_engine": "workspace/tesla_sweep/topa_bathymetry_benchmark_v1_1.py",
        "metric_rule": "FROZEN_GENERIC_MORPHOLOGY_V1_1_AFTER_SYNTHETIC_QA",
        "classification_ceiling": "MORPHOLOGICAL_CANDIDATE_ONLY",
        "field_runs": [],
    }

    if not candidates:
        result["status"] = "ACQUISITION_BLOCKED__NO_NUMERIC_GRID_DISCOVERED"
        result["verdict"] = "NO_FIELD_ACCURACY_CLAIM"
    else:
        # First field pass: prioritize one Puy 1 m-like product. Coarse 128 px / 128 px
        # non-overlapping tiles reduce compute and avoid pretending fine chimney accuracy.
        chosen = candidates[0]
        result["selected_first_grid"] = {
            "path": str(chosen.relative_to(WORK)),
            "bytes": chosen.stat().st_size,
            "sha256": sha256_file(chosen),
        }
        out = WORK / "field_candidate_output.json"
        cmd = [
            sys.executable,
            str(ROOT / "workspace/tesla_sweep/topa_bathymetry_benchmark_v1_1.py"),
            str(chosen),
            "--tile-size", "128",
            "--stride", "128",
            "--top-k", "30",
            "--output", str(out),
        ]
        import subprocess
        cp = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=720)
        result["field_runs"].append({
            "command": cmd,
            "returncode": cp.returncode,
            "stdout_tail": cp.stdout[-4000:],
            "stderr_tail": cp.stderr[-4000:],
        })
        if cp.returncode == 0 and out.exists():
            field = json.loads(out.read_text(encoding="utf-8"))
            result["field_candidate_output"] = field
            result["status"] = "NUMERIC_FIELD_GRID_SCORED__GROUND_TRUTH_COMPARISON_NOT_PERFORMED"
            result["verdict"] = "REAL_FIELD_MORPHOLOGY_CANDIDATES_ONLY"
        else:
            result["status"] = "GRID_ACQUIRED__BENCHMARK_EXECUTION_FAILED"
            result["verdict"] = "NO_FIELD_ACCURACY_CLAIM"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({"status": result["status"], "candidate_files": len(candidates), "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
