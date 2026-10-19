#!/usr/bin/env python3
"""Rehydrate V5-pinned official Fermi inputs for V6 (opt-in, byte-verified).

The source JSON defines the exclusive URL and SHA allowlist. The script does
not overwrite mismatching existing files or fetch unregistered URLs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/cosmos/COSMOS-ALI-V5-FERMI-LAT-GBM-CALDB-MANIFEST-2026-10-19.json"


def sha_file(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def safe_archive_url(url: str, is_gbm: bool) -> bool:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https" or parsed.username or parsed.password or parsed.query or parsed.fragment:
        return False
    if is_gbm:
        return parsed.hostname == "heasarc.gsfc.nasa.gov" and parsed.path.startswith("/FTP/fermi/data/gbm/triggers/")
    return parsed.hostname == "api.anaconda.org" and parsed.path == "/download/fermi/fermitools/2.5.8/linux-64/fermitools-2.5.8-py312h410524c_0.conda"


def download(url: str, dest: Path, expected_bytes: int | None, expected_sha: str, max_bytes: int) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        if (expected_bytes is None or dest.stat().st_size == expected_bytes) and sha_file(dest) == expected_sha:
            return
        raise RuntimeError(f"refusing to overwrite a mismatching pinned file: {dest.name}")
    with tempfile.NamedTemporaryFile(dir=dest.parent, prefix=".fermi-v6-", suffix=".part", delete=False) as stream:
        temp = Path(stream.name)
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "JANUS-Cosmos-V6-science-rehydration"}), timeout=120) as response:
                downloaded = 0
                h = hashlib.sha256()
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    downloaded += len(chunk)
                    if downloaded > max_bytes:
                        raise RuntimeError("external source exceeds pinned size cap")
                    stream.write(chunk)
                    h.update(chunk)
            if (expected_bytes is not None and downloaded != expected_bytes) or h.hexdigest() != expected_sha:
                raise RuntimeError(f"download SHA256/size mismatch: {dest.name}")
            stream.flush()
            os.replace(temp, dest)
        finally:
            temp.unlink(missing_ok=True)


def run(caldb: Path, gbm: Path, package: Path, only: str) -> dict:
    m = json.loads(MANIFEST.read_text())
    n_gbm = 0
    if only in ("all", "gbm"):
        for b in m["gbm_bursts"]:
            for detector in b["detectors"]:
                for meta in detector["event_and_response_files"].values():
                    if not safe_archive_url(meta["origin_url"], is_gbm=True):
                        raise ValueError("nonallowlisted GBM source URL")
                    name = meta["name"]
                    if name != Path(name).name:
                        raise ValueError("unsafe GBM manifest name")
                    download(meta["origin_url"], gbm / name, meta["bytes"], meta["sha256"], 60_000_000)
                    n_gbm += 1
    n_irf = 0
    if only in ("all", "caldb"):
        cm = m["lat_caldb"]
        if not safe_archive_url(cm["package_url"], is_gbm=False):
            raise ValueError("nonallowlisted CALDB package URL")
        download(cm["package_url"], package, None, cm["package_sha256"], 300_000_000)
        # Read only the exact allowlisted inner paths from the package.
        # The externally sourced *.tar.zst is hash-pinned *before* tar sees it.
        with zipfile.ZipFile(package) as arc:
            chunks = [n for n in arc.namelist() if n.startswith("pkg-") and n.endswith(".tar.zst")]
            if len(chunks) != 1:
                raise ValueError("unexpected Fermitools .conda structure")
            with tempfile.TemporaryDirectory(prefix="janus-v6-conda-") as temp:
                tarball = Path(temp) / "package.tar.zst"
                with arc.open(chunks[0]) as incoming, tarball.open("wb") as outgoing:
                    shutil.copyfileobj(incoming, outgoing)
                dest = caldb
                dest.mkdir(parents=True, exist_ok=True)
                paths = [entry["name"] for entry in cm["caldb_files"]]
                if any(not x.startswith("share/fermitools/data/caldb/data/glast/lat/bcf/") or ".." in Path(x).parts for x in paths):
                    raise ValueError("unexpected CALDB archive extraction path")
                subprocess.run(["tar", "--zstd", "--no-same-owner", "-xf", str(tarball), "-C", str(dest), *paths], check=True)
        for item in cm["caldb_files"]:
            p = caldb / item["name"]
            if not p.is_file() or p.stat().st_size != item["bytes"] or sha_file(p) != item["sha256"]:
                raise ValueError("extracted CALDB file not identical to pin")
            n_irf += 1
    return {"verified_gbm_files": n_gbm, "verified_caldb_irf_files": n_irf, "manifest": str(MANIFEST)}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--caldb-dir", type=Path, required=True)
    p.add_argument("--gbm-dir", type=Path, required=True)
    p.add_argument("--package", type=Path, required=True)
    p.add_argument("--only", choices=("all", "gbm", "caldb"), default="all")
    args = p.parse_args()
    print(json.dumps(run(args.caldb_dir, args.gbm_dir, args.package, args.only)))


if __name__ == "__main__":
    main()
