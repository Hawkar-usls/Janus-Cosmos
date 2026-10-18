#!/usr/bin/env python3
"""Rehydrate JANUS Cosmos V4 NASA FSSC query FITS; verify exact SHA256 bytes.

Query-result URLs may expire. If a resource is unavailable, repeat the NASA
Extended query and review the new bytes+hashes; NEVER silently replace data.
"""

import argparse
import hashlib
import json
import os
import re
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/"data/cosmos/COSMOS-ALI-V4-FERMI-RAW-FITS-MANIFEST-2026-10-18.json"


def trusted_url(url:str) -> bool:
    p=urllib.parse.urlparse(url)
    return (p.scheme=="https" and p.netloc=="fermi.gsfc.nasa.gov"
            and p.params==p.query==p.fragment==""
            and re.fullmatch(r"/FTP/fermi/data/lat/queries/L[A-Z0-9]+_(EV00|SC00)\.fits",p.path) is not None)


def verified(path:Path,entry:dict) -> bool:
    if not path.is_file() or path.stat().st_size!=entry["bytes"]:return False
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(65536),b""):h.update(chunk)
    return h.hexdigest()==entry["sha256"]


def retrieve(entry:dict,dest:Path) -> str:
    if not trusted_url(entry["url"]):raise ValueError("refusing non-NASA FSSC URL")
    if not re.fullmatch(r"\d{6}[A-Z]?_(EV00|SC00)\.fits",entry["filename"]):
        raise ValueError("unsafe filename")
    target=dest/entry["filename"]
    if verified(target,entry):return "already_verified"
    if target.exists():raise ValueError("existing binary SHA mismatch; refusing overwrite")
    fd,tmp=tempfile.mkstemp(prefix=".janus_v4_",dir=str(dest))
    try:
        h=hashlib.sha256();size=0
        with os.fdopen(fd,"wb") as f, urllib.request.urlopen(entry["url"],timeout=40) as remote:
            if not trusted_url(remote.geturl()):raise ValueError("unsafe NASA redirect")
            while True:
                block=remote.read(65536)
                if not block:break
                size+=len(block)
                if size>entry["bytes"]:raise ValueError("FITS download larger than manifest")
                h.update(block);f.write(block)
        if size!=entry["bytes"] or h.hexdigest()!=entry["sha256"]:
            raise ValueError("FITS download SHA256/size mismatch")
        os.replace(tmp,target)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)
    return "downloaded_verified"


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir",type=Path,required=True)
    parser.add_argument("--verify-only",action="store_true")
    args=parser.parse_args()
    manifest=json.loads(MANIFEST.read_text())
    if manifest["schema"]!="janus.cosmos.fermi-lat-fits-manifest.v4":
        raise ValueError("unknown manifest schema")
    if not args.verify_only:args.raw_dir.mkdir(parents=True,exist_ok=True)
    try:
        for burst in manifest["bursts"]:
            for role in ("event","spacecraft"):
                entry=burst[role]
                if args.verify_only:
                    if not verified(args.raw_dir/entry["filename"],entry):
                        raise ValueError("missing or altered FITS: "+entry["filename"])
                    status="verified"
                else:status=retrieve(entry,args.raw_dir)
                print(f"{burst['grb']} {role}: {status} SHA256={entry['sha256']}")
    except (OSError,ValueError,urllib.error.URLError) as exc:
        print(f"JANUS-FERMI-V4-FETCH-ERROR: {exc}",file=sys.stderr)
        return 2
    return 0


if __name__=="__main__":raise SystemExit(main())
