#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import requests
from rasterio.io import MemoryFile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
SRC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EMR0C-EXACT-SOURCE-SELECTION-FREEZE-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-ElizabethMiddleton-EMR0D-headers/1.0"}

def load_member(archive_url, expected_sha, expected_bytes, member, expected_member_bytes, expected_crc):
    r=requests.get(archive_url,headers=UA,timeout=360,allow_redirects=True);r.raise_for_status()
    blob=r.content
    if hashlib.sha256(blob).hexdigest()!=expected_sha:raise RuntimeError("archive SHA mismatch")
    if len(blob)!=expected_bytes:raise RuntimeError("archive byte length mismatch")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        zi=zf.getinfo(member)
        if zi.file_size!=expected_member_bytes or zi.CRC!=expected_crc:raise RuntimeError("member identity mismatch")
        tif=zf.read(member)
    with MemoryFile(tif) as mf:
        with mf.open() as ds:
            return {
              "member":member,
              "width":ds.width,"height":ds.height,"count":ds.count,
              "dtype":ds.dtypes[0],"crs":str(ds.crs) if ds.crs else None,
              "transform":list(ds.transform)[:6],
              "resolution":[abs(float(ds.transform.a)),abs(float(ds.transform.e))],
              "bounds":[float(ds.bounds.left),float(ds.bounds.bottom),float(ds.bounds.right),float(ds.bounds.top)],
              "nodata":ds.nodata,
              "member_sha256":hashlib.sha256(tif).hexdigest(),
              "pixel_values_read":False
            }

d=SRC["discovery_source"]
bathy=load_member(d["url"],d["archive_sha256"],d["archive_bytes"],d["member"],d["member_uncompressed_bytes"],d["member_crc32"])
h=SRC["heldout_backscatter_source"]
backs=[load_member(h["url"],h["archive_sha256"],h["archive_bytes"],m["name"],m["uncompressed_bytes"],m["crc32"]) for m in h["members"]]

out={
 "artifact_id":"JANUS-KUSTO-ELIZABETH-MIDDLETON-EMR0D-RASTER-HEADER-AUDIT-2026-09-24-v1.0",
 "source_freeze":SRC["artifact_id"],
 "bathymetry_header":bathy,
 "backscatter_headers":backs,
 "all_pixel_values_read":False,
 "groundtruth_assets_read":False,
 "claim_ceiling":"RASTER_HEADER_GEOMETRY_ONLY"
}
p=OUT/"JANUS-KUSTO-ELIZABETH-MIDDLETON-EMR0D-RASTER-HEADER-AUDIT-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(raw)
