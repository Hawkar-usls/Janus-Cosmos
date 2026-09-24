#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import requests
from rasterio.io import MemoryFile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
SRC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0E-EXACT-MEMBER-SELECTION-FREEZE-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-Arafura-A0F-headers/1.0"}

def header(spec,role):
    r=requests.get(spec["archive_url"],headers=UA,timeout=360,allow_redirects=True);r.raise_for_status()
    blob=r.content
    if hashlib.sha256(blob).hexdigest()!=spec["archive_sha256"]:raise RuntimeError(f"{role}: archive SHA mismatch")
    if len(blob)!=int(spec["archive_bytes"]):raise RuntimeError(f"{role}: archive bytes mismatch")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        zi=zf.getinfo(spec["member"])
        if zi.file_size!=int(spec["member_uncompressed_bytes"]) or zi.CRC!=int(spec["member_crc32"]):
            raise RuntimeError(f"{role}: member identity mismatch")
        tif=zf.read(spec["member"])
    with MemoryFile(tif) as mf:
        with mf.open() as ds:
            return {
              "role":role,"member":spec["member"],"member_sha256":hashlib.sha256(tif).hexdigest(),
              "width":ds.width,"height":ds.height,"count":ds.count,"dtype":ds.dtypes[0],
              "crs":str(ds.crs) if ds.crs else None,"transform":list(ds.transform)[:6],
              "resolution_native_units":[abs(float(ds.transform.a)),abs(float(ds.transform.e))],
              "bounds":[float(ds.bounds.left),float(ds.bounds.bottom),float(ds.bounds.right),float(ds.bounds.top)],
              "nodata":ds.nodata,"pixel_values_read":False
            }

b=header(SRC["bathymetry"],"bathymetry")
h=header(SRC["heldout_backscatter"],"heldout_backscatter_money_shoal")
out={
 "artifact_id":"JANUS-KUSTO-ARAFURA-A0F-RASTER-HEADER-AUDIT-2026-09-24-v1.0",
 "source_freeze":SRC["artifact_id"],"selected_subarea":SRC["selected_subarea"],
 "bathymetry_header":b,"backscatter_header":h,
 "all_pixel_values_read":False,"groundtruth_assets_read":False,
 "next_gate":"A0G_FREEZE_MONEY_SHOAL_PAIRED_DOMAIN_AND_REPRESENTATION_BRIDGE",
 "claim_ceiling":"RASTER_HEADER_GEOMETRY_ONLY"
}
p=OUT/"JANUS-KUSTO-ARAFURA-A0F-RASTER-HEADER-AUDIT-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(raw)
