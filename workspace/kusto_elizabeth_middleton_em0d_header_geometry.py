#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,os,zipfile
from pathlib import Path
import requests
import rasterio

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0B-ARCHIVE-MEMBER-RECEIPT-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-ElizabethMiddleton-EM0D-header/1.0"}
TMP=OUT/"em0d_tmp";TMP.mkdir(parents=True,exist_ok=True)

def fetch_to(role,spec):
    p=TMP/f"{role}.zip"
    with requests.get(spec["url"],headers=UA,timeout=300,allow_redirects=True,stream=True) as r:
        r.raise_for_status()
        h=hashlib.sha256();n=0
        with open(p,"wb") as f:
            for chunk in r.iter_content(1024*1024):
                if chunk:
                    f.write(chunk);h.update(chunk);n+=len(chunk)
    if n!=int(spec["bytes"]):raise RuntimeError(f"{role}: archive bytes mismatch")
    if h.hexdigest()!=spec["sha256"]:raise RuntimeError(f"{role}: archive SHA mismatch")
    return p

def header_from_zip(zpath,member,spec):
    with zipfile.ZipFile(zpath) as zf:
        zi=zf.getinfo(member)
        if int(zi.file_size)!=int(spec["uncompressed_bytes"]):raise RuntimeError("member size mismatch")
        if int(zi.CRC)!=int(spec["crc32"]):raise RuntimeError("member CRC mismatch")
    uri=f"zip://{zpath}!{member}"
    with rasterio.open(uri) as src:
        return {
          "name":member,
          "width":src.width,"height":src.height,"count":src.count,
          "dtype":str(src.dtypes[0]),"nodata":src.nodata,
          "crs":str(src.crs) if src.crs is not None else None,
          "transform":list(src.transform)[:6],
          "resolution_m":[abs(float(src.transform.a)),abs(float(src.transform.e))],
          "bounds":[float(src.bounds.left),float(src.bounds.bottom),float(src.bounds.right),float(src.bounds.top)],
          "pixel_values_read":False
        }

bz=fetch_to("bathymetry",REC["bathymetry_archive"])
kz=fetch_to("backscatter",REC["backscatter_archive"])
bath_spec=REC["bathymetry_archive"]["tif_members"][0]
bath=header_from_zip(bz,bath_spec["name"],bath_spec)
backs=[]
for spec in REC["backscatter_archive"]["tif_members"]:
    backs.append(header_from_zip(kz,spec["name"],spec))

out={
 "artifact_id":"JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0D-RASTER-HEADER-GEOMETRY-2026-09-24-v1.0",
 "member_receipt":REC["artifact_id"],
 "bathymetry_header":bath,
 "backscatter_headers":backs,
 "pixel_values_read":False,
 "raster_band_read_calls":0,
 "groundtruth_assets_read":False,
 "next_gate":"EM0E_FREEZE_PAIRED_GEOMETRY_AND_BACKSCATTER_UNION",
 "claim_ceiling":"RASTER_HEADER_GEOMETRY_ONLY"
}
p=OUT/"JANUS-KUSTO-ELIZABETH-MIDDLETON-EM0D-RASTER-HEADER-GEOMETRY-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(raw)
