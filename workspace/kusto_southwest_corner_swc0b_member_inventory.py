#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC0-SOURCE-PAIR-RECEIPT-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-SouthwestCorner-SWC0B-members/1.0"}

def inventory(role,url):
    r=requests.get(url,headers=UA,timeout=240,allow_redirects=True);r.raise_for_status()
    blob=r.content
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        members=[{
          "name":zi.filename,"file_size":zi.file_size,"compress_size":zi.compress_size,
          "crc32":zi.CRC,"is_dir":zi.is_dir()
        } for zi in zf.infolist()]
    return {"role":role,"url":url,"zip_sha256":hashlib.sha256(blob).hexdigest(),"zip_bytes":len(blob),"members":members}

rows=[
 inventory("bathymetry",REC["physical_resources"]["bathymetry_geotiff_zip"]),
 inventory("backscatter",REC["physical_resources"]["backscatter_geotiff_zip"])
]
out={
 "artifact_id":"JANUS-KUSTO-SOUTHWEST-CORNER-SWC0B-ARCHIVE-MEMBER-INVENTORY-2026-09-24-v1.0",
 "source_pair_receipt":REC["artifact_id"],
 "archives":rows,
 "raster_values_read":False,
 "members_decompressed_for_pixel_decode":False,
 "underwater_imagery_read":False,
 "claim_ceiling":"ARCHIVE_MEMBER_AND_HASH_INVENTORY_ONLY"
}
p=OUT/"JANUS-KUSTO-SOUTHWEST-CORNER-SWC0B-ARCHIVE-MEMBER-INVENTORY-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps(out,indent=2))
