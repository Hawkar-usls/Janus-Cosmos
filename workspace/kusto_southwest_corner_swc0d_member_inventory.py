#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC0C-CURRENT-ROUTE-RECOVERY-RECEIPT-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-SouthwestCorner-SWC0D-members/1.0"}

def inventory(role,url,expected_bytes):
    r=requests.get(url,headers=UA,timeout=300,allow_redirects=True)
    r.raise_for_status()
    blob=r.content
    if expected_bytes is not None and len(blob)!=int(expected_bytes):
        raise RuntimeError(f"{role}: byte length mismatch {len(blob)} != {expected_bytes}")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        members=[{
          "name":zi.filename,
          "file_size":zi.file_size,
          "compress_size":zi.compress_size,
          "crc32":zi.CRC,
          "is_dir":zi.is_dir()
        } for zi in zf.infolist()]
    return {
      "role":role,
      "url":url,
      "zip_sha256":hashlib.sha256(blob).hexdigest(),
      "zip_bytes":len(blob),
      "members":members
    }

rows=[
 inventory("bathymetry",REC["bathymetry"]["current_url"],REC["bathymetry"]["content_length"]),
 inventory("backscatter",REC["backscatter"]["current_url"],REC["backscatter"]["content_length"])
]
out={
 "artifact_id":"JANUS-KUSTO-SOUTHWEST-CORNER-SWC0D-ARCHIVE-MEMBER-INVENTORY-2026-09-24-v1.0",
 "route_receipt":REC["artifact_id"],
 "archives":rows,
 "raster_values_read":False,
 "raster_members_decompressed":False,
 "underwater_imagery_read":False,
 "claim_ceiling":"ARCHIVE_HASH_AND_MEMBER_INVENTORY_ONLY"
}
p=OUT/"JANUS-KUSTO-SOUTHWEST-CORNER-SWC0D-ARCHIVE-MEMBER-INVENTORY-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(raw)
