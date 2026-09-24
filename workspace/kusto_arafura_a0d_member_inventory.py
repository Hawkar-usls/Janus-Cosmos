#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0C-CURRENT-ROUTE-RECEIPT-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-Arafura-A0D-members/1.0"}

archives=[]
for role in ("bathymetry","backscatter"):
    spec=REC["current_resources"][role]
    r=requests.get(spec["url"],headers=UA,timeout=360,allow_redirects=True);r.raise_for_status()
    blob=r.content
    if len(blob)!=int(spec["content_range_total_bytes"]):raise RuntimeError(f"{role}: byte length mismatch")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        members=[{
          "name":zi.filename,"file_size":zi.file_size,"compress_size":zi.compress_size,
          "crc32":zi.CRC,"is_dir":zi.is_dir()
        } for zi in zf.infolist()]
    archives.append({
      "role":role,"url":spec["url"],"zip_sha256":hashlib.sha256(blob).hexdigest(),
      "zip_bytes":len(blob),"members":members,
      "tiff_members":[m for m in members if m["name"].lower().endswith((".tif",".tiff")) and not m["is_dir"]]
    })

out={
 "artifact_id":"JANUS-KUSTO-ARAFURA-A0D-ARCHIVE-MEMBER-INVENTORY-2026-09-24-v1.0",
 "source_receipt":REC["artifact_id"],
 "selected_subarea":REC["selected_subarea"],
 "archives":archives,
 "raster_values_read":False,
 "tiff_members_decompressed_for_pixel_decode":False,
 "seabed_samples_or_imagery_read":False,
 "claim_ceiling":"ARCHIVE_MEMBER_HASH_CRC_PROVENANCE_ONLY"
}
p=OUT/"JANUS-KUSTO-ARAFURA-A0D-ARCHIVE-MEMBER-INVENTORY-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "archives":[{"role":a["role"],"zip_sha256":a["zip_sha256"],"zip_bytes":a["zip_bytes"],"tiff_members":a["tiff_members"]} for a in archives],
 "raster_values_read":False
},indent=2,ensure_ascii=False))
