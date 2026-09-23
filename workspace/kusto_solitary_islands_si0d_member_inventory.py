#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI0C-CURRENT-TRANSPORT-RECEIPT-2026-09-24-v1.0.json").read_text())
URL=REC["current_resource"]["current_url"]
UA={"User-Agent":"JANUS-KUSTO-SolitaryIslands-SI0D-members/1.0"}

r=requests.get(URL,headers=UA,timeout=240,allow_redirects=True);r.raise_for_status()
blob=r.content
sha=hashlib.sha256(blob).hexdigest()
if len(blob)!=int(REC["current_resource"]["declared_total_bytes"]):
    raise RuntimeError(f"declared/current byte mismatch: {len(blob)}")
with zipfile.ZipFile(io.BytesIO(blob)) as zf:
    members=[{
      "name":zi.filename,"file_size":zi.file_size,"compress_size":zi.compress_size,
      "crc32":zi.CRC,"is_dir":zi.is_dir()
    } for zi in zf.infolist()]
    tif_members=[x for x in members if x["name"].lower().endswith((".tif",".tiff")) and not x["is_dir"]]

out={
 "artifact_id":"JANUS-KUSTO-SOLITARY-ISLANDS-SI0D-ARCHIVE-MEMBER-HASH-INVENTORY-2026-09-24-v1.0",
 "transport_receipt":REC["artifact_id"],
 "archive_url":URL,
 "archive_sha256":sha,
 "archive_bytes":len(blob),
 "members":members,
 "tif_members":tif_members,
 "raster_values_read":False,
 "members_decompressed_for_pixel_decode":False,
 "truth_layers_read":False,
 "next_gate":"SI0E_FREEZE_EXACT_BATHYMETRY_AND_BACKSCATTER_MEMBERS",
 "claim_ceiling":"ARCHIVE_HASH_AND_MEMBER_INVENTORY_ONLY"
}
p=OUT/"JANUS-KUSTO-SOLITARY-ISLANDS-SI0D-ARCHIVE-MEMBER-HASH-INVENTORY-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"archive_sha256":sha,"archive_bytes":len(blob),
 "tif_members":tif_members,"raster_values_read":False
},indent=2,ensure_ascii=False))
