#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI0-SOURCE-BINDING-RECEIPT-2026-09-24-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-SolitaryIslands-SI0B-members/1.0"}

def inventory(role,url):
    r=requests.get(url,headers=UA,timeout=240,allow_redirects=True)
    rec={"role":role,"requested_url":url,"status":r.status_code,"resolved_url":r.url,
         "content_type":r.headers.get("content-type"),"content_length":r.headers.get("content-length")}
    if not r.ok:
        rec["archive_opened"]=False
        rec["members"]=[]
        return rec
    blob=r.content
    rec.update({"archive_opened":True,"zip_sha256":hashlib.sha256(blob).hexdigest(),"zip_bytes":len(blob)})
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        rec["members"]=[{
          "name":zi.filename,"file_size":zi.file_size,"compress_size":zi.compress_size,
          "crc32":zi.CRC,"is_dir":zi.is_dir()
        } for zi in zf.infolist()]
    return rec

rows=[
 inventory("bathymetry",REC["pair_identity"]["bathymetry"]["data_url"]),
 inventory("backscatter",REC["pair_identity"]["backscatter"]["data_url"])
]
ok=all(x.get("archive_opened") for x in rows)
out={
 "artifact_id":"JANUS-KUSTO-SOLITARY-ISLANDS-SI0B-ARCHIVE-MEMBER-INVENTORY-2026-09-24-v1.0",
 "source_receipt":REC["artifact_id"],
 "archives":rows,
 "all_archives_opened":ok,
 "raster_values_read":False,
 "members_decompressed_for_pixel_decode":False,
 "towed_video_or_stills_read":False,
 "sediment_labels_read":False,
 "landform_or_substrate_classification_read":False,
 "next_gate":"SI0C_FREEZE_EXACT_BATHYMETRY_AND_BACKSCATTER_MEMBERS" if ok else "SOURCE_TRANSPORT_REPAIR_ONLY",
 "claim_ceiling":"ARCHIVE_ROUTE_MEMBER_AND_HASH_INVENTORY_ONLY"
}
p=OUT/"JANUS-KUSTO-SOLITARY-ISLANDS-SI0B-ARCHIVE-MEMBER-INVENTORY-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(raw)
