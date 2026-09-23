#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHELLHARBOUR-S0-SOURCE-BINDING-RECEIPT-2026-09-23-v1.0.json").read_text())
URL=REC["pair_identity"]["bathymetry_product"]["data_url"]
UA={"User-Agent":"JANUS-KUSTO-Shellharbour-S0C-bathy-members/1.0"}
r=requests.get(URL,headers=UA,timeout=180,allow_redirects=True);r.raise_for_status()
blob=r.content
with zipfile.ZipFile(io.BytesIO(blob)) as zf:
    rows=[]
    for zi in zf.infolist():
        rows.append({
          "name":zi.filename,
          "file_size":zi.file_size,
          "compress_size":zi.compress_size,
          "crc":zi.CRC,
          "is_dir":zi.is_dir()
        })
out={
 "artifact_id":"JANUS-KUSTO-SHELLHARBOUR-S0C-BATHYMETRY-MEMBER-INVENTORY-2026-09-23-v1.0",
 "source_receipt":REC["artifact_id"],
 "bathymetry_zip_url":URL,
 "zip_sha256":hashlib.sha256(blob).hexdigest(),
 "zip_bytes":len(blob),
 "members":rows,
 "raster_values_read":False,
 "backscatter_values_read":False,
 "claim_ceiling":"ARCHIVE_MEMBER_INVENTORY_ONLY"
}
p=OUT/"JANUS-KUSTO-SHELLHARBOUR-S0C-BATHYMETRY-MEMBER-INVENTORY-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
