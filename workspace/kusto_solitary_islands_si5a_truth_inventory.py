#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,zipfile
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI5-POSTRESULT-GROUNDTRUTH-CHARACTERIZATION-PREREG-2026-09-24-v1.0.json").read_text())
URL=PRE["truth_source"]["data_package_url"]
r=requests.get(URL,headers={"User-Agent":"JANUS-KUSTO-SolitaryIslands-SI5-truth-inventory/1.0"},timeout=240,allow_redirects=True);r.raise_for_status()
blob=r.content
with zipfile.ZipFile(io.BytesIO(blob)) as zf:
    members=[{
      "name":zi.filename,"file_size":zi.file_size,"compress_size":zi.compress_size,
      "crc32":zi.CRC,"is_dir":zi.is_dir()
    } for zi in zf.infolist()]

out={
 "artifact_id":"JANUS-KUSTO-SOLITARY-ISLANDS-SI5A-TRUTH-PACKAGE-INVENTORY-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "source_url":URL,
 "zip_sha256":hashlib.sha256(blob).hexdigest(),
 "zip_bytes":len(blob),
 "members":members,
 "candidate_file_names":[],
 "truth_file_contents_parsed":False,
 "formal_result_changed":False,
 "claim_ceiling":"POSTRESULT_TRUTH_PACKAGE_MEMBER_INVENTORY_ONLY"
}
keywords=("video","transect","sediment","sample","imagery","squidle","landform","substrate","ecosystem","csv","shp","gpkg","geojson")
out["candidate_file_names"]=[m["name"] for m in members if any(k in m["name"].lower() for k in keywords) and not m["is_dir"]]
p=OUT/"JANUS-KUSTO-SOLITARY-ISLANDS-SI5A-TRUTH-PACKAGE-INVENTORY-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"zip_sha256":out["zip_sha256"],"zip_bytes":out["zip_bytes"],
 "candidate_file_names":out["candidate_file_names"],"truth_file_contents_parsed":False
},indent=2,ensure_ascii=False))
