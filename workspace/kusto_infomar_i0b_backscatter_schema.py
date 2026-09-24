#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
BASE="https://image.marine.ie/arcgis/rest/services/INFOMAR/AllSurvey_Backscatter/ImageServer"
Q=BASE+"/query"
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I0B/1.0"})

tests=[
 ("return_ids",{"where":"1=1","returnIdsOnly":"true","f":"json"}),
 ("return_count",{"where":"1=1","returnCountOnly":"true","f":"json"}),
 ("features_1eq1",{"where":"1=1","outFields":"OBJECTID,Name,MinPS,MaxPS,Category,Tag,GroupName,ProductName,CenterX,CenterY","returnGeometry":"true","resultRecordCount":"10","f":"json"}),
 ("features_oid",{"where":"OBJECTID>0","outFields":"OBJECTID,Name,MinPS,MaxPS,Category,Tag,GroupName,ProductName,CenterX,CenterY","returnGeometry":"true","resultRecordCount":"10","f":"json"}),
 ("primary",{"where":"Category=1","outFields":"OBJECTID,Name,MinPS,MaxPS,Category,Tag,GroupName,ProductName,CenterX,CenterY","returnGeometry":"true","resultRecordCount":"10","f":"json"}),
]
rows=[]
for name,params in tests:
    r=S.get(Q,params=params,timeout=120,allow_redirects=True)
    text=r.text
    try:j=r.json()
    except Exception:j=None
    rows.append({
      "test":name,"status":r.status_code,"url":r.url,"bytes":len(r.content),
      "sha256":hashlib.sha256(r.content).hexdigest(),
      "json":j if isinstance(j,dict) else None,
      "text_prefix":text[:1000] if j is None else None
    })

out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I0B-BACKSCATTER-QUERY-SCHEMA-DIAGNOSTIC-2026-09-24-v1.0",
 "service":BASE,
 "tests":rows,
 "pixel_values_read":False,
 "export_image_called":False,
 "get_samples_called":False,
 "claim_ceiling":"ARCGIS_QUERY_SCHEMA_DIAGNOSTIC_ONLY"
}
p=OUT/"JANUS-KUSTO-INFOMAR-I0B-BACKSCATTER-QUERY-SCHEMA-DIAGNOSTIC-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "summary":[{"test":x["test"],"status":x["status"],
             "count":(x["json"] or {}).get("count"),
             "objectIds_len":len((x["json"] or {}).get("objectIds") or []),
             "features_len":len((x["json"] or {}).get("features") or []),
             "error":(x["json"] or {}).get("error")} for x in rows],
 "pixel_values_read":False
},indent=2,ensure_ascii=False))
