#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
LAYER="https://gsi.geodata.gov.ie/server/rest/services/Marine/IE_GSI_MI_Marine_Download_Seabed_Survey_Leg_Data_IE_Waters_WGS84/FeatureServer/0"
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I0E/1.0"})

meta=S.get(LAYER,params={"f":"json"},timeout=120);meta.raise_for_status();mj=meta.json()
fields=[{"name":f.get("name"),"alias":f.get("alias"),"type":f.get("type")} for f in mj.get("fields",[]) or []]
oid=mj.get("objectIdField") or mj.get("objectIdFieldName")
# Metadata rows only, no raster/file download.
q=LAYER+"/query"
params={
 "where":"1=1",
 "outFields":"*",
 "returnGeometry":"false",
 "orderByFields":(oid+" ASC") if oid else "",
 "resultRecordCount":"2000",
 "f":"json"
}
r=S.get(q,params=params,timeout=180);r.raise_for_status();j=r.json()
if "error" in j: raise RuntimeError(j["error"])
features=j.get("features",[]) or []
rows=[x.get("attributes") or {} for x in features]

# Keep fields potentially relevant to survey identity/download/product pairing prominent.
interesting=[]
for row in rows:
    low={str(k).lower():v for k,v in row.items()}
    interesting.append(row)

out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I0E-SURVEY-LEG-FEATURE-METADATA-INVENTORY-2026-09-24-v1.0",
 "layer":LAYER,
 "layer_name":mj.get("name"),
 "object_id_field":oid,
 "max_record_count":mj.get("maxRecordCount"),
 "fields":fields,
 "row_count":len(rows),
 "rows":rows,
 "query_url":r.url,
 "query_sha256":hashlib.sha256(r.content).hexdigest(),
 "geometry_read":False,
 "raster_pixels_read":False,
 "download_files_opened":False,
 "claim_ceiling":"INFOMAR_SURVEY_LEG_FEATURE_METADATA_ONLY"
}
p=OUT/"JANUS-KUSTO-INFOMAR-I0E-SURVEY-LEG-FEATURE-METADATA-INVENTORY-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)

# Print only field names and compact rows with values in url/download/product/survey/name-related fields.
keys=[f["name"] for f in fields]
focus=[k for k in keys if any(t in k.lower() for t in ["survey","name","title","download","url","bathy","back","grid","resol","year","leg","ship","vessel","cruise","product","id"])]
compact=[]
for row in rows[:100]:
    compact.append({k:row.get(k) for k in focus if row.get(k) not in (None,"")})
print(json.dumps({
 "artifact_id":out["artifact_id"],"row_count":len(rows),
 "field_names":keys,"focus_fields":focus,
 "first_rows_focus":compact[:30],
 "raster_pixels_read":False,"download_files_opened":False
},indent=2,ensure_ascii=False))
