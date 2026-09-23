#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-EMODNET-CDI-ALL8-SOURCE-INVENTORY-PREREG-2026-09-22-v1.0.json").read_text())
WFS="https://ows.emodnet-bathymetry.eu/wfs"
TYPE="emodnet:source_references"
UA={"User-Agent":"JANUS-KUSTO-EMODnet-E1-source-references/1.0"}
S=requests.Session();S.headers.update(UA)

rows={}
for t in PRE["candidates"]:
    pad=0.02
    params={
      "SERVICE":"WFS","REQUEST":"GetFeature","VERSION":"2.0.0",
      "TYPENAMES":TYPE,
      "SRSNAME":"EPSG:4326",
      "BBOX":f'{t["lon"]-pad},{t["lat"]-pad},{t["lon"]+pad},{t["lat"]+pad},EPSG:4326',
      "COUNT":"100",
      "OUTPUTFORMAT":"application/json"
    }
    r=S.get(WFS,params=params,timeout=120,allow_redirects=True)
    rec={"target":t,"status_code":r.status_code,"url":r.url,"bytes":len(r.content),"sha256":hashlib.sha256(r.content).hexdigest()}
    feats=[]
    parse_error=None
    if r.ok:
        try:
            j=r.json()
            for f in j.get("features",[]) or []:
                feats.append({
                  "id":f.get("id"),
                  "properties":f.get("properties") or {},
                  "geometry_intentionally_omitted":True
                })
        except Exception as e:
            parse_error=type(e).__name__+": "+str(e)
    rec["parse_error"]=parse_error
    rec["feature_count"]=len(feats)
    rec["features"]=feats
    rows[t["id"]]=rec

hits={k:v for k,v in rows.items() if v["feature_count"]>0}
out={
 "artifact_id":"JANUS-KUSTO-EMODNET-ALL8-E1-DIRECT-SOURCE-REFERENCES-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "correction":"JANUS-KUSTO-EMODNET-ALL8-WMS-INFORMATIVE-FLAG-CORRECTION-2026-09-23-v1.0",
 "wfs":WFS,
 "feature_type":TYPE,
 "bbox_halfwidth_deg":0.02,
 "targets":rows,
 "hit_target_ids":sorted(hits),
 "hit_target_count":len(hits),
 "depth_values_read":False,
 "composite_dtm_counted_as_independent_survey":False,
 "interpretation":"EXPLICIT_SOURCE_REFERENCE_FEATURES_FOUND" if hits else "NO_SOURCE_REFERENCE_FEATURES_IN_FROZEN_LOCAL_WINDOWS",
 "claim_ceiling":"EMODNET_SOURCE_REFERENCE_METADATA_ONLY"
}
p=OUT/"JANUS-KUSTO-EMODNET-ALL8-E1-DIRECT-SOURCE-REFERENCES-RUN-2026-09-23-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "hit_target_count":len(hits),
 "hit_target_ids":sorted(hits),
 "hits":{k:v["features"] for k,v in hits.items()},
 "depth_values_read":False,
 "interpretation":out["interpretation"]
},indent=2,ensure_ascii=False))
