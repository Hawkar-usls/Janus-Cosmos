#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,time
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-EMODNET-FOUR-BASIN-SOURCE-FIRST-PREREG-2026-09-23-v1.0.json").read_text())
REP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-EMODNET-FOUR-BASIN-WFS-SOURCE-QUERY-REPAIR-2026-09-23-v1.0.json").read_text())
WFS=PRE["services"]["wfs"]
UA={"User-Agent":"JANUS-KUSTO-EMODnet-four-basin-WFS-source/1.0","Accept-Encoding":"identity"}
S=requests.Session();S.headers.update(UA)

def get(params,timeout=120):
    last=None
    for a in range(5):
        try:
            r=S.get(WFS,params=params,timeout=timeout,allow_redirects=True)
            r.raise_for_status();return r
        except Exception as e:
            last=e;time.sleep(min(12,2**a))
    raise last

layers=["emodnet:source_references","emodnet:quality_index"]
probes={}
for p in PRE["fixed_probe_points"]:
    rows=[]
    pad=.005
    bbox=f'{p["lon"]-pad},{p["lat"]-pad},{p["lon"]+pad},{p["lat"]+pad},EPSG:4326'
    for layer in layers:
        r=get({"SERVICE":"WFS","REQUEST":"GetFeature","VERSION":"2.0.0",
               "TYPENAMES":layer,"BBOX":bbox,"COUNT":"50","OUTPUTFORMAT":"application/json"})
        raw=r.content
        try:j=r.json()
        except Exception:j={}
        feats=j.get("features") or []
        rows.append({
          "layer":layer,"url":r.url,"http_status":r.status_code,
          "bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),
          "feature_count":len(feats),
          "properties":[x.get("properties") or {} for x in feats]
        })
    probes[p["id"]]={"point":p,"bbox":bbox,"queries":rows}

# Explicit acquisition-level candidates are rows whose properties expose a CDI/survey/data-set identifier;
# do not treat composite DTM identities as independent acquisition evidence.
tokens=("cdi","survey","dataset","data_set","source")
summary={}
for basin,v in probes.items():
    cand=[]
    for q in v["queries"]:
        if q["layer"]!="emodnet:source_references":continue
        for props in q["properties"]:
            keys={str(k).lower():val for k,val in props.items()}
            looks=any(any(t in k for t in tokens) for k in keys)
            cand.append({"properties":props,"acquisition_metadata_fields_present":looks})
    summary[basin]=cand

out={
 "artifact_id":"JANUS-KUSTO-EMODNET-FOUR-BASIN-SOURCE-FIRST-RUN-2026-09-23-v1.1",
 "prereg":PRE["artifact_id"],"repair":REP["artifact_id"],
 "probes":probes,"source_reference_summary":summary,
 "bathymetric_depth_values_intentionally_not_requested":True,
 "composite_dtm_counted_as_independent_survey":False,
 "next_gate":"PRIMARY_SOURCE_AUDIT_OF_EXPLICIT_SURVEY_OR_CDI_IDENTITIES",
 "claim_ceiling":"EUROPEAN_SURVEY_LINEAGE_DISCOVERY_ONLY"
}
p=OUT/"JANUS-KUSTO-EMODNET-FOUR-BASIN-SOURCE-FIRST-RUN-2026-09-23-v1.1.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "counts":{k:{q["layer"]:q["feature_count"] for q in v["queries"]} for k,v in probes.items()},
 "source_properties":{k:[x["properties"] for x in vals] for k,vals in summary.items()},
 "depth_values_requested":False
},indent=2,ensure_ascii=False))
