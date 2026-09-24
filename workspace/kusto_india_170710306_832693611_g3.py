#!/usr/bin/env python3
import json, hashlib
from pathlib import Path
import requests

LAT=17.0710306; LON=83.2693611
OUT=Path("workspace/kusto_global_groundtruth_out/india_g3"); OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-INDIA-G3/1.0"}
TRACK="https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/multibeam_dynamic/MapServer/0/query"
MOSAIC="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_mosaic/ImageServer"
S=requests.Session(); S.headers.update(UA)

def get_json(url,params):
    r=S.get(url,params=params,timeout=120)
    rec={"status":r.status_code,"url":r.url,"content_type":r.headers.get("content-type"),"text_head":r.text[:1000]}
    try:
        rec["json"]=r.json()
    except Exception:
        rec["json"]=None
    return rec

point=json.dumps({"x":LON,"y":LAT,"spatialReference":{"wkid":4326}},separators=(",",":"))

mosaic_query=get_json(MOSAIC+"/query",{
    "where":"1=1","geometry":point,"geometryType":"esriGeometryPoint","inSR":"4326",
    "spatialRel":"esriSpatialRelIntersects","outFields":"OBJECTID,Name,MinPS,MaxPS,LowPS,HighPS,CenterX,CenterY,ProductName,Tag",
    "returnGeometry":"false","f":"json"
})
mosaic_samples=get_json(MOSAIC+"/getSamples",{
    "geometry":point,"geometryType":"esriGeometryPoint","returnFirstValueOnly":"false",
    "outFields":"OBJECTID,Name,MinPS,MaxPS,CenterX,CenterY,ProductName,Tag","f":"json"
})

boxes=[]
for h in [0.01,0.05,0.1,0.25,0.5]:
    env=json.dumps({"xmin":LON-h,"ymin":LAT-h,"xmax":LON+h,"ymax":LAT+h,
                    "spatialReference":{"wkid":4326}},separators=(",",":"))
    q=get_json(TRACK,{
      "where":"1=1","geometry":env,"geometryType":"esriGeometryEnvelope","inSR":"4326",
      "spatialRel":"esriSpatialRelIntersects",
      "outFields":"SURVEY_ID,PLATFORM,SURVEY_YEAR,SOURCE,NGDC_ID,INSTRUMENT,TRACK_LENGTH,DOWNLOAD_URL",
      "returnGeometry":"false","returnDistinctValues":"true","f":"json"
    })
    js=q.get("json") or {}
    feats=js.get("features",[]) if isinstance(js,dict) else []
    attrs=[x.get("attributes",{}) for x in feats]
    boxes.append({"half_span_deg":h,"query":{k:v for k,v in q.items() if k!="json"},
                  "feature_count":len(attrs),"features":attrs})

mq=(mosaic_query.get("json") or {})
ms=(mosaic_samples.get("json") or {})
mosaic_features=[x.get("attributes",{}) for x in mq.get("features",[])] if isinstance(mq,dict) else []
samples=ms.get("samples",[]) if isinstance(ms,dict) else []

out={
  "artifact_id":"JANUS-KUSTO-INDIA-170710306N-832693611E-G3-NCEI-MULTIBEAM-COVERAGE-RUN-2026-09-25-v1.0",
  "target":{"lat":LAT,"lon":LON},
  "ncei_exact_target":{
    "mosaic_footprint_feature_count":len(mosaic_features),
    "mosaic_footprint_features":mosaic_features,
    "sample_count":len(samples),
    "samples":samples,
    "query_diagnostics":{k:v for k,v in mosaic_query.items() if k!="json"},
    "sample_diagnostics":{k:v for k,v in mosaic_samples.items() if k!="json"}
  },
  "trackline_boxes":boxes,
  "interpretation_rule":{
    "exact_multibeam_at_target": bool(mosaic_features or samples),
    "nearest_box_with_track": next((b["half_span_deg"] for b in boxes if b["feature_count"]>0),None),
    "scope":"NCEI_ARCHIVE_ONLY"
  }
}
raw=json.dumps(out,indent=2,allow_nan=False)
(OUT/"JANUS-KUSTO-INDIA-170710306N-832693611E-G3-NCEI-MULTIBEAM-COVERAGE-RUN-2026-09-25-v1.0.json").write_text(raw)
print(json.dumps({
  "exact_mosaic_feature_count":len(mosaic_features),
  "exact_sample_count":len(samples),
  "boxes":[{"half_span_deg":b["half_span_deg"],"feature_count":b["feature_count"],
            "survey_ids":[f.get("SURVEY_ID") for f in b["features"][:20]]} for b in boxes],
  "exact_multibeam_at_target":out["interpretation_rule"]["exact_multibeam_at_target"],
  "nearest_box_with_track":out["interpretation_rule"]["nearest_box_with_track"],
  "sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2))
