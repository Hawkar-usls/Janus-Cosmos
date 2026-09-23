#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHELLHARBOUR-S0C-WMS-ROUTE-FREEZE-2026-09-23-v1.0.json").read_text())
END=PRE["endpoints"]["wfs"]
UA={"User-Agent":"JANUS-KUSTO-Shellharbour-S0C-geometry/1.0"}
S=requests.Session();S.headers.update(UA)

def fetch(role,layer):
    params={"SERVICE":"WFS","REQUEST":"GetFeature","VERSION":"2.0.0","TYPENAMES":layer,
            "COUNT":"1","OUTPUTFORMAT":"application/json","SRSNAME":"EPSG:32756"}
    r=S.get(END,params=params,timeout=120,allow_redirects=True);r.raise_for_status()
    raw=r.content;j=r.json()
    feats=j.get("features") or []
    if len(feats)!=1:raise RuntimeError(f"{role}: expected one coverage polygon feature, got {len(feats)}")
    geom=feats[0].get("geometry")
    if not geom:raise RuntimeError(f"{role}: geometry missing")
    coords=[]
    def walk(x):
        if isinstance(x,list) and len(x)>=2 and all(isinstance(v,(int,float)) for v in x[:2]):
            coords.append((float(x[0]),float(x[1])))
        elif isinstance(x,list):
            for q in x:walk(q)
    walk(geom.get("coordinates"))
    if not coords:raise RuntimeError(f"{role}: no coordinates")
    xs=[x for x,y in coords];ys=[y for x,y in coords]
    return {
      "layer":layer,"url":r.url,"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw),
      "geometry_type":geom.get("type"),"bbox_epsg32756":[min(xs),min(ys),max(xs),max(ys)],
      "properties":feats[0].get("properties") or {}
    }

bindings={role:fetch(role,layer) for role,layer in PRE["frozen_layers"].items()}
# Freeze the intersection envelope so both channels are evaluated only in common survey support.
b=bindings["BATHYMETRY_2M"]["bbox_epsg32756"];k=bindings["BACKSCATTER_5M"]["bbox_epsg32756"]
inter=[max(b[0],k[0]),max(b[1],k[1]),min(b[2],k[2]),min(b[3],k[3])]
if not (inter[0]<inter[2] and inter[1]<inter[3]):raise RuntimeError("no common projected bounding envelope")
width_m=inter[2]-inter[0];height_m=inter[3]-inter[1]
back_px=[int(round(width_m/5.0)),int(round(height_m/5.0))]
out={
 "artifact_id":"JANUS-KUSTO-SHELLHARBOUR-S0C-WMS-ROUTE-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],"bindings":bindings,
 "common_envelope_epsg32756":inter,
 "common_envelope_width_m":width_m,"common_envelope_height_m":height_m,
 "frozen_backscatter_wms_request":{
   "service":"WMS","version":"1.1.1","request":"GetMap",
   "layers":PRE["frozen_layers"]["BACKSCATTER_5M"],"srs":"EPSG:32756","bbox":inter,
   "width":back_px[0],"height":back_px[1],"format":"image/geotiff"
 },
 "wms_getmap_called":False,"raster_values_read":False,
 "next_gate":"S1_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE",
 "claim_ceiling":"TRANSPORT_ROUTE_AND_GEOMETRY_FREEZE_ONLY"
}
p=OUT/"JANUS-KUSTO-SHELLHARBOUR-S0C-WMS-ROUTE-RUN-2026-09-23-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({"artifact_id":out["artifact_id"],"common_envelope":inter,"backscatter_pixels":back_px,
 "wms_getmap_called":False,"raster_values_read":False},indent=2))
