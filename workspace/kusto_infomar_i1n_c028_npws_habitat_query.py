#!/usr/bin/env python3
from __future__ import annotations
import json, math, hashlib
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1N-C028-NPWS-REEF-AND-COMMUNITY-POINT-QUERY-PREREG-2026-09-25-v1.0.json").read_text())
T=PRE["target"]; lon=float(T["lon"]); lat=float(T["lat"])
BASE=PRE["service"]
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-I1N-C028/1.0"})

def q(layer,geometry,geometry_type,spatial_rel="esriSpatialRelIntersects"):
    params={
      "where":"1=1","geometry":json.dumps(geometry,separators=(",",":")),
      "geometryType":geometry_type,"inSR":"4326","spatialRel":spatial_rel,
      "outFields":"*","returnGeometry":"true","outSR":"4326","f":"json"
    }
    r=S.get(f"{BASE}/{layer}/query",params=params,timeout=180)
    r.raise_for_status(); j=r.json()
    if "error" in j: raise RuntimeError(j["error"])
    return j,r.url

pt={"x":lon,"y":lat,"spatialReference":{"wkid":4326}}
results={}
for name,layer in PRE["layers"].items():
    exact,url=q(layer,pt,"esriGeometryPoint")
    exact_features=exact.get("features",[]) or []
    # Secondary nearby context only. 50m approximate WGS84 envelope.
    m=50.0
    dlat=m/110540.0; dlon=m/(111320.0*math.cos(math.radians(lat)))
    env={"xmin":lon-dlon,"ymin":lat-dlat,"xmax":lon+dlon,"ymax":lat+dlat,"spatialReference":{"wkid":4326}}
    near,nurl=q(layer,env,"esriGeometryEnvelope")
    results[name]={
      "layer_id":layer,
      "exact_query_url":url,
      "exact_feature_count":len(exact_features),
      "exact_features":exact_features,
      "nearby_50m_query_url":nurl,
      "nearby_50m_feature_count":len(near.get("features",[]) or []),
      "nearby_50m_features":near.get("features",[]) or []
    }

rh=results["1170_Reefs"]["exact_feature_count"]>0
ch=results["Marine_Community_Types"]["exact_feature_count"]>0
cls="BOTH_HIT" if rh and ch else "REEF_ONLY" if rh else "COMMUNITY_ONLY" if ch else "NONE"

out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I1N-C028-NPWS-REEF-AND-COMMUNITY-POINT-QUERY-RUN-2026-09-25-v1.0",
 "prereg":PRE["artifact_id"],"target":T,
 "source":{"authority":"National Parks & Wildlife Service Ireland","service":BASE},
 "results":results,"primary_result_class":cls,
 "candidate_coordinate_changed":False,"candidate_radius_changed":False,
 "claim_ceiling":"EXACT_OFFICIAL_NPWS_POINT_IN_POLYGON_CHARACTERIZATION_ONLY"
}
raw=json.dumps(out,indent=2,ensure_ascii=False)
p=OUT/"JANUS-KUSTO-INFOMAR-I1N-C028-NPWS-REEF-AND-COMMUNITY-POINT-QUERY-RUN-2026-09-25-v1.0.json"
p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"primary_result_class":cls,
 "reef_exact_count":results["1170_Reefs"]["exact_feature_count"],
 "reef_attributes":[f.get("attributes") for f in results["1170_Reefs"]["exact_features"]],
 "community_exact_count":results["Marine_Community_Types"]["exact_feature_count"],
 "community_attributes":[f.get("attributes") for f in results["Marine_Community_Types"]["exact_features"]],
 "reef_nearby_50m_count":results["1170_Reefs"]["nearby_50m_feature_count"],
 "community_nearby_50m_count":results["Marine_Community_Types"]["nearby_50m_feature_count"],
 "output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2,ensure_ascii=False))
