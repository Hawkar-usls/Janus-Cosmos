#!/usr/bin/env python3
import requests,re,json,hashlib
from pathlib import Path

BASE="https://geodataindia.gov.in/"
LAT=17.0710306; LON=83.2693611
OUT=Path("workspace/kusto_global_groundtruth_out/ngdr_sm222_wfs"); OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-NGDR-SM222/5.0","Referer":BASE+"guestuser"})
h=S.get(BASE+"guestuser",timeout=30); h.raise_for_status()

queries=[]
for typename in ["baseline_data_gis_view","cite:baseline_data_gis_view"]:
  params={"service":"WFS","version":"1.0.0","request":"GetFeature","typeName":typename,
          "outputFormat":"application/json","CQL_FILTER":"baseid=5192"}
  try:
    r=S.get(BASE+"guestuser/wmsurl128/",params=params,timeout=30)
    rec={"typename":typename,"status":r.status_code,"url":r.url,"content_type":r.headers.get("content-type"),
         "bytes":len(r.content),"text_head":r.text[:3000]}
    try: rec["json"]=r.json()
    except: rec["json"]=None
    queries.append(rec)
  except Exception as e:
    queries.append({"typename":typename,"error":repr(e)})

# If WFS proxy doesn't expose, use WMS GetFeatureInfo exactly at target with target centered in tiny bbox.
half=.001
for layer in ["baseline_data_gis_view","cite:baseline_data_gis_view"]:
  params={"SERVICE":"WMS","VERSION":"1.1.1","REQUEST":"GetFeatureInfo","LAYERS":layer,"QUERY_LAYERS":layer,
          "SRS":"EPSG:4326","BBOX":f"{LON-half},{LAT-half},{LON+half},{LAT+half}",
          "WIDTH":"101","HEIGHT":"101","X":"50","Y":"50","INFO_FORMAT":"application/json",
          "CQL_FILTER":"baseid=5192","FEATURE_COUNT":"10","FORMAT":"image/png"}
  try:
    r=S.get(BASE+"guestuser/wmsurl128/",params=params,timeout=30)
    rec={"layer":layer,"status":r.status_code,"url":r.url,"content_type":r.headers.get("content-type"),
         "bytes":len(r.content),"text_head":r.text[:3000]}
    try: rec["json"]=r.json()
    except: rec["json"]=None
    queries.append(rec)
  except Exception as e:
    queries.append({"layer":layer,"error":repr(e)})

# Basic GeoJSON point-in-polygon if any WFS query succeeded.
def pip_ring(x,y,ring):
    inside=False
    n=len(ring)
    j=n-1
    for i in range(n):
        xi,yi=ring[i][0],ring[i][1]; xj,yj=ring[j][0],ring[j][1]
        if ((yi>y)!=(yj>y)):
            xx=(xj-xi)*(y-yi)/(yj-yi)+xi
            if x < xx: inside=not inside
        j=i
    return inside

def geom_contains(g):
    if not g:return False
    typ=g.get("type"); c=g.get("coordinates")
    if typ=="Polygon":
        return bool(c and pip_ring(LON,LAT,c[0]) and not any(pip_ring(LON,LAT,h) for h in c[1:]))
    if typ=="MultiPolygon":
        return any(bool(p and pip_ring(LON,LAT,p[0]) and not any(pip_ring(LON,LAT,h) for h in p[1:])) for p in c or [])
    if typ=="GeometryCollection":
        return any(geom_contains(x) for x in g.get("geometries",[]))
    return False

features=[]
for q in queries:
    js=q.get("json")
    if isinstance(js,dict) and isinstance(js.get("features"),list):
        for f in js["features"]:
            features.append({"properties":f.get("properties"),"geometry":f.get("geometry"),
                             "contains_target":geom_contains(f.get("geometry"))})

out={"artifact_id":"JANUS-KUSTO-NGDR-SM222-WFS-POINT-IN-POLYGON-2026-09-25-v1.0",
     "target":{"lat":LAT,"lon":LON},"baseid":5192,"queries":queries,"features":features,
     "any_geometry_contains_target":any(f["contains_target"] for f in features)}
raw=json.dumps(out,indent=2,ensure_ascii=False)
(OUT/"wfs_point_in_polygon.json").write_text(raw,encoding="utf-8")
print(json.dumps({"status":"PASS","query_summary":[{k:q.get(k) for k in ["typename","layer","status","content_type","bytes","text_head","error"] if k in q} for q in queries],
 "feature_count":len(features),"contains":[f["contains_target"] for f in features],
 "any_geometry_contains_target":out["any_geometry_contains_target"],
 "feature_properties":[f["properties"] for f in features[:10]],
 "sha256":hashlib.sha256(raw.encode()).hexdigest()},indent=2,ensure_ascii=False))
