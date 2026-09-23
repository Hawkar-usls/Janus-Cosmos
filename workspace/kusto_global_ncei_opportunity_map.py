#!/usr/bin/env python3
import json, math, time, requests
from pathlib import Path
from collections import defaultdict
from shapely.geometry import shape, Point
from shapely.ops import unary_union
from shapely.strtree import STRtree
from pyproj import Geod

OUT=Path("workspace/kusto_global_anomaly_out"); OUT.mkdir(parents=True,exist_ok=True)
LAYER="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0"
QUERY=LAYER+"/query"
UA={"User-Agent":"JANUS-KUSTO-global-opportunity-map/1.0"}
GEOD=Geod(ellps="WGS84")
FIELDS="OBJECTID,NCEI_ID,SURVEY_ID,PLATFORM,SOURCE,INSTRUMENT,START_TIME,END_TIME,SURVEY_YEAR,DOWNLOAD_URL,SURVEY_AND_VERSION"

def get(url,params=None):
    last=None
    for a in range(7):
        try:
            r=requests.get(url,params=params,headers=UA,timeout=120)
            if r.status_code==429:
                time.sleep(min(30,2**a)); continue
            r.raise_for_status()
            d=r.json()
            if isinstance(d,dict) and "error" in d: raise RuntimeError(d["error"])
            return d
        except Exception as e:
            last=e
            if a==6: raise
            time.sleep(min(30,2**a))
    raise last

meta=get(LAYER,{"f":"json"})
oid=meta.get("objectIdField") or meta.get("objectIdFieldName") or "OBJECTID"
ids=get(QUERY,{"where":"1=1","returnIdsOnly":"true","f":"json"}).get("objectIds",[])
ids=sorted(ids)
print("object ids",len(ids),flush=True)

features=[]
chunk=500
for start in range(0,len(ids),chunk):
    batch=ids[start:start+chunk]
    d=get(QUERY,{
      "objectIds":",".join(map(str,batch)),
      "outFields":FIELDS,
      "returnGeometry":"true",
      "outSR":"4326",
      "f":"geojson"
    })
    features.extend(d.get("features",[]))
    if start%5000==0: print("fetched",min(start+chunk,len(ids)),"/",len(ids),flush=True)

by_survey=defaultdict(list)
attrs={}
for f in features:
    p=f.get("properties") or {}
    sid=p.get("SURVEY_ID")
    if not sid or not f.get("geometry"): continue
    try:
        g=shape(f["geometry"])
        if not g.is_valid: g=g.buffer(0)
        if g.is_empty: continue
    except Exception:
        continue
    by_survey[str(sid)].append(g)
    attrs.setdefault(str(sid),{
      "survey_id":str(sid),
      "platform":p.get("PLATFORM"),
      "instrument":p.get("INSTRUMENT"),
      "survey_year":p.get("SURVEY_YEAR"),
      "start_time":p.get("START_TIME"),
      "end_time":p.get("END_TIME"),
      "source":p.get("SOURCE"),
      "download_url":p.get("DOWNLOAD_URL"),
      "survey_and_version":p.get("SURVEY_AND_VERSION"),
      "ncei_id":p.get("NCEI_ID")
    })

surveys=[]
for sid,gs in by_survey.items():
    try:
        g=unary_union(gs)
        if not g.is_valid:g=g.buffer(0)
        if g.is_empty:continue
    except Exception:
        continue
    surveys.append({"sid":sid,"geom":g,"attr":attrs[sid]})
surveys.sort(key=lambda x:x["sid"])
print("surveys dissolved",len(surveys),flush=True)

geoms=[x["geom"] for x in surveys]
tree=STRtree(geoms)
geom_to_idx={id(g):i for i,g in enumerate(geoms)}

def norm(s):
    return "" if s is None else str(s).strip().upper()

def year(v):
    try:return int(v)
    except:return None

def geod_area_km2(g):
    try:
        area,_=GEOD.geometry_area_perimeter(g)
        return abs(area)/1e6
    except Exception:
        return 0.0

pairs=[]
seen=set()
for i,a in enumerate(surveys):
    cand=tree.query(a["geom"])
    # shapely version may return geometries or integer indexes
    idxs=[]
    for q in cand:
        if isinstance(q,(int,)):
            idxs.append(int(q))
        else:
            try:
                import numpy as np
                if isinstance(q,np.integer): idxs.append(int(q)); continue
            except Exception: pass
            j=geom_to_idx.get(id(q))
            if j is not None: idxs.append(j)
    for j in idxs:
        if j<=i: continue
        b=surveys[j]
        key=(a["sid"],b["sid"])
        if key in seen:continue
        seen.add(key)
        aa=a["attr"]; bb=b["attr"]
        if not aa.get("download_url") or not bb.get("download_url"): continue
        try:
            inter=a["geom"].intersection(b["geom"])
        except Exception:
            continue
        if inter.is_empty:continue
        area=geod_area_km2(inter)
        if area<0.25:continue
        c=inter.representative_point()
        ya,yb=year(aa.get("survey_year")),year(bb.get("survey_year"))
        ysep=abs(ya-yb) if ya is not None and yb is not None else 0
        pairs.append({
          "survey_a":aa,"survey_b":bb,
          "overlap_area_km2":area,
          "overlap_point":{"lon":float(c.x),"lat":float(c.y)},
          "different_instrument":norm(aa.get("instrument"))!=norm(bb.get("instrument")),
          "different_platform":norm(aa.get("platform"))!=norm(bb.get("platform")),
          "year_separation":ysep,
          "intersection_wkt":inter.wkt
        })
    if i%500==0: print("pairs scan",i,"/",len(surveys),"eligible",len(pairs),flush=True)

# Count surveys covering each pair's frozen representative point.
for p in pairs:
    pt=Point(p["overlap_point"]["lon"],p["overlap_point"]["lat"])
    cover=[]
    for k in tree.query(pt):
        if isinstance(k,(int,)):
            idx=int(k)
        else:
            try:
                import numpy as np
                if isinstance(k,np.integer): idx=int(k)
                else: idx=geom_to_idx.get(id(k),-1)
            except Exception:
                idx=geom_to_idx.get(id(k),-1)
        if idx>=0:
            try:
                if surveys[idx]["geom"].intersects(pt): cover.append(surveys[idx]["sid"])
            except: pass
    p["independent_survey_count_at_overlap_point"]=len(sorted(set(cover)))
    p["surveys_at_overlap_point"]=sorted(set(cover))

pairs.sort(key=lambda p:(
  -p["independent_survey_count_at_overlap_point"],
  -int(p["different_instrument"]),
  -int(p["different_platform"]),
  -p["year_separation"],
  -p["overlap_area_km2"],
  p["survey_a"]["survey_id"],
  p["survey_b"]["survey_id"]
))

# Spatial deduplication of top regions by representative-point distance >=5 km.
def hav(a,b):
    R=6371.0088
    p1,p2=math.radians(a["lat"]),math.radians(b["lat"])
    dp=p2-p1; dl=math.radians(b["lon"]-a["lon"])
    q=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(min(1,math.sqrt(q)))

top_pairs=pairs[:500]
regions=[]
for p in top_pairs:
    pt=p["overlap_point"]
    if any(hav(pt,r["overlap_point"])<5.0 for r in regions): continue
    rr={k:v for k,v in p.items() if k!="intersection_wkt"}
    regions.append(rr)
    if len(regions)>=100:break

out={
 "artifact_id":"JANUS-KUSTO-GLOBAL-NCEI-MULTISURVEY-OPPORTUNITY-MAP-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-GLOBAL-BLIND-ANOMALY-HUNT-PREREG-2026-09-23-v1.0.json",
 "ranking_addendum":"data/cousteau/JANUS-KUSTO-GLOBAL-OPPORTUNITY-MAP-RANKING-ADDENDUM-2026-09-23-v1.0.json",
 "depth_values_read":False,
 "layer_metadata":{"name":meta.get("name"),"maxRecordCount":meta.get("maxRecordCount"),"objectIdField":oid},
 "raw_feature_count":len(features),
 "dissolved_survey_count":len(surveys),
 "eligible_overlap_pair_count":len(pairs),
 "top_500_pairs":[{k:v for k,v in p.items() if k!="intersection_wkt"} for p in top_pairs],
 "stage1_regions":regions,
 "stage1_region_count":len(regions),
 "claim_ceiling":"DEPTH_BLIND_GLOBAL_MULTI_SURVEY_OPPORTUNITY_RANKING_ONLY"
}
p=OUT/"JANUS-KUSTO-GLOBAL-NCEI-MULTISURVEY-OPPORTUNITY-MAP-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "raw_feature_count":out["raw_feature_count"],
 "dissolved_survey_count":out["dissolved_survey_count"],
 "eligible_overlap_pair_count":out["eligible_overlap_pair_count"],
 "stage1_region_count":out["stage1_region_count"],
 "top10":regions[:10]
},indent=2))
