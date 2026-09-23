#!/usr/bin/env python3
from __future__ import annotations
import json,re,time
from pathlib import Path
from urllib.parse import urljoin,urlparse
import requests
from shapely.geometry import shape
from shapely.ops import unary_union

OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
ARC="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
UA={"User-Agent":"JANUS-KUSTO-unseen-reference-replacement/1.0"}
EXCLUDED={"PF0501","EX2105","CD169","KN192-07","KNOX15RR"}

def get(u,params=None):
    last=None
    for i in range(7):
        try:
            r=requests.get(u,params=params,headers=UA,timeout=90)
            if r.status_code==429:
                time.sleep(min(30,2**i)); continue
            r.raise_for_status(); return r
        except Exception as e:
            last=e
            if i==6: raise
            time.sleep(min(30,2**i))
    raise last

def query_geojson(params):
    p=dict(params); p["f"]="geojson"
    d=get(ARC,p).json()
    if "error" in d: raise RuntimeError(d["error"])
    return d

pf=query_geojson({
 "where":"SURVEY_ID='PF0501'",
 "outFields":"SURVEY_ID,PLATFORM,INSTRUMENT,SURVEY_YEAR,START_TIME,END_TIME,SOURCE,DOWNLOAD_URL",
 "returnGeometry":"true","outSR":"4326"
})
pfg=[shape(f["geometry"]) for f in pf.get("features",[]) if f.get("geometry")]
if not pfg: raise RuntimeError("PF0501 footprint missing")
pfgeom=unary_union(pfg)
west,south,east,north=pfgeom.bounds

# Paginated candidate inventory over PF0501 bounding envelope.
features=[]
offset=0
while True:
    d=query_geojson({
      "where":"1=1",
      "geometry":f"{west},{south},{east},{north}",
      "geometryType":"esriGeometryEnvelope","inSR":"4326",
      "spatialRel":"esriSpatialRelIntersects",
      "outFields":"SURVEY_ID,PLATFORM,INSTRUMENT,SURVEY_YEAR,START_TIME,END_TIME,SOURCE,DOWNLOAD_URL",
      "returnGeometry":"true","outSR":"4326",
      "resultOffset":offset,"resultRecordCount":1000
    })
    fs=d.get("features",[])
    features.extend(fs)
    if len(fs)<1000: break
    offset += len(fs)
    if offset>20000: raise RuntimeError("unexpected candidate pagination >20000")

# Merge repeated features by survey.
by={}
for f in features:
    pr=f.get("properties",{}); sid=pr.get("SURVEY_ID")
    if not sid or not f.get("geometry"): continue
    by.setdefault(sid,{"properties":pr,"geoms":[]})["geoms"].append(shape(f["geometry"]))

ranked=[]
for sid,item in by.items():
    if sid in EXCLUDED: continue
    year=item["properties"].get("SURVEY_YEAR")
    try: year=int(year)
    except: continue
    if year<2006: continue
    g=unary_union(item["geoms"])
    inter=pfgeom.intersection(g)
    if inter.is_empty or inter.area<=0: continue
    ranked.append({
      "survey_id":sid,
      "survey_year":year,
      "platform":item["properties"].get("PLATFORM"),
      "instrument":item["properties"].get("INSTRUMENT"),
      "source":item["properties"].get("SOURCE"),
      "download_url":item["properties"].get("DOWNLOAD_URL"),
      "intersection_area_deg2":float(inter.area),
      "intersection_bounds":list(inter.bounds),
      "feature_count":len(item["geoms"])
    })
ranked.sort(key=lambda x:(-x["intersection_area_deg2"],x["survey_id"]))

def links(base):
    try:r=get(base)
    except:return []
    if "text/html" not in r.headers.get("Content-Type","").lower() and "<a " not in r.text.lower():return []
    return sorted(set(urljoin(base,h) for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I)))

def route_probe(rec,max_pages=80):
    roots=[]
    report=rec.get("download_url")
    if report:
        try:
            for u in links(report):
                if "data.ngdc.noaa.gov/" in u and (u.endswith("/") or "/multibeam/" in u):
                    roots.append(u)
        except:pass
    # Deterministic URL candidates derived from report path when discoverable.
    if report:
        parsed=urlparse(report)
        m=re.search(r'/ships/([^/]+)/([^/_]+)_mb\.html',parsed.path,re.I)
        if m:
            ship=m.group(1); sid=m.group(2)
            roots += [
              f"https://data.ngdc.noaa.gov/platforms/ocean/ships/{ship}/{sid}/multibeam/data/version1/MB/",
              f"https://data.ngdc.noaa.gov/platforms/ocean/ships/{ship}/{sid}/multibeam/data/version2/MB/"
            ]
    roots=sorted(set(roots))
    seen=set();queue=[(u,0) for u in roots]; pages=0
    fnv=[];fbt=[];dirs=[]
    while queue and pages<max_pages and not (fnv and fbt):
        u,depth=queue.pop(0)
        if u in seen or depth>4: continue
        seen.add(u);pages+=1
        try:ls=links(u)
        except:continue
        for x in ls:
            if not x.startswith("https://data.ngdc.noaa.gov/"): continue
            low=x.lower()
            if low.endswith(".fnv"):fnv.append(x)
            elif low.endswith(".fbt"):fbt.append(x)
            elif x.endswith("/") and x not in seen:
                b=urlparse(x).path.lower()
                if any(k in b for k in ["/mb/","generated","em","multibeam","data/version"]):
                    queue.append((x,depth+1));dirs.append(x)
    return {
      "roots_probed":roots,
      "pages_probed":pages,
      "fnv_found":len(fnv),
      "fbt_found":len(fbt),
      "sample_fnv":sorted(set(fnv))[:3],
      "sample_fbt":sorted(set(fbt))[:3],
      "public_fnv_fbt_pass":bool(fnv and fbt)
    }

selected=None
route_results=[]
# Rank order is frozen. Probe enough candidates until first public FNV+FBT route.
for idx,rec in enumerate(ranked,1):
    rr=route_probe(rec)
    row={"rank":idx,**rec,"route":rr}
    route_results.append(row)
    print("rank",idx,rec["survey_id"],"area",rec["intersection_area_deg2"],"route",rr["public_fnv_fbt_pass"],flush=True)
    if rr["public_fnv_fbt_pass"]:
        selected=row
        break
    if idx>=30:
        break

out={
 "artifact_id":"JANUS-KUSTO-UNSEEN-PF0501-REFERENCE-REPLACEMENT-SEARCH-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-UNSEEN-PF0501-REFERENCE-REPLACEMENT-SEARCH-PREREG-2026-09-23-v1.0.json",
 "depth_values_read":False,
 "PF0501_bounds":list(pfgeom.bounds),
 "candidate_footprint_features_in_envelope":len(features),
 "eligible_exact_intersection_surveys":len(ranked),
 "top_exact_intersections":ranked[:50],
 "route_probe_results":route_results,
 "selected_reference":selected,
 "candidate_found":selected is not None,
 "next_gate":"FREEZE_EXACT_PF0501_X_SELECTED_REFERENCE_FNV_OVERLAP" if selected else "NO_PUBLIC_FNV_FBT_REFERENCE_IN_FIRST_30_EXACT_OVERLAPS",
 "claim_ceiling":"UNSEEN_REFERENCE_DISCOVERY_BY_COVERAGE_ONLY"
}
p=OUT/"JANUS-KUSTO-UNSEEN-PF0501-REFERENCE-REPLACEMENT-SEARCH-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "candidate_footprint_features_in_envelope":out["candidate_footprint_features_in_envelope"],
 "eligible_exact_intersection_surveys":out["eligible_exact_intersection_surveys"],
 "selected_reference":selected,
 "top10":[{k:x[k] for k in ["survey_id","survey_year","platform","instrument","intersection_area_deg2","intersection_bounds"]} for x in ranked[:10]]
},indent=2))
