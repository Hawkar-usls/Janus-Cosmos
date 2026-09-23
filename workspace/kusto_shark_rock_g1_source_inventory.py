#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures, json, re, time
from pathlib import Path
from urllib.parse import urljoin, urlparse
import requests
from shapely.geometry import LineString, box

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHARK-ROCK-BLIND-TILE-PREREG-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-shark-rock-G1/1.0"}
LAYER="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
T=PRE["blind_tile"]; TILE=box(T["lon_min"],T["lat_min"],T["lon_max"],T["lat_max"])

def get(url,params=None):
    last=None
    for a in range(7):
        try:
            r=requests.get(url,params=params,headers=UA,timeout=120)
            if r.status_code==429:
                time.sleep(min(30,2**a)); continue
            r.raise_for_status(); return r
        except Exception as e:
            last=e
            if a==6: raise
            time.sleep(min(30,2**a))
    raise last

def hrefs(u):
    return sorted(set(urljoin(u,h) for h in re.findall(r'href=["\']([^"\']+)["\']',get(u).text,re.I)))

def derive_base(download_url):
    p=urlparse(download_url); parts=[x for x in p.path.split("/") if x]
    if len(parts)<2: raise RuntimeError(f"cannot derive base from {download_url}")
    ship=parts[-2]; survey=parts[-1]
    if survey.endswith("_mb.html"): survey=survey[:-8]
    else: survey=survey.split(".")[0]
    return f"https://data.ngdc.noaa.gov/platforms/ocean/ships/{ship}/{survey}/multibeam/data/"

def list_fnv(base,maxdepth=4):
    seen=set(); out=[]
    def walk(u,d):
        if u in seen or d>maxdepth:return
        seen.add(u)
        for x in hrefs(u):
            if not x.startswith(base) or x.rstrip("/")==u.rstrip("/"): continue
            if x.endswith("/"): walk(x,d+1)
            elif x.lower().endswith(".fnv"): out.append(x)
    walk(base,0)
    return sorted(set(out))

def parse_fnv(text):
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19: continue
        try:
            epoch=float(p[6]); plon=float(p[15]); plat=float(p[16]); slon=float(p[17]); slat=float(p[18])
        except: continue
        rows.append((epoch,plon,plat,slon,slat))
    return rows

def intersects_tile(url):
    rows=parse_fnv(get(url).text)
    hits=0; t0=None; t1=None
    for epoch,plon,plat,slon,slat in rows:
        seg=LineString([(plon,plat),(slon,slat)])
        if seg.intersects(TILE):
            hits+=1
            t0=epoch if t0 is None else min(t0,epoch)
            t1=epoch if t1 is None else max(t1,epoch)
    return {"fnv":url,"rows":len(rows),"tile_intersecting_segments":hits,"epoch_min":t0,"epoch_max":t1} if hits else None

def survey_features(sid):
    d=get(LAYER,params={
      "where":f"SURVEY_ID='{sid}'",
      "outFields":"SURVEY_ID,PLATFORM,INSTRUMENT,SURVEY_YEAR,START_TIME,END_TIME,DOWNLOAD_URL,SOURCE,NCEI_ID",
      "returnGeometry":"false","f":"json"
    }).json()
    return [f.get("attributes") or {} for f in d.get("features",[])]

surveys={}
for sid in PRE["survey_ids"]:
    feats=survey_features(sid)
    urls=sorted({f.get("DOWNLOAD_URL") for f in feats if f.get("DOWNLOAD_URL")})
    bases=[]
    for u in urls:
        try:bases.append(derive_base(u))
        except Exception:pass
    fnvs=[]
    for b in sorted(set(bases)):
        try:fnvs.extend(list_fnv(b))
        except Exception as e:
            fnvs.append({"_inventory_error":type(e).__name__+": "+str(e),"_base":b})
    actual=[u for u in fnvs if isinstance(u,str)]
    hit=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        for x in ex.map(intersects_tile,actual):
            if x: hit.append(x)
    surveys[sid]={
      "metadata_features":feats,
      "archive_bases":sorted(set(bases)),
      "fnv_total_discovered":len(actual),
      "fnv_tile_hits":sorted(hit,key=lambda x:x["fnv"]),
      "fnv_tile_hit_count":len(hit),
      "inventory_errors":[x for x in fnvs if isinstance(x,dict)]
    }
    print(sid, "fnv",len(actual),"tile_hits",len(hit),flush=True)

out={
 "artifact_id":"JANUS-KUSTO-SHARK-ROCK-G1-SOURCE-INVENTORY-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "tile":PRE["blind_tile"],
 "truth_coordinate_read":False,
 "bathymetric_depth_read":False,
 "video_or_rov_truth_read":False,
 "surveys":surveys,
 "eligible_surveys_for_blind_depth_stage":[sid for sid,v in surveys.items() if v["fnv_tile_hit_count"]>0],
 "hard_rules_preserved":["NO_TRUTH_COORDINATE_READ","NO_DEPTH_READ","NO_ROV_TRUTH_READ","NAV_GEOMETRY_ONLY"],
 "next_gate":"G2_FREEZE_NATIVE_SAMPLING_AND_BLIND_MORPHOLOGY_INPUTS",
 "claim_ceiling":"SOURCE_AND_NAVIGATION_SUPPORT_ONLY"
}
p=OUT/"JANUS-KUSTO-SHARK-ROCK-G1-SOURCE-INVENTORY-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "eligible":out["eligible_surveys_for_blind_depth_stage"],
 "counts":{k:{"fnv_total":v["fnv_total_discovered"],"tile_hits":v["fnv_tile_hit_count"]} for k,v in surveys.items()},
 "truth_coordinate_read":False,
 "bathymetric_depth_read":False
},indent=2))
