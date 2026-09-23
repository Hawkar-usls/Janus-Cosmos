#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
import hashlib
import json
import math
import re
import time
from pathlib import Path
from urllib.parse import urljoin

import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"
OUT.mkdir(parents=True,exist_ok=True)
PRE=ROOT/"data/cousteau/JANUS-KUSTO-SHARK-ROCK-MEASURED-SUPPORT-PREREG-2026-09-23-v1.0.json"
pre=json.loads(PRE.read_text())
UA={"User-Agent":"JANUS-KUSTO-SharkRock-FNV-support/1.0"}
LAT=float(pre["target"]["lat"]); LON=float(pre["target"]["lon"])
SURVEYS=[x["survey_id"] for x in pre["survey_roles"]]
M_PER_DEG=111320.0
COS0=math.cos(math.radians(LAT))

def get(url,timeout=90):
    last=None
    for a in range(6):
        try:
            r=requests.get(url,headers=UA,timeout=timeout)
            if r.status_code==429:
                time.sleep(min(20,2**a));continue
            r.raise_for_status()
            return r
        except Exception as e:
            last=e
            if a==5: raise
            time.sleep(min(20,2**a))
    raise last

def hrefs(url):
    r=get(url)
    return sorted(set(urljoin(url,h) for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I)))

def list_fnv(base,maxdepth=6):
    seen=set();out=[]
    def walk(u,d):
        if u in seen or d>maxdepth:return
        seen.add(u)
        for x in hrefs(u):
            if not x.startswith(base) or x.rstrip("/")==u.rstrip("/"):continue
            if x.endswith("/"):walk(x,d+1)
            elif x.lower().endswith(".fnv"):out.append(x)
    walk(base,0)
    return sorted(set(out))

def xy(lon,lat):
    return (lon-LON)*M_PER_DEG*COS0,(lat-LAT)*M_PER_DEG

def point_segment_distance(ax,ay,bx,by):
    vx=bx-ax;vy=by-ay
    vv=vx*vx+vy*vy
    if vv<=1e-12:return math.hypot(ax,ay)
    t=max(0.0,min(1.0,-(ax*vx+ay*vy)/vv))
    px=ax+t*vx;py=ay+t*vy
    return math.hypot(px,py)

def parse_scan(text):
    best=None;best_row=None;counts={str(x):0 for x in pre["report_only_counts_m"]}
    rows=0
    epochs=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:
            epoch=float(p[6]); plon=float(p[15]);plat=float(p[16]);slon=float(p[17]);slat=float(p[18])
        except Exception:continue
        rows+=1
        ax,ay=xy(plon,plat);bx,by=xy(slon,slat)
        d=point_segment_distance(ax,ay,bx,by)
        for th in pre["report_only_counts_m"]:
            if d<=float(th): counts[str(th)]+=1
        if best is None or d<best:
            best=d
            best_row={"epoch":epoch,"port_lon":plon,"port_lat":plat,"stbd_lon":slon,"stbd_lat":slat}
        epochs.append(epoch)
    return rows,best,best_row,counts

def scan_url(url):
    r=get(url)
    text=r.text
    rows,best,best_row,counts=parse_scan(text)
    return {
      "url":url,
      "sha256":hashlib.sha256(r.content).hexdigest(),
      "bytes":len(r.content),
      "parseable_rows":rows,
      "min_swath_line_distance_m":best,
      "nearest_row":best_row,
      "rows_within_m":counts
    }

def tier(d):
    if d is None:return "NO_LOCAL_FNV_SUPPORT"
    if d<=100:return "DIRECT_FNV_SWATH_LINE_SUPPORT"
    if d<=250:return "NEAR_FNV_SWATH_LINE_SUPPORT"
    if d<=1000:return "REGIONAL_FNV_SUPPORT_ONLY"
    return "NO_LOCAL_FNV_SUPPORT"

results=[]
for sid in SURVEYS:
    base=f"https://data.ngdc.noaa.gov/platforms/ocean/ships/okeanos_explorer/{sid}/multibeam/data/"
    rec={"survey_id":sid,"base_url":base}
    try:
        files=list_fnv(base)
        rec["fnv_file_count"]=len(files)
        scans=[]
        with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
            futs=[ex.submit(scan_url,u) for u in files]
            for f in concurrent.futures.as_completed(futs):
                try: scans.append(f.result())
                except Exception as e: scans.append({"error":type(e).__name__+": "+str(e)})
        good=[x for x in scans if x.get("min_swath_line_distance_m") is not None]
        good.sort(key=lambda x:x["min_swath_line_distance_m"])
        d=good[0]["min_swath_line_distance_m"] if good else None
        rec["parseable_fnv_file_count"]=len(good)
        rec["min_swath_line_distance_m"]=d
        rec["support_tier"]=tier(d)
        rec["nearest_fnv_records"]=good[:10]
        rec["failed_fnv_count"]=sum(1 for x in scans if x.get("error"))
        agg={str(x):0 for x in pre["report_only_counts_m"]}
        for x in good:
            for k,v in (x.get("rows_within_m") or {}).items(): agg[k]+=int(v)
        rec["all_rows_within_m"]=agg
    except Exception as e:
        rec["error"]=type(e).__name__+": "+str(e)
        rec["support_tier"]="NO_LOCAL_FNV_SUPPORT"
    results.append(rec)

out={
 "artifact_id":"JANUS-KUSTO-SHARK-ROCK-MEASURED-SUPPORT-RUN-2026-09-23-v1.0",
 "prereg":pre["artifact_id"],
 "target":{"lat":LAT,"lon":LON},
 "depth_values_read":False,
 "detector_run":False,
 "survey_results":results,
 "all_frozen_surveys_direct_support":all(r.get("support_tier")=="DIRECT_FNV_SWATH_LINE_SUPPORT" for r in results),
 "hard_rule_note":"FNV support is navigation/swath-edge geometry only and is not a morphology result.",
 "next_gate":"TARGET_HIDDEN_MULTISCALE_DISCOVERY_IF_SOURCE_BINDING_SUFFICIENT",
 "claim_ceiling":"FNV_MEASURED_NAVIGATION_SUPPORT_ONLY"
}
path=OUT/"JANUS-KUSTO-SHARK-ROCK-MEASURED-SUPPORT-RUN-2026-09-23-v1.0.json"
path.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "surveys":[{"survey_id":r["survey_id"],"fnv_file_count":r.get("fnv_file_count"),"parseable_fnv_file_count":r.get("parseable_fnv_file_count"),"min_swath_line_distance_m":r.get("min_swath_line_distance_m"),"support_tier":r.get("support_tier"),"failed_fnv_count":r.get("failed_fnv_count")} for r in results],
 "all_frozen_surveys_direct_support":out["all_frozen_surveys_direct_support"],
 "depth_values_read":False
},indent=2))
