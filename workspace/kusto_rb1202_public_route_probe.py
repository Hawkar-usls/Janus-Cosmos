#!/usr/bin/env python3
import json,re,requests,time
from urllib.parse import urljoin
from pathlib import Path
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-RB1202-route-probe/1.0"}
ROOTS=[
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/ronald_h_brown/RB1202/multibeam/data/version1/MB/",
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/ronald_h_brown/RB1202/multibeam/data/version1/MB/generated/",
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/ronald_h_brown/RB1202/multibeam/data/version2/MB/"
]
def get(u):
    for i in range(4):
        try:
            r=requests.get(u,headers=UA,timeout=30)
            if r.status_code==429:time.sleep(2**i);continue
            return r
        except Exception:
            if i==3:raise
            time.sleep(2**i)
def ls(u):
    r=get(u)
    if not r.ok:return {"url":u,"status":r.status_code,"links":[]}
    links=sorted(set(urljoin(u,h) for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I)))
    return {"url":u,"status":r.status_code,"links":links}
pages=[];queue=[(u,0) for u in ROOTS];seen=set();fnv=[];fbt=[];dirs=[]
while queue and len(seen)<40 and not(fnv and fbt):
    u,d=queue.pop(0)
    if u in seen or d>3:continue
    seen.add(u);p=ls(u);pages.append({"url":u,"status":p["status"],"link_count":len(p["links"])})
    for x in p["links"]:
        if "data.ngdc.noaa.gov/" not in x:continue
        low=x.lower()
        if low.endswith(".fnv"):fnv.append(x)
        elif low.endswith(".fbt"):fbt.append(x)
        elif x.endswith("/") and "RB1202" in x and x not in seen:
            queue.append((x,d+1));dirs.append(x)
out={"artifact_id":"JANUS-KUSTO-RB1202-PUBLIC-FNV-FBT-ROUTE-PROBE-2026-09-23-v1.0",
     "depth_values_read":False,"survey":"RB1202","pages":pages,
     "fnv_count":len(set(fnv)),"fbt_count":len(set(fbt)),
     "sample_fnv":sorted(set(fnv))[:5],"sample_fbt":sorted(set(fbt))[:5],
     "route_pass":bool(fnv and fbt),"claim_ceiling":"PUBLIC_ROUTE_AVAILABILITY_ONLY"}
(OUT/"JANUS-KUSTO-RB1202-PUBLIC-FNV-FBT-ROUTE-PROBE-2026-09-23-v1.0.json").write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
