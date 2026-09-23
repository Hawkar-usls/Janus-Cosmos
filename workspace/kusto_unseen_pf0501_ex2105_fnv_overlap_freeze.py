#!/usr/bin/env python3
import json,re,requests,time,concurrent.futures
from urllib.parse import urljoin
from pathlib import Path
from shapely.geometry import shape, LineString
from shapely.ops import unary_union

OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-unseen-fnv-overlap-freeze/1.0"}
FOOT="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
ROOTS={
 "PF0501":"https://data.ngdc.noaa.gov/platforms/ocean/ships/pathfinder/PF0501/multibeam/data/version1/MB/generated/",
 "EX2105":"https://data.ngdc.noaa.gov/platforms/ocean/ships/okeanos_explorer/EX2105/multibeam/data/version1/MB/em304/"
}

def get(u):
    last=None
    for i in range(7):
        try:
            r=requests.get(u,headers=UA,timeout=90)
            if r.status_code==429:
                time.sleep(min(30,2**i)); continue
            r.raise_for_status(); return r
        except Exception as e:
            last=e
            if i==6: raise
            time.sleep(min(30,2**i))
    raise last

def hrefs(u):
    r=get(u)
    return sorted(set(urljoin(u,h) for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I)))

def collect_fnv(root):
    out=[]
    for u in hrefs(root):
        if u.lower().endswith(".fnv"): out.append(u)
        elif u.endswith("/") and u.startswith(root) and u.rstrip("/")!=root.rstrip("/"):
            # only recurse one generated-like level
            for v in hrefs(u):
                if v.lower().endswith(".fnv"): out.append(v)
    return sorted(set(out))

def footprint(sid):
    p={"where":f"SURVEY_ID='{sid}'","outFields":"SURVEY_ID","returnGeometry":"true","outSR":"4326","f":"geojson"}
    d=get(FOOT+"?"+requests.compat.urlencode(p)).json()
    gs=[shape(f["geometry"]) for f in d.get("features",[]) if f.get("geometry")]
    return unary_union(gs)

inter=footprint("PF0501").intersection(footprint("EX2105"))
if inter.is_empty: raise RuntimeError("exact footprint intersection unexpectedly empty")

def parse_fnv(url):
    rows=[]; r=get(url)
    for line in r.text.splitlines():
        p=line.split()
        if len(p)<19: continue
        try:
            epoch=float(p[6]); portlon=float(p[15]); portlat=float(p[16]); stbdlon=float(p[17]); stbdlat=float(p[18])
        except: continue
        seg=LineString([(portlon,portlat),(stbdlon,stbdlat)])
        if seg.intersects(inter):
            rows.append({"epoch":epoch,"portlon":portlon,"portlat":portlat,"stbdlon":stbdlon,"stbdlat":stbdlat})
    return rows

survey={}
for sid,root in ROOTS.items():
    files=collect_fnv(root)
    hits=[]; file_hits=[]
    def one(u):
        rr=parse_fnv(u)
        return u,rr
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        futs=[ex.submit(one,u) for u in files]
        for k,fut in enumerate(concurrent.futures.as_completed(futs),1):
            u,rr=fut.result()
            if rr:
                hits.extend(rr)
                file_hits.append({
                  "fnv_url":u,
                  "fnv_file":u.rsplit("/",1)[-1],
                  "hit_rows":len(rr),
                  "epoch_min":min(x["epoch"] for x in rr),
                  "epoch_max":max(x["epoch"] for x in rr)
                })
            if k%50==0: print(sid,"scanned",k,"/",len(files),flush=True)
    file_hits.sort(key=lambda x:x["fnv_file"])
    survey[sid]={
      "fnv_files_discovered":len(files),
      "files_with_overlap_rows":len(file_hits),
      "overlap_row_count":len(hits),
      "file_hits":file_hits,
      "epoch_min":min((x["epoch"] for x in hits),default=None),
      "epoch_max":max((x["epoch"] for x in hits),default=None),
      "minimum_overlap_samples_pass":len(hits)>=100
    }

out={
 "artifact_id":"JANUS-KUSTO-UNSEEN-PF0501-X-EX2105-FNV-OVERLAP-FREEZE-2026-09-23-v1.0",
 "implementation":"v1.1_parallel_http_same_frozen_science",
 "parent_prereg":"data/cousteau/JANUS-KUSTO-UNSEEN-MULTIBEAM-PF0501-X-EX2105-SOURCE-INVENTORY-PREREG-2026-09-23-v1.0.json",
 "depth_values_read":False,
 "intersection_bounds":list(inter.bounds),
 "intersection_area_deg2":float(inter.area),
 "survey":survey,
 "next_gate_authorized":all(v["minimum_overlap_samples_pass"] for v in survey.values()),
 "next_gate":"PREREGISTER_UNSEEN_RESIDUAL_MODEL_VALIDATION" if all(v["minimum_overlap_samples_pass"] for v in survey.values()) else "UNSEEN_PAIR_INSUFFICIENT_FNV_OVERLAP",
 "claim_ceiling":"UNSEEN_FNV_COVERAGE_SUBSET_FREEZE_ONLY"
}
p=OUT/"JANUS-KUSTO-UNSEEN-PF0501-X-EX2105-FNV-OVERLAP-FREEZE-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
