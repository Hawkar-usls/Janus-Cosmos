#!/usr/bin/env python3
import concurrent.futures,json,re,requests,time
from urllib.parse import urljoin
from pathlib import Path
from shapely.geometry import shape,LineString
from shapely.ops import unary_union

OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-unseen-pf0501-rb1202-fnv-overlap/1.0"}
ARC="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
ROOTS={
 "PF0501":"https://data.ngdc.noaa.gov/platforms/ocean/ships/pathfinder/PF0501/multibeam/data/version1/MB/generated/",
 "RB1202":"https://data.ngdc.noaa.gov/platforms/ocean/ships/ronald_h._brown/RB1202/multibeam/data/version1/MB/generated/"
}

def get(u,params=None):
    last=None
    for i in range(7):
        try:
            r=requests.get(u,params=params,headers=UA,timeout=90)
            if r.status_code==429:
                time.sleep(min(30,2**i));continue
            r.raise_for_status();return r
        except Exception as e:
            last=e
            if i==6:raise
            time.sleep(min(30,2**i))
    raise last

def footprint(sid):
    d=get(ARC,{
      "where":f"SURVEY_ID='{sid}'","outFields":"SURVEY_ID",
      "returnGeometry":"true","outSR":"4326","f":"geojson"
    }).json()
    gs=[shape(f["geometry"]) for f in d.get("features",[]) if f.get("geometry")]
    if not gs:raise RuntimeError("missing footprint "+sid)
    return unary_union(gs)

inter=footprint("PF0501").intersection(footprint("RB1202"))
if inter.is_empty or inter.area<=0:raise RuntimeError("empty exact footprint intersection")

def list_fnv(root):
    r=get(root)
    return sorted(set(urljoin(root,h) for h in re.findall(r'href=["\']([^"\']+\.fnv)["\']',r.text,re.I)))

def one_fnv(u):
    r=get(u);rows=[]
    for line in r.text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:
            epoch=float(p[6]);plon=float(p[15]);plat=float(p[16]);slon=float(p[17]);slat=float(p[18])
        except:continue
        seg=LineString([(plon,plat),(slon,slat)])
        if seg.intersects(inter):
            rows.append({"epoch":epoch,"portlon":plon,"portlat":plat,"stbdlon":slon,"stbdlat":slat})
    return u,rows

survey={}
for sid,root in ROOTS.items():
    files=list_fnv(root);hits=[];fh=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        futs=[ex.submit(one_fnv,u) for u in files]
        for k,f in enumerate(concurrent.futures.as_completed(futs),1):
            u,rr=f.result()
            if rr:
                hits.extend(rr)
                fh.append({
                  "fnv_url":u,"fnv_file":u.rsplit("/",1)[-1],"hit_rows":len(rr),
                  "epoch_min":min(x["epoch"] for x in rr),"epoch_max":max(x["epoch"] for x in rr)
                })
            if k%50==0:print(sid,k,"/",len(files),flush=True)
    fh.sort(key=lambda x:x["fnv_file"])
    survey[sid]={
      "fnv_files_discovered":len(files),"files_with_overlap_rows":len(fh),
      "overlap_row_count":len(hits),"file_hits":fh,
      "epoch_min":min((x["epoch"] for x in hits),default=None),
      "epoch_max":max((x["epoch"] for x in hits),default=None),
      "minimum_overlap_samples_pass":len(hits)>=100
    }

out={
 "artifact_id":"JANUS-KUSTO-UNSEEN-PF0501-X-RB1202-FNV-OVERLAP-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-UNSEEN-PF0501-X-RB1202-FNV-OVERLAP-PREREG-2026-09-23-v1.0.json",
 "depth_values_read":False,
 "intersection_bounds":list(inter.bounds),"intersection_area_deg2":float(inter.area),
 "survey":survey,
 "next_gate_authorized":all(x["minimum_overlap_samples_pass"] for x in survey.values()),
 "next_gate":"PREREGISTER_FIRST_UNSEEN_ACROSS_TRACK_RESIDUAL_MODEL_VALIDATION_BEFORE_DEPTH" if all(x["minimum_overlap_samples_pass"] for x in survey.values()) else "UNSEEN_PAIR_INSUFFICIENT_FNV_OVERLAP",
 "claim_ceiling":"UNSEEN_FNV_COVERAGE_SUBSET_FREEZE_ONLY"
}
p=OUT/"JANUS-KUSTO-UNSEEN-PF0501-X-RB1202-FNV-OVERLAP-RUN-2026-09-23-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
