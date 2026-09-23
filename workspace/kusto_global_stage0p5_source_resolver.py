#!/usr/bin/env python3
import json,re,time,requests
from pathlib import Path
from urllib.parse import urljoin,urlparse

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"data/cousteau/JANUS-KUSTO-GLOBAL-NCEI-MULTISURVEY-OPPORTUNITY-MAP-RECEIPT-2026-09-23-v1.0.json"
OUT=ROOT/"workspace/kusto_global_anomaly_out"
OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-global-stage0p5-source-resolver/1.0"}

def get(u,stream=False):
    last=None
    for i in range(7):
        try:
            r=requests.get(u,headers=UA,timeout=90,stream=stream)
            if r.status_code==429:
                time.sleep(min(30,2**i)); continue
            return r
        except Exception as e:
            last=e
            if i==6: raise
            time.sleep(min(30,2**i))
    raise last

def hrefs(u):
    r=get(u)
    if not r.ok:return []
    return sorted(set(urljoin(u,h) for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I)))

def derive_base(download_url):
    # https://www.ngdc.noaa.gov/ships/<ship>/<survey>_mb.html
    p=urlparse(download_url)
    parts=[x for x in p.path.split("/") if x]
    if len(parts)>=3 and parts[-1].endswith("_mb.html"):
        ship=parts[-2]
        survey=parts[-1][:-8]
        return f"https://data.ngdc.noaa.gov/platforms/ocean/ships/{ship}/{survey}/multibeam/data/"
    return None

def inventory(base,maxdepth=4,max_dirs=250):
    seen=set(); files=[]; dirs=[]; errors=[]
    def walk(u,depth):
        if u in seen or depth>maxdepth or len(dirs)>=max_dirs:return
        seen.add(u)
        r=get(u)
        if not r.ok:
            errors.append({"url":u,"status":r.status_code});return
        dirs.append(u)
        for x in hrefs(u):
            if not x.startswith(base):continue
            if x.rstrip("/")==u.rstrip("/"):continue
            name=urlparse(x).path.rstrip("/").split("/")[-1]
            if not name:continue
            if x.endswith("/"):
                walk(x,depth+1)
            else:
                ext=Path(name).suffix.lower()
                files.append({"url":x,"name":name,"ext":ext,"depth":depth})
    walk(base,0)
    return files,dirs,errors

def classify(files):
    names=[f["name"].lower() for f in files]
    ext_counts={}
    for f in files:ext_counts[f["ext"]]=ext_counts.get(f["ext"],0)+1
    fbt=[f for f in files if f["name"].lower().endswith(".fbt")]
    fnv=[f for f in files if f["name"].lower().endswith(".fnv")]
    kmall=[f for f in files if ".kmall" in f["name"].lower()]
    allf=[f for f in files if f["name"].lower().endswith(".all") or f["name"].lower().endswith(".all.gz")]
    mb=[f for f in files if re.search(r'\.mb\d+(?:\.gz)?$',f["name"].lower())]
    grids=[f for f in files if f["ext"] in (".grd",".nc",".tif",".tiff",".xyz")]
    if fbt and fnv:
        state="READY_FNV_FBT"
    elif fbt or mb:
        state="READY_GENERATED_BATHY_OTHER"
    elif kmall:
        state="RAW_KMALL"
    elif allf:
        state="RAW_ALL_OR_EMOLDRAW"
    elif grids:
        state="GRID_ONLY"
    elif files:
        state="DOWNLOAD_ROUTE_PRESENT_UNRESOLVED"
    else:
        state="UNAVAILABLE"
    return {
      "state":state,
      "ext_counts":dict(sorted(ext_counts.items())),
      "fnv_count":len(fnv),"fbt_count":len(fbt),"kmall_count":len(kmall),
      "all_count":len(allf),"mb_count":len(mb),"grid_count":len(grids),
      "fnv_sample":[x["url"] for x in fnv[:5]],
      "fbt_sample":[x["url"] for x in fbt[:5]],
      "kmall_sample":[x["url"] for x in kmall[:5]],
      "all_sample":[x["url"] for x in allf[:5]],
      "mb_sample":[x["url"] for x in mb[:5]],
      "grid_sample":[x["url"] for x in grids[:5]]
    }

receipt=json.loads(RECEIPT.read_text())
regions=[]
for reg in receipt["stage1_regions"]:
    rr={"rank":reg["rank"],"overlap_point":reg["overlap_point"],"surveys":[]}
    for key in ["survey_a","survey_b"]:
        s=reg[key]
        base=derive_base(s["download_url"])
        if base is None:
            rec={**s,"base":None,"resolution":{"state":"UNAVAILABLE"}}
        else:
            files,dirs,errors=inventory(base)
            rec={**s,"base":base,"directory_count":len(dirs),"file_count":len(files),
                 "errors":errors[:20],"resolution":classify(files)}
        rr["surveys"].append(rec)
    states=[x["resolution"]["state"] for x in rr["surveys"]]
    ready_states={"READY_FNV_FBT","READY_GENERATED_BATHY_OTHER","RAW_KMALL","RAW_ALL_OR_EMOLDRAW","GRID_ONLY"}
    rr["stage1_ready"]=all(x in ready_states for x in states)
    rr["states"]=states
    regions.append(rr)

out={
 "artifact_id":"JANUS-KUSTO-GLOBAL-STAGE0P5-SOURCE-RESOLVER-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-GLOBAL-STAGE0P5-SOURCE-RESOLVER-PREREG-2026-09-23-v1.0.json",
 "stage0_receipt":str(RECEIPT.relative_to(ROOT)),
 "depth_values_read":False,
 "regions":regions,
 "stage1_ready_region_count":sum(r["stage1_ready"] for r in regions),
 "claim_ceiling":"GLOBAL_STAGE0P5_SOURCE_ROUTING_ONLY"
}
p=OUT/"JANUS-KUSTO-GLOBAL-STAGE0P5-SOURCE-RESOLVER-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
