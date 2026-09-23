#!/usr/bin/env python3
import json,re,requests,time
from urllib.parse import urljoin,urlparse
from pathlib import Path

OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-unseen-ex2105-em304-inventory-v3/1.0"}
ROOTS=[
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/okeanos_explorer/EX2105/multibeam/data/version1/MB/em304/",
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/okeanos_explorer/EX2105/multibeam/data/version2/MB/em304/"
]

def get(u,stream=False):
    last=None
    for i in range(6):
        try:
            r=requests.get(u,headers=UA,timeout=90,stream=stream)
            if r.status_code==429:
                time.sleep(min(20,2**i)); continue
            return r
        except Exception as e:
            last=e; time.sleep(min(20,2**i))
    raise last

def hrefs(r,base):
    if not r.ok:return []
    return sorted(set(urljoin(base,h) for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I)))

def basename(u): return urlparse(u).path.rstrip("/").split("/")[-1]
def is_data_candidate(u):
    b=basename(u).lower()
    return any(b.endswith(ext) for ext in [".kmall",".all",".raw",".kmwcd",".wcd",".txt",".nav",".csv",".json",".xml",".nc",".xyz",".gz",".bz2",".zip"])

def recurse(root,maxdepth=3):
    seen=set(); files=[]; dirs=[]; errors=[]
    def walk(u,depth):
        if u in seen or depth>maxdepth:return
        seen.add(u)
        r=get(u)
        if not r.ok:
            errors.append({"url":u,"status":r.status_code}); return
        for x in hrefs(r,u):
            if x.rstrip("/")==u.rstrip("/"):continue
            if not x.startswith(root):continue
            b=basename(x)
            if not b or b in ("em304","MB","version1","version2"):continue
            if x.endswith("/"):
                dirs.append({"url":x,"depth":depth})
                walk(x,depth+1)
            else:
                rec={"url":x,"basename":b,"extension":Path(b).suffix.lower(),"depth":depth}
                if is_data_candidate(x):
                    try:
                        h=get(x,stream=True)
                        rec["status"]=h.status_code
                        cl=h.headers.get("Content-Length")
                        rec["bytes"]=int(cl) if cl and cl.isdigit() else None
                        h.close()
                    except Exception as e:
                        rec["error"]=repr(e)
                files.append(rec)
    walk(root,0)
    return {"root":root,"directories":dirs,"files":files,"errors":errors}

inventories=[recurse(r,4) for r in ROOTS]
for inv in inventories:
    ext={}
    for f in inv["files"]:
        ext[f["extension"]]=ext.get(f["extension"],0)+1
    inv["summary"]={
      "directory_count":len(inv["directories"]),
      "file_count":len(inv["files"]),
      "extension_counts":dict(sorted(ext.items())),
      "data_candidate_count":sum(is_data_candidate(f["url"]) for f in inv["files"])
    }

out={
 "artifact_id":"JANUS-KUSTO-UNSEEN-EX2105-EM304-RAW-INVENTORY-V3-2026-09-23-v1.0",
 "parent_prereg":"data/cousteau/JANUS-KUSTO-UNSEEN-MULTIBEAM-PF0501-X-EX2105-SOURCE-INVENTORY-PREREG-2026-09-23-v1.0.json",
 "depth_values_read":False,
 "inventories":inventories,
 "claim_ceiling":"UNSEEN_EX2105_RAW_SOURCE_ROUTE_INVENTORY_ONLY"
}
p=OUT/"JANUS-KUSTO-UNSEEN-EX2105-EM304-RAW-INVENTORY-V3-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "depth_values_read":False,
 "summaries":[{"root":x["root"],**x["summary"]} for x in inventories],
 "sample_files":[f for x in inventories for f in x["files"][:20]]
},indent=2))
