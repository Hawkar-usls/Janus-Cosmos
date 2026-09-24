#!/usr/bin/env python3
from __future__ import annotations
import json, os
from pathlib import Path
from urllib.parse import urlparse
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-MONEY-SHOAL-UNSEEN-V2-PREREG-2026-09-23-v1.0.json").read_text())
COR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ARAFURA-A0-PREREG-METHOD-INVARIANT-CORRECTION-2026-09-24-v1.0.json").read_text())

API="https://www.data.gov.au/data/api/3/action/package_search"
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-Arafura-A0C-current-route/1.0"})
queries=['"Arafura Marine Park survey" SOL7491 GA0366 bathymetry','"Arafura Marine Park survey" SOL7491 GA0366 backscatter']
packages={}
for q in queries:
    r=S.get(API,params={"q":q,"rows":100},timeout=120);r.raise_for_status()
    for p in (r.json().get("result") or {}).get("results",[]) or []:
        packages[p.get("id")]=p

def candidates(role):
    rows=[]
    for p in packages.values():
        title=(p.get("title") or "").lower()
        if "arafura marine park" not in title: continue
        if role=="bathymetry" and "backscatter" in title: continue
        if role=="backscatter" and "backscatter" not in title: continue
        for res in p.get("resources",[]) or []:
            rn=(res.get("name") or "").lower()
            fm=(res.get("format") or "").lower()
            if "metadata" in rn: continue
            if not res.get("url"): continue
            # numeric raster/archive resources only
            if role=="bathymetry" and not any(k in (rn+" "+fm) for k in ["tif","geotiff","zip"]): continue
            if role=="backscatter" and not any(k in (rn+" "+fm) for k in ["geotif","geotiff","zip","tif"]): continue
            rows.append({
              "package_id":p.get("id"),"package_name":p.get("name"),"package_title":p.get("title"),
              "metadata_modified":p.get("metadata_modified"),"resource_id":res.get("id"),
              "resource_name":res.get("name"),"format":res.get("format"),"url":res.get("url")
            })
    return rows

cand={"bathymetry":candidates("bathymetry"),"backscatter":candidates("backscatter")}
probes=[]
for role,rows in cand.items():
    for x in rows[:40]:
        pkg=x["package_name"];rid=x["resource_id"];decl=x["url"]
        routes=[("declared_url",decl)]
        if pkg and rid:
            page=f"https://www.data.gov.au/data/dataset/{pkg}/resource/{rid}"
            routes.append(("ckan_download",page+"/download"))
            bn=os.path.basename(urlparse(decl).path)
            if bn:routes.append(("ckan_download_with_basename",page+"/download/"+bn))
        for kind,u in routes:
            try:
                rr=S.get(u,headers={"Range":"bytes=0-0","Referer":f"https://www.data.gov.au/data/dataset/{pkg}"},timeout=90,allow_redirects=True,stream=True)
                probes.append({
                  "role":role,"route_kind":kind,"requested":u,"status":rr.status_code,"resolved":rr.url,
                  "content_type":rr.headers.get("content-type"),"content_length":rr.headers.get("content-length"),
                  "content_range":rr.headers.get("content-range"),"resource":x,"body_consumed":False
                })
            except Exception as e:
                probes.append({"role":role,"route_kind":kind,"requested":u,"error":type(e).__name__+": "+str(e),"resource":x,"body_consumed":False})

usable=[]
for p in probes:
    st=int(p.get("status",0) or 0);ct=(p.get("content_type") or "").lower()
    if st in (200,206) and "text/html" not in ct:
        usable.append(p)
# collapse catalog revisions/aliases by role + resolved physical URL
ded={}
for p in usable:ded[(p["role"],p.get("resolved") or p["requested"])]=p
usable=list(ded.values())
by_role={r:[x for x in usable if x["role"]==r] for r in ("bathymetry","backscatter")}

out={
 "artifact_id":"JANUS-KUSTO-ARAFURA-A0C-CURRENT-ROUTE-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],"method_invariant_correction":COR["artifact_id"],
 "selected_subarea":PRE["dataset"]["selected_subarea"],"selection_rule":PRE["dataset"]["selection_rule"],
 "queries":queries,"package_count":len(packages),"semantic_candidates":cand,
 "transport_probes":probes,"usable_routes_by_role":by_role,
 "status":"PASS_CURRENT_ROUTES_RESOLVED" if all(len(by_role[r])>=1 for r in by_role) else "CURRENT_ROUTE_BLOCKED",
 "raster_values_read":False,"response_bodies_consumed":False,"groundtruth_assets_read":False,
 "next_gate":"A0D_ARCHIVE_MEMBER_AND_RASTER_HEADER_AUDIT" if all(len(by_role[r])>=1 for r in by_role) else "SOURCE_ROUTE_REPAIR_ONLY",
 "claim_ceiling":"CURRENT_SOURCE_ROUTE_RESOLUTION_ONLY"
}
p=OUT/"JANUS-KUSTO-ARAFURA-A0C-CURRENT-ROUTE-RUN-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"status":out["status"],
 "usable_routes_by_role":{k:[{"route_kind":x["route_kind"],"resolved":x.get("resolved"),"status":x.get("status"),"content_range":x.get("content_range"),"resource":x["resource"]} for x in v] for k,v in by_role.items()},
 "raster_values_read":False,"groundtruth_assets_read":False
},indent=2,ensure_ascii=False))
