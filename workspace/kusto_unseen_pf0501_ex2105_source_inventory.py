#!/usr/bin/env python3
from __future__ import annotations
import json,re,time
from pathlib import Path
from urllib.parse import urljoin
import requests

OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-unseen-source-inventory/1.0"}
SURVEYS={
 "PF0501":{
   "report":"https://www.ngdc.noaa.gov/ships/pathfinder/PF0501_mb.html",
   "generated_candidates":[
     "https://data.ngdc.noaa.gov/platforms/ocean/ships/pathfinder/PF0501/multibeam/data/version1/MB/generated/",
     "https://www.ngdc.noaa.gov/ships/pathfinder/PF0501_mb.html"
   ]
 },
 "EX2105":{
   "report":"https://www.ngdc.noaa.gov/ships/okeanos_explorer/EX2105_mb.html",
   "generated_candidates":[
     "https://data.ngdc.noaa.gov/platforms/ocean/ships/okeanos_explorer/EX2105/multibeam/data/version1/MB/generated/",
     "https://www.ngdc.noaa.gov/ships/okeanos_explorer/EX2105_mb.html"
   ]
 }
}

def get(url,timeout=90):
    last=None
    for i in range(6):
        try:
            r=requests.get(url,headers=UA,timeout=timeout)
            if r.status_code==429:
                time.sleep(min(20,2**i));continue
            return r
        except Exception as e:
            last=e;time.sleep(min(20,2**i))
    raise last

def links(html,base):
    return sorted(set(urljoin(base,h) for h in re.findall(r'href=["\']([^"\']+)["\']',html,re.I)))

def parse_fnv_text(text):
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:
            portlon=float(p[15]);portlat=float(p[16]);stbdlon=float(p[17]);stbdlat=float(p[18])
        except:continue
        if all(-180<=x<=180 for x in (portlon,stbdlon)) and all(-90<=x<=90 for x in (portlat,stbdlat)):
            rows.append((portlon,portlat,stbdlon,stbdlat))
    return rows

def summarize_rows(rows):
    if not rows:return None
    lons=[x for r in rows for x in (r[0],r[2])]
    lats=[x for r in rows for x in (r[1],r[3])]
    return {"rows":len(rows),"west":min(lons),"east":max(lons),"south":min(lats),"north":max(lats)}

result={}
for sid,spec in SURVEYS.items():
    item={"report":spec["report"],"candidates":[],"selected_generated_base":None,"fnv_files":[],"fnv_total_rows":0,"fnv_envelope":None}
    all_rows=[]
    for base in spec["generated_candidates"]:
        r=get(base)
        rec={"url":base,"status":r.status_code,"content_type":r.headers.get("content-type"),"bytes":len(r.content)}
        if r.ok and "text" in (r.headers.get("content-type") or "").lower() or (r.ok and b"href" in r.content[:20000].lower()):
            html=r.text
            ls=links(html,base)
            rec["link_count"]=len(ls)
            rec["fnv_links"]=[u for u in ls if u.lower().endswith(".fnv")]
            rec["fbt_links"]=[u for u in ls if u.lower().endswith(".fbt")]
            rec["raw_like_links"]=[u for u in ls if any(u.lower().endswith(ext) for ext in (".mb51.gz",".kmall.gz",".all.gz"))][:500]
            rec["xyz_links"]=[u for u in ls if ".xyz" in u.lower()][:100]
            if rec["fnv_links"] and item["selected_generated_base"] is None:
                item["selected_generated_base"]=base
                for u in rec["fnv_links"]:
                    fr=get(u)
                    if not fr.ok:continue
                    rows=parse_fnv_text(fr.text)
                    item["fnv_files"].append({"url":u,"rows":len(rows),"bytes":len(fr.content)})
                    all_rows.extend(rows)
        item["candidates"].append(rec)
    item["fnv_total_rows"]=len(all_rows)
    item["fnv_envelope"]=summarize_rows(all_rows)
    result[sid]=item

# Bounding-envelope intersection is discovery-only, not the exact overlap gate.
a=result["PF0501"].get("fnv_envelope");b=result["EX2105"].get("fnv_envelope")
if a and b:
    inter={
      "west":max(a["west"],b["west"]),"east":min(a["east"],b["east"]),
      "south":max(a["south"],b["south"]),"north":min(a["north"],b["north"])
    }
    inter["positive_area_envelope"]=bool(inter["west"]<inter["east"] and inter["south"]<inter["north"])
else:inter=None

out={
 "artifact_id":"JANUS-KUSTO-UNSEEN-MULTIBEAM-PF0501-X-EX2105-SOURCE-INVENTORY-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-UNSEEN-MULTIBEAM-PF0501-X-EX2105-SOURCE-INVENTORY-PREREG-2026-09-23-v1.0.json",
 "depth_values_read":False,
 "surveys":result,
 "coordinate_envelope_intersection_discovery_only":inter,
 "next_gate_authorized":bool(
    result["PF0501"]["fnv_total_rows"]>=100 and result["EX2105"]["fnv_total_rows"]>=100 and inter and inter["positive_area_envelope"]
 ),
 "next_gate":"EXACT_FNV_SWATH_INTERSECTION_ONLY" if (
    result["PF0501"]["fnv_total_rows"]>=100 and result["EX2105"]["fnv_total_rows"]>=100 and inter and inter["positive_area_envelope"]
 ) else "SOURCE_ROUTE_REPAIR_OR_RAW_NAV_EXTRACTION",
 "claim_ceiling":"UNSEEN_VALIDATION_SOURCE_AND_COVERAGE_INVENTORY_ONLY"
}
p=OUT/"JANUS-KUSTO-UNSEEN-MULTIBEAM-PF0501-X-EX2105-SOURCE-INVENTORY-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "depth_values_read":False,
 "PF0501":{"selected_generated_base":result["PF0501"]["selected_generated_base"],"fnv_files":len(result["PF0501"]["fnv_files"]),"fnv_total_rows":result["PF0501"]["fnv_total_rows"],"fnv_envelope":result["PF0501"]["fnv_envelope"]},
 "EX2105":{"selected_generated_base":result["EX2105"]["selected_generated_base"],"fnv_files":len(result["EX2105"]["fnv_files"]),"fnv_total_rows":result["EX2105"]["fnv_total_rows"],"fnv_envelope":result["EX2105"]["fnv_envelope"]},
 "coordinate_envelope_intersection_discovery_only":inter,
 "next_gate_authorized":out["next_gate_authorized"],
 "next_gate":out["next_gate"]
},indent=2))
