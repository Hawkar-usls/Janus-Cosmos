#!/usr/bin/env python3
import json,re,requests
from urllib.parse import urljoin
from pathlib import Path
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-unseen-route-repair/1.0"}

def get(u):
    try:return requests.get(u,headers=UA,timeout=45)
    except Exception as e:return None

def hrefs(r,base):
    if r is None or not r.ok:return []
    return sorted(set(urljoin(base,h) for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I)))

pf="https://data.ngdc.noaa.gov/platforms/ocean/ships/pathfinder/PF0501/multibeam/data/version1/MB/generated/"
pr=get(pf);pls=hrefs(pr,pf)
pfnv=[u for u in pls if u.lower().endswith(".fnv")]
pf_probe={"base_status":None if pr is None else pr.status_code,"fnv_count":len(pfnv),"first_fnv":pfnv[0] if pfnv else None,"first_data_lines":[]}
if pfnv:
    rr=get(pfnv[0])
    if rr and rr.ok:
        pf_probe["first_data_lines"]=[x for x in rr.text.splitlines() if x.strip() and not x.lstrip().startswith("#")][:12]

bases=[]
for version in ["version1","version2"]:
    root=f"https://data.ngdc.noaa.gov/platforms/ocean/ships/okeanos_explorer/EX2105/multibeam/data/{version}/MB/"
    for suffix in ["","generated/","raw/","processed/","products/"]:
        bases.append(root+suffix)

ex=[]
for base in bases:
    r=get(base);ls=hrefs(r,base)
    ex.append({
      "url":base,
      "status":None if r is None else r.status_code,
      "bytes":None if r is None else len(r.content),
      "href_count":len(ls),
      "fnv":[u for u in ls if u.lower().endswith(".fnv")][:20],
      "kmall":[u for u in ls if ".kmall" in u.lower()][:20],
      "gsf":[u for u in ls if ".gsf" in u.lower()][:20],
      "xyz":[u for u in ls if ".xyz" in u.lower()][:20]
    })

out={
 "artifact_id":"JANUS-KUSTO-UNSEEN-PF0501-EX2105-SOURCE-ROUTE-REPAIR-PROBE-2026-09-23-v1.0",
 "depth_values_read":False,
 "PF0501_fnv_format_probe":pf_probe,
 "EX2105_route_probe":ex,
 "claim_ceiling":"SOURCE_ROUTE_AND_TEXT_FORMAT_PROBE_ONLY"
}
p=OUT/"JANUS-KUSTO-UNSEEN-PF0501-EX2105-SOURCE-ROUTE-REPAIR-PROBE-2026-09-23-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
