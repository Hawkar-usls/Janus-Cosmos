#!/usr/bin/env python3
from __future__ import annotations
import json,re
from pathlib import Path
from urllib.parse import urljoin
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-LAKE-CHELAN-UNSEEN-DUALCHANNEL-PREREG-2026-09-23-v1.0.json").read_text())
PAGE=PRE["dataset"]["source_page"]
UA={"User-Agent":"Mozilla/5.0 (compatible; JANUS-KUSTO-LakeChelan-source-resolver/1.0; research source binding)"}

expected={}
for tile,v in PRE["tiles"].items():
    expected[v["bathy"]]={"tile":tile,"channel":"BATHYMETRY","md5":v["bathy_md5"]}
    expected[v["backscatter"]]={"tile":tile,"channel":"BACKSCATTER","md5":v["backscatter_md5"]}

r=requests.get(PAGE,headers=UA,timeout=90,allow_redirects=True)
page_meta={"status":r.status_code,"resolved_url":r.url,"content_type":r.headers.get("content-type"),"bytes":len(r.content)}
html=r.text if r.status_code==200 else ""
hrefs=[urljoin(r.url,h) for h in re.findall(r'href=["\']([^"\']+)["\']',html,re.I)]
links={}
for name,meta in expected.items():
    matches=[u for u in hrefs if name in u or name in requests.utils.unquote(u)]
    links[name]=sorted(set(matches))

head_results={}
for name,urls in links.items():
    rows=[]
    for u in urls:
        try:
            h=requests.head(u,headers=UA,timeout=60,allow_redirects=True)
            rows.append({"url":u,"status":h.status_code,"resolved_url":h.url,"content_length":h.headers.get("content-length"),"content_type":h.headers.get("content-type"),"etag":h.headers.get("etag")})
        except Exception as e:
            rows.append({"url":u,"error":type(e).__name__+": "+str(e)})
    head_results[name]=rows

resolved_count=sum(1 for name,urls in links.items() if urls)
out={
 "artifact_id":"JANUS-KUSTO-LAKE-CHELAN-LC0-SOURCE-RESOLVER-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "raster_values_read":False,
 "page":page_meta,
 "expected_files":expected,
 "resolved_links":links,
 "head_results":head_results,
 "resolved_expected_file_count":resolved_count,
 "expected_file_count":len(expected),
 "status":"PASS_ALL_FILE_LINKS_RESOLVED" if resolved_count==len(expected) else ("PARTIAL_FILE_LINK_RESOLUTION" if resolved_count else "SOURCE_PAGE_LINK_RESOLUTION_BLOCKED"),
 "next_gate":"LC1_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if resolved_count==len(expected) else "ALTERNATE_TRACEABLE_DOWNLOAD_ROUTE_DISCOVERY",
 "claim_ceiling":"SOURCE_RESOLUTION_ONLY__NO_RASTER_VALUES_READ"
}
p=OUT/"JANUS-KUSTO-LAKE-CHELAN-LC0-SOURCE-RESOLVER-RUN-2026-09-23-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "page":page_meta,
 "resolved_expected_file_count":resolved_count,
 "expected_file_count":len(expected),
 "status":out["status"],
 "resolved_names":[n for n,u in links.items() if u],
 "raster_values_read":False
},indent=2))
