#!/usr/bin/env python3
from __future__ import annotations
import json,re
from pathlib import Path
from urllib.parse import urljoin
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-LAKE-SAMMAMISH-UNSEEN-DUALCHANNEL-PREREG-2026-09-23-v1.0.json").read_text())
PAGE=PRE["dataset"]["source_page"]
UA={"User-Agent":"Mozilla/5.0 (compatible; JANUS-KUSTO-LakeSammamish-source-resolver/1.0; research source binding)"}

expected={v["name"]:{**v,"channel":k.upper()} for k,v in PRE["dataset"]["files"].items()}
r=requests.get(PAGE,headers=UA,timeout=90,allow_redirects=True);r.raise_for_status()
html=r.text
hrefs=[urljoin(r.url,h) for h in re.findall(r'href=["\']([^"\']+)["\']',html,re.I)]
resolved={}
for name,meta in expected.items():
    urls=sorted(set(u for u in hrefs if name in u or name in requests.utils.unquote(u)))
    rows=[]
    for u in urls:
        h=requests.head(u,headers=UA,timeout=60,allow_redirects=True)
        rows.append({"url":u,"status":h.status_code,"resolved_url":h.url,"content_length":h.headers.get("content-length"),"content_type":h.headers.get("content-type"),"etag":h.headers.get("etag")})
    resolved[name]={"expected":meta,"links":urls,"head":rows}

ok=sum(1 for v in resolved.values() if v["links"])
out={
 "artifact_id":"JANUS-KUSTO-LAKE-SAMMAMISH-LS0-SOURCE-BINDING-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "raster_values_read":False,
 "source_page":{"url":PAGE,"resolved_url":r.url,"status":r.status_code,"bytes":len(r.content)},
 "files":resolved,
 "resolved_expected_file_count":ok,
 "expected_file_count":len(expected),
 "status":"PASS_ALL_FILE_LINKS_RESOLVED" if ok==len(expected) else "SOURCE_BINDING_INCOMPLETE",
 "next_gate":"LS1_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if ok==len(expected) else "ALTERNATE_TRACEABLE_DOWNLOAD_ROUTE_DISCOVERY",
 "claim_ceiling":"SOURCE_RESOLUTION_ONLY__NO_RASTER_VALUES_READ"
}
p=OUT/"JANUS-KUSTO-LAKE-SAMMAMISH-LS0-SOURCE-BINDING-RUN-2026-09-23-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps({"artifact_id":out["artifact_id"],"resolved":ok,"expected":len(expected),"status":out["status"],"files":{k:v["links"] for k,v in resolved.items()},"raster_values_read":False},indent=2))
