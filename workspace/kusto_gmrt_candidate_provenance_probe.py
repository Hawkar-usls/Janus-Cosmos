#!/usr/bin/env python3
import json, requests
from pathlib import Path
TARGETS=[
 {"id":"KN19207_CAND_002","lat":-3.9727527956056825,"lon":-12.272824298723462},
 {"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432}
]
BASE="https://www.gmrt.org/services/GridServer"
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
s=requests.Session();s.headers.update({"User-Agent":"JANUS-KUSTO-GMRT-provenance-probe/1.0"})
rows=[]
for t in TARGETS:
    pad=.03
    p={"west":t["lon"]-pad,"east":t["lon"]+pad,"south":t["lat"]-pad,"north":t["lat"]+pad,
       "resolution":"max","format":"netcdf","mformat":"json"}
    rec={**t,"query":p}
    try:
        r=s.get(BASE+"/metadata",params=p,timeout=120)
        rec["http_status"]=r.status_code
        rec["content_type"]=r.headers.get("content-type")
        rec["response_text"]=r.text[:200000]
        try:rec["json"]=r.json()
        except:pass
    except Exception as e:rec["error"]=repr(e)
    rows.append(rec)
out={"artifact_id":"JANUS-KUSTO-GMRT-CAND002003-PROVENANCE-PROBE-2026-09-22-v1.0",
     "targets":rows,
     "claim_ceiling":"OPEN_METADATA_SOURCE_DISCOVERY_ONLY"}
p=OUT/"JANUS-KUSTO-GMRT-CAND002003-PROVENANCE-PROBE-2026-09-22-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
