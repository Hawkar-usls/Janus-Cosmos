#!/usr/bin/env python3
import requests,re,json,hashlib
from pathlib import Path

BASE="https://geodataindia.gov.in/"
OUT=Path("workspace/kusto_global_groundtruth_out/ngdr_sm222_catalog")
OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-NGDR-SM222/3.0","Referer":BASE+"guestuser"})
h=S.get(BASE+"guestuser",timeout=120); h.raise_for_status()
html=h.text
csrf=re.search(r'id=["\']csrfvalue["\']\s+value=["\']([^"\']+)',html)
if not csrf: raise RuntimeError("CSRF not found")
tok=csrf.group(1)
url=BASE+"guestuser/getbasereport?_csrf="+tok
r=S.post(url,data={"geophysical_report":"Geological Report"},timeout=180)
rec={"status":r.status_code,"url":r.url,"content_type":r.headers.get("content-type"),"bytes":len(r.content)}
try: data=r.json()
except Exception:
    data=None
    rec["text_head"]=r.text[:4000]
if not isinstance(data,list):
    raise RuntimeError("Unexpected catalog response: "+json.dumps(rec,indent=2))

keys=["sm-222","sm222","pudimadaka","godavari","multibeam","continental slope","kakinada","marine"]
hits=[]
for row in data:
    text=json.dumps(row,ensure_ascii=False).lower()
    matched=[k for k in keys if k in text]
    if matched:
        hits.append({"matched":matched,"row":row})

out={"artifact_id":"JANUS-KUSTO-NGDR-SM222-CATALOG-QUERY-2026-09-25-v1.0",
     "endpoint":url,"catalog_count":len(data),"keys":keys,"hit_count":len(hits),"hits":hits,
     "transport":rec}
raw=json.dumps(out,indent=2,ensure_ascii=False)
(OUT/"catalog_hits.json").write_text(raw,encoding="utf-8")
print(json.dumps({
  "status":"PASS","catalog_count":len(data),"hit_count":len(hits),
  "hits":[{"matched":x["matched"],"row":x["row"]} for x in hits[:100]],
  "sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2,ensure_ascii=False))
