#!/usr/bin/env python3
import json,re,hashlib
from pathlib import Path
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup

OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-JC169-independence-audit/1.0"}
URLS=[
 "https://www.marine-geo.org/tools/entry/JC169",
 "https://www.marine-geo.org/tools/search/entry.php?id=JC169",
 "https://gmrt.org/services/GmrtCruises.php"
]
sess=requests.Session();sess.headers.update(UA)
pages=[]
links=[]
for url in URLS:
    rec={"requested_url":url}
    try:
        r=sess.get(url,timeout=120,allow_redirects=True)
        rec.update({"status":r.status_code,"final_url":r.url,"content_type":r.headers.get("content-type"),"bytes":len(r.content),"sha256":hashlib.sha256(r.content).hexdigest()})
        txt=r.text
        if "GmrtCruises.php" in url:
            try:
                data=r.json()
                j=[x for x in data if str(x.get("entry_id","")).upper()=="JC169" or str(x.get("survey_id","")).upper()=="JC169"]
                rec["jc169_records"]=j
            except Exception as e:rec["json_error"]=repr(e)
        else:
            rec["text_snippet"]=txt[:50000]
            soup=BeautifulSoup(txt,"html.parser")
            for a in soup.find_all("a",href=True):
                href=urljoin(r.url,a["href"])
                text=" ".join(a.get_text(" ",strip=True).split())
                if any(k in (href+" "+text).lower() for k in ["jc169","download","multibeam","swath","file","mbsystem","gsf","all","mb"]):
                    links.append({"source_page":r.url,"text":text,"href":href})
    except Exception as e:
        rec["error"]=repr(e)
    pages.append(rec)

# Probe same-origin candidate links metadata-only; don't download large/binary data.
probes=[]
seen=set()
for x in links:
    u=x["href"]
    if u in seen:continue
    seen.add(u)
    if len(probes)>=100:break
    try:
        rr=sess.head(u,timeout=30,allow_redirects=True)
        probes.append({"url":u,"status":rr.status_code,"final_url":rr.url,"content_type":rr.headers.get("content-type"),"content_length":rr.headers.get("content-length"),"last_modified":rr.headers.get("last-modified"),"etag":rr.headers.get("etag")})
    except Exception as e:probes.append({"url":u,"error":repr(e)})

out={
 "artifact_id":"JANUS-KUSTO-GMRT-JC169-INDEPENDENCE-AUDIT-PROBE-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-GMRT-JC169-INDEPENDENCE-AUDIT-PREREG-2026-09-22-v1.0.json",
 "pages":pages,
 "candidate_links":links,
 "link_head_probes":probes,
 "depth_values_read":False,
 "claim_ceiling":"ACQUISITION_IDENTITY_AND_INDEPENDENCE_ONLY"
}
p=OUT/"JANUS-KUSTO-GMRT-JC169-INDEPENDENCE-AUDIT-PROBE-2026-09-22-v1.0.json";p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "pages":[{k:v for k,v in x.items() if k not in ["text_snippet"]} for x in pages],
 "candidate_links":links[:100],
 "link_head_probes":probes[:100]
},indent=2,ensure_ascii=False))
