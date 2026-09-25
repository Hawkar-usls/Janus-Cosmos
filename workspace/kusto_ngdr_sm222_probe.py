#!/usr/bin/env python3
import json,re,hashlib,urllib.parse
from pathlib import Path
import requests

BASE="https://geodataindia.gov.in"
START=BASE+"/guestuser"
OUT=Path("workspace/kusto_global_groundtruth_out/ngdr_sm222_probe")
OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session()
S.headers.update({"User-Agent":"JANUS-KUSTO-NGDR-SM222/1.0","Accept":"text/html,application/xhtml+xml,application/json,*/*"})

r=S.get(START,timeout=120)
r.raise_for_status()
html=r.text
(OUT/"guestuser.html").write_text(html,encoding="utf-8")

srcs=[]
for m in re.finditer(r'''<script[^>]+src=["']([^"']+)["']''',html,re.I):
    src=urllib.parse.urljoin(r.url,m.group(1))
    if src not in srcs: srcs.append(src)

hrefs=[]
for m in re.finditer(r'''href=["']([^"']+)["']''',html,re.I):
    href=urllib.parse.urljoin(r.url,m.group(1))
    if href not in hrefs: hrefs.append(href)

keywords=["sm-222","pudimadaka","keyword","report","search","guestuser","nuid","ajax","api/","mapservice","baseline","bathym"]
assets=[]
all_snips=[]
for i,url in enumerate(srcs):
    try:
        rr=S.get(url,timeout=120)
        ct=rr.headers.get("content-type","")
        txt=rr.text if ("javascript" in ct or url.lower().endswith(".js") or "text" in ct) else ""
        sha=hashlib.sha256(rr.content).hexdigest()
        rec={"url":url,"status":rr.status_code,"bytes":len(rr.content),"content_type":ct,"sha256":sha}
        if txt:
            low=txt.lower()
            hits=[]
            for kw in keywords:
                pos=0
                while True:
                    p=low.find(kw,pos)
                    if p<0: break
                    a=max(0,p-350); b=min(len(txt),p+650)
                    snippet=txt[a:b]
                    hits.append({"keyword":kw,"offset":p,"snippet":snippet})
                    all_snips.append({"asset":url,"keyword":kw,"offset":p,"snippet":snippet})
                    pos=p+len(kw)
                    if len(hits)>=80: break
            rec["hit_count"]=len(hits)
        assets.append(rec)
    except Exception as e:
        assets.append({"url":url,"error":repr(e)})

# Candidate endpoint/string mining from hit snippets only.
cand=set()
url_pat=re.compile(r'''(?:"|')((?:https?://[^"'\\\s]+)|(?:/[A-Za-z0-9_./?=&:%{}-]{4,}))(?:"|')''')
for x in all_snips:
    for mm in url_pat.finditer(x["snippet"]):
        s=mm.group(1)
        if any(k in s.lower() for k in ["search","report","api","map","guest","nuid","keyword","layer"]):
            cand.add(s)

# Inline HTML endpoint clues around search forms / JS.
inline=[]
low=html.lower()
for kw in ["by keyword","search result","nuid","search","report"]:
    pos=0
    while True:
        p=low.find(kw,pos)
        if p<0: break
        inline.append({"keyword":kw,"offset":p,"snippet":html[max(0,p-500):min(len(html),p+900)]})
        pos=p+len(kw)
        if len(inline)>100: break

out={
 "artifact_id":"JANUS-KUSTO-NGDR-SM222-ENDPOINT-DISCOVERY-2026-09-25-v1.0",
 "start_url":r.url,
 "html_bytes":len(r.content),
 "script_count":len(srcs),
 "href_count":len(hrefs),
 "scripts":srcs,
 "assets":assets,
 "candidate_endpoint_strings":sorted(cand),
 "matching_snippets":all_snips[:1000],
 "inline_html_snippets":inline[:100]
}
raw=json.dumps(out,indent=2,ensure_ascii=False)
(OUT/"JANUS-KUSTO-NGDR-SM222-ENDPOINT-DISCOVERY-2026-09-25-v1.0.json").write_text(raw,encoding="utf-8")

print(json.dumps({
 "status":"PASS",
 "start_url":r.url,
 "script_count":len(srcs),
 "scripts_with_hits":[{"url":a.get("url"),"hit_count":a.get("hit_count",0)} for a in assets if a.get("hit_count",0)>0],
 "candidate_endpoint_strings":sorted(cand)[:200],
 "total_matching_snippets":len(all_snips),
 "output_sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2,ensure_ascii=False))
