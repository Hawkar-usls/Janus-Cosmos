#!/usr/bin/env python3
import requests,re,json,hashlib,urllib.parse
from pathlib import Path

BASE="https://geodataindia.gov.in"
OUT=Path("workspace/kusto_global_groundtruth_out/ngdr_sm222_endpoints")
OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-NGDR-SM222/2.0"})
assets=[
"/customol/js/keywordsearch/keywordsearch.js",
"/mis/js/customscripts/Baselineexplorationdatasearch.js",
"/mis/js/customscripts/explorationdatasearch.js",
"/customol/js/get_explo_data/get_explo_data.js",
"/mis/js/customscripts/multiCriteriaSearch.js",
"/mis/js/customscripts/Baselinemulticriteriasearch.js",
"/mis/js/services/misservices.js",
"/customol/js/services/services.js",
"/customol/js/map_guestuser.js"
]
results=[]
for path in assets:
    url=BASE+path
    r=S.get(url,timeout=120)
    rec={"url":url,"status":r.status_code,"bytes":len(r.content),"sha256":hashlib.sha256(r.content).hexdigest()}
    txt=r.text
    name=Path(path).name
    (OUT/name).write_text(txt,encoding="utf-8")
    # pull ajax blocks and any endpoint-looking string literals
    blocks=[]
    for m in re.finditer(r'\$\.ajax\s*\(\s*\{',txt):
        p=m.start(); q=txt.find('});',p)
        if q<0: q=min(len(txt),p+2500)
        blocks.append(txt[p:min(len(txt),q+3)])
    eps=[]
    for m in re.finditer(r'''url\s*:\s*([^,\n]+)''',txt,re.I):
        eps.append(m.group(1).strip())
    funcs=[]
    for m in re.finditer(r'function\s+([A-Za-z0-9_]+)\s*\(([^)]*)\)',txt):
        fn=m.group(1)
        if any(k in fn.lower() for k in ['search','explor','keyword','report','data','baseline','query','geom']):
            funcs.append({"name":fn,"args":m.group(2)})
    relevant=[b for b in blocks if any(k in b.lower() for k in ['search','keyword','baseline','explor','report','nuid','geom','data'])]
    rec.update({"url_expressions":eps[:500],"functions":funcs[:300],"ajax_blocks":relevant[:300]})
    results.append(rec)

guest=S.get(BASE+"/guestuser",timeout=120).text
hidden={}
for m in re.finditer(r'<input[^>]+(?:id|name)=["\']([^"\']+)["\'][^>]*>',guest,re.I):
    tag=m.group(0); key=m.group(1)
    vm=re.search(r'value=["\']([^"\']*)["\']',tag,re.I)
    if vm and any(k in key.lower() for k in ['context','csrf','token','user']):
        hidden[key]=vm.group(1)
out={"artifact_id":"JANUS-KUSTO-NGDR-SM222-ENDPOINTS-2026-09-25-v1.0","base":BASE,"hidden":hidden,"assets":results}
raw=json.dumps(out,indent=2,ensure_ascii=False)
(OUT/"endpoints.json").write_text(raw,encoding="utf-8")
print(json.dumps({
 "status":"PASS",
 "hidden":hidden,
 "assets":[{"url":x["url"],"status":x["status"],"functions":x["functions"][:20],
            "url_expressions":x["url_expressions"][:50]} for x in results],
 "sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2,ensure_ascii=False))
