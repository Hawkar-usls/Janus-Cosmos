#!/usr/bin/env python3
"""R13B access-only repair. Keep original candidate tests; never infer coverage from HTTP status."""
import hashlib,json,pathlib,sys
from urllib.parse import urlparse,urlunparse
import requests
import r13_complete_archive as r13
import r9_expanded_search as eng
import r4_chunk_replay as base
ROOT=pathlib.Path("runtime/r13b");ROOT.mkdir(parents=True,exist_ok=True)
eng.OUT=ROOT
def emit(k,x):print(k,json.dumps(x),flush=True)
def save(p,x):p.write_text(json.dumps(x,indent=2)+"\n")
def probe(url,expected):
    result={"url":url}
    try:
        with requests.get(url,headers={"Range":"bytes=0-0","Accept-Encoding":"identity"},stream=True,timeout=(8,15)) as r:
            result.update(status=r.status_code,content_range=r.headers.get("Content-Range"),content_type=r.headers.get("Content-Type"),final_url=r.url)
            cr=r.headers.get("Content-Range","")
            result["eligible"]=r.status_code==206 and cr.startswith("bytes 0-0/") and cr.split("/")[-1].isdigit() and int(cr.split("/")[-1])==expected
    except Exception as e:result["error"]=str(e)
    return result
def alternatives(f):
    u=f["url"];p=urlparse(u)
    # Same official archive host and exact original basename, only test public mirror path or HTTPS.
    urls=[u,urlunparse(p._replace(scheme="https"))]
    if p.path.startswith("/dl2/"):
        for scheme in ("https","http"):
            urls.append(urlunparse(p._replace(scheme=scheme,path=p.path.replace("/dl2/","/dl/",1))))
    return list(dict.fromkeys(urls))
def candidate_valid(h,f):
    # R9 broad selection required unrelated windows at 2800 etc.
    # Only require R13 candidate and its fixed model windows, checked independently during inspect.
    return h["source"].upper()==r13.TARGET+"_"+f["position"] and len(h["shape"])==3 and 12<=h["shape"][0]<=64 and 0<abs(h["df"])*1e6<=10 and 0<h["ts"]<=60 and abs((h["mjd"]-f["day"])*86400-f["seconds"])<=3 and h["bounds"][0]<r13.F-.5 and h["bounds"][1]>r13.F+.5
r13.valid=candidate_valid
def main():
    fs,q=eng.query(r13.TARGET)
    # Same exact archive index as R13, no changing sample after outcomes.
    assert q["sha256"]=="5ee07a6143715b1c0c078bee40f5a5d1a1ff2e5c1215ca16240bd4a13e65a73f"
    selected=sorted((f for f in fs if f["day"]==57708 or "/dl2/" in urlparse(f["url"]).path),key=lambda f:(f["day"],f["seconds"],f["url"]))
    m={"source_index":q,"selected":[{"name":f["name"],"url":f["url"],"day":f["day"],"position":f["position"]} for f in selected],"selection":"R13 access failures: old header band and dl2 URLs; frozen archive index","limits":"No full-file download; validate range and exact indexed length."}
    save(ROOT/"manifest.json",m);emit("R13B_MANIFEST",m)
    rows=[]
    for f in selected:
        expected=int(f["metadata"]["size"]);probes=[];result=None
        for url in alternatives(f):
            p=probe(url,expected);probes.append(p)
            if not p.get("eligible"):continue
            v=r13.inspect({**f,"url":url},ROOT)
            if v["status"]=="COMPLETE":
                # Coordinates absent in old header are explicit limitation, not assumed pointing verification.
                v["original_url"]=f["url"];v["size_verified"]=expected
                v["pointing_header_missing"]=v["header"]["src_raj"]=="0.0" and v["header"]["src_dej"]=="0.0"
                result=v;break
            probes[-1]["read_error"]=v.get("error")
        if result is None:result={"status":"HOLD_NO_VERIFIED_RANGE_ACCESS","name":f["name"],"day":f["day"],"position":f["position"],"url":f["url"],"tests":{}}
        result["access_probes"]=probes;rows.append(result);emit("R13B_SCAN",result)
        save(ROOT/"result.json",{"manifest":m,"rows":rows})
    summary={"selected":len(selected),"completed":sum(r["status"]=="COMPLETE" for r in rows),"held":sum(r["status"]!="COMPLETE" for r in rows),"candidate_strong_scans":sum(any(v.get("strong") for k,v in r["tests"].items() if k!="reverse_control_fixed") for r in rows),"pointing_missing_scans":sum(r.get("pointing_header_missing",False) for r in rows),"confirmed_technosignatures":0,"independently_confirmed":0}
    save(ROOT/"result.json",{"manifest":m,"rows":rows,"summary":summary});emit("R13B_SUMMARY",summary)
if __name__=="__main__":main()
