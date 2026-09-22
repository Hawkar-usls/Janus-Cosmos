#!/usr/bin/env python3
import ftplib, hashlib, json, os
from pathlib import Path

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/EM12"
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
TEXT_EXT={".txt",".log",".lst",".list",".csv",".asc",".ascii",".info",".inf",".cfg",".conf",".dat",".prn",".nav",".md"}

def conn():
    f=ftplib.FTP(timeout=120); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I"); return f

def list_dir(f,path):
    rows=[]
    try:
        for name,facts in f.mlsd(path):
            if name in (".",".."): continue
            rows.append((name,facts))
        return rows
    except Exception:
        out=[]
        for p in f.nlst(path):
            out.append((p.rstrip("/").split("/")[-1],{}))
        return out

def walk(f,path,depth=0,maxdepth=3):
    items=[]
    for name,facts in list_dir(f,path):
        full=path.rstrip("/")+"/"+name
        typ=facts.get("type")
        if typ=="dir" and depth<maxdepth:
            items.append({"path":full,"type":"dir","children":walk(f,full,depth+1,maxdepth)})
        else:
            size=None
            try:size=f.size(full)
            except:pass
            items.append({"path":full,"type":typ or "file","size":size})
    return items

def flatten(items):
    for x in items:
        yield x
        for y in flatten(x.get("children",[])): yield y

f=conn()
try:
    tree=walk(f,ROOT,0,4)
    files=[x for x in flatten(tree) if x.get("type")!="dir"]
    extracts=[]
    for x in files:
        p=x["path"]; ext=Path(p).suffix.lower(); sz=x.get("size")
        base=Path(p).name.lower()
        wanted=(ext in TEXT_EXT or any(k in base for k in ["readme","process","header","meta","line","svp","sound","mermaid"]))
        if not wanted or sz is None or sz>2_000_000: continue
        data=bytearray()
        try:
            f.retrbinary("RETR "+p,data.extend)
            b=bytes(data); h=hashlib.sha256(b).hexdigest()
            txt=b.decode("utf-8","replace")
            if txt.count("\x00")>10: txt=""
            extracts.append({"path":p,"size":len(b),"sha256":h,"text":txt[:50000]})
        except Exception as e:
            extracts.append({"path":p,"error":repr(e)})
finally:
    try:f.quit()
    except: f.close()

out={
 "artifact_id":"JANUS-KUSTO-CD169-EM12-OPEN-ARCHIVE-PROCESSING-INVENTORY-2026-09-22-v1.0",
 "root":ROOT,
 "file_count":len(files),
 "tree":tree,
 "small_textlike_extracts":extracts,
 "claim_ceiling":"SOURCE_INVENTORY_ONLY__NO_PROCESSING_CAUSALITY_CLAIM"
}
p=OUT/"JANUS-KUSTO-CD169-EM12-OPEN-ARCHIVE-PROCESSING-INVENTORY-2026-09-22-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "file_count":len(files),
 "paths":[x["path"] for x in files],
 "extracts":[{"path":x["path"],"size":x.get("size"),"sha256":x.get("sha256"),"text":x.get("text","")[:4000]} for x in extracts]
},indent=2,ensure_ascii=False))
