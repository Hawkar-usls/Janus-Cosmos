#!/usr/bin/env python3
import ftplib, json
from pathlib import Path

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281"
MAX_DEPTH=4
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)

def conn():
    f=ftplib.FTP(timeout=120);f.connect(HOST,21);f.login("anonymous","janus-kusto@example.invalid");f.voidcmd("TYPE I");return f

def mlsd(path):
    f=conn()
    try:
        rows=[]
        try:
            for name,facts in f.mlsd(path):
                if name in (".",".."):continue
                rows.append((name,facts))
            return rows
        except Exception:
            out=[]
            for p in f.nlst(path):
                name=p.rstrip("/").split("/")[-1]
                try:
                    size=f.size(p);typ="file"
                except:
                    size=None;typ="unknown"
                out.append((name,{"type":typ,"size":size}))
            return out
    finally:
        try:f.quit()
        except:f.close()

seen=[]
queue=[(ROOT,0)]
while queue:
    path,depth=queue.pop(0)
    try: rows=mlsd(path)
    except Exception as e:
        seen.append({"path":path,"depth":depth,"error":repr(e)});continue
    for name,facts in sorted(rows):
        full=path.rstrip("/")+"/"+name
        typ=facts.get("type","unknown")
        rec={"path":full,"name":name,"depth":depth+1,"type":typ,
             "size":facts.get("size"),"modify":facts.get("modify")}
        seen.append(rec)
        if depth+1<MAX_DEPTH and typ=="dir":
            queue.append((full,depth+1))

interesting=[]
for r in seen:
    n=r.get("name","").lower()
    if any(k in n for k in ["em12","raw","line","neptune","simrad","accept","proc","xyz","all","81","28","29","30","31","32","33"]):
        interesting.append(r)

out={
 "artifact_id":"JANUS-KUSTO-CD169-SD11281-ARCHIVE-TREE-INVENTORY-2026-09-22-v1.0",
 "root":ROOT,"max_depth":MAX_DEPTH,
 "entry_count":len(seen),
 "entries":seen,
 "interesting_entries":interesting,
 "claim_ceiling":"ARCHIVE_PATH_INVENTORY_ONLY"
}
p=OUT/"JANUS-KUSTO-CD169-SD11281-ARCHIVE-TREE-INVENTORY-2026-09-22-v1.0.json"
p.write_text(json.dumps(out,indent=2),encoding="utf-8")
print(json.dumps({"entry_count":len(seen),"interesting_entries":interesting},indent=2))
