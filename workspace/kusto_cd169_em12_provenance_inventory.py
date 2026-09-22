#!/usr/bin/env python3
import ftplib, hashlib, json, os, posixpath, re
from pathlib import Path

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/EM12"
FOCUS={"B1-81-1_Acceptl28-33.xyz.ascii","B1-81-1_Acceptl15-22.xyz.ascii"}
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)

def conn():
    f=ftplib.FTP(timeout=120);f.connect(HOST,21);f.login("anonymous","janus-kusto@example.invalid");f.voidcmd("TYPE I");return f

def list_dir(path):
    f=conn()
    try:
        try:
            rows=[]
            for name,facts in f.mlsd(path):
                if name in (".",".."):continue
                rows.append((name,facts.get("type"),facts))
            return rows
        except Exception:
            names=f.nlst(path);out=[]
            for p in names:
                name=p.rstrip("/").split("/")[-1]
                try:
                    size=f.size(p)
                    typ="file"
                except Exception:
                    size=None;typ="unknown"
                out.append((name,typ,{"size":size}))
            return out
    finally:
        try:f.quit()
        except: f.close()

def fetch_bytes(path):
    f=conn();chunks=[]
    try:
        f.retrbinary("RETR "+path,chunks.append)
        return b"".join(chunks)
    finally:
        try:f.quit()
        except:f.close()

rows=list_dir(ROOT)
inventory=[]
for name,typ,facts in rows:
    full=ROOT+"/"+name
    size=None
    try:size=int(facts.get("size")) if facts.get("size") is not None else None
    except:size=None
    inventory.append({"name":name,"path":full,"type":typ,"size":size})

focus={}
for item in inventory:
    if item["name"] not in FOCUS: continue
    b=fetch_bytes(item["path"])
    text=b.decode("ascii","replace")
    parsed=[]
    for ln,line in enumerate(text.splitlines(),1):
        p=line.replace(","," ").split()
        if len(p)<3: continue
        try: lon=float(p[0]);lat=float(p[1]);z=float(p[2])
        except: continue
        parsed.append((ln,lon,lat,z))
    focus[item["name"]]={
      "sha256":hashlib.sha256(b).hexdigest(),
      "bytes":len(b),
      "parsed_xyz_count":len(parsed),
      "first_5_rows":parsed[:5],
      "last_5_rows":parsed[-5:],
      "lon_min_max":[min(x[1] for x in parsed),max(x[1] for x in parsed)] if parsed else None,
      "lat_min_max":[min(x[2] for x in parsed),max(x[2] for x in parsed)] if parsed else None,
      "depth_min_max":[min(x[3] for x in parsed),max(x[3] for x in parsed)] if parsed else None
    }

small_text=[]
for item in inventory:
    n=item["name"].lower()
    if item["name"] in FOCUS: continue
    if item["size"] is not None and item["size"]<=2_000_000 and (n.endswith((".txt",".log",".cfg",".readme",".asc",".ascii")) or "read" in n or "line" in n or "proc" in n):
        try:
            b=fetch_bytes(item["path"])
            small_text.append({
              "name":item["name"],"sha256":hashlib.sha256(b).hexdigest(),
              "preview":b.decode("utf-8","replace")[:12000]
            })
        except Exception as e:
            small_text.append({"name":item["name"],"error":repr(e)})

out={
 "artifact_id":"JANUS-KUSTO-CD169-EM12-PRODUCT-PROVENANCE-INVENTORY-2026-09-22-v1.0",
 "root":ROOT,
 "inventory_count":len(inventory),
 "inventory":sorted(inventory,key=lambda x:x["name"]),
 "focus_files":focus,
 "small_text_metadata":small_text,
 "external_report_binding":{
   "report":"RRS Charles Darwin Cruise 169 report",
   "appendix_D_note":"CD169 logged and processed 85 EM12 swath lines; line 29 was active in the ~4deg04S/12deg16W sector and lines 28-33 span the relevant Leg-1 survey interval.",
   "processing_warning":"Cruise report states newly acquired EM12 data were processed and merged with existing data for bathymetry plots. File-level provenance must therefore not be assumed solely from directory placement."
 },
 "claim_ceiling":"ARCHIVE_PRODUCT_INVENTORY_AND_FILENAME/EXTENT_PROVENANCE_ONLY"
}
p=OUT/"JANUS-KUSTO-CD169-EM12-PRODUCT-PROVENANCE-INVENTORY-2026-09-22-v1.0.json"
p.write_text(json.dumps(out,indent=2),encoding="utf-8")
print(json.dumps({
 "inventory_count":len(inventory),
 "names":[x["name"] for x in sorted(inventory,key=lambda x:x["name"])],
 "focus_files":focus,
 "small_text_metadata":small_text
},indent=2))
