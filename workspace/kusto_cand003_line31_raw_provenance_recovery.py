#!/usr/bin/env python3
import ftplib, hashlib, json, os, re
from pathlib import Path

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI"
MAX_DEPTH=5
MAX_TEXT=3_000_000
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)

NAME_TERMS=["em12","raw","all","dat","svp","neptune","mermaid","line31","line_31","acceptl28","acceptl31","processing","process","readme"]
TEXT_TERMS=["em12","raw.all","line 31","line31","acceptl28-33","acceptl28","acceptl31","neptune","mermaid","svp","sound velocity","clock","processing"]
TEXT_EXT={".txt",".log",".lst",".list",".csv",".asc",".ascii",".info",".inf",".cfg",".conf",".dat",".nav",".svp",".md"}

def conn():
    f=ftplib.FTP(timeout=120); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I"); return f

def classify(f,path):
    cur=f.pwd()
    try:
        f.cwd(path); f.cwd(cur); return "dir",None
    except Exception:
        try:f.cwd(cur)
        except:pass
        try:return "file",f.size(path)
        except:return "file",None

def list_dir(path):
    f=conn()
    try:
        try:
            out=[]
            for name,facts in f.mlsd(path):
                if name in (".",".."):continue
                typ=facts.get("type")
                full=path.rstrip("/")+"/"+name
                if typ not in ("dir","file"):
                    typ,size=classify(f,full)
                else:
                    size=facts.get("size")
                    if size is not None:
                        try:size=int(size)
                        except:size=None
                out.append((name,typ,size,facts.get("modify")))
            return out
        except Exception:
            out=[]
            for p in f.nlst(path):
                name=p.rstrip("/").split("/")[-1]
                typ,size=classify(f,p)
                out.append((name,typ,size,None))
            return out
    finally:
        try:f.quit()
        except:f.close()

entries=[]
queue=[(ROOT,0)]
while queue:
    path,depth=queue.pop(0)
    try: rows=list_dir(path)
    except Exception as e:
        entries.append({"path":path,"depth":depth,"error":repr(e)}); continue
    for name,typ,size,modify in sorted(rows):
        full=path.rstrip("/")+"/"+name
        rec={"path":full,"name":name,"depth":depth+1,"type":typ,"size":size,"modify":modify}
        entries.append(rec)
        if typ=="dir" and depth+1<MAX_DEPTH:
            queue.append((full,depth+1))

interesting=[]
for r in entries:
    if r.get("type")!="file": continue
    n=r.get("name","").lower()
    ext=Path(n).suffix.lower()
    reasons=[]
    if any(t in n for t in NAME_TERMS): reasons.append("filename_term")
    if ext in {".all",".raw",".bin",".mb",".svp",".dat"}: reasons.append("rawlike_extension")
    if reasons:
        interesting.append({**r,"reasons":reasons})

# Text-content scan only for small textlike files; do not treat large XYZ bathymetry as text metadata.
text_hits=[]
text_scanned=0
for r in entries:
    if r.get("type")!="file" or r.get("size") is None: continue
    n=r["name"].lower(); ext=Path(n).suffix.lower()
    likely_text=(ext in TEXT_EXT or "readme" in n or "process" in n or "header" in n)
    if not likely_text or r["size"]>MAX_TEXT: continue
    # skip giant processed XYZ/ascii datasets unless name itself is metadata-like
    if (".xyz" in n or ext==".ascii") and not any(k in n for k in ["readme","process","log","info"]):
        continue
    f=conn(); data=bytearray()
    try:
        f.retrbinary("RETR "+r["path"],data.extend)
    except Exception as e:
        try:f.quit()
        except:pass
        text_hits.append({"path":r["path"],"error":repr(e)}); continue
    finally:
        try:f.quit()
        except:pass
    b=bytes(data); text_scanned+=1
    txt=b.decode("utf-8","replace")
    low=txt.lower()
    terms=[t for t in TEXT_TERMS if t in low]
    if terms:
        snippets=[]
        lines=txt.splitlines()
        for i,line in enumerate(lines):
            ll=line.lower()
            if any(t in ll for t in terms):
                snippets.append({"line":i+1,"text":line[:1000]})
                if len(snippets)>=80: break
        text_hits.append({
          "path":r["path"],"size":len(b),"sha256":hashlib.sha256(b).hexdigest(),
          "matched_terms":terms,"snippets":snippets
        })

raw_candidates=[x for x in interesting if "rawlike_extension" in x["reasons"] or any(k in x["name"].lower() for k in ["raw","em12","mermaid","neptune"])]
out={
 "artifact_id":"JANUS-KUSTO-CAND003-LINE31-DOMINANT-CELLS-RAW-PROVENANCE-RECOVERY-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-LINE31-DOMINANT-CELLS-RAW-PROVENANCE-RECOVERY-PREREG-2026-09-23-v1.0.json",
 "root":ROOT,
 "max_depth":MAX_DEPTH,
 "entry_count":len(entries),
 "file_count":sum(1 for x in entries if x.get("type")=="file"),
 "directory_count":sum(1 for x in entries if x.get("type")=="dir"),
 "interesting_filename_or_extension_count":len(interesting),
 "raw_candidate_count":len(raw_candidates),
 "text_files_scanned":text_scanned,
 "text_hit_count":len([x for x in text_hits if "matched_terms" in x]),
 "raw_candidates":raw_candidates,
 "text_hits":text_hits,
 "all_entries":entries,
 "provisional_verdict":"RAW_OR_LINE31_PROVENANCE_FOUND" if (raw_candidates or any("matched_terms" in x for x in text_hits)) else "NO_RAW_OR_LINE31_SPECIFIC_PROVENANCE_IN_ACCESSIBLE_BODC_ARCHIVE",
 "claim_ceiling":"ARCHIVAL_PROVENANCE_RECOVERY_ONLY"
}
p=OUT/"JANUS-KUSTO-CAND003-LINE31-DOMINANT-CELLS-RAW-PROVENANCE-RECOVERY-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "entry_count":out["entry_count"],"file_count":out["file_count"],"directory_count":out["directory_count"],
 "raw_candidate_count":out["raw_candidate_count"],"text_files_scanned":out["text_files_scanned"],
 "text_hit_count":out["text_hit_count"],"raw_candidates":raw_candidates[:200],"text_hits":text_hits[:200],
 "provisional_verdict":out["provisional_verdict"]
},indent=2,ensure_ascii=False))
