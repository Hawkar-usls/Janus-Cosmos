#!/usr/bin/env python3
from __future__ import annotations
import ftplib, hashlib, json, math, os, posixpath, re
from datetime import datetime, timedelta, timezone
from pathlib import Path

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281"
NAV=ROOT+"/Nav/cd169leg1_1min.listit"
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
M=111320.0; LAT0=-4.0; COS0=math.cos(math.radians(LAT0))
CELLS=[
 {"id":"DOMINANT_12","lat":-4.0165853,"lon":-12.2985657},
 {"id":"DOMINANT_10","lat":-4.0152462,"lon":-12.299318}
]
LINE31_START="2005-02-26T19:11:00Z";LINE31_END="2005-02-27T14:40:00Z"
TOKENS=["raw","mermaid","merlin","opu","em12","line31","line_31","cd169_31","svp","soundvelocity","sound_velocity","neptune","beam","swath","proc","process","history"]
EXTS={".all",".raw",".dat",".bin",".log",".lst",".list",".txt",".cfg",".conf",".svp",".vel",".nav",".ascii",".xyz"}
PATTERNS=["line 31","line31","cd169_31","mermaid","opu","neptune","svp","sound velocity","raw","em12"]
MAX_DEPTH=8;MAX_TEXT=2_000_000

def xy(lon,lat):return lon*M*COS0,lat*M
def ts(s):return datetime.fromisoformat(s.replace("Z","+00:00")).timestamp()
def iso(t):return datetime.fromtimestamp(t,tz=timezone.utc).isoformat().replace("+00:00","Z")

def conn():
    f=ftplib.FTP(timeout=120);f.connect(HOST,21);f.login("anonymous","janus-kusto@example.invalid");f.voidcmd("TYPE I");return f

def list_dir(f,path):
    try:
        rows=[]
        for name,facts in f.mlsd(path):
            if name in (".",".."):continue
            rows.append((name,facts.get("type","unknown"),facts))
        return rows
    except Exception:
        rows=[]
        for p in f.nlst(path):
            name=p.rstrip("/").split("/")[-1]
            typ="unknown";size=None
            try:
                size=f.size(p);typ="file"
            except:pass
            rows.append((name,typ,{"size":size}))
        return rows

def walk(f,path,depth=0):
    out=[]
    if depth>MAX_DEPTH:return out
    for name,typ,facts in list_dir(f,path):
        full=path.rstrip("/")+"/"+name
        size=None
        try:size=int(facts.get("size")) if facts.get("size") not in (None,"") else None
        except:pass
        rec={"path":full,"name":name,"type":typ,"size":size,"depth":depth}
        out.append(rec)
        if typ=="dir" and depth<MAX_DEPTH:
            try:out.extend(walk(f,full,depth+1))
            except Exception as e:rec["walk_error"]=repr(e)
    return out

def fetch(f,path):
    chunks=[];f.retrbinary("RETR "+path,chunks.append);return b"".join(chunks)

f=conn()
try:
    inventory=walk(f,ROOT,0)
    navb=fetch(f,NAV)
    name_matches=[]
    text_matches=[]
    for rec in inventory:
        if rec.get("type")=="dir":continue
        low=rec["name"].lower();ext=Path(low).suffix
        hits=sorted({t for t in TOKENS if t in low})
        if hits or ext in EXTS:
            name_matches.append({**rec,"token_hits":hits,"extension":ext})
        if rec.get("size") is not None and rec["size"]<=MAX_TEXT and (hits or ext in {".txt",".log",".lst",".list",".cfg",".conf",".svp",".vel",".nav"}):
            try:
                b=fetch(f,rec["path"])
                txt=b.decode("utf-8","replace")
                found=[]
                lines=txt.splitlines()
                for k,p in enumerate(PATTERNS):
                    pl=p.lower()
                    idx=[i for i,line in enumerate(lines) if pl in line.lower()]
                    if idx:
                        found.append({"pattern":p,"line_numbers":[i+1 for i in idx[:50]],"contexts":[lines[max(0,i-1):min(len(lines),i+2)] for i in idx[:10]]})
                if found:text_matches.append({"path":rec["path"],"size":len(b),"sha256":hashlib.sha256(b).hexdigest(),"matches":found})
            except Exception as e:
                text_matches.append({"path":rec["path"],"error":repr(e)})
finally:
    try:f.quit()
    except:f.close()

# Parse frozen 1-minute nav.
base=datetime(2005,1,1,tzinfo=timezone.utc);rows=[]
for raw in navb.decode("ascii","ignore").splitlines():
    p=raw.split()
    if len(p)<7:continue
    try:
        jd=int(p[1]);hh,mm,ss=map(int,p[2].split(":"));lat=float(p[3]);lon=float(p[5])
    except:continue
    dt=base+timedelta(days=jd-1,hours=hh,minutes=mm,seconds=ss)
    x,y=xy(lon,lat);rows.append({"t":dt.timestamp(),"lat":lat,"lon":lon,"x":x,"y":y})
rows.sort(key=lambda r:r["t"])
a,b=ts(LINE31_START),ts(LINE31_END);line=[r for r in rows if a<=r["t"]<=b]

def project_time(cell):
    px,py=xy(cell["lon"],cell["lat"])
    best=None
    for i,(u,v) in enumerate(zip(line,line[1:])):
        ax,ay=u["x"],u["y"];bx,by=v["x"],v["y"]
        vx,vy=bx-ax,by-ay;den=vx*vx+vy*vy
        q=0.0 if den==0 else ((px-ax)*vx+(py-ay)*vy)/den
        q=max(0.0,min(1.0,q))
        qx,qy=ax+q*vx,ay+q*vy
        d=math.hypot(px-qx,py-qy)
        if best is None or d<best["cross_track_m"]:
            t=u["t"]+q*(v["t"]-u["t"])
            lat=u["lat"]+q*(v["lat"]-u["lat"])
            lon=u["lon"]+q*(v["lon"]-u["lon"])
            best={"cross_track_m":d,"projected_ship_track_utc":iso(t),"segment_index":i,
                  "segment_start_utc":iso(u["t"]),"segment_end_utc":iso(v["t"]),
                  "segment_fraction":q,"projected_ship_lat":lat,"projected_ship_lon":lon}
    return {**cell,**best}

bindings=[project_time(c) for c in CELLS]

out={
 "artifact_id":"JANUS-KUSTO-CD169-EM12-DOMINANT-CELLS-RAW-PROVENANCE-HUNT-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CD169-EM12-DOMINANT-CELLS-RAW-PROVENANCE-HUNT-PREREG-2026-09-23-v1.0.json",
 "archive_root":ROOT,
 "inventory_count":len(inventory),
 "directory_count":sum(r.get("type")=="dir" for r in inventory),
 "file_count":sum(r.get("type")!="dir" for r in inventory),
 "inventory":inventory,
 "name_or_extension_matches":name_matches,
 "small_text_pattern_matches":text_matches,
 "navigation":{"path":NAV,"sha256":hashlib.sha256(navb).hexdigest(),"parsed_rows":len(rows),"line31_rows":len(line)},
 "dominant_cell_line31_time_bindings":bindings,
 "warnings":[
   "Projected ship-track times are geometric nearest-track interpolations, not raw EM12 ping times.",
   "Absence from this accessible tree does not establish absence from NOC/BODC offline or originator holdings.",
   "TOBI-originator processed/unprocessed holdings are not assumed to contain shipborne EM12 raw data."
 ],
 "claim_ceiling":"ACCESSIBLE_ARCHIVE_RAW_PROVENANCE_DISCOVERY_AND_APPROXIMATE_LINE31_TIME_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-CD169-EM12-DOMINANT-CELLS-RAW-PROVENANCE-HUNT-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "inventory_count":out["inventory_count"],"directory_count":out["directory_count"],"file_count":out["file_count"],
 "matched_paths":[{"path":r["path"],"size":r.get("size"),"hits":r.get("token_hits"),"ext":r.get("extension")} for r in name_matches],
 "text_matches":text_matches,
 "time_bindings":bindings
},indent=2,ensure_ascii=False))
