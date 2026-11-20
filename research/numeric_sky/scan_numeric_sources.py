#!/usr/bin/env python3
"""JANUS NUMERIC SKY AUDIT R0: code-only numeric corpus with provenance.
Never interprets arbitrary numeric literals as real celestial directions.
"""
import argparse
import ast
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import token
import tokenize

ROOT=Path("runtime/numeric_sky")
PROTO=json.loads(Path("research/numeric_sky/audit_protocol.json").read_text())
EXS=set(PROTO["scope"]["source_extensions"])
DENY=set(PROTO["scope"]["excluded_path_segments"])
LIMIT=int(PROTO["scope"]["max_file_size_bytes"])
MAX_TOTAL=250*1024*1024
RENUM=re.compile(r"(?<![\w.])(?:0[xX][0-9a-fA-F]+|0[bB][01]+|0[oO][0-7]+|(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)(?![\w.])")
STRINGS=re.compile(r'(?:"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\')')
RARE=re.compile(r"(?:^|_)(?:ra|right_ascension)(?:_|$)",re.I)
DECRE=re.compile(r"(?:^|_)(?:dec|declination)(?:_|$)",re.I)

def sha(b):return hashlib.sha256(b).hexdigest()
def run(cmd, **kw):
    a=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,**kw)
    if a.returncode:raise RuntimeError("exit="+str(a.returncode)+"; stderr="+a.stderr[-500:])
    return a.stdout
def numeric(v):
    if isinstance(v,bool) or not isinstance(v,(int,float)):return None
    try:
        f=float(v)
        if math.isfinite(f):return f
    except (ValueError,OverflowError,TypeError):pass
    return None
def emit_num(out,repo,path,rev,src_sha,kind,label,value,lineno=None,context=None):
    v=numeric(value)
    if v is None:return False
    o={"repo":repo,"path":path,"revision":rev,"source_sha256":src_sha,"kind":kind,
       "field":str(label)[:220],"number":v,"line":lineno,"context":(context or "")[:180]}
    out.write((json.dumps(o,ensure_ascii=False,separators=(",",":"))+"\n").encode())
    return True
def make_coord(repo,path,rev,src_sha,fields,loc,kind):
    a=b=None
    for k,v in fields.items():
        if RARE.search(str(k)) and a is None:a=numeric(v)
        if DECRE.search(str(k)) and b is None:b=numeric(v)
    if a is None or b is None:return None
    if not (0<=a<360 and -90<=b<=90):return None
    # Frame and unit verification is required for scientific use. Units stay visible.
    lc=str(loc).lower()+" "+str(path).lower()
    context=(str(fields.get("frame",""))+" "+str(fields.get("epoch",""))+" "+str(loc)).lower()
    strong=bool(re.search(r"icrs|j2000|sky|galact|astronom|cosmos|target|orion|osiris|love|edem|abydos|sirius|beacon|star|catalog|wise|gaia",lc+" "+context))
    deg=bool(any("deg" in str(k).lower() or "j2000" in str(k).lower() for k in fields) or fields.get("unit") in ("degree","degrees","deg"))
    return {"ra_deg":a,"dec_deg":b,"repo":repo,"path":path,"revision":rev,"source_sha256":src_sha,
            "origin":str(loc)[:250],"kind":kind,"frame":str(fields.get("frame",fields.get("epoch","UNKNOWN"))),
            "unit_evidence":deg,"semantic_evidence":strong,
            "scientific_qualified":bool((deg or "icrs" in context or "j2000" in context) and strong)}
def constnum(n):
    if isinstance(n,ast.Constant):return numeric(n.value)
    if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.USub,ast.UAdd)):
        v=constnum(n.operand);return (-v if isinstance(n.op,ast.USub) else v) if v is not None else None
    if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Mult)):
        l,r=n.left,n.right
        if isinstance(l,ast.Attribute) and l.attr=="deg":return constnum(r)
        if isinstance(r,ast.Attribute) and r.attr=="deg":return constnum(l)
    return None
def json_scan(value,prefix,emit,cands,repo,path,rev,src_sha):
    if isinstance(value,dict):
        for k,v in value.items():
            if numeric(v) is not None:emit("json_number",prefix+"."+str(k) if prefix else str(k),v)
        if any(RARE.search(str(x)) for x in value):
            c=make_coord(repo,path,rev,src_sha,value,prefix,"json_named_pair")
            if c:cands.append(c)
        if "unit_vector_icrs" in value and isinstance(value["unit_vector_icrs"],dict):
            v=value["unit_vector_icrs"];x,y,z=(numeric(v.get(k)) for k in ("x","y","z"))
            if None not in (x,y,z):
                norm=math.sqrt(x*x+y*y+z*z)
                if norm>0 and abs(norm-1)<.0001:
                    c={"ra_deg":math.degrees(math.atan2(y,x))%360,"dec_deg":math.degrees(math.asin(max(-1,min(1,z/norm)))),"repo":repo,"path":path,
                       "revision":rev,"source_sha256":src_sha,"origin":prefix+".unit_vector_icrs","kind":"icrs_unit_vector",
                       "frame":"ICRS","unit_evidence":True,"semantic_evidence":True,"scientific_qualified":True}
                    cands.append(c)
        for k,v in value.items():json_scan(v,prefix+"."+str(k) if prefix else str(k),emit,cands,repo,path,rev,src_sha)
    elif isinstance(value,list):
        for i,v in enumerate(value):
            if numeric(v) is not None:emit("json_number",prefix+"["+str(i)+"]",v)
            else:json_scan(v,prefix+"["+str(i)+"]",emit,cands,repo,path,rev,src_sha)
def py_scan(s,emit,cands,repo,path,rev,src_sha):
    n=0;lines=s.splitlines()
    try:
        toks=tokenize.generate_tokens(io.StringIO(s).readline)
        for t in toks:
            if t.type!=token.NUMBER:continue
            try:v=ast.literal_eval(t.string)
            except (ValueError,SyntaxError):continue
            if numeric(v) is not None:
                context=lines[t.start[0]-1] if t.start[0]<=len(lines) else ""
                n+=bool(emit("python_number",str(t.string),v,t.start[0],context))
    except (tokenize.TokenError,IndentationError) as e:
        pass
    try:tree=ast.parse(s)
    except (SyntaxError,ValueError,MemoryError):return n
    globals_={}
    for node in tree.body:
        if isinstance(node,(ast.Assign,ast.AnnAssign)):
            names=[x.id for x in node.targets if isinstance(x,ast.Name)] if isinstance(node,ast.Assign) else [node.target.id] if isinstance(node.target,ast.Name) else []
            value=constnum(node.value) if node.value is not None else None
            if value is not None:
                for k in names:globals_[k]=value
    for node in ast.walk(tree):
        if isinstance(node,ast.Dict):
            f={}
            for k,v in zip(node.keys,node.values):
                if isinstance(k,ast.Constant) and isinstance(k.value,str):
                    z=constnum(v)
                    if z is not None:f[k.value]=z
                    elif isinstance(v,ast.Constant) and isinstance(v.value,str):f[k.value]=v.value
            if f:
                c=make_coord(repo,path,rev,src_sha,f,"line:"+str(node.lineno),"python_dictionary")
                if c:cands.append(c)
        if isinstance(node,ast.Call):
            fname=node.func.id if isinstance(node.func,ast.Name) else node.func.attr if isinstance(node.func,ast.Attribute) else ""
            if fname=="SkyCoord":
                va=constnum(node.args[0]) if len(node.args)>0 else None
                vb=constnum(node.args[1]) if len(node.args)>1 else None
                kw={x.arg:x.value for x in node.keywords}
                if va is None and "ra" in kw:va=constnum(kw["ra"])
                if vb is None and "dec" in kw:vb=constnum(kw["dec"])
                if va is not None and vb is not None:
                    cc=make_coord(repo,path,rev,src_sha,{"ra_deg":va,"dec_deg":vb,"frame":"ICRS"},"SkyCoord:line:"+str(node.lineno),"python_skycoord")
                    if cc:cands.append(cc)
    ra=next((v for k,v in globals_.items() if RARE.search(k) and ("deg" in k.lower() or k.lower()=="ra")),None)
    dec=next((v for k,v in globals_.items() if DECRE.search(k) and ("deg" in k.lower() or k.lower()=="dec")),None)
    if ra is not None and dec is not None:
        c=make_coord(repo,path,rev,src_sha,{"ra_deg":ra,"dec_deg":dec,"frame":"ICRS"},"python_module_constants","python_named_constants")
        if c:cands.append(c)
    return n
def lexical_scan(s,emit):
    n=0
    for i,line in enumerate(s.splitlines(),1):
        if len(line)>20000:continue
        clean=STRINGS.sub("",line)
        for m in RENUM.finditer(clean):
            try:
                raw=m.group(0); v=int(raw,0) if raw.lower().startswith(("0x","0b","0o")) else float(raw)
            except (ValueError,OverflowError):continue
            n+=bool(emit("lexical_number",raw,v,i,line))
    return n
def source_filter(path,size):
    parts=Path(path).parts
    if any(p in DENY for p in parts):return "excluded_vendor"
    if any(path.startswith(p) for p in PROTO["scope"].get("excluded_generated_paths",[])):return "excluded_generated"
    if Path(path).suffix not in EXS:return "not_code_extension"
    if size>LIMIT or size<0:return "excluded_size"
    return "selected"
def scan_repository(repo,base,number_sink):
    record={"repo":repo,"status":"NOT_STARTED","selected_files":0,"numbers":0,"candidate_pairs":0,"head_sha":None,
            "skipped":{"excluded_vendor":0,"excluded_generated":0,"not_code_extension":0,"excluded_size":0,"size_budget":0},"errors":[]}
    dest=base/("clone_"+str(abs(hash(repo)))[:10])
    cands=[]
    try:
        run(["git","clone","-q","--depth","1","--filter=blob:none","--no-checkout","--single-branch","https://github.com/"+repo+".git",str(dest)],timeout=420)
        record["head_sha"]=run(["git","-C",str(dest),"rev-parse","HEAD"],timeout=15).strip()
        listing=run(["git","-C",str(dest),"ls-tree","-r","-l","HEAD"],timeout=60)
        selected=[];total=0
        for line in listing.splitlines():
            try:
                meta,path=line.split("\t",1)
                a=meta.split();size=int(a[-1]) if a[-1]!="-" else -1
            except (ValueError,IndexError):continue
            st=source_filter(path,size)
            if st!="selected":record["skipped"][st]+=1;continue
            if total+size>MAX_TOTAL:record["skipped"]["size_budget"]+=1;continue
            selected.append(path);total+=size
        record["eligible_bytes"]=total;record["selected_files"]=len(selected)
        if selected:
            run(["git","-C",str(dest),"sparse-checkout","init","--no-cone"],timeout=60)
            run(["git","-C",str(dest),"sparse-checkout","set","--no-cone","--stdin"],input="\n".join(selected)+"\n",timeout=900)
            run(["git","-C",str(dest),"checkout","-f","HEAD"],timeout=900)
        for path in selected:
            p=dest/path
            try:
                b=p.read_bytes()
                if len(b)>LIMIT:raise ValueError("unexpected_size_after_checkout")
                s=b.decode("utf-8")
                source_hash=sha(b)
                def emit(kind,label,val,line=None,context=None):
                    ok=emit_num(number_sink,repo,path,record["head_sha"],source_hash,kind,label,val,line,context)
                    if ok:record["numbers"]+=1
                    return ok
                if path.endswith(".json"):
                    try:json_scan(json.loads(s),"",emit,cands,repo,path,record["head_sha"],source_hash)
                    except (json.JSONDecodeError,RecursionError) as e:record["errors"].append({"path":path,"error":"json_parse:"+str(e)[:100]})
                elif path.endswith(".ipynb"):
                    try:
                        nb=json.loads(s)
                        for ci,cell in enumerate(nb.get("cells",[])):
                            if cell.get("cell_type")=="code":py_scan("".join(cell.get("source",[])),emit,cands,repo,path,record["head_sha"],source_hash)
                    except (json.JSONDecodeError,RecursionError) as e:record["errors"].append({"path":path,"error":"notebook_parse:"+str(e)[:100]})
                elif path.endswith(".py"):
                    py_scan(s,emit,cands,repo,path,record["head_sha"],source_hash)
                else:lexical_scan(s,emit)
            except (UnicodeError,ValueError,OSError) as e:
                record["errors"].append({"path":path,"error":str(e)[:130]})
        record["candidate_pairs"]=len(cands)
        record["status"]="PARTIAL" if record["errors"] or record["skipped"]["size_budget"] else "COMPLETE_SELECTED_CODE"
    except Exception as e:
        record["status"]="FAILED"
        record["errors"].append({"stage":"clone_or_checkout","error":str(e)[:500]})
    finally:
        shutil.rmtree(dest,ignore_errors=True)
    return record,cands
def main():
    pa=argparse.ArgumentParser();pa.add_argument("--shard",type=int,required=True);pa.add_argument("--shards",type=int,default=8)
    args=pa.parse_args();assert 0<=args.shard<args.shards
    allrepos=PROTO["provenance"]["public_repositories"]
    targets=[v for i,v in enumerate(allrepos) if i%args.shards==args.shard]
    outdir=ROOT/f"shard-{args.shard:02d}";outdir.mkdir(parents=True,exist_ok=True)
    inventory=[];candidates=[]
    with gzip.open(outdir/"numbers.jsonl.gz","wb",compresslevel=6) as sink:
        for repo in targets:
            entry,cand=scan_repository(repo,outdir,sink)
            inventory.append(entry);candidates.extend(cand)
            print("NUMSKY_REPO",json.dumps({k:entry[k] for k in ("repo","status","head_sha","selected_files","numbers","candidate_pairs")}),flush=True)
    (outdir/"inventory.json").write_text(json.dumps({"shard":args.shard,"shards":args.shards,"repos":inventory,"private_excluded":PROTO["provenance"]["not_scannable_private"]},indent=2))
    (outdir/"explicit_icrs_candidates.json").write_text(json.dumps({"shard":args.shard,"candidates":candidates},indent=2,allow_nan=False))
    verdict={"shard":args.shard,"repositories":len(inventory),"complete":sum(r["status"]=="COMPLETE_SELECTED_CODE" for r in inventory),
             "failed":sum(r["status"]=="FAILED" for r in inventory),"selected_files":sum(r["selected_files"] for r in inventory),
             "numbers":sum(r["numbers"] for r in inventory),"coordinate_candidates":len(candidates)}
    (outdir/"audit_summary.json").write_text(json.dumps(verdict,indent=2))
    print("NUMSKY_SHARD",json.dumps(verdict),flush=True)
    if verdict["failed"]:print("NUMSKY_SHARD_REPO_ERRORS_RECORDED",flush=True)
if __name__=="__main__":main()
