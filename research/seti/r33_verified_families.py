#!/usr/bin/env python3
"""R33: externally documented RFI occurrences, not inferred origin labels."""
import hashlib,json,math,pathlib,sys,re,os,datetime
from urllib.parse import unquote
import numpy as np
import r4_chunk_replay as base
import r9_expanded_search as eng
import r15c_r16c_spectral as fil
import r30c_hannah_memory as mem
P=pathlib.Path("research/seti/r33_verified_families_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r33");OUT.mkdir(parents=True,exist_ok=True);mem.OUT=OUT
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,v):mem.save(n,v)
def emit(k,v):mem.emit(k,v)
def parse_name(s):
    m=re.fullmatch(r"(blc\d+)_(.+)_([SR])_(2019-04-29T\d\d:\d\d:\d\d(?:\.\d+)?)\.fil",s,re.I)
    if not m:return None
    return {"backend":m[1],"target":m[2],"position":m[3],"utc":m[4],"time":datetime.datetime.fromisoformat(m[4]).timestamp()}
def covers(h,c):
    lo,hi=sorted((h["fch1"],h["fch1"]+(h["nchans"]-1)*h["foff"]))
    return lo<c-.0004 and hi>c+.0004
def freeze():
    sel=S["archive_selection"];rows=[];indexes=[];blocked=[]
    for q in sel["query_index"]:
        raw,url=base.fetch(eng.API,{"target":q["target"],"telescope":"Parkes","file-type":"HDF5"},limit=32*base.MIB)
        actual=sha(raw);assert actual==q["sha256"],"Archive index digest changed"
        indexes.append({"target":q["target"],"url":url,"sha256":actual})
        d=json.loads(raw);rs=d if isinstance(d,list) else next(d[k] for k in ("data","files","results") if isinstance(d.get(k),list))
        for v in rs:
            u=v.get("url","")
            if not isinstance(u,str) or not eng.allowed(u):continue
            if int(float(v.get("mjd",0)))!=sel["mjd"]:continue
            n=unquote(base.basename(u));p=parse_name(n)
            if p is None:continue
            if q["target"]=="ProxCen" and not n.endswith(sel["anchor_suffix"]):continue
            if q["target"]=="1421-490" and p["target"]!="1421-490":continue
            rows.append({"name":n,"url":u,"metadata":v,**p,"role":"ON" if q["target"]=="ProxCen" else "OFF"})
    rows=list({r["url"]:r for r in rows}.values())
    assert 0<len(rows)<=512
    on=[];cache={}
    def header(f):
        if f["url"] not in cache:
            with base.RangeFile(f["url"],1024*1024) as r:
                h=fil.parse_header(r);assert int(h["tstart"])==sel["mjd"] and int(f["metadata"]["size"])==r.size
                cache[f["url"]]=h
        return cache[f["url"]]
    for f in sorted((r for r in rows if r["role"]=="ON"),key=lambda r:r["name"]):
        try:on.append({**f,"header":header(f)})
        except Exception as e:blocked.append({"name":f["name"],"error":str(e)})
    selected=[]
    for target in S["targets"]:
        c=target["frequency_mhz"];eligible=[f for f in on if covers(f["header"],c)]
        if len(eligible)!=1:
            selected.append({"target":target,"status":"BLOCKED_NO_UNIQUE_ON_BAND","on_matches":[f["name"] for f in eligible]});continue
        f=eligible[0];off=[]
        pool=[r for r in rows if r["role"]=="OFF" and r["backend"]==f["backend"] and 0<abs(r["time"]-f["time"])<=sel["max_off_gap_seconds"]]
        for direction in (-1,1):
            pool2=sorted((r for r in pool if (r["time"]-f["time"])*direction>0),key=lambda r:abs(r["time"]-f["time"]))
            # Predeclared nearest metadata-eligible scan only, no intensity-based selection.
            for g in pool2[:4]:
                try:
                    gh=header(g);h=f["header"]
                    assert covers(gh,c) and (gh.get("src_raj"),gh.get("src_dej"))!=(h.get("src_raj"),h.get("src_dej"))
                    assert math.isclose(gh["foff"],h["foff"],rel_tol=1e-12) and math.isclose(gh["tsamp"],h["tsamp"],rel_tol=1e-8)
                    assert all(k in gh and k in h for k in ("src_raj","src_dej"))
                    off.append({**g,"header":gh});break
                except Exception as e:blocked.append({"name":g["name"],"error":str(e)})
        selected.append({"target":target,"status":"READY" if off else "BLOCKED_NO_OFF","on":f,"off":off})
    result={"schema":"janus.hannah.r33.frozen_source_manifest.v1","protocol_sha256":sha(P.read_bytes()),"indexes":indexes,"selected":selected,"blocked_headers":blocked,"no_intensity_read":True}
    save("manifest.json",result)
    emit("R33_MANIFEST",{"selected":[{"target":v["target"]["id"],"status":v["status"],"on":v.get("on",{}).get("name"),"off":[r["name"] for r in v.get("off",[])]} for v in selected],"blocked_headers":blocked,"indexes":indexes})
def peak(raw,t,c,tol):
    y=mem.norm(raw.astype(float));freq=t["channel0_mhz"]*1e6+np.arange(y.shape[1])*t["df_hz"]
    z=np.median(y,axis=0);occ=np.mean(y>3,axis=0)
    # A distinct spectral local maximum rather than all neighboring bins of one peak.
    ids=[i for i in range(4,len(freq)-4) if abs(freq[i]-c*1e6)<=tol+1e-6 and z[i]>=6 and occ[i]>=.75 and z[i]>z[i-1] and z[i]>z[i+1]]
    if len(ids)!=1:return None,{"eligible_peaks":len(ids),"max_median":float(max(z)),"reason":"NO_UNIQUE_FIXED_ZERO_DRIFT_PEAK"}
    i=ids[0];path=np.full(raw.shape[0],i,int)
    return path,{"eligible_peaks":1,"frequency_hz":float(freq[i]),"median_sigma":float(z[i]),"occupancy":float(occ[i]),"drift_hz_s":0.,"bin":i}
def diagnostic(raw,t,c):
    y=mem.norm(raw.astype(float));n,w=y.shape;freq=t["channel0_mhz"]*1e6+np.arange(w)*t["df_hz"];basebins=np.flatnonzero(abs(freq-c*1e6)<=100)
    found=[]
    # Frozen frequency prediction; both drift signs tested and retained, never admission.
    for drift in (-.06,-.04,-.028,-.02,-.01,0,.01,.02,.028,.04,.06):
        d=np.rint(np.arange(n)*t["ts"]*drift/t["df_hz"]).astype(int)
        for j in basebins:
            p=j+d
            if p.min()<4 or p.max()>=w-4:continue
            v=y[np.arange(n),p];found.append((float(np.median(v)),float(np.mean(v>3)),int(j),float(drift)))
    if not found:return None
    score,occ,j,dr=max(found)
    return {"frequency_hz":float(freq[j]),"median_sigma":score,"occupancy":occ,"trial_drift_hz_s":dr,"scope":"frequency prediction only; not published individual identity"}
def scan():
    m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==sha(P.read_bytes())
    parent=pathlib.Path(S["parent"]["path"]).read_bytes();assert sha(parent)==S["parent"]["sha256"]
    old,_=mem.append(json.loads(parent),[]);assert len(old["entries"])==169
    a=os.environ["R33_MANIFEST_ARTIFACT"];b=os.environ["R33_MANIFEST_DIGEST"].removeprefix("sha256:")
    assert a.isdecimal() and mem.ishash(b);artifact={"id":int(a),"digest":b}
    exact=[];diagnostics=[];blocked=[];receipts=[]
    for v in m["selected"]:
        t0=v["target"];c=t0["frequency_mhz"]
        if v["status"]!="READY":blocked.append({"target":t0["id"],"reason":v["status"]});continue
        pairs=[]
        for f in [v["on"]]+v["off"]:
            try:
                with base.RangeFile(f["url"],8*base.MIB) as r:
                    h=fil.parse_header(r);assert h==f["header"]
                    x,t=fil.fil_window(r,h,c,S["archive_selection"]["window_half_hz"])
                    assert x.shape[0]>=16 and sha(x.tobytes())==t["slice_sha256"]
                    name=t0["id"]+"-"+sha(f["url"].encode())[:16]+".npy";np.save(OUT/name,x,allow_pickle=False)
                    rec={"target":t0["id"],"scan":f["name"],"role":f["role"],"header":h,"url":f["url"],"tile":t,"artifact":name,"bytes_read":r.total}
                if t0["gate"]=="EXACT_PUBLISHED_ZERO_DRIFT":
                    p,info=peak(x,t,c,t0["tolerance_hz"]);rec["peak"]=info
                    if p is not None:
                        meta=mem.meta(t,f["name"],artifact);meta.update({"source_artifact_role":"FROZEN_METADATA_MANIFEST_NOT_RAW_TILE","archive_url":f["url"],"role":f["role"],"day":S["archive_selection"]["mjd"],"source_telescope":"Parkes/Murriyang"})
                        passport=mem.make_record(x,p,meta);pairs.append((rec,passport))
                else:rec["diagnostic"]=diagnostic(x,t,c);diagnostics.append({k:rec[k] for k in ("target","scan","role","diagnostic")})
                receipts.append(rec);emit("R33_SCAN",rec)
            except Exception as e:blocked.append({"target":t0["id"],"scan":f["name"],"reason":str(e)})
        if t0["gate"]!="EXACT_PUBLISHED_ZERO_DRIFT":continue
        on=[(r,p) for r,p in pairs if r["role"]=="ON"];off=[(r,p) for r,p in pairs if r["role"]=="OFF"]
        if len(on)!=1 or not off:
            blocked.append({"target":t0["id"],"reason":"NO_PAIRED_PUBLISHED_PEAK"});continue
        rr,pp=on[0];good=[(r,p) for r,p in off if abs(p["frequency_hz"]-pp["frequency_hz"])<=4.000001 and p["provenance"]["source_tile_sha256"]!=pp["provenance"]["source_tile_sha256"]]
        if not good:blocked.append({"target":t0["id"],"reason":"NO_EXACT_DISTINCT_OFF_MATCH"});continue
        for r,p in [on[0]]+good:
            evidence={"schema":"janus.hannah.external_label_wrapper.v1","label":"VERIFIED_RFI_CANDIDATE_REVIEW_REQUIRED","family":t0["family"],"occurrence_group":t0["id"]+"_20190429","source_doi":S["publication"]["doi"],"source_section":t0["source_section"],"published_frequency_mhz":c,"scope":"Published central interferer occurrence only; no blanket lookalike membership or hardware attribution","paired_scans":[rr["scan"]]+[x["scan"] for x,_ in good],"paired_tile_sha256":[pp["provenance"]["source_tile_sha256"]]+[q["provenance"]["source_tile_sha256"] for _,q in good],"passport":p,"measured_peak":r["peak"]}
            evidence["evidence_sha256"]=mem.digest(evidence);exact.append(evidence)
    assert len({v["passport"]["id"] for v in exact})==len(exact)
    save("raw_receipts.json",receipts);save("external_label_candidates.json",{"schema":"janus.hannah.r33.external_label_candidates.v1","entries":exact,"parent_sha256":S["parent"]["sha256"],"authority":"REVIEW_REQUIRED_NO_VETO"})
    save("diagnostics.json",diagnostics)
    d={"stage":"R33A","status":"PASS_SOURCE_JOIN_CANDIDATES" if exact else "HOLD_NO_EXACT_SOURCE_JOIN","candidate_occurrences":len(exact),"candidate_family_groups":len({e["family"] for e in exact}),"families":sorted({e["family"] for e in exact}),"raw_tiles":len(receipts),"diagnostics":len(diagnostics),"blocked":blocked,"parent_entries":169,"verified_in_parent":1,"verified_promoted":0,"candidate_sha256":sha((OUT/"external_label_candidates.json").read_bytes()),"manifest_sha256":sha((OUT/"manifest.json").read_bytes()),"production_adoption":False}
    save("result.json",d);emit("R33_SUMMARY",d)
def selftest():
    rng=np.random.default_rng(33001);x=rng.normal(size=(32,201)).astype(np.float32);x[:,100]+=20.;t={"channel0_mhz":1399.9996,"df_hz":4.,"ts":18.}
    p,d=peak(x,t,1400.,4.);assert p is not None and np.all(p==100)
    p,d=peak(rng.normal(size=x.shape).astype(np.float32),t,1400.,4.);assert p is None
    return {"exact_peak":"PASS","noise_abstain":"PASS"}
if __name__=="__main__":
    phase=sys.argv[1];{"freeze":freeze,"scan":scan,"selftest":lambda:emit("R33_SELFTEST",selftest())}[phase]()
