#!/usr/bin/env python3
"""Bounded figure-guided drift replay. Raw tracks are not automatic origin labels."""
import hashlib,json,math,os,pathlib,sys
from urllib.parse import unquote
import numpy as np
import r33_verified_families as old
import r4_chunk_replay as base
import r9_expanded_search as eng
import r15c_r16c_spectral as fil
import r30c_hannah_memory as mem
P=pathlib.Path("research/seti/r33b2_drift_replay_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r33b2");OUT.mkdir(parents=True,exist_ok=True);mem.OUT=OUT
A=S["acquire"];G=S["search"]
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,v):mem.save(n,v)
def emit(k,v):mem.emit(k,v)
def freeze():
    rows=[];indices=[];errors=[]
    for q in S["index"]:
        raw,url=base.fetch(eng.API,{"target":q["target"],"telescope":"Parkes","file-type":"HDF5"},limit=32*base.MIB)
        assert sha(raw)==q["sha256"],"Archive index drift"
        indices.append({"url":url,"sha256":sha(raw)})
        d=json.loads(raw);rs=d if isinstance(d,list) else next(d[k] for k in ("data","files","results") if isinstance(d.get(k),list))
        for v in rs:
            u=v.get("url","")
            if not isinstance(u,str) or not eng.allowed(u) or int(float(v.get("mjd",0)))!=58602:continue
            n=unquote(base.basename(u));p=old.parse_name(n)
            if not p:continue
            if q["target"]=="ProxCen" and not n.endswith(S["first_on"]):continue
            if q["target"]=="1421-490" and p["target"]!="1421-490":continue
            rows.append({"name":n,"url":u,"metadata":v,**p,"role":"ON" if q["target"]=="ProxCen" else "OFF"})
    rows=list({r["url"]:r for r in rows}.values());assert 0<len(rows)<=512
    cached={}
    def read_header(f):
        if f["url"] not in cached:
            with base.RangeFile(f["url"],1024*1024) as r:
                h=fil.parse_header(r);assert int(h["tstart"])==58602 and int(f["metadata"]["size"])==r.size
                cached[f["url"]]=h
        return cached[f["url"]]
    ons=[]
    for f in sorted((f for f in rows if f["role"]=="ON"),key=lambda f:f["name"]):
        try:ons.append({**f,"header":read_header(f)})
        except Exception as e:errors.append({"scan":f["name"],"error":str(e)})
    def covers(h,c):
        lo,hi=sorted((h["fch1"],h["fch1"]+(h["nchans"]-1)*h["foff"]))
        return lo<c-A["half_hz"]/1e6 and hi>c+A["half_hz"]/1e6
    chosen=[]
    for q in S["targets"]:
        candidates=[f for f in ons if covers(f["header"],q["rounded_mhz"])]
        if len(candidates)!=1:
            chosen.append({"target":q,"status":"HOLD_NO_UNIQUE_BAND","matches":[f["name"] for f in candidates]});continue
        f=candidates[0];off=[]
        for sign in (-1,1):
            pool=sorted((g for g in rows if g["role"]=="OFF" and g["backend"]==f["backend"] and 0<(g["time"]-f["time"])*sign<=A["off_max_gap_seconds"]),key=lambda g:abs(g["time"]-f["time"]))
            for g in pool[:4]:
                try:
                    h=read_header(g);fh=f["header"]
                    assert covers(h,q["rounded_mhz"]) and all(k in h and k in fh for k in ("src_raj","src_dej"))
                    assert (h["src_raj"],h["src_dej"])!=(fh["src_raj"],fh["src_dej"])
                    assert math.isclose(h["foff"],fh["foff"],rel_tol=1e-10) and math.isclose(h["tsamp"],fh["tsamp"],rel_tol=1e-8)
                    off.append({**g,"header":h});break
                except Exception as e:errors.append({"scan":g["name"],"error":str(e)})
        chosen.append({"target":q,"on":f,"off":off,"status":"READY" if off else "HOLD_NO_OFF"})
    d={"selected":chosen,"index_receipts":indices,"header_errors":errors,"protocol_sha256":sha(P.read_bytes()),"no_intensity_read":True}
    save("manifest.json",d);emit("R33B2_MANIFEST",{"targets":len(chosen),"ready":sum(q["status"]=="READY" for q in chosen),"blocked":[q["target"]["id"] for q in chosen if q["status"]!="READY"],"header_errors":errors})
def search(x,t,rate,tolerance,center,rad):
    y=mem.norm(x.astype(float));n,w=y.shape;df=t["df_hz"];ts=t["ts"];f0=t["channel0_mhz"]*1e6
    hits=[];trials=0
    # All selected paths are scored on measured samples; no synthetic passports.
    for length in sorted(set(min(n,k) for k in G["windows"])):
        if length<G["rows_min"]:continue
        starts=sorted(set(list(range(0,n-length+1,max(1,length//2)))+[n-length,min(n-length,round(900/ts))]))
        for st in starts:
            offsets=set()
            for drift in np.linspace(rate-tolerance,rate+tolerance,G["drift_grid_points"]):
                if drift*rate<=0:continue
                ds=np.rint(np.arange(length)*ts*drift/df).astype(int)
                key=tuple(ds)
                if key in offsets:continue
                offsets.add(key);trials+=1
                lo=max(4,4-int(ds.min()));hi=min(w-4,w-4-int(ds.max()))
                if hi<=lo:continue
                vals=sum(y[st+i,lo+k:hi+k] for i,k in enumerate(ds));med=np.median(vals);sigma=1.4826*np.median(abs(vals-med))
                if sigma<=0:continue
                actual_d=(ds[-1]-ds[0])*df/((length-1)*ts)
                if actual_d*rate<=0:continue
                freqs=f0+np.arange(lo,hi)*df
                # Publication frequency precision is not a unique hit identifier.
                at0=freqs-actual_d*st*ts
                ix=np.flatnonzero(abs(at0-center)<=rad)
                if not len(ix):continue
                z=(vals-med)/sigma
                take=ix[np.argsort(z[ix])[-min(12,len(ix)):]]
                for j in take:
                    if z[j]<G["integrated_robust_score_min"]:continue
                    p=int(lo+j)+ds;occ=float(np.mean(y[st+np.arange(length),p]>3))
                    if occ<G["occupancy_min"]:continue
                    hits.append({"start":st,"rows":length,"path":p.tolist(),"score":float(z[j]),"occupancy":occ,"frequency_hz":float(freqs[j]),"scan_start_frequency_hz":float(at0[j]),"drift_hz_s":float(actual_d),"end_frequency_hz":float(f0+p[-1]*df)})
    hits.sort(key=lambda h:(-h["score"],-h["rows"],h["start"],h["frequency_hz"]))
    distinct=[]
    # Collapse only nearby trajectories at the common scan reference time.
    for h in hits:
        if all(abs(h["scan_start_frequency_hz"]-v["scan_start_frequency_hz"])>G["distinct_track_separation_hz"] for v in distinct):
            distinct.append(h)
        if len(distinct)==8:break
    best=distinct[0] if distinct else None
    unique=bool(best and (len(distinct)==1 or distinct[1]["score"]<best["score"]*G["ambiguity_score_ratio"]))
    return {"best":best,"unique":unique,"alternatives":distinct,"trials":trials}
def acquire(f,c,target,artifact):
    with base.RangeFile(f["url"],A["range_budget_mib_per_file"]*base.MIB) as r:
        h=fil.parse_header(r);assert h==f["header"]
        x,t=fil.fil_window(r,h,c,A["half_hz"]);assert sha(x.tobytes())==t["slice_sha256"]
        nm=target+"-"+sha(f["url"].encode())[:16]+".npy";np.save(OUT/nm,x,allow_pickle=False)
        return x,t,{"scan":f["name"],"url":f["url"],"role":f["role"],"header":h,"tile":t,"raw_artifact":nm,"bytes_read":r.total,"range_requests":r.requests_count,"manifest_artifact":artifact}
def passport(x,t,meta,hit,f,artifact):
    st=hit["start"];nr=hit["rows"];raw=x[st:st+nr];p=np.array(hit["path"],int)
    assert sha(x.tobytes())==t["slice_sha256"]
    m=mem.meta(t,f["name"],artifact,st);m.update({"source_artifact_role":"FROZEN_METADATA_MANIFEST_NOT_RAW_TILE","archive_url":f["url"],"day":58602,"role":f["role"],"row_slice_sha256":sha(raw.tobytes()),"source_telescope":"Parkes/Murriyang"})
    r=mem.make_record(raw,p,m,verify_full=False);mem.validate(r);return r
def scan():
    parent=pathlib.Path(S["parent"]["path"]);assert sha(parent.read_bytes())==S["parent"]["sha256"]
    s=pathlib.Path(S["stationary_sidecar"]["path"]);assert sha(s.read_bytes())==S["stationary_sidecar"]["sha256"]
    m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==sha(P.read_bytes())
    i=os.environ["R33B2_MANIFEST_ID"];d=os.environ["R33B2_MANIFEST_SHA"].removeprefix("sha256:");assert i.isdecimal() and mem.ishash(d);art={"id":int(i),"digest":d}
    rows=[];out=[]
    for q in m["selected"]:
        label=q["target"];c=label["rounded_mhz"];rate=label["sign"]*S["first_panel_rate"]*c/S["first_panel_reference_mhz"]
        row={"published":label,"expected_drift_hz_s":rate,"status":q["status"],"measurements":[],"off_matches":[]}
        try:
            if q["status"]=="READY":
                on=q["on"];x,t,r=acquire(on,c,label["id"],art)
                h=search(x,t,rate,abs(rate)*S["relative_rate_tolerance"],c*1e6,A["start_frequency_within_hz"]);row["measurements"].append({**r,"search":h})
                row["status"]="HOLD_NO_UNIQUE_ON_TRACK"
                if h["unique"]:
                    best=h["best"];pr=passport(x,t,r,best,on,art)
                    row["status"]="HOLD_NO_UNIQUE_OFF_CORROBORATION"
                    for f in q["off"]:
                        try:
                            X,T,R=acquire(f,c,label["id"],art);dt=(f["header"]["tstart"]-on["header"]["tstart"])*86400
                            pred=best["scan_start_frequency_hz"]+dt*best["drift_hz_s"]
                            # Frozen extrapolation, do not re-center on test peak.
                            j=search(X,T,best["drift_hz_s"],abs(best["drift_hz_s"])*G["off_drift_tolerance_fraction"],pred,G["off_prediction_tolerance_hz"])
                            row["measurements"].append({**R,"predicted_frequency_hz":pred,"search":j})
                            if j["unique"] and T["slice_sha256"]!=t["slice_sha256"]:
                                other=passport(X,T,R,j["best"],f,art)
                                row["off_matches"].append({"passport":other,"predicted_frequency_hz":pred,"frequency_residual_hz":j["best"]["scan_start_frequency_hz"]-pred,"measured_track":j["best"]})
                        except Exception as e:row["measurements"].append({"scan":f["name"],"status":"HOLD_READ_OR_TRACK_ERROR","reason":str(e)[:300]})
                    if row["off_matches"]:
                        row["status"]="SOURCE_LINKED_DRIFT_CANDIDATE_REVIEW_REQUIRED"
                        candidate={"published":label,"label":"SOURCE_LINKED_DRIFT_CANDIDATE_REVIEW_REQUIRED","origin_scope":"Rounded published illustration plus unique ON/OFF replay; exact original turboSETI event identity not established","campaign":"Parkes-20190429","passport":pr,"on_track":best,"off_matches":row["off_matches"]}
                        candidate["evidence_sha256"]=mem.digest(candidate);out.append(candidate)
        except Exception as e:row["status"]="HOLD_READ_OR_TRACK_ERROR";row["error"]=str(e)[:300]
        rows.append(row);save("partial.json",rows)
        emit("R33B2_TARGET",{"id":label["id"],"direction":label["sign"],"status":row["status"],"measurements":len(row["measurements"]),"off_matches":len(row["off_matches"]),"on_best":next((v["search"]["best"] and {k:v["search"]["best"][k] for k in ("score","drift_hz_s","frequency_hz","occupancy")} for v in row["measurements"] if "search" in v),None),"error":row.get("error")})
    assert len(rows)==36
    for e in out:
        mem.validate(e["passport"]);assert e["passport"]["label"]=="UNKNOWN"
        assert e["evidence_sha256"]==mem.digest({k:v for k,v in e.items() if k!="evidence_sha256"})
    save("raw_receipts.json",rows);save("drift_candidates.json",{"schema":"janus.hannah.r33b2.drift_candidates.v1","parent_sha256":S["parent"]["sha256"],"entries":out,"authority":"REVIEW_ONLY_NO_VETO"})
    d={"status":"PASS_DRIFT_CANDIDATE_REPLAY" if out else "HOLD_NO_DRIFT_JOIN","published_targets":len(rows),"unique_on":sum(bool(r["measurements"] and r["measurements"][0].get("search",{}).get("unique")) for r in rows),"direct":sum(e["published"]["sign"]==1 for e in out),"mirrored":sum(e["published"]["sign"]==-1 for e in out),"candidate_count":len(out),"off_corroborations":sum(len(e["off_matches"]) for e in out),"holds":{k:sum(r["status"]==k for r in rows) for k in sorted({r["status"] for r in rows})},"candidate_sha256":sha((OUT/"drift_candidates.json").read_bytes()),"manifest_sha256":sha((OUT/"manifest.json").read_bytes()),"verified_promoted":0,"independent_family_generalization":False,"r34":"HOLD_INDEPENDENT_LABELS_AND_PROTECTED_GRID_REQUIRED","production_adoption":False}
    save("result.json",d);emit("R33B2_SUMMARY",d)
def selftest():
    rng=np.random.default_rng(33);x=rng.normal(size=(32,241));t={"channel0_mhz":1000.,"df_hz":4.,"ts":16.}
    for sign in (-1,1):
        z=x.copy();p=120+np.rint(np.arange(32)*16*.04*sign/4).astype(int);z[np.arange(32),p]+=20
        v=search(z,t,.04*sign,.01,1e9+480,20);assert v["best"] and v["best"]["drift_hz_s"]*sign>0
    assert not search(x,t,.04,.01,1e9+480,20)["unique"]
    return {"positive_negative_drift":"PASS","noise_abstention":"PASS"}
if __name__=="__main__":
    {"freeze":freeze,"scan":scan,"selftest":lambda:emit("R33B2_SELFTEST",selftest())}[sys.argv[1]]()
