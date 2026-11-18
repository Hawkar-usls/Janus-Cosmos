#!/usr/bin/env python3
"""R30C: provenance-bound radio fingerprints and append-only research memory."""
import hashlib,json,math,os,pathlib,sys,time
import numpy as np
P=pathlib.Path("research/seti/r30c_hannah_protocol.json")
S=json.loads(P.read_text());OUT=pathlib.Path("runtime/r30c");OUT.mkdir(parents=True,exist_ok=True)
SCHEMA=S["ledger"]["schema"]; VERSION=S["feature_contract"]["version"]
def canon(v):return json.dumps(v,sort_keys=True,separators=(",",":"),ensure_ascii=False,allow_nan=False).encode()
def digest(v):return hashlib.sha256(canon(v)).hexdigest()
def filehash(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,v):
    p=OUT/name;t=p.with_suffix(p.suffix+".tmp");t.write_text(json.dumps(v,indent=2,allow_nan=False)+"\n");t.replace(p)
def emit(k,v):print(k,json.dumps(v,allow_nan=False),flush=True)
def ishash(s):return isinstance(s,str) and len(s)==64 and all(x in "0123456789abcdef" for x in s)
def norm(X):
    m=np.median(X,axis=1,keepdims=True);s=1.4826*np.median(abs(X-m),axis=1,keepdims=True)
    if not np.isfinite(X).all() or not np.all(s>0):raise ValueError("Invalid measured power or zero robust noise")
    return (X-m)/s
def features(raw,path,df,ts):
    x=np.asarray(raw);p=np.asarray(path)
    if x.ndim!=2 or len(p)!=len(x) or len(p)<16 or not np.issubdtype(p.dtype,np.integer):raise ValueError("Need >=16 actual rows and integer track")
    if p.min()<4 or p.max()>=x.shape[1]-4 or df==0 or ts<=0 or not np.isfinite([df,ts]).all():raise ValueError("Incomplete measured support")
    y=norm(x.astype(float));i=np.arange(len(p))[:,None];v=y[i,p[:,None]+np.arange(-4,5)[None,:]]
    # Orientation follows physical frequency, not file bin direction.
    if df<0:v=v[:,::-1]
    v=np.clip(v,-8,32);core=v[:,4];profile=np.median(v,axis=0)
    rel=profile-profile[[0,8]].mean();spectral=np.tanh(rel/8)
    series=np.interp(np.linspace(0,len(core)-1,16),np.arange(len(core)),core)
    dt=np.diff(p)*df/ts
    vals=np.r_[spectral,np.tanh(series/8),np.tanh(np.mean(core)/8),np.mean(core>3),
               np.tanh(np.std(core)/8),np.tanh(np.median(abs(dt))/.1),
               np.tanh(np.mean(abs(np.diff(dt)))/.1)]
    assert len(vals)==30 and np.isfinite(vals).all()
    return {"version":VERSION,"vector":[float(round(v,8)) for v in vals],
            "profile_9_bins":[float(round(v,8)) for v in profile],
            "median_track_robust_sigma":float(np.median(core)),
            "occupancy_gt_3":float(np.mean(core>3)),
            "drift_hz_s":float((int(p[-1])-int(p[0]))*df/((len(p)-1)*ts)),
            "channel_width_hz":float(abs(df)),"sample_interval_s":float(ts),"rows":len(p)}
def make_record(raw,path,meta,*,label="UNKNOWN",evidence=None,verify_full=True):
    fields=("source_tile_sha256","source_scan","source_artifact","artifact_sha256","channel0_mhz","df_hz","ts","start_row","kind")
    if raw is None or any(k not in meta for k in fields):raise ValueError("BLOCKED_MISSING_MEASUREMENTS_OR_PROVENANCE")
    if meta["kind"]!="MEASURED_ARCHIVE":raise ValueError("Synthetic samples cannot enter measured memory")
    if not ishash(meta["source_tile_sha256"]) or not ishash(meta["artifact_sha256"]):raise ValueError("Invalid source hashes")
    if not isinstance(meta["source_scan"],str) or not meta["source_scan"] or not isinstance(meta["source_artifact"],int):raise ValueError("Invalid source identity")
    if verify_full and hashlib.sha256(np.ascontiguousarray(raw).tobytes()).hexdigest()!=meta["source_tile_sha256"]:raise ValueError("Raw hash mismatch")
    p=np.asarray(path);f=features(raw,p,meta["df_hz"],meta["ts"])
    origin=float(meta["channel0_mhz"]*1e6+p[0]*meta["df_hz"])
    key={"version":VERSION,"source_scan":meta["source_scan"],"source_tile_sha256":meta["source_tile_sha256"],"start_row":int(meta["start_row"]),"track":p.astype(int).tolist()}
    rec={"id":digest(key),"identity":key,"fingerprint_sha256":digest(f),"features":f,"frequency_hz":origin,
         "provenance":dict(meta),"label":label,"evidence":evidence,"authority":"RETRIEVAL_ONLY"}
    if label=="VERIFIED_RFI":check_verified(rec)
    elif label!="UNKNOWN" or evidence is not None:raise ValueError("Unreviewed labels forbidden")
    rec["record_sha256"]=digest(rec);return rec
def check_verified(rec):
    v=S["verified_seed"];m=rec["provenance"];e=rec["evidence"] or {}
    if not (m["source_scan"]==v["scan"] and m["source_artifact"]==S["sources"]["blc_artifact"]["id"]
            and m["artifact_sha256"]==S["sources"]["blc_artifact"]["digest"] and e.get("id")==v["id"]
            and e.get("source")==v["source"] and e.get("receipt")==v["receipt"]
            and abs(rec["frequency_hz"]-v["frequency_mhz"]*1e6)<=10
            and abs(e.get("measured_search_drift_hz_s",1e6)-v["drift_hz_s"])<=.012
            and e.get("measured_search_score",0)>=10):
        raise ValueError("Verified RFI admission rejected")
def validate(rec):
    if rec["id"]!=digest(rec["identity"]) or rec["fingerprint_sha256"]!=digest(rec["features"]):raise ValueError("Fingerprint tampering")
    if rec["record_sha256"]!=digest({k:v for k,v in rec.items() if k!="record_sha256"}):raise ValueError("Record tampering")
    m=rec["provenance"];key=rec["identity"]
    assert m["kind"]=="MEASURED_ARCHIVE" and ishash(m["source_tile_sha256"]) and ishash(m["artifact_sha256"])
    assert key["source_tile_sha256"]==m["source_tile_sha256"] and key["source_scan"]==m["source_scan"] and key["start_row"]==m["start_row"]
    assert key["version"]==VERSION==rec["features"]["version"] and len(rec["features"]["vector"])==30 and rec["authority"]=="RETRIEVAL_ONLY"
    if rec["label"]=="VERIFIED_RFI":check_verified(rec)
    elif rec["label"]!="UNKNOWN" or rec["evidence"] is not None:raise ValueError("Unreviewed labels forbidden")
def empty():return {"schema":SCHEMA,"entries":[]}
def append(ledger,records):
    assert ledger["schema"]==SCHEMA
    by={}
    for rec in ledger["entries"]:
        validate(rec)
        if rec["id"] in by:raise ValueError("Duplicate parent IDs")
        by[rec["id"]]=rec
    added=duplicate=0
    for rec in records:
        validate(rec)
        if rec["id"] in by:
            if by[rec["id"]]!=rec:raise ValueError("Conflicting immutable record")
            duplicate+=1;continue
        by[rec["id"]]=rec;added+=1
    return {"schema":SCHEMA,"entries":list(by.values())},{"added":added,"duplicates":duplicate,"total":len(by)}
def compare(a,b):
    u=a["features"];v=b["features"];x=np.asarray(u["vector"]);y=np.asarray(v["vector"])
    compatible=math.isclose(u["channel_width_hz"],v["channel_width_hz"],rel_tol=.02) and math.isclose(u["sample_interval_s"],v["sample_interval_s"],rel_tol=.02)
    return {"shape_similarity":float(np.exp(-np.sqrt(np.mean((x-y)**2)))),"grid_compatible":compatible,
       "frequency_delta_hz":abs(a["frequency_hz"]-b["frequency_hz"]),"drift_delta_hz_s":abs(u["drift_hz_s"]-v["drift_hz_s"]),
       "distinct_scan":a["provenance"]["source_scan"]!=b["provenance"]["source_scan"],"distinct_tile":a["provenance"]["source_tile_sha256"]!=b["provenance"]["source_tile_sha256"],"authority":"RETRIEVAL_ONLY"}
def review(ledger):
    entries=ledger["entries"];known=[r for r in entries if r["label"]=="VERIFIED_RFI"]
    unknown=[r for r in entries if r["label"]=="UNKNOWN"];recurrence=set();pairs=[]
    # Different scans and compatible grids; never count same tile twice or grant a physical label.
    for i,a in enumerate(unknown):
        for b in unknown[i+1:]:
            d=compare(a,b)
            if d["distinct_scan"] and d["distinct_tile"] and d["grid_compatible"] and d["frequency_delta_hz"]<=8.4 and d["drift_delta_hz_s"]<=.011:
                recurrence.update([a["id"],b["id"]]);pairs.append({"a":a["id"],"b":b["id"],**d})
    matches=[]
    for a in unknown:
        q=sorted([{"reference":b["id"],**compare(a,b)} for b in known if a["provenance"]["source_scan"]!=b["provenance"]["source_scan"]],key=lambda r:r["shape_similarity"],reverse=True)
        if q:matches.append({"query":a["id"],"nearest":q[:3]})
    return {"verified_rfi":len(known),"rfi_candidates":len(recurrence),"unknown_without_recurrence":len(unknown)-len(recurrence),"recurrence_pairs":pairs,"derived_candidate_ids":sorted(recurrence),"nearest_verified":matches,
            "notice":"Candidate tier is derived review evidence only. No automatic promotion, veto or origin decision."}
def search_tile(raw):
    y=norm(raw.astype(float));n,w=y.shape;assert n==16;best=[]
    for end in range(-32,33):
        d=np.rint(np.arange(16)*end/15).astype(int);lo=4-int(d.min());hi=w-4-int(d.max())
        if hi<=lo:continue
        z=sum(y[i,lo+k:hi+k] for i,k in enumerate(d))/4
        count=min(2,len(z))
        for j in np.argpartition(z,-count)[-count:]:
            best.append((float(z[j]),lo+int(j),end))
    keep=[]
    for score,j,end in sorted(best,reverse=True):
        if score<S["observations"]["score_min"]:break
        if any(abs(j-v[1])<8 for v in keep):continue
        keep.append((score,j,end))
        if len(keep)==2:break
    return [(s,j+np.rint(np.arange(16)*end/15).astype(int)) for s,j,end in keep]
def meta(t,scan,artifact,start=0):
    return {"source_tile_sha256":t["slice_sha256"],"source_scan":scan,"source_artifact":artifact["id"],"artifact_sha256":artifact["digest"],
      "channel0_mhz":t["channel0_mhz"],"df_hz":t["df_hz"],"ts":t["ts"],"start_row":start,"kind":"MEASURED_ARCHIVE"}
def sources():
    # Import existing verified artifact reader only after focused memory selftests.
    import r26c_reference_alignment as r
    r.BLC_OUT=OUT
    cfg=S["sources"];fp=pathlib.Path(cfg["family_protocol"])
    assert hashlib.sha1(("blob "+str(len(fp.read_bytes()))+"\0").encode()+fp.read_bytes()).hexdigest()==cfg["family_blob"]
    family=json.loads(fp.read_text());assert family["source"]==cfg["family_artifact"] and family["manifest"]==cfg["family_manifest"]
    r.BLC_P["artifacts"].update({"hannah_family":cfg["family_artifact"],"hannah_manifest":cfg["family_manifest"]})
    z,n,zp=r.archive("hannah_family");m,mn,mp=r.archive("hannah_manifest")
    manifest=r.get(m,mn,"manifest.json");assert z.read(n["manifest.json"])==m.read(mn["manifest.json"])
    by={f["name"]:f for f in family["files"]};assert len(by)==24
    records=[];observations=[];blocked=[]
    for q in manifest["selected"]:
        for f in q["files"]:
            scan=f["name"];ff=by[scan]
            for t in ff["tiles"]:
                a=r.tile(z,n,t) # exact hash of actual measured source tile
                hits=search_tile(a);source=meta(t,scan,cfg["family_artifact"])
                ids=[]
                for sc,p in hits:
                    rec=make_record(a,p,source);records.append(rec);ids.append(rec["id"])
                observations.append({"scan":scan,"role":f["role"],"target":q["target"],"tile_sha256":t["slice_sha256"],"candidate_ids":ids,"searched":True})
    z.close();m.close();zp.unlink();mp.unlink()
    assert len(observations)==96
    assert r.BLC_P["artifacts"]["blc"]==cfg["blc_artifact"] and r.BLC_P["artifacts"]["blc_manifest"]==cfg["blc_manifest"]
    z,n,zp=r.archive("blc");m,mn,mp=r.archive("blc_manifest")
    v=S["verified_seed"];d=r.get(z,n,"result.json");rows=[x for x in d["rows"] if x["name"]==v["scan"] and "tile" in x];assert len(rows)==1
    row=rows[0];t=row["tile"];a=r.tile(z,n,t,row["artifact"]);from r9_expanded_search import normalized
    best=r.blc_search(normalized(a).astype(float),t)
    match=[x for x in best if x is not None and x["score"]>=10 and abs(x["frequency_mhz"]-v["frequency_mhz"])*1e6<=10 and abs(x["drift_hz_s"]-v["drift_hz_s"])<=.012]
    if match:
        hit=match[0];start=hit["start"];nr=hit["rows"];tt=np.arange(nr)*t["ts"]
        shifts=np.rint((hit["drift_hz_s"]*tt+.5*hit["drift_change_hz_s"]*tt**2/max((a.shape[0]-1)*t["ts"],1))/t["df_hz"]).astype(int)
        q0=int(round((hit["frequency_mhz"]-t["channel0_mhz"])*1e6/t["df_hz"]));p=q0+shifts
        evidence={"id":v["id"],"source":v["source"],"receipt":v["receipt"],"measured_search_score":hit["score"],"measured_search_drift_hz_s":hit["drift_hz_s"],"measured_search":hit}
        # r.tile verified the full tile before slicing; record the exact row slice.
        record=make_record(a[start:start+nr],p,meta(t,v["scan"],cfg["blc_artifact"],start),label="VERIFIED_RFI",evidence=evidence,verify_full=False)
        records.append(record);observations.append({"scan":v["scan"],"tile_sha256":t["slice_sha256"],"candidate_ids":[record["id"]],"searched":True,"known_anchor":True})
    else:blocked.append({"source":v["scan"],"reason":"KNOWN_ANCHOR_NOT_REPRODUCED_NO_ADMISSION"})
    z.close();m.close();zp.unlink();mp.unlink()
    return records,observations,blocked
def load_parent():
    aid=os.environ.get("HANNAH_PARENT_ARTIFACT","").strip();sha=os.environ.get("HANNAH_PARENT_SHA256","").strip()
    if not aid and not sha:return empty(),None
    if not aid.isdecimal() or not ishash(sha):raise ValueError("Require explicit numeric parent ID and pinned digest")
    import r26c_reference_alignment as r
    r.BLC_OUT=OUT;r.BLC_P["artifacts"]["hannah_parent"]={"id":int(aid),"digest":sha}
    z,n,p=r.archive("hannah_parent")
    ledger=r.get(z,n,"library.json");z.close();p.unlink()
    validated,_=append(ledger,[])
    return validated,{"artifact_id":int(aid),"sha256":sha}
def selftest():
    rng=np.random.default_rng(30);x=rng.normal(size=(16,128));p=np.full(16,64,int);x[:,64]+=12
    # Synthetic fixture exercises code only, never persisted to measured archive ledger.
    f=features(x,p,2.8,18.);g=features(x.copy(),p,2.8,18.)
    assert f==g and len(f["vector"])==30
    m={"source_tile_sha256":hashlib.sha256(x.tobytes()).hexdigest(),"source_scan":"TEST_FIXTURE_ONLY","source_artifact":0,"artifact_sha256":"0"*64,"channel0_mhz":1000.,"df_hz":2.8,"ts":18.,"start_row":0,"kind":"SYNTHETIC_TEST"}
    def blocked(fn):
        try:fn()
        except (ValueError,AssertionError):return
        raise AssertionError("Expected blocked")
    blocked(lambda:make_record(x,p,m));blocked(lambda:make_record(None,p,m))
    # Mock measured metadata in ephemeral selftest to test ledger contract only; no exported records.
    mm={**m,"kind":"MEASURED_ARCHIVE"}
    a=make_record(x,p,mm);b=make_record(x,p,{**mm,"source_scan":"TEST_FIXTURE_2"})
    assert a["fingerprint_sha256"]==b["fingerprint_sha256"] and a["id"]!=b["id"]
    assert compare(a,a)["shape_similarity"]==1 and not compare(a,a)["distinct_scan"]
    ledger,c=append(empty(),[a,a]);assert c=={"added":1,"duplicates":1,"total":1}
    ledger,c=append(ledger,[b]);assert c["total"]==2 and review(ledger)["rfi_candidates"]==0 and review(ledger)["verified_rfi"]==0
    x2=x.copy();x2[0,0]+=0.01;mm2={**mm,"source_scan":"TEST_FIXTURE_3","source_tile_sha256":hashlib.sha256(x2.tobytes()).hexdigest()};cc=make_record(x2,p,mm2)
    ledger,c=append(ledger,[cc]);assert c["total"]==3 and review(ledger)["rfi_candidates"]==3 and review(ledger)["verified_rfi"]==0
    bad={**a,"frequency_hz":0};blocked(lambda:append(ledger,[bad]))
    blocked(lambda:make_record(x,p,mm,label="VERIFIED_RFI",evidence={}))
    blocked(lambda:make_record(x,p,{**mm,"source_tile_sha256":"1"*64}))
    y=rng.normal(size=(16,128));y[:,60]+=12;u=features(y,np.full(16,60),2.8,18.)
    v=features(y,np.arange(16)+30,2.8,18.)
    sim=lambda k:float(np.exp(-np.sqrt(np.mean((np.asarray(f["vector"])-np.asarray(k["vector"]))**2))))
    assert sim(u)>sim(v)
    return {k:"PASS" for k in ("determinism","label_firewall","source_identity","duplicate_replay","tamper_detection","verified_admission","missing_raw","synthetic_rejected","distinct_scan","lookalike_retrieval")}
def main():
    checks=selftest();emit("R30C_SELFTEST",checks);t=time.perf_counter()
    ledger,parent=load_parent();records,observations,blocked=sources()
    result,delta=append(ledger,records);repeat,again=append(result,records)
    assert result==repeat and again["added"]==0 and again["duplicates"]==len(records)
    summary=review(result);save("library.json",result);save("observations.json",observations);save("review.json",summary)
    d={"status":"DEVELOPMENT_ARCHIVE_MEMORY","selftests":checks,"parent":parent,"delta":delta,"replay":again,"observations":len(observations),"blocked":blocked,
      "verified_rfi":summary["verified_rfi"],"rfi_candidates":summary["rfi_candidates"],"unknown_without_recurrence":summary["unknown_without_recurrence"],
      "recurrence_pairs":len(summary["recurrence_pairs"]),"library_sha256":filehash(OUT/"library.json"),"protocol_sha256":filehash(P),
      "seconds":time.perf_counter()-t,"gate":"PASS_MEMORY_DEVELOPMENT" if summary["verified_rfi"]>=1 and len(records)>1 and not blocked else "HOLD",
      "production_adoption":False,"scope":"Actual archive catalogue + published BLC1 anchor; synthetic fixture tests only. No fresh classification exam. No automatic RFI labeling or veto. Parent artifact chain is explicit and pinned, not autonomous repository promotion."}
    save("result.json",d);emit("R30C_SUMMARY",d)
if __name__=="__main__":
    if len(sys.argv)>1 and sys.argv[1]=="selftest":emit("R30C_SELFTEST",selftest())
    else:main()
