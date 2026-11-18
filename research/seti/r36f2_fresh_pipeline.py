#!/usr/bin/env python3
"""R36F2: source-disjoint development benchmark, full patch pipeline and exact twin R28."""
import concurrent.futures,hashlib,io,json,os,pathlib,platform,resource,sys,time,zipfile
import numpy as np
import requests
import r35h_integrated as h
import r36_galaxy_swarm as g
import r36f1_hybrid as f
import r35b_grand_benchmark as b
P=pathlib.Path("research/seti/r36f2_fresh_pipeline_protocol.json");S=json.loads(P.read_text());H=S["r36f2"]
OUT=pathlib.Path("runtime/r36f2");OUT.mkdir(parents=True,exist_ok=True)
IN=OUT/"inputs";IN.mkdir(exist_ok=True)
h.S=S;h.P=P;h.H=S["r35h"];h.OUT=OUT;h.mem.OUT=OUT;h.prev.S=S;h.prev.OUT=OUT
b.P=P;b.S=S;b.OUT=OUT;b.sh.P=P;b.sh.S=S;b.sh.C=S["shield"];b.sh.OUT=OUT;b.h.S=S;b.rb.S=S;b.ref.S=S;b.od.S=S;b.ref.OUT=OUT;b.sh.eng.OUT=OUT
b.sh.hold.P=P;b.sh.hold.OUT=OUT
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,v):(OUT/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def emit(k,v):print(k,json.dumps(v,allow_nan=False),flush=True)
def timed(fn):
    w=time.perf_counter();c=time.process_time();x=fn()
    return x,{"wall":time.perf_counter()-w,"cpu":time.process_time()-c}
def add(parts):
    return {k:sum(p[k] for p in parts) for k in ("wall","cpu")}
def exclusion():
    p=H["prior_manifest"];u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['id']}/zip"
    resp=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=60);resp.raise_for_status()
    raw=resp.content;assert len(raw)<16*1024**2 and sha(raw)==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        ns=[n for n in z.namelist() if pathlib.PurePosixPath(n).name=="manifest.json"];assert len(ns)==1
        m=json.loads(z.read(ns[0]))
    ss=m["selected"];assert len(ss)==8 and len({q["day"] for q in ss})==8
    return {"artifact":p,"selected":[{"target":q["target"],"day":q["day"],"files":[f["name"] for f in q["files"]]} for q in ss],"targets":sorted(set(S["exclude_extra"])|{q["target"] for q in ss}),"days":sorted(set(S["shield"]["exclude_days"])|{q["day"] for q in ss})}
def freeze():
    x=exclusion();save("exclusion.json",x)
    S["shield"]["exclude_days"]=x["days"]
    b.sh.hold.S={**b.sh.hold.S,"seed":S["seed"],"exclude":sorted(set(b.sh.hold.S["exclude"]+x["targets"]))}
    b.sh.freeze()
    m=json.loads((OUT/"manifest.json").read_text());assert not {q["target"] for q in m["selected"]}&set(x["targets"])
    assert not {q["day"] for q in m["selected"]}&set(x["days"])
    assert not {f["name"] for q in m["selected"] for f in q["files"]}&{f for q in x["selected"] for f in q["files"]}
    emit("R36F2_EXCLUSIONS",{"targets":len(x["targets"]),"days":len(x["days"]),"prior_r35h":x["selected"],"sha256":sha((OUT/"exclusion.json").read_bytes())})
def prep(rows,R):
    raw=np.stack([v[0] for v in rows]);X=b.ref.norm(raw);Y,_,_,_=b.ref.register(X,b.ref.norm(R))
    return X.astype(np.float32),Y.astype(np.float32)
def pathscore(V):
    return g.scan(np.ascontiguousarray(V,dtype=np.float64),g.D,g.R,g.L,False)
def od(V):return b.od.sliced(V,0,None)[:2]
def original(rows,R):
    (X,Y),prep_cost=timed(lambda:prep(rows,R));V=Y.astype(np.float64)
    scores,cost=timed(lambda:(b.taylor.straight(V),od(V)))
    return scores,X,Y,prep_cost,cost
def twin(rows,R,pool,raw_needed=True):
    (X,Y),prep_cost=timed(lambda:prep(rows,R));V=Y.astype(np.float64)
    def work():
        a=pool.submit(pathscore,V);z=pool.submit(od,V)
        lin=f.get(a.result(),0);o=z.result()
        r=f.get(g.scan(np.ascontiguousarray(X,dtype=np.float64),g.D,g.R,g.T,False),1) if raw_needed else None
        return lin,o,r
    scores,cost=timed(work)
    return scores,X,Y,prep_cost,cost
def eq(a,b):
    assert all(np.array_equal(x,y) for x,y in zip(a,b)),"score/path mismatch"
def raw_equal(X,r):
    eq(b.template(X.astype(np.float64),"raw"),r)
def cut(s,rate=.05):
    z=np.sort(s)[::-1];assert len(z)>=64
    return float(np.nextafter(z[int(len(z)*rate)],np.inf))
def margin_fit(fit):
    y=np.concatenate([d["y"] for d in fit]);L={}
    for i,k in enumerate(("linear","odonto")):
        z=np.concatenate([d["o"][i][0] for d in fit])[y==0];th=cut(z);L[k+"_threshold"]=th;L[k+"_scale"]=max(.1,th-float(np.quantile(z,.8)))
    return L
def labels(rows):return np.asarray([r[1] for r in rows])
def make_part(rows,R,pool):
    a,X,Y,ap,ac=original(rows,R);o,x,v,gp,gc=twin(rows,R,pool)
    assert np.array_equal(X,x) and np.array_equal(Y,v)
    eq(a[0],o[0]);eq(a[1],o[1]);raw_equal(X,o[2])
    return {"rows":rows,"R":R,"y":labels(rows),"o":a,"g":o,"Y":Y,"original_prep":ap,"original_search":ac,"galaxy_prep":gp,"galaxy_search":gc}
def read_group(q,guard,pool):
    got,cost=timed(lambda:h.prev.gather(q,guard));rows,R,ch,ids=got;assert len(ids)==16
    d=make_part(rows,R,pool);shared=make_part([v[0] for v in ch],np.stack([v[1] for v in ch]),pool)
    return d,shared,cost
def evidence(d,z,p,alarm,known,unknown):
    return [h.evidence(row,R,path,known,unknown) if a else {"shape":0.,"reference":None,"review_shape":0.,"review_reference":None,"off":-1e6,"compatible":0,"valid":False} for row,R,path,a in zip(d["rows"],d["R"],p,alarm)]
def extract(d,L,th,known,unknown):
    (z,p),ot=timed(lambda:f.fuse(*d["o"],L))
    (zz,pp),gt=timed(lambda:f.fuse(*d["g"][:2],L));eq((z,p),(zz,pp))
    r,rp=d["g"][2];alarm=z>=th
    ev,et=timed(lambda:evidence(d,z,p,alarm,known,unknown))
    rows=d["rows"];truth=np.stack([v[7] for v in rows]);act=np.stack([v[8] for v in rows]);amp=np.asarray([v[4] for v in rows])
    loc=lambda paths:((abs(paths-truth)<=2)&act).sum(1)>=.75*act.sum(1)
    x={"y":d["y"],"weak":amp<=1.25,"very_weak":amp==.75,"alarm":alarm,"loc":loc(p),"raw_loc":loc(rp),"raw":r,"p":p,"rp":rp,"z":z,"ev":ev}
    return x,{"original_prep":d["original_prep"],"original_search_fuse":add([d["original_search"],ot]),"galaxy_prep":d["galaxy_prep"],"galaxy_search_fuse":add([d["galaxy_search"],gt]),"hannah":et}
def export(d,x,label,name,policy,costs):
    rows=d["rows"];(keep,loc,hold,rescue),policy_apply=timed(lambda:h.apply(x,policy));base=x["alarm"]&x["loc"]&(x["y"]==1);joint=keep&loc&(x["y"]==1)
    a={"V":d["Y"],"y":x["y"],"truth":np.stack([r[7] for r in rows]),"active":np.stack([r[8] for r in rows]),"amp":np.asarray([r[4] for r in rows]),"model":np.asarray([r[2] for r in rows],dtype="U32"),"tile_sha":np.asarray([r[6] for r in rows],dtype="U64"),"r28_alert":x["alarm"],"r28_localized":base,"galaxy_alert":keep,"galaxy_localized":joint,"galaxy_hold":hold,"galaxy_rescue":rescue,"r28_score":x["z"],"r28_path":x["p"],"raw_score":x["raw"],"raw_path":x["rp"]}
    path=IN/(label+"_"+name+".npz")
    _,dt=timed(lambda:np.savez_compressed(path,**a))
    # Full decision/output records include advisory evidence; not saved in operational memory.
    ev=OUT/(label+"_"+name+"_evidence.json")
    def write():
        ev.write_text(json.dumps({"evidence":x["ev"],"held":hold.tolist(),"rescued":rescue.tolist()},allow_nan=False)+"\n")
    _,et=timed(write)
    return {"name":path.name,"sha256":sha(path.read_bytes()),"array_sha256":sha(a["V"].tobytes()),"examples":len(rows),"timing":{**costs,"common_export":dt,"galaxy_evidence_write":et,"galaxy_policy_apply":policy_apply},"baseline":h.metric(x,{}),"galaxy":h.metric(x,policy)}
def main():
    m=json.loads((OUT/"manifest.json").read_text());ex=json.loads((OUT/"exclusion.json").read_text())
    assert m["protocol_sha256"]==sha(P.read_bytes());ss=m["selected"];assert len(ss)==8 and len({q["day"] for q in ss})==8
    assert [q["split"] for q in ss]==["component_fit"]*2+["null_calibration"]*2+["evaluation"]*4
    assert not {q["target"] for q in ss}&set(ex["targets"]) and not {q["day"] for q in ss}&set(ex["days"])
    assert not {f["name"] for q in ss for f in q["files"]}&{f for q in ex["selected"] for f in q["files"]}
    (known,unknown,old),memcost=timed(h.loadmemory);guard=h.Guard(old)
    x=np.random.default_rng(36062).normal(size=(4,16,128)).astype(np.float32).astype(np.float64)
    _,compile=timed(lambda:pathscore(x));eq(b.taylor.straight(x),f.get(pathscore(x),0))
    raw_equal(x,f.get(g.scan(x,g.D,g.R,g.T,False),1))
    records=[];acq=[];evals=[];shared=[];fit_phase=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        fit=[]
        for q in ss[:2]:
            d,l,t=read_group(q,guard,pool);fit.append(d);acq.append({"target":q["target"],"split":q["split"],"cost":t})
        L,margin_cost=timed(lambda:margin_fit(fit))
        for d in fit:fit_phase.append({"original":add([d["original_prep"],d["original_search"]]),"galaxy":add([d["galaxy_prep"],d["galaxy_search"]])})
        del fit
        pairs=[]
        for q in ss[2:4]:
            d,l,t=read_group(q,guard,pool);pairs.append((d,l));acq.append({"target":q["target"],"split":q["split"],"cost":t})
        def thresholds():
            y=np.concatenate([d["y"] for d,l in pairs]);z=np.concatenate([f.fuse(*d["o"],L)[0] for d,l in pairs])
            r=np.concatenate([d["g"][2][0] for d,l in pairs]);return cut(z[y==0]),[cut(r[y==0],v) for v in h.H["policy"]["raw_rates"]]
        (th,rawcuts),threshold_cost=timed(thresholds)
        cal=[extract(d,L,th,known,unknown) for d,l in pairs];cl=[extract(l,L,th,known,unknown) for d,l in pairs]
        (policy,trials),policy_cost=timed(lambda:h.fit_policy([x for x,t in cal],[x for x,t in cl],rawcuts))
        assert h.sum_metrics([x for x,t in cal],policy)["control_alerts"]<=12
        lock={"margins":L,"threshold":th,"rawcuts":rawcuts,"policy":policy,"trials":trials,"manifest_sha256":sha((OUT/"manifest.json").read_bytes()),"exclusion_sha256":sha((OUT/"exclusion.json").read_bytes()),"memory":h.H["memory"]}
        save("lock.json",lock);lh=sha((OUT/"lock.json").read_bytes());emit("R36F2_LOCK",{"sha256":lh,"policy":policy,"trials":len(trials),"admissible":sum(t["admissible"] for t in trials),"baseline_controls":h.sum_metrics([x for x,t in cal],{})["control_alerts"],"joint_controls":h.sum_metrics([x for x,t in cal],policy)["control_alerts"]})
        for q,(d,l),(x,t),(v,u) in zip(ss[2:4],pairs,cal,cl):
            records.append(export(d,x,"cal",q["target"],policy,t))
            # Calibration shared challenge used for safety only, not native null thresholds.
            fit_phase.append({"original":add([t["original_prep"],t["original_search_fuse"],u["original_prep"],u["original_search_fuse"]]),"galaxy":add([t["galaxy_prep"],t["galaxy_search_fuse"],t["hannah"],u["galaxy_prep"],u["galaxy_search_fuse"],u["hannah"]])})
        del pairs,cal,cl
        for q in ss[4:]:
            assert sha((OUT/"lock.json").read_bytes())==lh
            d,l,t=read_group(q,guard,pool);acq.append({"target":q["target"],"split":q["split"],"cost":t})
            x,u=extract(d,L,th,known,unknown);v,w=extract(l,L,th,known,unknown)
            evals.append(x);shared.append(v)
            records.append(export(d,x,"eval",q["target"],policy,u));records.append(export(l,v,"shared",q["target"],policy,w))
            save("partial.json",records);emit("R36F2_TARGET",{"target":q["target"],"eval":records[-2],"shared":records[-1]})
    assert len(guard.seen)==128 and len(records)==10 and sha((OUT/"lock.json").read_bytes())==lh
    E={k:h.sum_metrics(evals,p) for k,p in (("baseline",{}),("galaxy",policy))};J={k:h.sum_metrics(shared,p) for k,p in (("baseline",{}),("galaxy",policy))}
    assert all(v["positives"]==v["controls"]==512 for v in E.values()) and all(v["positives"]==128 for v in J.values())
    a=E["baseline"];z=E["galaxy"];j=J["galaxy"]
    gate={"overall":z["localized"]>a["localized"],"weak":z["weak"]>a["weak"],"fewer_controls":z["control_alerts"]<a["control_alerts"],"no_paired_losses":z["lost"]==z["held_positives"]==j["lost"]==j["held_positives"]==0}
    t={k:add([r["timing"][k] for r in records if not r["name"].startswith("cal_")]) for k in records[0]["timing"]}
    result={"status":"PASS_EXACT_COMMON_INPUT_REPLAY","science_gate":"PASS_DEVELOPMENT" if all(gate.values()) else "HOLD","gate_checks":gate,"eval":E,"shared":J,"policy":policy,"records":records,"acquisition_and_patch_cost":acq,"evaluation_cost":t,"fit_and_null_cost":{k:add([r[k] for r in fit_phase]) for k in ("original","galaxy")},"margin_cost":margin_cost,"threshold_cost":threshold_cost,"policy_cost":policy_cost,"memory_load_cost":memcost,"jit_compile":compile,"lock_sha256":lh,"manifest_sha256":lock["manifest_sha256"],"exclusion_sha256":lock["exclusion_sha256"],"groups":[{"target":q["target"],"day":q["day"],"split":q["split"]} for q in ss],"source_tile_hashes":sorted(guard.seen),"verified_reference_occurrences":len(known),"review_only":len(unknown),"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"host":platform.node(),"cpu_count":os.cpu_count(),"native_scope":"Same registered float32 arrays; shared acquisition and normalization cost separately measured. Synthetic injections only, no origin labels; actual compatible-grid useful source still HOLD.","production_adoption":False}
    save("receipt.json",result);emit("R36F2_SUMMARY",{k:v for k,v in result.items() if k not in ("records","source_tile_hashes","acquisition_and_patch_cost")})
if __name__=="__main__":
    if sys.argv[1]=="freeze":freeze()
    else:main()
