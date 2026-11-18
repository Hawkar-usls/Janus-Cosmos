#!/usr/bin/env python3
"""R32: adaptive family-retrieval ablations with calibrated no-loss advisory HOLD."""
import json,pathlib,sys,time,math,itertools
import numpy as np
import r31_hannah_paired as prev
import r30_false_alarm_shield as shield
import r30c_hannah_memory as mem
import r32_g3_recover as recover
import hannah_observation_hook as hook
import r26c_reference_alignment as ref
P=pathlib.Path("research/seti/r32_hannah_adaptive_protocol.json");S=json.loads(P.read_text());H=S["hannah"]
OUT=pathlib.Path("runtime/r32");OUT.mkdir(parents=True,exist_ok=True)
shield.P=P;shield.S=S;shield.C=S["shield"];shield.OUT=OUT
shield.h.S=S;shield.rb.S=S;ref.S=S;shield.od.S=S;ref.OUT=OUT;shield.eng.OUT=OUT
shield.hold.S={**shield.hold.S,"seed":S["seed"],"exclude":sorted(set(shield.hold.S["exclude"]+S["exclude_extra"]))}
shield.hold.P=P;shield.hold.OUT=OUT;mem.OUT=OUT;prev.S=S;prev.OUT=OUT
def save(n,v):mem.save(n,v)
def emit(k,v):mem.emit(k,v)
def sha(p):return prev.sha(p)
def load():
    path=pathlib.Path(H["parent_path"])
    if path.exists():
        assert sha(path)==H["parent_sha256"]
        lib,_=mem.append(json.loads(path.read_text()),[])
    else:lib,raw=recover.restore()
    assert len(lib["entries"])==H["parent_entries"]
    return lib
def shapevec(f):
    prof=np.asarray(f["profile_9_bins"],float);prof=prof-(prof[0]+prof[-1])/2
    prof=prof/max(1e-8,float(np.max(abs(prof))))
    t=np.asarray(f["vector"][9:25],float);t=t/max(1e-8,float(np.max(abs(t))))
    occ=float(f["occupancy_gt_3"])
    return prof,t,occ,float(f["vector"][27])
def family(u,v):
    if not (math.isclose(u["channel_width_hz"],v["channel_width_hz"],rel_tol=.02) and math.isclose(u["sample_interval_s"],v["sample_interval_s"],rel_tol=.02)):return 0.
    a,t,o,st=shapevec(u);b,w,p,su=shapevec(v)
    pd=min(float(np.mean((a-b)**2)),float(np.mean((a-b[::-1])**2)))
    td=float(np.mean((t-w)**2))
    od=(o-p)**2;sd=(st-su)**2
    drift=min(abs(u["drift_hz_s"]-v["drift_hz_s"]),abs(u["drift_hz_s"]+v["drift_hz_s"]))
    return float(np.exp(-math.sqrt(.45*pd+.35*td+.15*od+.05*sd)-.10*min(drift/.15,3.)))
def off_support(raw,path):
    try:y=mem.norm(raw.astype(float))
    except ValueError:return -1e6
    t=np.arange(len(path));best=-1e6
    for start,end in itertools.product(range(-2,3),(-2,0,2)):
        p=path+start+np.rint(t*end/max(1,len(t)-1)).astype(int)
        if p.min()<0 or p.max()>=y.shape[1]:continue
        v=y[t,p];score=float(np.sort(v.reshape(4,4).mean(1))[1])
        best=max(best,score)
    return best
def evidence(row,R,path,entries,recurrent):
    m=row[9];sup=off_support(R,path)
    out={"valid":False,"exact":False,"all":0.,"recurrent":0.,"off":sup,"best_all":None,"best_recurrent":None}
    try:f=mem.features(row[0],path,m["df_hz"],m["ts"])
    except ValueError:return out
    out["valid"]=True;freq=m["channel0_mhz"]*1e6+int(path[0])*m["df_hz"]
    q={"features":f,"frequency_hz":freq,"provenance":{"source_scan":"QUERY_ONLY","source_tile_sha256":row[6]}}
    for e in entries:
        assert row[6]!=e["provenance"]["source_tile_sha256"]
        sc=family(f,e["features"])
        if sc>out["all"]:out["all"]=sc;out["best_all"]=e["id"]
        if e["id"] in recurrent and sc>out["recurrent"]:out["recurrent"]=sc;out["best_recurrent"]=e["id"]
        if abs(freq-e["frequency_hz"])<=H["exact"]["frequency_hz"] and abs(f["drift_hz_s"]-e["features"]["drift_hz_s"])<=H["exact"]["drift_hz_s"]:
            d=mem.compare(q,e)
            if d["grid_compatible"] and d["shape_similarity"]>=H["exact"]["shape_similarity"]:out["exact"]=True
    return out
def extract(rows,R,L,entries,recurrent):
    pre=prev.panel(rows,R);q=pre[0];z,p,_,_=shield.h.fuse(q,L,"registered_hybrid")
    alarm=z>=L["thresholds"]["registered_hybrid"];y=shield.h.labels(rows);loc=shield.rb.loc(p,rows);weak=np.array([v[4]<=1.25 for v in rows])
    ev=[evidence(row,rr,path,entries,recurrent) if al else {"valid":False,"exact":False,"all":0.,"recurrent":0.,"off":-1e6,"best_all":None,"best_recurrent":None} for row,rr,path,al in zip(rows,R,p,alarm)]
    return {"alarm":alarm,"y":y,"loc":loc,"weak":weak,"evidence":ev,"z":z,"p":p,"source_tile":[v[6] for v in rows]}
def held(d,arm,policy):
    a=d["alarm"];v=d["evidence"]
    if arm=="baseline_full":return np.zeros(len(a),bool)
    if arm=="exact_g3":return a&np.array([x["exact"] for x in v])
    if policy is None:return np.zeros(len(a),bool)
    shape_key="recurrent" if "recurrent" in arm else "all"
    sh=policy.get("shape");off=policy.get("off")
    return a&np.asarray([(sh is None or (x["valid"] and x[shape_key]>=sh)) and (off is None or x["off"]>=off) for x in v],bool)
def metrics(d,arm,policy):
    a=d["alarm"];y=d["y"];loc=d["loc"];weak=d["weak"];hh=held(d,arm,policy);keep=a&~hh
    return {"localized":int((keep&loc&(y==1)).sum()),"weak_localized":int((keep&loc&(y==1)&weak).sum()),"control_alerts":int((keep&(y==0)).sum()),"controls":int((y==0).sum()),"positives":int((y==1).sum()),"held_controls":int((hh&(y==0)).sum()),"held_positives":int((hh&(y==1)).sum()),"lost_localized":int((hh&loc&(y==1)).sum()),"lost_weak":int((hh&loc&weak&(y==1)).sum())}
def fit_policy(cal,leak):
    policies={};audit={}
    for arm in H["arms"][2:]:
        shape=[None] if arm=="off_only" else H["grid"]["shape"]
        off=H["grid"]["off"] if arm.endswith("_off") or arm=="off_only" else [None]
        trials=[]
        for sh,of in itertools.product(shape,off):
            p={"shape":sh,"off":of}
            mm=[metrics(d,arm,p) for d in cal];ll=[metrics(d,arm,p) for d in leak]
            removed=sum(x["held_controls"] for x in mm);lost=sum(x["held_positives"] for x in mm+ll)
            trials.append({"policy":p,"removed":removed,"held_positive":lost})
        valid=[t for t in trials if t["held_positive"]==0 and t["removed"]>0]
        winner=max(valid,key=lambda t:(t["removed"],t["policy"]["shape"] or 0,t["policy"]["off"] or 0)) if valid else None
        policies[arm]=winner["policy"] if winner else None
        audit[arm]={"grid_trials":len(trials),"usable_trials":len(valid),"chosen":winner,"trials":trials}
    return policies,audit
def selftest():
    rng=np.random.default_rng(32001);x=rng.normal(size=(16,128));x[:,64]+=12;p=np.full(16,64,int)
    f=mem.features(x,p,2.8,18.);assert abs(family(f,f)-1)<1e-9
    z=np.roll(x,9,axis=1);g=mem.features(z,p+9,2.8,18.);assert abs(family(f,g)-1)<1e-7
    assert family(f,{**f,"channel_width_hz":10})==0
    assert off_support(x,p)>off_support(rng.normal(size=x.shape),p)
    # No safe control-removing policy can be chosen if an identical positive would be held.
    d={"alarm":np.array([True,True]),"y":np.array([0,1]),"loc":np.array([False,True]),"weak":np.array([False,True]),"evidence":[{"valid":True,"exact":True,"all":1.,"recurrent":1.,"off":20.}]*2}
    pp,_=fit_policy([d],[]);assert all(v is None for v in pp.values())
    return {k:"PASS" for k in ("frequency_translation","identity_shape","grid_incompatibility","off_support","unsafe_policy_abstains")}
def main():
    t0=time.perf_counter();checks={**mem.selftest(),**selftest()};shield.a.init_tree();checks.update(shield.h.selftest())
    lib=load();entries=lib["entries"];recurrent=set(mem.review(lib)["derived_candidate_ids"]);assert len(recurrent)==12
    m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==sha(P)
    selected=m["selected"];assert len(selected)==8 and len({q["day"] for q in selected})==8
    assert [q["split"] for q in selected]==["component_fit"]*2+["null_calibration"]*2+["evaluation"]*4
    assert not set(q["day"] for q in selected)&set(S["shield"]["exclude_days"])
    acc=hook.ObservationAccumulator(lib,hook.manifest_artifact());seen=set()
    def fetch(q):
        rows,R,ch,hs=prev.gather(q,acc);assert len(hs)==16 and not seen&hs;seen.update(hs)
        emit("R32_ACQUIRED",{"target":q["target"],"day":q["day"],"split":q["split"],"tiles":len(hs)})
        return rows,R,ch
    fit=[]
    for q in selected[:2]:
        rows,R,ch=fetch(q);fit.append((rows,prev.panel(rows,R)))
    L=shield.h.fit(fit);cal=[];calleak=[]
    def eval_pair(q):
        rows,R,ch=fetch(q)
        return extract(rows,R,L,entries,recurrent),extract([a[0] for a in ch],np.stack([a[1] for a in ch]),L,entries,recurrent)
    for q in selected[2:4]:
        d,dd=eval_pair(q);cal.append(d);calleak.append(dd)
    policy,audit=fit_policy(cal,calleak)
    save("lock.json",{"baseline":L,"policy":policy,"calibration_trials":audit,"parent_sha256":H["parent_sha256"],"protocol_sha256":sha(P),"manifest_sha256":sha(OUT/"manifest.json")})
    lock=sha(OUT/"lock.json");emit("R32_LOCK",{"sha256":lock,"policies":policy,"calibration":{k:{"grid_trials":v["grid_trials"],"usable_trials":v["usable_trials"],"chosen":v["chosen"]} for k,v in audit.items()}})
    out={};leak={};records={}
    for q in selected[4:]:
        d,dd=eval_pair(q);name=q["target"];out[name]={k:metrics(d,k,policy.get(k)) for k in H["arms"]};leak[name]={k:metrics(dd,k,policy.get(k)) for k in H["arms"]}
        records[name]={"evidence":d["evidence"],"alarm":d["alarm"].tolist(),"y":d["y"].tolist(),"localized":d["loc"].tolist(),"weak":d["weak"].tolist(),"path":d["p"].tolist(),"scores":d["z"].tolist(),"source_tile":d["source_tile"],"held":{k:held(d,k,policy.get(k)).tolist() for k in H["arms"]},"coverage":{"exact":sum(x["exact"] for x in d["evidence"]),"shape_all_above_0_9":sum(x["all"]>=.9 for x in d["evidence"]),"shape_recurrent_above_0_9":sum(x["recurrent"]>=.9 for x in d["evidence"])}}
        save("partial.json",out);emit("R32_TARGET",{"target":name,"metrics":out[name],"coverage":records[name]["coverage"]})
    assert sha(OUT/"lock.json")==lock
    fields=list(next(iter(out.values()))["baseline_full"]);total={k:{f:sum(v[k][f] for v in out.values()) for f in fields} for k in H["arms"]};lt={k:{f:sum(v[k][f] for v in leak.values()) for f in fields} for k in H["arms"]}
    assert all(v["controls"]==512 and v["positives"]==512 for v in total.values()) and all(v["positives"]==128 for v in lt.values())
    b=total["baseline_full"];v=total["shape_recurrent_off"];gate={"less_controls":v["control_alerts"]<b["control_alerts"],"increment_over_off_only":v["control_alerts"]<total["off_only"]["control_alerts"],"no_held_positives":v["held_positives"]==0,"no_lost_localized":v["lost_localized"]==0,"no_lost_weak":v["lost_weak"]==0,"no_shared_loss":lt["shape_recurrent_off"]["held_positives"]==0}
    res={"status":"FRESH_ARCHIVE_PAIRED_DEVELOPMENT","gate":"PASS_DEVELOPMENT" if all(gate.values()) else "HOLD","gate_checks":gate,"modes":total,"leak_challenge":lt,"per_target":out,"coverage":{k:sum(v["coverage"][k] for v in records.values()) for k in next(iter(records.values()))["coverage"]},"policies":policy,"selftests":checks,"targets":[{k:q[k] for k in ("target","day","split")} for q in selected],"lock_sha256":lock,"memory_sha256":H["parent_sha256"],"scope":S["scope"],"production_adoption":False}
    save("scores.json",records);save("result.json",res);emit("R32_SUMMARY",{k:v for k,v in res.items() if k!="per_target"})
    assert len(acc.observations)==128;new,delta,replay=acc.seal();save("library.json",new);save("observations.json",acc.observations);rev=mem.review(new);save("review.json",rev)
    hook_result={"status":"PASS_HOOK_DEVELOPMENT","delta":delta,"replay":replay,"observations":128,"library_sha256":sha(OUT/"library.json"),"parent_sha256":H["parent_sha256"],"verified_rfi":rev["verified_rfi"],"rfi_candidates":rev["rfi_candidates"],"unknown_without_recurrence":rev["unknown_without_recurrence"],"seconds":time.perf_counter()-t0,"production_adoption":False}
    save("hook_result.json",hook_result);emit("R32_HOOK_SUMMARY",hook_result)
if __name__=="__main__":
    if sys.argv[1]=="freeze":shield.freeze()
    elif sys.argv[1]=="selftest":emit("R32_SELFTEST",selftest())
    elif sys.argv[1]=="scan":main()
    else:raise ValueError("phase")
