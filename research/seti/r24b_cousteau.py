#!/usr/bin/env python3
"""R24B: Cousteau-inspired background shape/scale correction and stratified calibration."""
import hashlib,json,pathlib,sys,time,statistics
import numpy as np
import r24_odontocrystal as od
import r20d_weak_joint as reg
import r17c_frozen_holdout as hold
import juxtapose_undina_r18b2_vector as vec
import juxtapose_undina_r18bc as old
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r24b_cousteau_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r24b");OUT.mkdir(parents=True,exist_ok=True)
od.S=S;od.OUT=OUT;hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))}
hold.OUT=OUT;hold.P=P;old.OUT=OUT;old.S={**old.S,"seed":S["seed"]};eng.OUT=OUT
def emit(k,d):print(k,json.dumps(d),flush=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def background(V):
    med=np.median(V,axis=1)
    scale=1.4826*np.median(np.abs(V-med[:,None,:]),axis=1)
    # Shrink unreliable finite-time column scales toward one local robust-noise unit.
    scale=np.clip(.5*scale+.5,1,4)
    X=(V-.75*np.maximum(med,0)[:,None,:])/scale[:,None,:]
    regime=np.quantile(med,.95,axis=1)
    assert np.isfinite(X).all() and np.isfinite(regime).all()
    return X,regime
def cutoff(x):
    a=np.sort(np.asarray(x,dtype=float))[::-1]
    assert len(a)>=64
    return float(np.nextafter(a[int(len(a)*.05)],np.inf))
def get(V,kind):
    if kind=="raw":
        a,p=reg.reg_sliced(V,.7,True);return a,p,np.zeros_like(p,bool)
    if kind=="fast":return od.sliced(V,0,None)
    X,_=background(V)
    return od.sliced(X,4 if kind=="gap" else 0,3)
def thresholds(a,y,g):
    return [cutoff(a[(y==0)&(g==i)]) for i in (0,1)]
def calibrate(rows):
    t=time.perf_counter();V=vec.normalize([r[0] for r in rows]).astype(np.float64);y=np.asarray([r[1] for r in rows])
    _,f=background(V);split=float(np.median(f[y==0]));g=(f>split).astype(int)
    s={k:get(V,k)[0] for k in ("raw","residual","gap")}
    pool={k:cutoff(a[y==0]) for k,a in s.items()};th={k:thresholds(a,y,g) for k,a in s.items()}
    scale={k:[max(.001,float(np.std(a[(y==0)&(g==i)]))) for i in (0,1)] for k,a in s.items()}
    def z(k):return (s[k]-np.asarray(th[k])[g])/np.asarray(scale[k])[g]
    u=np.maximum(z("raw"),z("residual"));ut=thresholds(u,y,g)
    return {"split":split,"pooled":pool,"stratified":th,"scale":scale,"union":ut,"calibration_controls":[int(np.sum((g==i)&(y==0))) for i in (0,1)],"calibration_alerts":{"raw":[int(np.sum((s["raw"]>=th["raw"][i])&(g==i)&(y==0))) for i in (0,1)],"residual":[int(np.sum((s["residual"]>=th["residual"][i])&(g==i)&(y==0))) for i in (0,1)],"union":[int(np.sum((u>=ut[i])&(g==i)&(y==0))) for i in (0,1)]},"seconds":time.perf_counter()-t}
def route(V,k,L):
    _,f=background(V);g=(f>L["split"]).astype(int)
    if k=="raw_residual_union":
        a=get(V,"raw");b=get(V,"residual")
        def z(q,name):return (q[0]-np.asarray(L["stratified"][name])[g])/np.asarray(L["scale"][name])[g]
        za=z(a,"raw");zb=z(b,"residual");pick=za>=zb;v=np.maximum(za,zb)
        return v>=np.asarray(L["union"])[g],np.where(pick[:,None],a[1],b[1]),np.where(pick[:,None],a[2],b[2]),g
    name="fast" if k=="raw_fast_exact" else "raw" if k.startswith("raw") else "gap" if "gap" in k else "residual"
    a,p,m=get(V,name);key="raw" if name=="fast" else name
    th=L["pooled"][key] if "pooled" in k or name=="fast" else np.asarray(L["stratified"][key])[g]
    return a>=th,p,m,g
def selftest():
    rng=np.random.default_rng(24201);V=rng.normal(size=(16,16,128));V[:,:,40]+=15
    X,g=background(V);assert np.all(np.median(X,axis=1)[:,40]<np.median(V,axis=1)[:,40])
    a,_,_=get(V,"raw");b,_,_=get(V,"fast");assert np.allclose(a,b,atol=1e-10,rtol=1e-12)
    return {"shape_scale_control":"PASS","raw_exact_equivalence_cases":16}
def bench(rows,L):
    logs={k:[] for k in S["modes"]};out={};results={}
    for rep in range(S["repeats"]):
        modes=list(S["modes"]);np.random.default_rng(24202+rep).shuffle(modes)
        for k in modes:
            t=time.perf_counter();cpu=time.process_time();V=vec.normalize([r[0] for r in rows]).astype(np.float64)
            q=route(V,k,L);w=time.perf_counter()-t;c=time.process_time()-cpu
            if k in results:assert np.array_equal(q[0],results[k][0])
            results[k]=q;logs[k].append((w,c))
    V=vec.normalize([r[0] for r in rows]).astype(np.float64)
    a=get(V,"raw")[0];b=get(V,"fast")[0];assert np.allclose(a,b,atol=1e-10,rtol=1e-12)
    assert np.array_equal(results["raw_pooled"][0],results["raw_fast_exact"][0])
    y=np.array([r[1] for r in rows]);truth=np.stack([r[7] for r in rows]);active=np.stack([r[8] for r in rows])
    ref=results["raw_pooled"][0]
    for k,(p,path,m,g) in results.items():
        loc=((np.abs(path-truth)<=2)&active).sum(1)>=.75*active.sum(1)
        out[k]={**od.measure(p,path,m,rows),"wall":statistics.median(a[0] for a in logs[k]),"cpu":statistics.median(a[1] for a in logs[k]),"regimes":{str(i):{"localized":int(np.sum(p&loc&(y==1)&(g==i))),"positives":int(np.sum((y==1)&(g==i))),"control_alerts":int(np.sum(p&(y==0)&(g==i))),"controls":int(np.sum((y==0)&(g==i)))} for i in (0,1)},"paired":{ "added_positive_alerts":int(np.sum(p&~ref&(y==1))),"lost_positive_alerts":int(np.sum(~p&ref&(y==1))),"added_controls":int(np.sum(p&~ref&(y==0))),"removed_controls":int(np.sum(~p&ref&(y==0)))}}
    return out
def main():
    checks=selftest();m=json.loads((OUT/"manifest.json").read_text())
    assert len(m["selected"])==6 and m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest()
    assert len({q["target"] for q in m["selected"]})==6 and all(q["target"] not in hold.S["exclude"] for q in m["selected"])
    t=time.perf_counter();cal=[]
    for q in m["selected"][:2]:cal+=od.gather(q)
    load=time.perf_counter()-t;L=calibrate(cal);L["load_seconds"]=load;L["targets"]=[q["target"] for q in m["selected"][:2]]
    save("lock.json",L);emit("R24B_LOCK",L)
    seen={r[6] for r in cal};out={};load=0
    for q in m["selected"][2:]:
        t=time.perf_counter();rows=od.gather(q);load+=time.perf_counter()-t;h={r[6] for r in rows};assert not seen&h;seen|=h
        out[q["target"]]=bench(rows,L);save("partial.json",out)
    keys=("detected","localized","positives","control_alerts","controls","wall","cpu")
    total={k:{a:sum(x[k][a] for x in out.values()) for a in keys} for k in S["modes"]}
    for k in total:
        for f,vals in (("amplitudes",S["amplitudes"]),("gaps",S["gap_levels"]),("models",S["injection_models"]),("regimes",[0,1])):
            total[k][f]={str(v):{a:sum(x[k][f][str(v)][a] for x in out.values()) for a in (("localized","positives","control_alerts","controls") if f=="regimes" else keys[:5])} for v in vals}
    d={"status":"FRESH_ARCHIVE_SYNTHETIC_EXAM","samples":1024,"targets":list(out),"selftests":checks,"lock":L,"modes":total,"per_target":out,"evaluation_load_seconds":load,"raw_exact_decision_mismatches":0,"scope":"Candidate-local robust shape correction can attenuate stationary sky signals. Stratified calibration is trained on distinct targets, not target labels. No real source inference or production adoption.","production_adoption":False}
    save("result.json",d);emit("R24B_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="scan":main()
    else:raise ValueError("phase")
