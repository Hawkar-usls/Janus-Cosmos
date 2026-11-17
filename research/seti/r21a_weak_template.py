#!/usr/bin/env python3
"""R21A: frozen smooth-template bank and robust clipped integration on fresh spectra."""
import hashlib,json,pathlib,sys,statistics,time,math
import numpy as np
import r20d_weak_joint as reg
import r20a_vector_dp as dp
import r19e_integrated_exam as prior
import r17c_frozen_holdout as hold
import juxtapose_undina_r18b2_vector as vec
import juxtapose_undina_r18bc as old
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r21a_weak_template_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r21a");OUT.mkdir(parents=True,exist_ok=True)
hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))}
hold.OUT=OUT;hold.P=P;eng.OUT=OUT
old.OUT=OUT;old.S={**old.S,"seed":S["seed"]}
def emit(k,d):print(k,json.dumps(d),flush=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def gather(q):
    data=old.gather(q);rng=np.random.default_rng(int.from_bytes(hashlib.sha256((S["seed"]+"|"+q["target"]).encode()).digest()[:8],"big"))
    plan=[(a,m,g) for a in S["amplitudes"] for m in S["injection_models"] for g in S["gap_levels"]]*4
    rng.shuffle(plan);rows=[]
    for k,(raw,label,hard,sha) in enumerate(v for v in data if v[1]==0):
        n,w=raw.shape;t=np.arange(n);u=t-(n-1)/2
        amp,model,gap=plan[k]
        if model=="straight":
            d=np.rint(t*int(rng.integers(-32,33))/(n-1)).astype(int)
        elif model=="off_grid_quadratic":
            slope=float(rng.choice([-1.05,-.85,.35,1.05]));a=float(rng.choice([-.137,-.067,.063,.117]))
            d=np.rint(slope*u+a*u*u).astype(int)
        elif model=="slope_break":
            a,b=rng.choice([-1.6,-.7,.6,1.7],2,replace=False);cut=int(rng.integers(5,11))
            d=np.rint(a*np.minimum(t,cut)+b*np.maximum(t-cut,0)).astype(int)
        elif model=="smooth_sinusoid":
            d=np.rint(float(rng.choice([-.6,.6]))*u+float(rng.choice([3.2,5.4]))*np.sin(2*np.pi*t/n+float(rng.uniform(0,2*np.pi)))).astype(int)
        else:raise ValueError(model)
        d-=d[0];lo=8+max(0,-int(d.min()));hi=w-8-max(0,int(d.max()))
        p=int(rng.integers(lo,hi))+d
        assert p.min()>=0 and p.max()<w and np.max(np.abs(np.diff(p)))<=3
        active=np.ones(n,bool)
        if gap:active[rng.choice(n,size=int(n*gap),replace=False)]=False
        scale=1.4826*np.median(np.abs(raw-np.median(raw,axis=1)[:,None]),axis=1)
        assert np.all(scale>0)
        for y in (0,1):
            x=raw.copy()
            if y:x[t[active],p[active]]+=amp*scale[active]
            rows.append((x,y,model,gap,amp,hard,sha))
    assert len(rows)==256 and len(plan)==128
    return rows
def bank(n):
    t=np.arange(n);u=t-(n-1)/2
    paths=[np.rint(t*s/(n-1)).astype(int) for s in range(-32,33)]
    for s in S["template_slopes"]:
        for a in S["template_curvatures"]:
            d=np.rint(s*u+a*u*u).astype(int);d-=d[0]
            if np.max(np.abs(np.diff(d)))<=3:paths.append(d)
    return np.asarray(sorted(set(tuple(p) for p in paths)),dtype=np.int32)
def template(V,clip=None):
    V=np.asarray(V,dtype=np.float64)
    if clip is not None:V=np.clip(V,-clip,clip)
    B,n,w=V.shape;D=bank(n);N=len(D)
    best=np.full(B,-np.inf);bestp=np.empty((B,n),dtype=np.int32);best_template=np.full(B,-1,dtype=np.int32)
    cols=np.arange(w);rows=np.arange(B)
    for i in range(0,N,12):
        ds=D[i:i+12];T=len(ds);acc=np.zeros((B,T,w),dtype=np.float64);valid=np.ones((T,w),bool)
        for k in range(n):
            ix=cols[None,:]+ds[:,k,None]
            valid&=(ix>=0)&(ix<w)
            acc+=V[:,k,np.clip(ix,0,w-1)]
        acc=np.where(valid[None],acc,-np.inf).reshape(B,-1)
        ix=np.argmax(acc,axis=1);v=acc[rows,ix];good=v>best
        best[good]=v[good];best_template[good]=i+ix[good]//w
        bestp[good]=(ix[good]%w)[:,None]+ds[ix[good]//w]
    assert np.all(best_template>=0)
    assert np.all((bestp>=0)&(bestp<w)) and np.all(np.abs(np.diff(bestp,axis=1))<=3)
    assert np.array_equal(bestp-bestp[:,0,None],D[best_template])
    replay=np.sum(V[rows[:,None],np.arange(n)[None,:],bestp],axis=1,dtype=np.float64)
    assert np.allclose(replay,best,rtol=1e-12,atol=1e-10)
    return replay/np.sqrt(n),bestp
def template_sliced(V,clip):
    a=[];p=[]
    for i in range(0,len(V),32):
        s,q=template(V[i:i+32],clip);a.append(s);p.append(q)
    return np.concatenate(a),np.concatenate(p)
def selftest():
    rng=np.random.default_rng(21001);V=rng.normal(size=(4,16,40));V[0]=0
    D=bank(16)
    for clip in [None,3,5]:
        a,p=template(V,clip);X=V if clip is None else np.clip(V,-clip,clip)
        for i,x in enumerate(X):
            best=-np.inf
            for d in D:
                for c in range(max(0,-int(d.min())),min(x.shape[1],x.shape[1]-int(d.max()))):
                    best=max(best,float(np.sum(x[np.arange(len(d)),c+d])/4))
            assert abs(best-a[i])<1e-10
    return {"template_count":len(D),"brute_force_cases":12,"status":"PASS"}
def score(V,k):
    if k=="straight":return np.asarray([hold.classical(v) for v in V]),0
    if k=="dp":
        a,p=dp.sliced(V,256,True);return a,len(p)
    if k=="regular":
        a,p=reg.reg_sliced(V,S["penalty"],True);return a,len(p)
    if k=="template_raw":a,p=template_sliced(V,None);return a,len(p)
    if k.startswith("template_clip"):a,p=template_sliced(V,int(k.split("clip")[-1]));return a,len(p)
    raise ValueError(k)
def cutoff(x):
    a=sorted(map(float,x),reverse=True)
    return float(np.nextafter(a[math.floor(len(a)*S["calibration_control_fraction"])],np.inf))
def met(p,y,m):
    return {"detected":int(np.sum(p&(y==1)&m)),"positives":int(np.sum((y==1)&m)),"control_alerts":int(np.sum(p&(y==0)&m)),"controls":int(np.sum((y==0)&m))}
BASE=["straight","dp","regular","template_raw","template_clip3","template_clip5"]
def calibrate(rows):
    t=time.perf_counter();V=vec.normalize([r[0] for r in rows]);y=np.array([r[1] for r in rows]);sc={k:score(V,k)[0] for k in BASE}
    th={k:cutoff(v[y==0]) for k,v in sc.items()}
    assert all(v>0 for v in th.values())
    def joint(k,l):
        v=np.maximum(sc[k]/th[k],sc[l]/th[l]);t=cutoff(v[y==0])
        return {k:th[k]*t,l:th[l]*t}
    jreg=joint("straight","regular")
    candidates=[]
    for k in BASE[3:]:
        jt=joint(k,"regular");p=(sc[k]>=jt[k])|(sc["regular"]>=jt["regular"])
        candidates.append({"template":k,"thresholds":jt,**met(p,y,np.ones(len(y),bool))})
    pick=min(candidates,key=lambda x:(-x["detected"],x["control_alerts"],x["template"]))
    return {"thresholds":th,"joint_regular":jreg,"joint_template_regular":pick["thresholds"],"selected_template":pick["template"],"joint_trials":candidates,"calibration_seconds":time.perf_counter()-t,"calibration_metrics":{k:met(sc[k]>=th[k],y,np.ones(len(y),bool)) for k in BASE}}
def bench(rows,lock):
    X=[r[0] for r in rows];y=np.array([r[1] for r in rows]);mod=np.array([r[2] for r in rows]);gap=np.array([r[3] for r in rows]);amp=np.array([r[4] for r in rows]);hard=np.array([r[5] for r in rows])
    result={};rec={k:[] for k in S["modes"]}
    for rep in range(S["repeats"]):
        modes=list(S["modes"]);np.random.default_rng(21040+rep).shuffle(modes)
        for k in modes:
            t=time.perf_counter();cpu=time.process_time();V=vec.normalize(X)
            if k in BASE:
                a,w=score(V,k);p=a>=lock["thresholds"][k]
            else:
                th=lock[k];keys=list(th);assert len(keys)==2
                a,w=score(V,keys[0]);b,w2=score(V,keys[1]);w+=w2
                p=(a>=th[keys[0]])|(b>=th[keys[1]])
            dt=time.perf_counter()-t;dc=time.process_time()-cpu
            if k in result:assert np.array_equal(result[k],p)
            result[k]=p;rec[k].append({"wall":dt,"cpu":dc,"path_replays":w})
    b=result["joint_regular"]
    return {k:{**met(p,y,np.ones(len(y),bool)),"amplitudes":{str(a):met(p,y,amp==a) for a in S["amplitudes"]},"models":{a:met(p,y,mod==a) for a in S["injection_models"]},"gaps":{str(a):met(p,y,gap==a) for a in S["gap_levels"]},"hard":met(p,y,hard),"paired_vs_joint_regular":{"added_positive":int(np.sum(p&~b&(y==1))),"lost_positive":int(np.sum(~p&b&(y==1))),"added_control":int(np.sum(p&~b&(y==0))),"removed_control":int(np.sum(~p&b&(y==0)))},"wall":statistics.median(r["wall"] for r in rec[k]),"cpu":statistics.median(r["cpu"] for r in rec[k]),"path_replays":rec[k][0]["path_replays"],"repeats":rec[k]} for k,p in result.items()}
def main():
    check=selftest();m=json.loads((OUT/"manifest.json").read_text())
    assert m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest() and len(m["selected"])==6
    assert len({q["target"] for q in m["selected"]})==6 and all(q["target"] not in hold.S["exclude"] for q in m["selected"])
    t=time.perf_counter();cal=[]
    for q in m["selected"][:2]:cal+=gather(q)
    load=time.perf_counter()-t;lock=calibrate(cal);lock["load_seconds"]=load;lock["targets"]=[q["target"] for q in m["selected"][:2]]
    save("lock.json",lock);emit("R21A_LOCK",lock)
    seen={r[6] for r in cal};out={};load=0
    for q in m["selected"][2:]:
        t=time.perf_counter();rows=gather(q);load+=time.perf_counter()-t
        assert not seen&{r[6] for r in rows};seen|={r[6] for r in rows}
        out[q["target"]]=bench(rows,lock);save("partial.json",out);emit("R21A_TARGET",{"target":q["target"],"metrics":out[q["target"]]})
    keys=("detected","positives","control_alerts","controls","wall","cpu","path_replays")
    total={k:{a:sum(x[k][a] for x in out.values()) for a in keys} for k in S["modes"]}
    for k in total:
        for field,values in (("amplitudes",[str(a) for a in S["amplitudes"]]),("models",S["injection_models"]),("gaps",[str(a) for a in S["gap_levels"]])):
            total[k][field]={v:{a:sum(x[k][field][v][a] for x in out.values()) for a in keys[:4]} for v in values}
        total[k]["paired_vs_joint_regular"]={a:sum(x[k]["paired_vs_joint_regular"][a] for x in out.values()) for a in ("added_positive","lost_positive","added_control","removed_control")}
    d={"status":"FRESH_ARCHIVE_SYNTHETIC_EXAM","samples":1024,"targets":list(out),"lock":lock,"selftest":check,"modes":total,"per_target":out,"evaluation_load_seconds":load,"scope":"Amplitude and morphology balanced weak injected signals on fresh archive OFF. Template bank finite, per-path exact numeric witness; no astronomical origin inference. Realized false alert differences must be considered, not merely calibration parity.","production_adoption":False}
    save("result.json",d);emit("R21A_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="scan":main()
    else:raise ValueError("phase")
