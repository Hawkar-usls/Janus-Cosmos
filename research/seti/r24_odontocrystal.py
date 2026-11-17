#!/usr/bin/env python3
"""R24: exact appositional growth of observed spectral evidence, inspired by OdontoForge."""
import hashlib,json,pathlib,sys,time,statistics,itertools,math
import numpy as np
import r21a_weak_template as prior
import r20d_weak_joint as reg
import r17c_frozen_holdout as hold
import juxtapose_undina_r18b2_vector as vec
import juxtapose_undina_r18bc as old
import r9_expanded_search as eng
import donor_odontoforge_pnp as donor
P=pathlib.Path("research/seti/r24_odonto_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r24");OUT.mkdir(parents=True,exist_ok=True)
hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))}
hold.OUT=OUT;hold.P=P;eng.OUT=OUT;old.OUT=OUT;old.S={**old.S,"seed":S["seed"]}
RAD=S["growth"]["radius"];PEN=S["growth"]["curvature_penalty"];COST=S["growth"]["skip_cost"]
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
            rows.append((x,y,model,gap,amp,hard,sha,p.copy(),active.copy()))
    assert len(rows)==256 and len(plan)==128
    return rows

def witness(V,p,missing,pen,cost,clip):
    B,n,w=V.shape
    assert p.shape==missing.shape==(B,n) and np.all((p>=0)&(p<w))
    assert np.all(np.abs(np.diff(p,axis=1))<=RAD)
    v=V[np.arange(B)[:,None],np.arange(n)[None,:],p]
    if clip is not None:v=np.minimum(v,clip)
    return (np.where(missing,-cost,v).sum(1,dtype=np.float64)-pen*np.abs(np.diff(p,n=2,axis=1)).sum(1))/np.sqrt(n)
def grow(V,G,pen=PEN,cost=COST,clip=None):
    """State is (gap count, incoming slope, frequency); no heuristic pruning.
    Pulp: immutable measurements. Dentin: additive evidence. Tubules: feasible states.
    Regenerate: bounded skip transitions, never imputed power. Enamel: replayed path.
    Dominance merges states only when future continuation and cost are identical.
    """
    V=np.asarray(V,dtype=np.float64);B,n,w=V.shape
    assert n>=2 and 0<=G<n and np.isfinite(V).all()
    Z=V if clip is None else np.minimum(V,clip)
    K=2*RAD+1;slopes=np.arange(-RAD,RAD+1);cols=np.arange(w)
    prev=np.full((B,G+1,K,w),-np.inf)
    prev[:,0,:,:]=Z[:,0,None,:]
    if G:prev[:,1,:,:]=-cost
    back=[]
    for t in range(1,n):
        # Max-plus L1 distance transform on slope axis, preserves exact maximum.
        best=prev.copy();par=np.broadcast_to(np.arange(K,dtype=np.int8)[None,None,:,None],best.shape).copy()
        for k in range(1,K):
            z=best[:,:,k-1]-pen;take=z>best[:,:,k]
            best[:,:,k]=np.where(take,z,best[:,:,k]);par[:,:,k]=np.where(take,par[:,:,k-1],par[:,:,k])
        for k in range(K-2,-1,-1):
            z=best[:,:,k+1]-pen;take=z>best[:,:,k]
            best[:,:,k]=np.where(take,z,best[:,:,k]);par[:,:,k]=np.where(take,par[:,:,k+1],par[:,:,k])
        cur=np.full_like(prev,-np.inf);parent=np.zeros(prev.shape,dtype=np.int8);skip=np.zeros(prev.shape,dtype=bool)
        for k,s in enumerate(slopes):
            lo=max(0,int(s));hi=min(w,w+int(s));src=slice(lo-int(s),hi-int(s));dst=slice(lo,hi)
            v=best[:,:,k,src];a=v+Z[:,t,None,dst]
            cur[:,:,k,dst]=a;parent[:,:,k,dst]=par[:,:,k,src]
            if G:
                alt=v[:,:-1]-cost;take=alt>a[:,1:]
                cur[:,1:,k,dst]=np.where(take,alt,a[:,1:])
                parent[:,1:,k,dst]=np.where(take,par[:,:-1,k,src],par[:,1:,k,src])
                skip[:,1:,k,dst]=take
        prev=cur;back.append((parent,skip))
    flat=prev.reshape(B,-1);r=np.arange(B);j=np.argmax(flat,axis=1);score=flat[r,j]/np.sqrt(n)
    g=j//(K*w);k=(j//w)%K;x=j%w
    p=np.empty((B,n),dtype=np.int32);missing=np.zeros((B,n),bool)
    for t in range(n-1,0,-1):
        p[:,t]=x;par,sk=back[t-1];m=sk[r,g,k,x];missing[:,t]=m
        nk=par[r,g,k,x];x=x-slopes[k];g=g-m.astype(int);k=nk
    p[:,0]=x;missing[:,0]=g==1
    assert np.all((g==0)|(g==1)) and np.all(missing.sum(1)<=G)
    check=witness(V,p,missing,pen,cost,clip)
    assert np.allclose(check,score,rtol=1e-12,atol=1e-10)
    return check,p,missing
def sliced(V,G,clip):
    a=[];p=[];m=[]
    for i in range(0,len(V),16):
        x,y,z=grow(V[i:i+16],G,clip=clip);a.append(x);p.append(y);m.append(z)
    return np.concatenate(a),np.concatenate(p),np.concatenate(m)
def straight(V):
    B,n,w=V.shape;t=np.arange(n);r=np.arange(B)
    best=np.full(B,-np.inf);path=np.zeros((B,n),int)
    for sh in range(-32,33):
        d=np.rint(t*sh/(n-1)).astype(int)
        if np.max(np.abs(np.diff(d)))>RAD:continue
        lo=max(0,-int(d.min()));hi=min(w,w-int(d.max()))
        if hi<=lo:continue
        a=sum(V[:,i,lo+k:hi+k] for i,k in enumerate(d))/np.sqrt(n)
        j=np.argmax(a,axis=1);z=a[r,j];take=z>best
        best[take]=z[take];path[take]=lo+j[take,None]+d
    assert np.allclose(best,V[r[:,None],t[None,:],path].sum(1)/np.sqrt(n),atol=1e-10)
    return best,path,np.zeros_like(path,bool)
def selftest():
    donor.selftest()
    rng=np.random.default_rng(24001);V=rng.normal(size=(6,4,5));V[0]=0
    cases=0
    for G,clip in [(0,None),(1,None),(2,3)]:
        z,p,m=grow(V,G,clip=clip);X=V if clip is None else np.minimum(V,clip)
        for i,v in enumerate(X):
            best=-np.inf
            for path in itertools.product(range(5),repeat=4):
                if max(abs(a-b) for a,b in zip(path,path[1:]))>RAD:continue
                q=np.array(path);power=v[np.arange(4),q]
                for mask in itertools.product((False,True),repeat=4):
                    if sum(mask)>G:continue
                    s=(np.where(mask,-COST,power).sum()-PEN*np.abs(np.diff(q,n=2)).sum())/2
                    best=max(best,s)
            assert abs(best-z[i])<1e-10
            cases+=1
    V=rng.normal(size=(32,16,128));a,p,m=sliced(V,0,None);b,_=reg.reg_sliced(V,PEN,True)
    assert np.allclose(a,b,atol=1e-10,rtol=1e-12)
    return {"donor_selftest":"PASS","brute_force_cases":cases,"zero_gap_regular_equivalence_cases":32}
def base(V,k):
    if k=="straight":return straight(V)
    if k=="regular":
        a,p=reg.reg_sliced(V,PEN,True);return a,p,np.zeros_like(p,bool)
    if k=="odonto_gap":return sliced(V,4,None)
    if k=="odonto_clip":return sliced(V,0,3)
    if k=="odonto_gap_clip":return sliced(V,4,3)
    raise ValueError(k)
BASE=["straight","regular","odonto_gap","odonto_clip","odonto_gap_clip"]
def cutoff(x):
    a=sorted(map(float,x),reverse=True)
    return float(np.nextafter(a[int(len(a)*S["calibration_control_fraction"])],np.inf))
def measure(pred,p,m,rows):
    y=np.array([r[1] for r in rows]);amp=np.array([r[4] for r in rows]);gap=np.array([r[3] for r in rows]);model=np.array([r[2] for r in rows])
    truth=np.stack([r[7] for r in rows]);active=np.stack([r[8] for r in rows])
    loc=((np.abs(p-truth)<=2)&active).sum(1)>=.75*active.sum(1)
    def met(mask):
        return {"detected":int(np.sum(pred&mask&(y==1))),"localized":int(np.sum(pred&loc&mask&(y==1))),"positives":int(np.sum(mask&(y==1))),"control_alerts":int(np.sum(pred&mask&(y==0))),"controls":int(np.sum(mask&(y==0)))}
    return {**met(np.ones(len(y),bool)),"amplitudes":{str(a):met(amp==a) for a in S["amplitudes"]},"gaps":{str(a):met(gap==a) for a in S["gap_levels"]},"models":{a:met(model==a) for a in S["injection_models"]},"skipped_layers_in_alerts":int(m[pred].sum()),"correct_actual_gap_skips":int((m&~active&pred[:,None]&(y==1)[:,None]).sum()),"skipped_active_layers_in_positive_alerts":int((m&active&pred[:,None]&(y==1)[:,None]).sum())}
def union(a,b,ta,tb):
    z=np.maximum(a[0]/ta,b[0]/tb);choose=a[0]/ta>=b[0]/tb
    return z,np.where(choose[:,None],a[1],b[1]),np.where(choose[:,None],a[2],b[2])
def calibrate(rows):
    t=time.perf_counter();V=vec.normalize([r[0] for r in rows]).astype(np.float64);y=np.asarray([r[1] for r in rows])
    a={k:base(V,k) for k in BASE};th={k:cutoff(q[0][y==0]) for k,q in a.items()};assert all(v>0 for v in th.values())
    j={};metrics={}
    for k in BASE:
        q=a[k];metrics[k]=measure(q[0]>=th[k],q[1],q[2],rows)
        if k=="straight":continue
        q=union(a["straight"],a[k],th["straight"],th[k]);t0=cutoff(q[0][y==0]);m=measure(q[0]>=t0,q[1],q[2],rows)
        j[k]={"threshold":t0,"metrics":m}
    # Freeze Odonto variant selection by localized weak recovery over calibration only.
    candidates=BASE[2:]
    pick=min(candidates,key=lambda k:(-sum(j[k]["metrics"]["amplitudes"][str(x)]["localized"] for x in (.75,1.25)),-j[k]["metrics"]["localized"],j[k]["metrics"]["control_alerts"],k))
    return {"thresholds":th,"joint_regular":j["regular"]["threshold"],"joint_odonto":j[pick]["threshold"],"selected_odonto":pick,"trials":j,"base_metrics":metrics,"seconds":time.perf_counter()-t}
def bench(rows,lock):
    a={};out={};logs={k:[] for k in S["modes"]};preds={};th=lock["thresholds"]
    for rep in range(S["repeats"]):
        order=list(S["modes"]);np.random.default_rng(24003+rep).shuffle(order)
        for k in order:
            t=time.perf_counter();cpu=time.process_time();V=vec.normalize([r[0] for r in rows]).astype(np.float64)
            if k in BASE:q=base(V,k);t0=th[k]
            else:
                z="regular" if k=="joint_regular" else lock["selected_odonto"]
                q=union(base(V,"straight"),base(V,z),th["straight"],th[z]);t0=lock[k]
            pred=q[0]>=t0;dt=time.perf_counter()-t;dc=time.process_time()-cpu
            if k in preds:assert np.array_equal(pred,preds[k])
            preds[k]=pred;a[k]=q;logs[k].append({"wall":dt,"cpu":dc})
    y=np.array([r[1] for r in rows]);b=preds["joint_regular"]
    for k,p in preds.items():
        q=a[k];out[k]={**measure(p,q[1],q[2],rows),"wall":statistics.median(x["wall"] for x in logs[k]),"cpu":statistics.median(x["cpu"] for x in logs[k]),"paired_vs_joint_regular":{"added_alerts_positive":int(np.sum(p&~b&(y==1))),"lost_alerts_positive":int(np.sum(~p&b&(y==1))),"added_controls":int(np.sum(p&~b&(y==0))),"removed_controls":int(np.sum(~p&b&(y==0)))}}
    return out
def main():
    checks=selftest();m=json.loads((OUT/"manifest.json").read_text())
    assert len(m["selected"])==6 and m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest()
    assert len({q["target"] for q in m["selected"]})==6 and all(q["target"] not in hold.S["exclude"] for q in m["selected"])
    t=time.perf_counter();cal=[]
    for q in m["selected"][:2]:cal+=gather(q)
    load=time.perf_counter()-t;lock=calibrate(cal);lock["load_seconds"]=load;lock["targets"]=[q["target"] for q in m["selected"][:2]]
    save("lock.json",lock);emit("R24_LOCK",lock)
    seen={r[6] for r in cal};out={};load=0
    for q in m["selected"][2:]:
        t=time.perf_counter();rows=gather(q);load+=time.perf_counter()-t
        assert not seen&{r[6] for r in rows};seen|={r[6] for r in rows}
        out[q["target"]]=bench(rows,lock);save("partial.json",out)
        emit("R24_TARGET",{"target":q["target"],"metrics":out[q["target"]]})
    keys=("detected","localized","positives","control_alerts","controls","wall","cpu","skipped_layers_in_alerts","correct_actual_gap_skips","skipped_active_layers_in_positive_alerts")
    total={k:{a:sum(x[k][a] for x in out.values()) for a in keys} for k in S["modes"]}
    for k in total:
        for f,vals in (("amplitudes",S["amplitudes"]),("gaps",S["gap_levels"]),("models",S["injection_models"])):
            total[k][f]={str(v):{a:sum(x[k][f][str(v)][a] for x in out.values()) for a in keys[:5]} for v in vals}
        total[k]["paired_vs_joint_regular"]={a:sum(x[k]["paired_vs_joint_regular"][a] for x in out.values()) for a in next(iter(out.values()))[k]["paired_vs_joint_regular"]}
    d={"status":"FRESH_ARCHIVE_SYNTHETIC_EXAM","samples":1024,"targets":list(out),"lock":lock,"selftests":checks,"modes":total,"per_target":out,"evaluation_load_seconds":load,"scope":"Adaptation of OdontoForge growth/regeneration, not direct SAT solver on spectra. Exact finite score model; skip cost is not evidence. Positive alerts not necessarily localized or astronomical; no production adoption.","production_adoption":False}
    save("result.json",d);emit("R24_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="scan":main()
    elif sys.argv[1]=="selftest":emit("R24_SELFTEST",selftest())
    else:raise ValueError("phase")
