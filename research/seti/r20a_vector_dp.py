#!/usr/bin/env python3
"""R20A: batched exact bounded-step DP with trace certificates and admissible bounds."""
import hashlib,json,pathlib,statistics,sys,time,math,resource
import numpy as np
import r19e_integrated_exam as prior
import r19a_physarum_curves as curved
import r17c_frozen_holdout as hold
import juxtapose_undina_r18b2_vector as vec
import juxtapose_undina_r18bc as old
import r9_expanded_search as eng

P=pathlib.Path("research/seti/r20a_dp_vector_protocol.json")
S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r20a");OUT.mkdir(parents=True,exist_ok=True)
hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))}
hold.OUT=OUT;hold.P=P;eng.OUT=OUT
old.OUT=OUT;old.S={**old.S,"seed":S["seed"]}
prior.S={**prior.S,**S};prior.OUT=OUT
R=tuple(S["rank"]);RAD=S["radius"]
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def emit(k,d):print(k,json.dumps(d),flush=True)
def scalar(v):
    x=v.astype(np.float64,copy=False);p=curved.dp(x)
    return float(np.sum(x[np.arange(len(x)),p],dtype=np.float64)/np.sqrt(len(x))),p
def batch(V,trace=True):
    V=np.asarray(V,dtype=np.float64)
    b,n,w=V.shape
    assert b and n>1 and w>RAD and np.isfinite(V).all()
    cur=V[:,0].copy();back=[]
    for i in range(1,n):
        pad=np.pad(cur,((0,0),(RAD,RAD)),constant_values=-np.inf)
        best=pad[:,:w].copy()
        if trace:step=np.full((b,w),-RAD,dtype=np.int8)
        for k in range(1,2*RAD+1):
            nxt=pad[:,k:k+w]
            if trace:
                good=nxt>best
                step[good]=k-RAD
            np.maximum(best,nxt,out=best)
        cur=V[:,i]+best
        if trace:back.append(step)
    idx=np.argmax(cur,axis=1);rows=np.arange(b)
    if not trace:return cur[rows,idx]/np.sqrt(n),None
    p=np.empty((b,n),dtype=np.int32);p[:,-1]=idx
    for i in range(n-1,0,-1):
        idx=idx+back[i-1][rows,idx];p[:,i-1]=idx
    assert np.all((p>=0)&(p<w)) and np.all(np.abs(np.diff(p,axis=1))<=RAD)
    vals=V[rows[:,None],np.arange(n)[None,:],p]
    score=np.sum(vals,axis=1,dtype=np.float64)/np.sqrt(n)
    assert np.allclose(score,cur[rows,np.argmax(cur,axis=1)]/np.sqrt(n),rtol=1e-12,atol=1e-10)
    return score,p
def sliced(V,size,trace):
    a=[];paths=[]
    for i in range(0,len(V),size):
        z,p=batch(V[i:i+size],trace);a.append(z)
        if trace:paths.append(p)
    return np.concatenate(a),np.concatenate(paths) if trace else None
def bounds(V,t,size):
    X=np.asarray(V,dtype=np.float64);b,n,w=X.shape
    factor=np.sqrt(n);cols=np.sum(X,axis=1,dtype=np.float64)/factor
    lo=np.max(cols,axis=1);idx=np.argmax(cols,axis=1)
    hi=np.sum(np.max(X,axis=2),axis=1,dtype=np.float64)/factor
    margin=1e-10*(1+np.sum(np.max(np.abs(X),axis=2),axis=1))
    yes=lo>=t+margin;no=hi<t-margin
    assert not np.any(yes&no)
    u=np.flatnonzero(~(yes|no))
    pred=yes.copy();paths=np.full((b,n),-1,dtype=np.int32)
    paths[yes]=idx[yes,None]
    if len(u):
        s,p=sliced(X[u],size,True)
        pred[u]=s>=t;paths[u]=p
    return pred,paths,{"lower_accept":int(yes.sum()),"upper_reject":int(no.sum()),"exact_dp":len(u)}
def verify_witness(V,p,pred,t):
    ix=np.flatnonzero(pred)
    if not len(ix):return 0
    a=p[ix];assert np.all((a>=0)&(a<V.shape[2])) and np.all(np.abs(np.diff(a,axis=1))<=RAD)
    s=np.sum(V[ix[:,None],np.arange(V.shape[1])[None,:],a],axis=1,dtype=np.float64)/np.sqrt(V.shape[1])
    assert np.all(s>=t)
    return len(ix)
def selftest():
    rng=np.random.default_rng(20001)
    X=rng.normal(size=(32,16,128));X[0]=0;X[1]=1;X[2,:,0]=30;X[3,:,-1]=30
    ref=np.asarray([scalar(x)[0] for x in X]);z,p=batch(X)
    assert np.allclose(z,ref,rtol=1e-12,atol=1e-10)
    for t in (0,1,4,30,ref[4],np.nextafter(ref[4],np.inf)):
        pred,paths,counts=bounds(X,t,16)
        assert np.array_equal(pred,ref>=t)
        verify_witness(X,paths,pred,t)
    return {"random_and_edge_cases":32,"threshold_cases":6,"status":"PASS"}
def calibration(rows):
    t=time.perf_counter();V=vec.normalize([r[0] for r in rows]);y=np.array([r[1] for r in rows])
    ds=np.asarray([scalar(x)[0] for x in V]);ss=np.asarray([hold.classical(x) for x in V])
    d=prior.cutoff(ds[y==0]);s=prior.cutoff(ss[y==0])
    return {"dp":d,"straight":s,"controls":int((y==0).sum()),"dp_alerts":int(np.sum(ds[y==0]>=d)),"straight_alerts":int(np.sum(ss[y==0]>=s)),"seconds":time.perf_counter()-t}
def bench(rows,lock):
    X=[r[0] for r in rows];y=np.array([r[1] for r in rows]);sub=np.array([r[2] for r in rows]);amp=np.array([r[3] for r in rows]);hard=np.array([r[4] for r in rows])
    V0=vec.normalize(X);ref=np.asarray([scalar(v)[0] for v in V0])
    legacy=np.asarray([prior.dp_score(v) for v in V0])
    assert np.allclose(ref,legacy,rtol=1e-6,atol=1e-5)
    assert np.array_equal(ref>=lock["dp"],legacy>=lock["dp"])
    refstraight=np.asarray([hold.classical(v)>=lock["straight"] for v in V0])
    preds={};records={k:[] for k in S["modes"]}
    for rep in range(S["repeats"]):
        order=list(S["modes"]);np.random.default_rng(20030+rep).shuffle(order)
        for k in order:
            start=time.perf_counter();cpu=time.process_time();V=vec.normalize(X);cert=0;counts={}
            if k=="scalar_dp":
                a=[scalar(v) for v in V];sc=np.asarray([v[0] for v in a]);p=np.stack([v[1] for v in a]);pred=sc>=lock["dp"]
                cert=verify_witness(V,p,pred,lock["dp"])
            elif k.startswith("vector_dp_trace"):
                size=int(k.split("trace")[-1]);sc,p=sliced(V,size,True);pred=sc>=lock["dp"]
                cert=verify_witness(V,p,pred,lock["dp"])
            elif k=="vector_dp_score64":
                sc,p=sliced(V,64,False);pred=sc>=lock["dp"]
                near=np.flatnonzero(np.abs(sc-lock["dp"])<1e-9)
                for i in near:pred[i]=scalar(V[i])[0]>=lock["dp"]
            elif k=="bounded_dp_trace64":
                pred,p,counts=bounds(V,lock["dp"],64);cert=verify_witness(V,p,pred,lock["dp"])
            elif k=="vector_straight":
                vec.TH=lock["straight"];pred=vec.vector_route(V,R,S["wave_width"])[0]
            else:raise ValueError(k)
            dt=time.perf_counter()-start;dc=time.process_time()-cpu
            assert np.array_equal(pred,refstraight if k=="vector_straight" else ref>=lock["dp"]),k
            if k in preds:assert np.array_equal(pred,preds[k])
            preds[k]=pred
            records[k].append({"wall":dt,"cpu":dc,"verified_witnesses":cert,**counts})
    def met(p,m):
        return {"detected":int(np.sum(p&m&(y==1))),"positives":int(np.sum(m&(y==1))),"control_alerts":int(np.sum(p&m&(y==0))),"controls":int(np.sum(m&(y==0)))}
    return {k:{**met(p,np.ones(len(y),bool)),"subtypes":{str(i):met(p,sub==i) for i in range(4)},"amplitudes":{str(a):met(p,amp==a) for a in S["amplitudes"]},"hard":met(p,hard),"wall":statistics.median(r["wall"] for r in records[k]),"cpu":statistics.median(r["cpu"] for r in records[k]),"verified_witnesses":records[k][0]["verified_witnesses"],"bounds":{a:records[k][0].get(a,0) for a in ("lower_accept","upper_reject","exact_dp")},"repeats":records[k]} for k,p in preds.items()}
def main():
    checks=selftest();m=json.loads((OUT/"manifest.json").read_text())
    assert len(m["selected"])==6 and m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest()
    assert len({q["target"] for q in m["selected"]})==6
    assert all(q["target"] not in hold.S["exclude"] for q in m["selected"])
    st=time.perf_counter();cal=[]
    for q in m["selected"][:2]:cal+=prior.gather(q)
    load=time.perf_counter()-st;lock=calibration(cal)
    lock.update({"targets":[q["target"] for q in m["selected"][:2]],"load_seconds":load,"manifest_sha256":hashlib.sha256((OUT/"manifest.json").read_bytes()).hexdigest()})
    save("lock.json",lock);emit("R20A_LOCK",lock)
    per={};load=0;seen={r[5] for r in cal}
    for q in m["selected"][2:]:
        st=time.perf_counter();rows=prior.gather(q);load+=time.perf_counter()-st
        assert not seen&{r[5] for r in rows};seen|={r[5] for r in rows}
        per[q["target"]]=bench(rows,lock);save("partial.json",per)
        emit("R20A_TARGET",{"target":q["target"],"metrics":per[q["target"]]})
    keys=("detected","positives","control_alerts","controls","wall","cpu","verified_witnesses")
    total={k:{a:sum(x[k][a] for x in per.values()) for a in keys} for k in S["modes"]}
    for k in total:
        for field,values in (("subtypes",[str(i) for i in range(4)]),("amplitudes",[str(a) for a in S["amplitudes"]])):
            total[k][field]={i:{a:sum(x[k][field][i][a] for x in per.values()) for a in keys[:4]} for i in values}
        total[k]["bounds"]={a:sum(x[k]["bounds"][a] for x in per.values()) for a in ("lower_accept","upper_reject","exact_dp")}
    d={"targets":list(per),"samples":1024,"dp_decision_mismatches":0,"legacy_decision_mismatches":0,"straight_decision_mismatches":0,"selftests":checks,"lock":lock,"modes":total,"per_target":per,"evaluation_load_seconds":load,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"scope":"Finite bounded-step arbitrary track. Float64 numerical replay, not exact real arithmetic. Weak synthetic injections in fresh OFF observations; score-only mode has no path certificate. No source inference.","production_adoption":False}
    save("result.json",d);emit("R20A_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="scan":main()
    else:raise ValueError("phase")
