#!/usr/bin/env python3
"""R20D: certified two-pass DP, matched-control joint and curvature-regularized DP."""
import hashlib,json,pathlib,statistics,time,itertools
import numpy as np
import r20c_verified_tile_cache as cache
import r20a_vector_dp as dp
import r19e_integrated_exam as prior
import r17c_frozen_holdout as hold
import juxtapose_undina_r18b2_vector as vec
import juxtapose_undina_r18bc as old
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r20d_weak_joint_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r20d");OUT.mkdir(parents=True,exist_ok=True)
cache.S=S;cache.OUT=OUT;cache.CACHE=OUT/"cache";cache.CACHE.mkdir(exist_ok=True)
cache.COLD=OUT/"cold"
R=tuple(dp.S["rank"]);RAD=dp.RAD
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def emit(k,d):print(k,json.dumps(d),flush=True)
def regular(V,penalty,trace=True):
    # State (frequency, incoming integer slope); exact max-sum finite recurrence.
    V=np.asarray(V,dtype=np.float64);b,n,w=V.shape
    slopes=np.arange(-RAD,RAD+1);K=len(slopes);cols=np.arange(w)
    prev=np.full((b,K,w),-np.inf);backs=[]
    for k,s in enumerate(slopes):
        ix=cols-s;valid=(ix>=0)&(ix<w)
        prev[:,k,valid]=V[:,0,ix[valid]]+V[:,1,valid]
    for i in range(2,n):
        cur=np.full_like(prev,-np.inf)
        back=np.zeros((b,K,w),dtype=np.int8) if trace else None
        for k,s in enumerate(slopes):
            ix=cols-s;valid=(ix>=0)&(ix<w);dest=np.flatnonzero(valid);src=ix[valid]
            z=prev[:,:,src]-penalty*np.abs(slopes-s)[None,:,None]
            cur[:,k,dest]=np.max(z,axis=1)+V[:,i,dest]
            if trace:back[:,k,dest]=np.argmax(z,axis=1)
        prev=cur
        if trace:backs.append(back)
    flat=prev.reshape(b,-1);j=np.argmax(flat,axis=1);rows=np.arange(b)
    score=flat[rows,j]/np.sqrt(n)
    if not trace:return score,None
    k=j//w;idx=j%w;p=np.zeros((b,n),dtype=np.int32);p[:,-1]=idx
    for i in range(n-1,1,-1):
        oldk=k.copy();k=backs[i-2][rows,oldk,idx]
        idx=idx-slopes[oldk];p[:,i-1]=idx
    p[:,0]=idx-slopes[k]
    val=regular_witness(V,p,penalty)
    assert np.allclose(val,score,atol=1e-10,rtol=1e-12)
    return val,p
def regular_witness(V,p,penalty):
    assert np.all((p>=0)&(p<V.shape[2])) and np.all(np.abs(np.diff(p,axis=1))<=RAD)
    v=V[np.arange(len(V))[:,None],np.arange(V.shape[1])[None,:],p]
    return (np.sum(v,axis=1,dtype=np.float64)-penalty*np.sum(np.abs(np.diff(p,n=2,axis=1)),axis=1))/np.sqrt(V.shape[1])
def reg_sliced(V,penalty,trace):
    a=[];ps=[]
    for i in range(0,len(V),32):
        z,p=regular(V[i:i+32],penalty,trace);a.append(z)
        if trace:ps.append(p)
    return np.concatenate(a),np.concatenate(ps) if trace else None
def selftest():
    rng=np.random.default_rng(20204);V=rng.normal(size=(8,4,6));V[0]=0
    for pen in S["regularized_penalties"]:
        a,p=regular(V,pen);b,_=regular(V,pen,False)
        assert np.allclose(a,b,atol=1e-10)
        paths=np.asarray([q for q in itertools.product(range(6),repeat=4) if max(abs(x-y) for x,y in zip(q,q[1:]))<=RAD])
        penv=pen*np.sum(np.abs(np.diff(paths,n=2,axis=1)),axis=1)
        for i,v in enumerate(V):
            truth=np.max((v[np.arange(4)[None,:],paths].sum(1)-penv)/2)
            assert abs(a[i]-truth)<1e-10
    return "PASS_24_BRUTE_FORCE_CASES"
def met(p,y,mask):
    return {"detected":int(np.sum(p&(y==1)&mask)),"positives":int(np.sum((y==1)&mask)),"control_alerts":int(np.sum(p&(y==0)&mask)),"controls":int(np.sum((y==0)&mask))}
def calibration(rows):
    t=time.perf_counter();V=vec.normalize([r[0] for r in rows]);y=np.array([r[1] for r in rows])
    ss=np.asarray([hold.classical(v) for v in V]);ds,_=dp.sliced(V,256,True)
    st=prior.cutoff(ss[y==0]);dt=prior.cutoff(ds[y==0]);jt=prior.cutoff(np.maximum(ss[y==0]/st,ds[y==0]/dt))
    trials=[];reg={}
    for pen in S["regularized_penalties"]:
        a,_=reg_sliced(V,pen,True);th=prior.cutoff(a[y==0])
        trials.append({"penalty":pen,"threshold":th,**met(a>=th,y,np.ones(len(y),bool))});reg[pen]=a
    best=min(trials,key=lambda a:(-a["detected"],a["control_alerts"],a["penalty"]))
    pen=best["penalty"];rt=best["threshold"];rr=reg[pen]
    jrt=prior.cutoff(np.maximum(ss[y==0]/st,rr[y==0]/rt))
    return {"straight":st,"dp":dt,"joint_straight":st*jt,"joint_dp":dt*jt,"penalty":pen,"regular":rt,"joint_reg_straight":st*jrt,"joint_reg":rt*jrt,"regular_trials":trials,"calibration_seconds":time.perf_counter()-t,"calibration_controls":int((y==0).sum()),"calibration_alerts":{"straight":int(np.sum(ss[y==0]>=st)),"dp":int(np.sum(ds[y==0]>=dt)),"joint":int(np.sum((ss[y==0]>=st*jt)|(ds[y==0]>=dt*jt))),"regular":int(np.sum(rr[y==0]>=rt)),"joint_regular":int(np.sum((ss[y==0]>=st*jrt)|(rr[y==0]>=rt*jrt)))}}
def exact_two(V,t):
    a,_=dp.sliced(V,64,False);pos=a>=t-1e-9
    p=np.zeros((len(V),V.shape[1]),dtype=np.int32);pred=np.zeros(len(V),bool)
    ix=np.flatnonzero(pos)
    if len(ix):
        s,pp=dp.sliced(V[ix],256,True);pred[ix]=s>=t;p[ix]=pp
    return pred,dp.verify_witness(V,p,pred,t),len(ix)
def straight(V,t):
    vec.TH=t;return vec.vector_route(V,R,2)[0]
def bench(rows,lock):
    X=[r[0] for r in rows];y=np.array([r[1] for r in rows]);sub=np.array([r[2] for r in rows]);amp=np.array([r[3] for r in rows])
    V=vec.normalize(X);reference=np.asarray([dp.scalar(v)[0] for v in V])>=lock["dp"]
    sr=np.asarray([hold.classical(v)>=lock["straight"] for v in V])
    # Direct exact-score references for every joint decision, outside timed paths.
    ss=np.asarray([hold.classical(v) for v in V])
    ds,_=dp.sliced(V,256,True);rs,_=reg_sliced(V,lock["penalty"],True)
    ref={"dp_trace":reference,"dp_two_pass":reference,"straight_scalar":sr,"joint":(ss>=lock["joint_straight"])|(ds>=lock["joint_dp"]),"regularized_dp":rs>=lock["regular"],"joint_regularized":(ss>=lock["joint_reg_straight"])|(rs>=lock["joint_reg"])}
    out={k:[] for k in S["modes"]};preds={}
    for rep in range(S["repeats"]):
        modes=list(S["modes"]);np.random.default_rng(20300+rep).shuffle(modes)
        for k in modes:
            t=time.perf_counter();cpu=time.process_time();V=vec.normalize(X);cert=0;traces=0
            if k=="dp_trace":
                s,p=dp.sliced(V,256,True);pred=s>=lock["dp"];cert=dp.verify_witness(V,p,pred,lock["dp"]);traces=len(V)
            elif k=="dp_two_pass":pred,cert,traces=exact_two(V,lock["dp"])
            elif k=="straight_scalar":pred=np.asarray([hold.classical(v)>=lock["straight"] for v in V])
            elif k=="joint":
                pred=straight(V,lock["joint_straight"]);ix=np.flatnonzero(~pred)
                if len(ix):
                    h,c,tr=exact_two(V[ix],lock["joint_dp"]);pred[ix]=h;cert=c;traces=tr
            elif k=="regularized_dp":
                s,p=reg_sliced(V,lock["penalty"],True);pred=s>=lock["regular"];traces=len(V)
                cert=int(pred.sum());assert np.all(regular_witness(V[pred],p[pred],lock["penalty"])>=lock["regular"])
            elif k=="joint_regularized":
                pred=straight(V,lock["joint_reg_straight"]);ix=np.flatnonzero(~pred);traces=len(ix)
                if len(ix):
                    s,p=reg_sliced(V[ix],lock["penalty"],True);h=s>=lock["joint_reg"];pred[ix]=h;cert=int(h.sum());assert np.all(regular_witness(V[ix][h],p[h],lock["penalty"])>=lock["joint_reg"])
            else:raise ValueError(k)
            dt=time.perf_counter()-t;dc=time.process_time()-cpu
            assert np.array_equal(pred,ref[k]),k
            preds[k]=pred;out[k].append({"wall":dt,"cpu":dc,"dp_witnesses":cert,"traces":traces})
    return {k:{**met(p,y,np.ones(len(y),bool)),"subtypes":{str(i):met(p,y,sub==i) for i in range(4)},"amplitudes":{str(a):met(p,y,amp==a) for a in dp.S["amplitudes"]},"wall":statistics.median(r["wall"] for r in out[k]),"cpu":statistics.median(r["cpu"] for r in out[k]),"dp_witnesses":out[k][0]["dp_witnesses"],"traces":out[k][0]["traces"],"repeats":out[k]} for k,p in preds.items()}
def main():
    check=selftest();sources,by,manifest,restore=cache.restore()
    eng.get_tiles=lambda f:cache.cached_get(f,sources,by);old.OUT=cache.CACHE
    # Reuse R20A frozen seed, amplitude and target manifest, not R20D protocol for injections.
    assert manifest["protocol_sha256"]==hashlib.sha256(dp.P.read_bytes()).hexdigest()
    t=time.perf_counter();cal=[]
    for q in manifest["selected"][:2]:cal+=prior.gather(q)
    cal_load=time.perf_counter()-t;lock=calibration(cal);lock["calibration_patch_seconds"]=cal_load
    save("lock.json",lock);emit("R20D_LOCK",lock)
    rows={};load=0
    for q in manifest["selected"][2:]:
        t=time.perf_counter();data=prior.gather(q);load+=time.perf_counter()-t
        rows[q["target"]]=bench(data,lock);save("partial.json",rows)
    keys=("detected","positives","control_alerts","controls","wall","cpu","dp_witnesses","traces")
    total={k:{a:sum(v[k][a] for v in rows.values()) for a in keys} for k in S["modes"]}
    for k in total:
        for fld,vals in (("subtypes",[str(i) for i in range(4)]),("amplitudes",[str(a) for a in dp.S["amplitudes"]])):
            total[k][fld]={i:{a:sum(v[k][fld][i][a] for v in rows.values()) for a in keys[:4]} for i in vals}
    # Exact recreation of the previously published R20A control result is mandatory.
    assert (total["dp_trace"]["detected"],total["dp_trace"]["control_alerts"])==(160,9)
    assert (total["straight_scalar"]["detected"],total["straight_scalar"]["control_alerts"])==(119,9)
    assert total["dp_two_pass"]["dp_witnesses"]==169
    d={"status":"DEVELOPMENT_ONLY","samples":1024,"exact_mismatches":0,"selftest":check,"lock":lock,"modes":total,"per_target":rows,"cache_restore":restore,"evaluation_patch_seconds":load,"scope":"Previously viewed R20A cohort. Thresholds jointly calibrated; realized alert counts need not match. Regularity penalty defines a new detector, not equivalent to arbitrary DP. No astronomical origin inference or production adoption.","production_adoption":False}
    save("result.json",d);emit("R20D_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":main()
