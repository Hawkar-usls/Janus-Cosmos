#!/usr/bin/env python3
"""R25A: exact complete boundary-state messages and checkpoint traceback.
Adaptation of Fundamentum factorization/boundary-join concepts, not a SAT carrier replay.
"""
import hashlib,json,pathlib,statistics,time,sys
import numpy as np
import r24b_cousteau as cb
import r24_odontocrystal as od
import r20c_verified_tile_cache as cache
import juxtapose_undina_r18bc as old
import juxtapose_undina_r18b2_vector as vec
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r25a_boundary_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r25a");OUT.mkdir(parents=True,exist_ok=True)
RAD=od.RAD;PEN=od.PEN;COST=od.COST;K=2*RAD+1;slopes=np.arange(-RAD,RAD+1)
def emit(k,d):print(k,json.dumps(d),flush=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def advance(prev,Z,t,G,pen,cost):
    B,_,w=Z.shape
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

    return cur,parent,skip
def checkpoint(V,G,block,pen=PEN,cost=COST,clip=None):
    V=np.asarray(V,dtype=np.float64);B,n,w=V.shape
    assert n>=2 and 0<=G<n and np.isfinite(V).all()
    Z=V if clip is None else np.minimum(V,clip)
    prev=np.full((B,G+1,K,w),-np.inf);prev[:,0]=Z[:,0,None,:]
    if G:prev[:,1]=-cost
    cuts={0:prev.copy()}
    for t in range(1,n):
        prev,_,_=advance(prev,Z,t,G,pen,cost)
        if t%block==0 and t<n-1:cuts[t]=prev.copy()
    j=np.argmax(prev.reshape(B,-1),axis=1);r=np.arange(B)
    optimum=prev.reshape(B,-1)[r,j]/np.sqrt(n)
    g=j//(K*w);k=j//w%K;x=j%w
    p=np.empty((B,n),dtype=np.int32);m=np.zeros((B,n),bool)
    end=n-1;peak=0;recomputed=0
    for start in reversed(sorted(cuts)):
        state=cuts[start];backs=[]
        for t in range(start+1,end+1):
            state,pa,sk=advance(state,Z,t,G,pen,cost);backs.append((pa,sk));recomputed+=1
        peak=max(peak,sum(a.nbytes+b.nbytes for a,b in backs))
        for t in range(end,start,-1):
            p[:,t]=x;pa,sk=backs[t-start-1];q=sk[r,g,k,x];m[:,t]=q
            nk=pa[r,g,k,x];x=x-slopes[k];g=g-q.astype(int);k=nk
        end=start
    p[:,0]=x;m[:,0]=g==1
    assert np.all((g==0)|(g==1)) and np.all(m.sum(1)<=G) and recomputed==n-1
    score=od.witness(V,p,m,pen,cost,clip);assert np.allclose(score,optimum,atol=1e-10,rtol=1e-12)
    return score,p,m,{"retained_trace_and_checkpoint_bytes":sum(a.nbytes for a in cuts.values())+peak,"full_trace_bytes":(n-1)*B*(G+1)*K*w*2,"recomputed_steps":recomputed}
def batched(V,G,block):
    a=[];p=[];m=[];mem=[];full=[]
    for i in range(0,len(V),S["batch"]):
        v=V[i:i+S["batch"]]
        if block:
            s,q,r,c=checkpoint(v,G,block);mem.append(c["retained_trace_and_checkpoint_bytes"]);full.append(c["full_trace_bytes"])
        else:
            s,q,r=od.grow(v,G);f=(v.shape[1]-1)*len(v)*(G+1)*K*v.shape[2]*2;mem.append(f);full.append(f)
        a.append(s);p.append(q);m.append(r)
    return np.concatenate(a),np.concatenate(p),np.concatenate(m),max(mem),max(full)
def selftest():
    rng=np.random.default_rng(25000);cases=0
    for n in (2,7,16,17,33):
        v=rng.normal(size=(6,n,21));v[0]=0;v[1,:,0]=5
        for g in (0,min(4,n-1)):
            a,p,m=od.grow(v,g)
            for b in S["block_sizes"]:
                z,q,r,c=checkpoint(v,g,b)
                assert np.allclose(a,z,atol=1e-10,rtol=1e-12)
                assert np.array_equal(p,q) and np.array_equal(m,r)
                cases+=len(v)
    return {"boundary_equivalence_cases":cases,"status":"PASS"}
def cutoff(a):
    z=np.sort(a)[::-1];return float(np.nextafter(z[int(.05*len(z))],np.inf))
def replay():
    cache.S=S;cache.OUT=OUT;cache.CACHE=OUT/"cache";cache.CACHE.mkdir(exist_ok=True)
    source,by,m,restore=cache.restore()
    assert m["protocol_sha256"]==hashlib.sha256(cb.P.read_bytes()).hexdigest()
    eng.get_tiles=lambda f:cache.cached_get(f,source,by);old.OUT=cache.CACHE
    od.S=cb.S;old.S={**old.S,"seed":cb.S["seed"]}
    # Recreate R24B data with its original target-specific injection seed.
    rows=[]
    for q in m["selected"][2:]:rows+=od.gather(q)
    V=vec.normalize([r[0] for r in rows]).astype(np.float64)
    a,p,m,_,_=batched(V,0,0)
    z,q,r,mem,full=batched(V,0,8)
    assert np.allclose(a,z,atol=1e-10,rtol=1e-12)
    assert np.array_equal(p,q) and np.array_equal(m,r)
    t=5.9291978865861905;pred=z>=t;assert np.array_equal(pred,a>=t)
    met=od.measure(pred,q,r,rows)
    assert (met["detected"],met["localized"],met["control_alerts"])==(141,110,22)
    return {"samples":len(rows),"mismatches":0,"metrics":met,"stored_bytes_checkpoint":mem,"stored_bytes_full":full,"cache_restore":restore}
def generate(n,pairs,seed):
    rng=np.random.default_rng(seed+n);w=S["width"]
    plan=[(a,r,g) for a in S["amplitudes"] for r in S["noise"] for g in S["gap_fraction"]]*(pairs//32);rng.shuffle(plan)
    V=[];Y=[];truth=[];active=[];meta=[]
    for amp,noise,gap in plan:
        x=rng.normal(size=(n,w))
        if noise=="ar1":
            for i in range(1,n):x[i]=.65*x[i-1]+np.sqrt(1-.65**2)*x[i]
        elif noise=="persistent_rfi":x[:,int(rng.integers(8,w-8))]+=float(rng.uniform(2,8))
        elif noise=="impulsive":
            for _ in range(n//2):x[int(rng.integers(n)),int(rng.integers(w))]+=float(rng.uniform(4,15))
        tt=np.arange(n);center=float(rng.uniform(36,92));ph=float(rng.uniform(0,2*np.pi))
        path=np.rint(center+rng.uniform(-.25,.25)*(tt-(n-1)/2)+rng.uniform(4,min(9,2*n/(2*np.pi)))*np.sin(2*np.pi*tt/n+ph)).astype(int)
        assert path.min()>=0 and path.max()<w and np.max(np.abs(np.diff(path)))<=RAD
        ok=np.ones(n,bool)
        if gap:ok[rng.choice(n,int(n*gap),replace=False)]=False
        p=x.copy();scale=1.4826*np.median(np.abs(x-np.median(x,axis=1)[:,None]),axis=1)
        p[tt[ok],path[ok]]+=amp*scale[ok]
        V.extend([x,p]);Y.extend([0,1]);truth.extend([path,path]);active.extend([ok,ok]);meta.extend([(amp,noise,gap)]*2)
    return vec.normalize(V).astype(np.float64),np.array(Y),np.stack(truth),np.stack(active),meta
def independent(V):
    a=[];p=[]
    for i in range(0,V.shape[1],16):
        z,q,_,_,_=batched(V[:,i:i+16],0,0);a.append(z*4);p.append(q)
    return np.sum(a,axis=0)/np.sqrt(V.shape[1]),np.concatenate(p,axis=1)
def synth(n):
    cal=generate(n,S["cal_pairs"],S["cal_seed"]);test=generate(n,S["test_pairs"],S["test_seed"])
    V,y,truth,active,meta=test;X,Y,*_=cal
    th={"continuous":cutoff(batched(X,0,0)[0][Y==0]),"independent_blocks":cutoff(independent(X)[0][Y==0])}
    def metrics(s,p,threshold):
        z=s>=threshold;loc=((np.abs(p-truth)<=2)&active).sum(1)>=.75*active.sum(1)
        def m(mask):return {"positive_alerts":int(np.sum(z&mask&(y==1))),"localized":int(np.sum(z&loc&mask&(y==1))),"control_alerts":int(np.sum(z&mask&(y==0))),"positives":int(np.sum(mask&(y==1))),"controls":int(np.sum(mask&(y==0)))}
        return {**m(np.ones(len(y),bool)),"amplitudes":{str(a):m(np.array([q[0]==a for q in meta])) for a in S["amplitudes"]},"noise":{a:m(np.array([q[1]==a for q in meta])) for a in S["noise"]}}
    ref=batched(V,0,0);out={}
    for k in ["full","checkpoint8","checkpoint16","checkpoint32","independent_blocks"]:
        ts=[];mem=None
        for rep in range(S["repeats"]):
            t=time.perf_counter()
            if k=="independent_blocks":s,p=independent(V)
            else:
                block=0 if k=="full" else int(k[10:]);s,p,m,mem,full=batched(V,0,block)
            ts.append(time.perf_counter()-t)
            if k!="independent_blocks":
                assert np.allclose(s,ref[0],atol=1e-10,rtol=1e-12) and np.array_equal(p,ref[1]) and np.array_equal(m,ref[2])
                assert np.array_equal(s>=th["continuous"],ref[0]>=th["continuous"])
        out[k]={**metrics(s,p,th["independent_blocks" if k=="independent_blocks" else "continuous"]),"wall_seconds":statistics.median(ts),"stored_trace_and_checkpoint_bytes":mem}
    return {"length":n,"thresholds":th,"modes":out}
def main():
    check=selftest();emit("R25A_SELFTEST",check)
    a=replay();save("archive.json",a);emit("R25A_ARCHIVE",a)
    res={}
    for n in S["lengths"]:
        res[str(n)]=synth(n);save("synthetic_partial.json",res);emit("R25A_LENGTH",res[str(n)])
    d={"status":"RESEARCH_ONLY","selftests":check,"archive":a,"synthetic":res,"scope":"Pinned replay, plus independently generated synthetic lengths, not new long telescope scans. Retained traceback/checkpoint bytes are not process RSS or total allocation; recomputation cost included. Independent-block path may jump at joins. Exact checkpoint algorithm gives same detector, sensitivity changes only due longer inputs.","production_adoption":False}
    save("result.json",d);emit("R25A_SUMMARY",d)
if __name__=="__main__":main()
