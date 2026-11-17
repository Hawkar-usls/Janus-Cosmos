#!/usr/bin/env python3
"""R21BC: timestamp-aligned noncoherent integration across two actual scans."""
import hashlib,json,pathlib,sys,time,statistics,math
import numpy as np
import r21a_weak_template as prior
import r17c_frozen_holdout as hold
import r9_expanded_search as eng
import juxtapose_undina_r18b2_vector as vec
P=pathlib.Path("research/seti/r21bc_multiscan_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r21bc");OUT.mkdir(parents=True,exist_ok=True)
hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))}
hold.OUT=OUT;hold.P=P;eng.OUT=OUT
GRID=np.linspace(-1,1,81);W=S["patch_width"]
def emit(k,d):print(k,json.dumps(d),flush=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def gather(q):
    fs=[x for x in q["files"] if x["role"].startswith("OFF")];assert len(fs)==2
    h,g=[f["header"] for f in fs]
    assert h["src_raj"]==g["src_raj"] and h["src_dej"]==g["src_dej"] and h["ts"]==g["ts"] and h["df"]==g["df"]
    delta=(g["mjd"]-h["mjd"])*86400/h["ts"];assert 16<delta<100
    ts=[np.arange(16,dtype=float),np.arange(16,dtype=float)+delta]
    tiles=[eng.get_tiles(f) for f in fs];assert len(tiles[0])==len(tiles[1])==4
    rng=np.random.default_rng(int.from_bytes(hashlib.sha256((S["seed"]+"|"+q["target"]).encode()).digest()[:8],"big"))
    plan=[(a,p,g) for a in S["amplitudes"] for p in S["presence"] for g in S["gap_levels"]]*8;rng.shuffle(plan)
    rows=[];challenges=[];k=0
    for a,b in zip(*tiles):
        assert a["center_mhz"]==b["center_mhz"] and a["df_hz"]==b["df_hz"] and a["channel0_mhz"]==b["channel0_mhz"] and a["shape"]==b["shape"]
        x=[np.load(OUT/t["artifact"],allow_pickle=False) for t in (a,b)]
        for v,t in zip(x,(a,b)):
            assert hashlib.sha256(v.tobytes()).hexdigest()==t["slice_sha256"] and v.shape[0]==16
        hard=np.argsort(-np.maximum(np.median(vec.normalize([x[0]])[0],axis=0),np.median(vec.normalize([x[1]])[0],axis=0)),kind="stable")
        hc=[]
        for j in hard:
            if W//2<=j<x[0].shape[1]-W//2 and all(abs(j-i)>=W for i in hc):hc.append(int(j))
            if len(hc)==8:break
        assert len(hc)==8
        for trial in range(32):
            start=hc[trial-24]-W//2 if trial>=24 else int(rng.integers(0,x[0].shape[1]-W+1))
            raw=[z[:,start:start+W].copy() for z in x]
            amp,presence,gap=plan[k];k+=1
            slope=float(rng.choice(S["injection_slopes"]));ds=[np.rint(slope*t).astype(int) for t in ts]
            lo=8+max(0,-min(int(d.min()) for d in ds));hi=W-8-max(0,max(int(d.max()) for d in ds))
            center=int(rng.integers(lo,hi));paths=np.stack([center+d for d in ds])
            scales=[1.4826*np.median(np.abs(v-np.median(v,axis=1)[:,None]),axis=1) for v in raw];assert all(np.all(v>0) for v in scales)
            active=np.ones((2,16),bool)
            if gap:
                for i in range(2):active[i,rng.choice(16,size=int(16*gap),replace=False)]=False
            if presence=="first_only":active[1]=False
            positive=[v.copy() for v in raw]
            for i in range(2):
                tt=np.flatnonzero(active[i]);positive[i][tt,paths[i,tt]]+=amp*scales[i][tt]
            meta={"amp":amp,"presence":presence,"gap":gap,"hard":trial>=24,"paths":paths.tolist(),"slope":slope,"source":[a["slice_sha256"],b["slice_sha256"]],"start":start}
            rows+=[(np.stack(raw),0,meta),(np.stack(positive),1,meta)]
            if trial<4:
                # Repeating stationary synthetic interference and incompatible moving tracks.
                for kind in ("stationary_rfi","discordant"):
                    z=[v.copy() for v in raw]
                    if kind=="stationary_rfi":p=[np.full(16,128),np.full(16,128)]
                    else:p=[np.rint(48+.7*np.arange(16)).astype(int),np.rint(208-.7*np.arange(16)).astype(int)]
                    for i in range(2):z[i][np.arange(16),p[i]]+=2*scales[i]
                    challenges.append((np.stack(z),kind))
    assert len(rows)==256 and len(challenges)==32
    return rows,challenges,delta
def maps(V,offset):
    # Each hypothesis (slope, starting global frequency bin) shares the same key.
    V=np.asarray(V,dtype=np.float64);B,n,w=V.shape
    d=np.rint(GRID[:,None]*(np.arange(n)[None,:]+offset)).astype(int)
    a=np.zeros((B,len(GRID),w));valid=np.ones((len(GRID),w),bool);cols=np.arange(w)
    for k in range(n):
        ix=cols[None,:]+d[:,k,None];valid&=(ix>=0)&(ix<w)
        a+=V[:,k,np.clip(ix,0,w-1)]
    return np.where(valid[None],a,-np.inf)
def maximum(a):
    flat=a.reshape(len(a),-1);i=np.argmax(flat,axis=1)
    return flat[np.arange(len(a)),i],i
def verify(V,delta,idx,scores,mode):
    B,_,n,w=V.shape;g=idx//w;c=idx%w
    p=np.stack([c[:,None]+np.rint(GRID[g,None]*(np.arange(n)[None,:]+off)).astype(int) for off in (0,delta)],axis=1)
    assert np.all((p>=0)&(p<w))
    v=V[np.arange(B)[:,None,None],np.arange(2)[None,:,None],np.arange(n)[None,None,:],p]
    r=v.sum((1,2),dtype=np.float64)/np.sqrt(2*n)
    assert np.allclose(r,scores,rtol=1e-12,atol=1e-10),mode
    return p
def join_semijoin(a,b,t):
    A=a.reshape(len(a),-1);B=b.reshape(len(b),-1);limit=t*np.sqrt(32)
    # Necessary conditions only. No individual minimum significance requirement.
    am=np.max(A,axis=1);bm=np.max(B,axis=1)
    margin=1e-10*(1+np.maximum(np.abs(am),np.abs(bm))+abs(limit))
    support=(A+bm[:,None]>=limit-margin[:,None])&(B+am[:,None]>=limit-margin[:,None])
    p=np.zeros(len(a),bool);indices=np.zeros(len(a),dtype=int);scores=np.full(len(a),-np.inf)
    for i in range(len(a)):
        k=np.flatnonzero(support[i])
        if len(k):
            v=A[i,k]+B[i,k];j=int(np.argmax(v));indices[i]=k[j];scores[i]=v[j]/np.sqrt(32);p[i]=scores[i]>=t
    return p,scores,indices,int(support.sum()),int(support.size)
def normalized(rows):
    x=np.stack([r[0] for r in rows]);b,s,n,w=x.shape
    return vec.normalize(list(x.reshape(b*s,n,w))).reshape(b,s,n,w).astype(np.float64)
def calculate(V,delta,k,t=None):
    a=maps(V[:,0],0);x,i=maximum(a);n=V.shape[2];w=V.shape[3]
    if k=="single1":
        path=i[:,None]%w+np.rint(GRID[i//w,None]*np.arange(n)[None,:]).astype(int)
        return x/np.sqrt(n),path[:,None],0,0
    b=maps(V[:,1],delta if k in ("timestamp_joint","timestamp_joint_semijoin") else 0)
    y,j=maximum(b)
    if k=="single_best":
        sel=x>=y;ind=np.where(sel,i,j);p=ind[:,None]%w+np.rint(GRID[ind//w,None]*np.arange(n)[None,:]).astype(int)
        # Record scan chosen for independent location check.
        return np.maximum(x,y)/np.sqrt(n),np.concatenate((p[:,None],np.broadcast_to(np.where(sel,0,1)[:,None,None],(len(p),1,n))),axis=1),0,0
    if k=="unconstrained_sum":return (x+y)/np.sqrt(2*n),None,0,0
    if k=="timestamp_joint_semijoin":
        pred,sc,idx,support,total=join_semijoin(a,b,t)
        p=np.full((len(V),2,n),-1,dtype=int);ii=np.flatnonzero(pred)
        if len(ii):p[ii]=verify(V[ii],delta,idx[ii],sc[ii],k)
        return sc,p,support,total
    z,idx=maximum(a+b);sc=z/np.sqrt(2*n)
    p=verify(V,delta if k=="timestamp_joint" else 0,idx,sc,k)
    return sc,p,0,0
def batched(V,delta,k,t=None):
    a=[];p=[];sup=total=0
    for j in range(0,len(V),16):
        z,q,s,n=calculate(V[j:j+16],delta,k,t);a.append(z)
        if q is not None:p.append(q)
        sup+=s;total+=n
    return np.concatenate(a),np.concatenate(p) if p else None,sup,total
BASE=S["modes"][:-1]
def cutoff(a):
    z=sorted(map(float,a),reverse=True);return float(np.nextafter(z[math.floor(len(z)*S["calibration_control_fraction"])],np.inf))
def calibrate(groups):
    t=time.perf_counter();sc={k:[] for k in BASE}
    for rows,_,delta in groups:
        V=normalized(rows);y=np.array([r[1] for r in rows])
        for k in BASE:sc[k].extend(batched(V,delta,k)[0][y==0])
    th={k:cutoff(v) for k,v in sc.items()};th["timestamp_joint_semijoin"]=th["timestamp_joint"]
    return {"thresholds":th,"controls":len(sc["single1"]),"calibration_alerts":{k:int(np.sum(np.asarray(v)>=th[k])) for k,v in sc.items()},"seconds":time.perf_counter()-t}
def localization(k,p,y,meta,pred):
    if p is None:return None
    hit=0;active=0
    for i in np.flatnonzero(pred&(y==1)):
        m=meta[i];truth=np.asarray(m["paths"])
        if k=="single1":ok=np.mean(np.abs(p[i,0]-truth[0])<=2)>=.75
        elif k=="single_best":
            s=int(p[i,1,0]);ok=(s==0 or m["presence"]=="both") and np.mean(np.abs(p[i,0]-truth[s])<=2)>=.75
        else:
            n=2 if m["presence"]=="both" else 1
            ok=all(np.mean(np.abs(p[i,j]-truth[j])<=2)>=.75 for j in range(n))
        hit+=int(ok);active+=1
    return {"localized_positive":hit,"detected_positive":active}
def metrics(p,y,m):return {"detected":int(np.sum(p&(y==1)&m)),"positives":int(np.sum((y==1)&m)),"control_alerts":int(np.sum(p&(y==0)&m)),"controls":int(np.sum((y==0)&m))}
def bench(rows,challenges,delta,lock):
    y=np.array([r[1] for r in rows]);meta=[r[2] for r in rows];amp=np.asarray([m["amp"] for m in meta]);presence=np.asarray([m["presence"] for m in meta]);gap=np.asarray([m["gap"] for m in meta])
    th=lock["thresholds"];first={};paths={};logs={k:[] for k in S["modes"]}
    for rep in range(S["repeats"]):
        modes=list(S["modes"]);np.random.default_rng(21200+rep).shuffle(modes)
        for k in modes:
            t=time.perf_counter();cpu=time.process_time();V=normalized(rows);s,p,sup,total=batched(V,delta,k,th[k]);z=s>=th[k]
            dt=time.perf_counter()-t;dc=time.process_time()-cpu
            if k in first:assert np.array_equal(first[k],z)
            first[k]=z;paths[k]=p;logs[k].append({"wall":dt,"cpu":dc,"support":sup,"hypotheses":total})
    assert np.array_equal(first["timestamp_joint"],first["timestamp_joint_semijoin"])
    # Same shared-hypothesis witness selected for each positive.
    ix=first["timestamp_joint"];assert np.array_equal(paths["timestamp_joint"][ix],paths["timestamp_joint_semijoin"][ix])
    CV=normalized([(r,0,{}) for r,_ in challenges]);labels=np.asarray([k for _,k in challenges]);out={}
    for k,p in first.items():
        cs=batched(CV,delta,k,th[k])[0]>=th[k]
        out[k]={**metrics(p,y,np.ones(len(y),bool)),"amplitudes":{str(a):metrics(p,y,amp==a) for a in S["amplitudes"]},"presence":{a:metrics(p,y,presence==a) for a in S["presence"]},"gaps":{str(a):metrics(p,y,gap==a) for a in S["gap_levels"]},"amplitude_presence":{str(a)+"_"+b:metrics(p,y,(amp==a)&(presence==b)) for a in S["amplitudes"] for b in S["presence"]},"localization":localization(k,paths[k],y,meta,p),"challenges":{a:int(np.sum(cs&(labels==a))) for a in ("stationary_rfi","discordant")},"wall":statistics.median(r["wall"] for r in logs[k]),"cpu":statistics.median(r["cpu"] for r in logs[k]),"support":logs[k][0]["support"],"hypotheses":logs[k][0]["hypotheses"],"repeats":logs[k]}
    return out
def selftest():
    rng=np.random.default_rng(21203);V=rng.normal(size=(8,2,16,256));V[0]=0;V[1,:, :,128]=5
    a=maps(V[:,0],0);b=maps(V[:,1],36.1);z,i=maximum(a+b);sc=z/np.sqrt(32)
    for t in (0,4,10,float(sc[2]),float(np.nextafter(sc[2],np.inf))):
        p,s,j,_,_=join_semijoin(a,b,t);assert np.array_equal(p,sc>=t);assert np.array_equal(i[p],j[p]);verify(V[p],36.1,j[p],s[p],"selftest")
    # Brute verify all hypotheses for two selected slopes and every start.
    for k in (5,41):
        for start in range(256):
            d=[start+np.rint(GRID[k]*(np.arange(16)+off)).astype(int) for off in (0,36.1)]
            if all(np.all((p>=0)&(p<256)) for p in d):
                assert abs(a[3,k,start]+b[3,k,start]-sum(V[3,j,np.arange(16),d[j]].sum() for j in range(2)))<1e-9
    return {"threshold_tests":5,"brute_slopes":2,"status":"PASS"}
def main():
    checks=selftest();m=json.loads((OUT/"manifest.json").read_text())
    assert m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest() and len(m["selected"])==6
    assert len({q["target"] for q in m["selected"]})==6 and all(q["target"] not in hold.S["exclude"] for q in m["selected"])
    t=time.perf_counter();cal=[gather(q) for q in m["selected"][:2]]
    load=time.perf_counter()-t;lock=calibrate(cal);lock["load_seconds"]=load;lock["targets"]=[q["target"] for q in m["selected"][:2]]
    save("lock.json",lock);emit("R21BC_LOCK",lock)
    seen={s for rows,_,_ in cal for _,_,meta in rows for s in meta["source"]};out={};deltas={};load=0
    for q in m["selected"][2:]:
        t=time.perf_counter();rows,ch,delta=gather(q);load+=time.perf_counter()-t
        hashes={s for _,_,meta in rows for s in meta["source"]};assert not seen&hashes;seen|=hashes
        out[q["target"]]=bench(rows,ch,delta,lock);deltas[q["target"]]=delta;save("partial.json",out)
        emit("R21BC_TARGET",{"target":q["target"],"delta_rows":delta,"metrics":out[q["target"]]})
    keys=("detected","positives","control_alerts","controls","wall","cpu","support","hypotheses")
    total={k:{a:sum(x[k][a] for x in out.values()) for a in keys} for k in S["modes"]}
    for k in total:
        for field,values in (("amplitudes",[str(a) for a in S["amplitudes"]]),("presence",S["presence"]),("gaps",[str(a) for a in S["gap_levels"]]),("amplitude_presence",[str(a)+"_"+b for a in S["amplitudes"] for b in S["presence"]])):
            total[k][field]={v:{a:sum(x[k][field][v][a] for x in out.values()) for a in keys[:4]} for v in values}
        total[k]["localization"]=None if k=="unconstrained_sum" else {a:sum(x[k]["localization"][a] for x in out.values()) for a in ("localized_positive","detected_positive")}
        total[k]["challenges"]={a:sum(x[k]["challenges"][a] for x in out.values()) for a in ("stationary_rfi","discordant")}
    d={"status":"FRESH_ARCHIVE_SYNTHETIC_EXAM","samples":1024,"challenge_cases":128,"targets":list(out),"delta_rows":deltas,"lock":lock,"selftests":checks,"semijoin_mismatches":0,"modes":total,"per_target":out,"evaluation_load_seconds":load,"scope":"Two real same-pointing OFF scans per target; injected constant drift and missing-scan cases. No phase coherence, no origin classification. Semijoin exact only on fixed finite shared hypothesis bank; all score maps computed first, no claim that support count equals total cost. Artificial RFI alerts are separate diagnostics.","production_adoption":False}
    save("result.json",d);emit("R21BC_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="scan":main()
    else:raise ValueError("phase")
