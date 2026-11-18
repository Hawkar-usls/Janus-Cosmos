#!/usr/bin/env python3
"""R36A-E: pinned Kenshi and exact two-bank twin search on historical common spectra."""
import hashlib,io,itertools,json,os,pathlib,resource,time,zipfile,concurrent.futures
import numpy as np
import requests
from numba import njit
from r36_kenshi_blackboard import KenshiSwarmBlackboard,self_test
P=pathlib.Path("research/seti/r36_galaxy_swarm_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r36");OUT.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def bank():
    b=json.loads(pathlib.Path(S["bank"]["source"]).read_text())["r35"]["bank"];t=np.arange(16);u=t-7.5
    linear=[tuple(np.rint(t*s/15).astype(int)) for s in range(-32,33)]
    a=[np.array(x) for x in linear]
    a += [np.rint(s*u+k*u*u).astype(int) for s,k in itertools.product(b["slopes"],b["curvature"])]
    a += [np.rint(s*u+v*np.sin(2*np.pi*t/16+k*2*np.pi/b["phases"])).astype(int) for s,v,k in itertools.product(b["sinusoid_slopes"],b["sinusoid_amplitude"],range(b["phases"]))]
    smooth=sorted(set(tuple(x-x[0]) for x in a if np.max(abs(np.diff(x)))<=b["max_step"]))
    full=sorted(set(linear+smooth));ids={d:i for i,d in enumerate(full)}
    rank=np.full((len(full),2),-1,np.int64)
    L=np.array([ids[d] for d in linear],np.int64);T=np.array([ids[d] for d in smooth],np.int64)
    rank[L,0]=np.arange(len(L));rank[T,1]=np.arange(len(T))
    return np.array(full,np.int64),rank,L,T
D,R,L,T=bank();ALL=np.arange(len(D),dtype=np.int64);COST=128-(D.max(1)-D.min(1))
@njit(nogil=True,fastmath=False,cache=False)
def scan(V,D,R,ids,prune):
    B,n,w=V.shape;best=np.full((B,2),-np.inf);bid=np.full((B,2),-1,np.int64);start=np.full((B,2),-1,np.int64)
    reads=0;paths=0;pruned=0;bound_reads=0
    for i in range(B):
        suffix=np.zeros(n+1);ma=0.
        if prune:
            for t in range(n-1,-1,-1):
                mx=-np.inf
                for k in range(w):
                    x=V[i,t,k];mx=max(mx,x);ma=max(ma,abs(x));bound_reads+=1
                suffix[t]=suffix[t+1]+mx
        eps=1e-10*(1+n*ma)
        for k in ids:
            d=D[k];lo=max(0,-d.min());hi=min(w,w-d.max())
            for j in range(lo,hi):
                paths+=1;s=0.;stop=False
                cutoff=np.inf
                for a in range(2):
                    if R[k,a]>=0:cutoff=min(cutoff,best[i,a])
                for t in range(n):
                    s+=V[i,t,j+d[t]];reads+=1
                    if prune and t<n-1 and s+suffix[t+1]+eps<cutoff*4:
                        pruned+=1;stop=True;break
                if stop:continue
                z=s/4.
                for a in range(2):
                    if R[k,a]<0:continue
                    old=bid[i,a]
                    if z>best[i,a] or (z==best[i,a] and (old<0 or R[k,a]<R[old,a] or (R[k,a]==R[old,a] and j<start[i,a]))):
                        best[i,a]=z;bid[i,a]=k;start[i,a]=j
    return best,bid,start,np.array([reads,paths,pruned,bound_reads],np.int64)
def exact(a,b):assert all(np.array_equal(x,y) for x,y in zip(a[:3],b[:3])),"score or canonical witness mismatch"
def separate(V):
    a=scan(V,D,R,L,False);b=scan(V,D,R,T,False)
    return tuple(np.stack([a[i][:,0],b[i][:,1]],axis=1) for i in range(3))+(a[3]+b[3],)
def merge(a,b):
    s,i,j=[x.copy() for x in a[:3]]
    for k in range(2):
        for n in range(len(s)):
            u=a[1][n,k];v=b[1][n,k]
            if v>=0 and (u<0 or b[0][n,k]>s[n,k] or (b[0][n,k]==s[n,k] and (R[v,k],b[2][n,k])<(R[u,k],j[n,k]))):
                s[n,k]=b[0][n,k];i[n,k]=v;j[n,k]=b[2][n,k]
    return s,i,j,a[3]+b[3]
def partition(order):
    v=[[],[]];load=[0,0]
    for k in order:
        j=0 if load[0]<=load[1] else 1;v[j].append(k);load[j]+=int(COST[k])
    assert set(v[0]).isdisjoint(v[1]) and sorted(v[0]+v[1])==ALL.tolist()
    return [np.array(x,np.int64) for x in v]
def twin(V,ids,prune,executor):
    fut=[executor.submit(scan,V,D,R,x,prune) for x in ids]
    return merge(*(f.result() for f in fut))
def tail_rank(V):
    # TailGEX-inspired measured exceedance density, no hash randomness assumption.
    t=time.perf_counter();vals=[]
    for k,d in enumerate(D):
        lo=max(0,-int(d.min()));hi=min(V.shape[-1],V.shape[-1]-int(d.max()))
        a=sum(V[:,i,lo+v:hi+v] for i,v in enumerate(d))/4
        score=sum(w*int(np.count_nonzero(a>=z)) for w,z in zip(S["calibration"]["tail_weights"],S["calibration"]["tail_thresholds"]))
        vals.append(score/a.size)
    order=np.array(sorted(range(len(D)),key=lambda k:(-vals[k],k)),np.int64)
    return order,vals,time.perf_counter()-t
def compare_numpy(V,base):
    for a,ids in enumerate((L,T)):
        best=np.full(len(V),-np.inf);bi=np.full(len(V),-1);bj=bi.copy()
        for k in ids:
            d=D[k];lo=max(0,-d.min());hi=min(V.shape[-1],V.shape[-1]-d.max())
            z=sum(V[:,i,lo+v:hi+v] for i,v in enumerate(d))/4
            j=z.argmax(1);s=z[np.arange(len(z)),j];mask=s>best
            best[mask]=s[mask];bi[mask]=k;bj[mask]=lo+j[mask]
        assert np.array_equal(best,base[0][:,a]) and np.array_equal(bi,base[1][:,a]) and np.array_equal(bj,base[2][:,a])
def fetch():
    s=S["source"];u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{s['artifact']}/zip"
    r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90);r.raise_for_status();raw=r.content
    assert len(raw)<64*1024**2 and sha(raw)==s["sha256"]
    z=zipfile.ZipFile(io.BytesIO(raw));assert len(z.namelist())==6 and sum(i.file_size for i in z.infolist())<64*1024**2
    data=[]
    for n in sorted(z.namelist()):
        with np.load(io.BytesIO(z.read(n)),allow_pickle=False) as q:
            d={k:q[k].copy() for k in q.files};V=d["V"];assert V.ndim==3 and V.shape[1:]==(16,128) and np.isfinite(V).all()
            data.append((n,d,sha(V.tobytes())))
    z.close();assert sum(n.startswith("cal_") for n,d,h in data)==2 and sum(n.startswith("eval_") for n,d,h in data)==4
    return data
def output(r):
    v=r[0];a=v.argmax(1);ix=np.arange(len(v));k=r[1][ix,a];j=r[2][ix,a];assert np.all(k>=0)
    return v[ix,a],D[k]+j[:,None]
def cut(v):
    u=np.sort(v)[::-1];return float(np.nextafter(u[int(len(u)*.05)],np.inf))
def metrics(q,r,th):
    z,p=output(r);y=q["y"];act=q["active"];loc=((abs(p-q["truth"])<=2)&act).sum(1)>=.75*act.sum(1)
    alert=z>=th;hit=alert&loc&(y==1);base=q["registered_hybrid_alert"]&q["registered_hybrid_localized"]&(y==1)
    return {"localized":int(hit.sum()),"weak":int((hit&(q["amp"]<=1.25)).sum()),"very_weak":int((hit&(q["amp"]==.75)).sum()),"controls":int((y==0).sum()),"positives":int((y==1).sum()),"control_alerts":int((alert&(y==0)).sum()),"added_vs_r28":int((hit&~base).sum()),"lost_vs_r28":int((base&~hit).sum()),"r28_localized":int(base.sum()),"r28_control_alerts":int((q["registered_hybrid_alert"]&(y==0)).sum())}
def timed(fn):
    s=time.perf_counter();c=time.process_time();v=fn();return v,{"wall":time.perf_counter()-s,"cpu":time.process_time()-c}
def measure(V,order,ids,executor):
    fs={"duplicate":lambda:separate(V),"kenshi_shared":lambda:scan(V,D,R,ALL,False),"twin_fixed":lambda:twin(V,ids,False,executor),"tail_order_prune":lambda:scan(V,D,R,order,True),"tail_twin":lambda:twin(V,partition(order),True,executor)}
    timings={k:[] for k in fs};results={}
    for rep in range(S["measurement"]["repeat"]):
        for k in list(fs)[rep:]+list(fs)[:rep]:
            a,t=timed(fs[k]);timings[k].append(t);results[k]=a
        for k,a in results.items():exact(results["duplicate"],a)
    compare_numpy(V,results["duplicate"])
    return results["duplicate"],{k:{"median_wall_seconds":float(np.median([v["wall"] for v in timings[k]])),"median_cpu_seconds":float(np.median([v["cpu"] for v in timings[k]])),"all_repeats":timings[k],"row_reads":int(a[3][0]),"candidate_paths":int(a[3][1]),"pruned_paths":int(a[3][2]),"bound_reads":int(a[3][3]),"exact":True} for k,a in results.items()}
def core_and_failure(V,order,expected,executor):
    st=time.perf_counter();bs=S["measurement"]["chunk"];chunks=[V[i:i+bs] for i in range(0,len(V),bs)]
    priorities=[float(np.max(np.abs(x))) for x in chunks];priority_order=sorted(range(len(chunks)),key=lambda i:(-priorities[i],i));prep=time.perf_counter()-st
    ids=partition(order)
    def schedule(seq):
        t=time.perf_counter();cpu=time.process_time();results={};done={}
        for k in seq:
            r=twin(chunks[k],ids,True,executor);results[k]=r;done[k]=time.perf_counter()-t
        got=tuple(np.concatenate([results[k][j] for k in range(len(chunks))]) for j in range(3))+(sum((r[3] for r in results.values()),np.zeros(4,np.int64)),)
        exact(expected,got)
        return {"wall_seconds":time.perf_counter()-t,"cpu_seconds":time.process_time()-cpu,"priority_top_quartile_completion_seconds":float(np.mean([done[k] for k in priority_order[:max(1,len(chunks)//4)]])),"completed":len(done),"row_reads":int(got[3][0]),"bound_reads":int(got[3][3])}
    fifo=schedule(list(range(len(chunks))));priority=schedule(priority_order)
    bb=KenshiSwarmBlackboard(active_window=16)
    for a,role in [("CORE2","CONTROLLER"),("ANCHOR","STABLE"),("GLADIUS","EXPLORER")]:bb.register_agent(a,role)
    x=chunks[0];task=sha(x.tobytes()+ids[1].tobytes()+P.read_bytes());commits={};turn=1
    def claim(agent,task,turn):
        return bb.claim_lease(agent,task,now_turn=turn,ttl_turns=2)
    def commit(agent,task,turn,result):
        lease=bb.leases[task]
        if lease["owner_agent_id"]!=agent or lease["expires_after_turn"]<turn or task in commits:raise RuntimeError("stale or duplicate commit")
        digest=sha(b"".join(r.tobytes() for r in result));commits[task]=digest
        bb.append_event(agent,"CAUSAL_ACK",{"task":task,"output_sha256":digest},parent_event_ids=[lease["event_id"]],source_bound=True)
    t=time.perf_counter();c=time.process_time();ev=claim("GLADIUS",task,1)
    interrupted=scan(x,D,R,ids[1],True)
    try:raise RuntimeError("injected failure before commit")
    except RuntimeError:bb.mark_stale("GLADIUS","injected precommit failure")
    denied=0
    try:claim("ANCHOR",task,2)
    except RuntimeError:denied+=1
    claim("ANCHOR",task,4)
    try:commit("GLADIUS",task,4,interrupted)
    except RuntimeError:denied+=1
    retry=scan(x,D,R,ids[1],True);exact(interrupted,retry);commit("ANCHOR",task,4,retry)
    try:commit("ANCHOR",task,4,retry)
    except RuntimeError:denied+=1
    a=scan(x,D,R,ids[0],True);got=merge(a,retry);exact(tuple(v[:len(x)] for v in expected[:3]),got)
    bb.recover("GLADIUS",{"status":"RECOVERED","last_task":task});bb.validate();bb.save(OUT/"kenshi_failure_ledger.json")
    assert denied==3 and len(commits)==1 and len(bb.agents)==3
    return {"priority_preprocessing_seconds":prep,"priority":priority,"fifo":fifo,"priority_order":priority_order,"scope":"Label-free peak proxy, full coverage; not confidence or recall improvement."},{"exact_recovered":True,"denied_conflict_stale_duplicate":denied,"unique_commits":len(commits),"wasted_row_reads":int(interrupted[3][0]),"wasted_bound_reads":int(interrupted[3][3]),"retry_row_reads":int(retry[3][0]),"events":len(bb.events),"wall_seconds":time.perf_counter()-t,"cpu_seconds":time.process_time()-c,"injected_precommit_failure":True,"source_bound_is_not_independent_replication":True}
def main():
    # Verify byte-for-byte donor ancestry, then actual donor unit test.
    raw=pathlib.Path(S["donors"]["kenshi"]["local"]).read_bytes();assert hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()==S["donors"]["kenshi"]["blob"];self_test()
    t=time.perf_counter();c=time.process_time();rng=np.random.default_rng(36001);v=rng.normal(size=(8,16,128));v[0]=0;v[1]=1;v[2]=-1;v[3]*=1e8
    a=separate(v);b=scan(v,D,R,ALL[::-1].copy(),True);exact(a,b);compare_numpy(v,a);compile={"wall_seconds":time.perf_counter()-t,"cpu_seconds":time.process_time()-c}
    data=fetch();cal=[d for n,d,h in data if n.startswith("cal_")];Vcal=np.ascontiguousarray(np.concatenate([d["V"] for d in cal]),dtype=np.float64)
    order,tails,setup=tail_rank(Vcal);ids=partition(ALL)
    scores=[]
    for d in cal:
        z,p=output(scan(np.ascontiguousarray(d["V"],dtype=np.float64),D,R,ALL,False));scores.append(z[d["y"]==0])
    threshold=cut(np.concatenate(scores));lock={"threshold":threshold,"order":order.tolist(),"tail_density":tails,"bank_sha256":sha(D.tobytes()),"membership_sha256":sha(R.tobytes()),"calibration_arrays":[h for n,d,h in data if n.startswith("cal_")],"calibration_alerts":int((np.concatenate(scores)>=threshold).sum())}
    save("lock.json",lock);lockhash=sha((OUT/"lock.json").read_bytes());emit("R36_LOCK",{"sha256":lockhash,"linear":len(L),"smooth":len(T),"unique":len(D),"overlap":len(L)+len(T)-len(D),"threshold":threshold})
    rows=[];core=None;fault=None
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for n,d,h in data:
            if not n.startswith("eval_"):continue
            V=np.ascontiguousarray(d["V"],dtype=np.float64);ref,m=measure(V,order,ids,pool);r={"target":n,"array_sha256":h,"metrics":metrics(d,ref,threshold),"arms":m}
            if core is None:core,fault=core_and_failure(V,order,ref,pool)
            rows.append(r);save("partial.json",rows);emit("R36_TARGET",r)
    assert len(rows)==4 and sha((OUT/"lock.json").read_bytes())==lockhash
    total={k:sum(x["metrics"][k] for x in rows) for k in rows[0]["metrics"]};assert total["positives"]==total["controls"]==512
    agg={k:{f:sum(x["arms"][k][f] for x in rows) for f in ("median_wall_seconds","median_cpu_seconds","row_reads","candidate_paths","pruned_paths","bound_reads")} for k in rows[0]["arms"]}
    base=agg["duplicate"];gates={k:{"exact":True,"fewer_row_reads_including_bound":v["row_reads"]+v["bound_reads"]<base["row_reads"],"lower_wall":v["median_wall_seconds"]<base["median_wall_seconds"],"lower_cpu":v["median_cpu_seconds"]<base["median_cpu_seconds"],"lower_wall_with_all_compile_and_tail_setup":v["median_wall_seconds"]+compile["wall_seconds"]+setup<base["median_wall_seconds"]} for k,v in agg.items() if k!="duplicate"}
    result={"status":"PASS_EXACT_DEVELOPMENT_NOT_FULL_DETECTOR","A_dedup":{"linear":len(L),"smooth":len(T),"unique":len(D),"duplicate_templates":len(L)+len(T)-len(D)},"B_twin":{"workers":2,"exclusive_template_counts":[len(x) for x in ids],"total_coverage":len(D)},"C_tail":{"setup_seconds":setup,"fixed_calibration_only":True},"D_core":core,"E_failure":fault,"F":"HOLD_FULL_NATIVE_COMPARISON_REQUIRES_QUALIFIED_INTEGRATED_INPUTS","metrics":total,"timing":agg,"gates":gates,"compile":compile,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"cpu_count":os.cpu_count(),"lock_sha256":lockhash,"source":S["source"],"production_adoption":False,"scope":"Exact isolated two-bank max, no new ODONTO/HANNAH behavior, historical registered arrays and synthetic injections. No end-to-end or independent science claim."}
    save("result.json",result);emit("R36_SUMMARY",result)
if __name__=="__main__":main()
