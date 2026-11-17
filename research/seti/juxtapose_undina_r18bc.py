#!/usr/bin/env python3
"""R18B/C pilot: native bounded admission and first-order conditional pivot ranking."""
import concurrent.futures as cf
import hashlib,json,math,pathlib,statistics,sys,time,resource
import numpy as np
import juxtapose_r18a as j
import r17c_frozen_holdout as hold
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r18bc_protocol.json")
S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r18bc");OUT.mkdir(parents=True,exist_ok=True)
hold.S={**hold.S,"seed":S["seed"],"exclude":hold.S["exclude"]+S["exclude_extra"]}
hold.OUT=OUT;hold.P=P;eng.OUT=OUT
R=tuple(S["rank"]);ORDER=tuple(s for g in R for s in j.GROUPS[g])
assert set(ORDER)==set(j.ASC) and len(ORDER)==65 and j.TH==S["threshold"]
def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2)+"\n")
def emit(k,x):print(k,json.dumps(x),flush=True)
def hashbytes(x):return hashlib.sha256(x).hexdigest()
def calibrate():
    t=time.perf_counter();X,y,_,targets,_=j.build(j.S["source_artifacts"][0])
    cache=j.cache_cost(X)
    hits=cache<=np.array([len(g) for g in j.GROUPS])[None,:]
    # Dirichlet-smoothed P(next group hits | last group failed), per cost.
    prob={}
    for prev in range(-1,8):
        mask=np.ones(len(X),bool) if prev==-1 else ~hits[:,prev]
        prob[str(prev)]=[(int(hits[mask,g].sum())+1)/(int(mask.sum())+2)/len(j.GROUPS[g]) for g in range(8)]
    remaining=set(range(8));prev=-1;rank=[]
    while remaining:
        g=min(remaining,key=lambda g:(-prob[str(prev)][g],R.index(g)))
        rank.append(g);remaining.remove(g);prev=g
    lock={"protocol_sha256":hashbytes(P.read_bytes()),"markov_order":rank,"conditional_hit_per_cost":prob,"training_targets":sorted(targets),"calibration_seconds":time.perf_counter()-t,"calibration_samples":len(X),"status":"FROZEN_BEFORE_NEW_SPECTRAL_POWER"}
    save("lock.json",lock);emit("R18BC_LOCK",lock);return lock
def gather(q):
    data=[]
    for f in q["files"]:
        if not f["role"].startswith("OFF"):continue
        for tile in eng.get_tiles(f):
            raw=np.load(OUT/tile["artifact"],allow_pickle=False)
            assert hashbytes(raw.tobytes())==tile["slice_sha256"] and raw.shape[0]==16
            rng=np.random.default_rng(int.from_bytes(hashlib.sha256((S["seed"]+"|"+f["name"]+"|"+str(tile["center_mhz"])).encode()).digest()[:8],"big"))
            centers=[]
            for k in np.argsort(-np.median(j.normalized(raw),axis=0),kind="stable"):
                if 64<=k<raw.shape[1]-64 and all(abs(int(k)-c)>=128 for c in centers):centers.append(int(k))
                if len(centers)==4:break
            assert len(centers)==4
            for trial in range(16):
                start=centers[trial-12]-64 if trial>=12 else int(rng.integers(0,raw.shape[1]-127))
                p=raw[:,start:start+128].copy();shift=int(rng.integers(-32,33));ds=np.rint(np.arange(16)*shift/15).astype(int)
                c=int(rng.integers(max(16,16-ds.min()),min(112,112-ds.max())))
                scale=1.4826*np.median(np.abs(p-np.median(p,axis=1)[:,None]),axis=1)
                for label in (0,1):
                    x=p.copy()
                    if label:x[np.arange(16),c+ds]+=(3,4,6,8)[trial%4]*scale
                    data.append((x,label,trial>=12,tile["slice_sha256"]))
    assert len(data)==256;return data
def check(v,s):
    a,b=j.score(v,s)
    assert math.isfinite(a)
    if a>=j.TH:
        w=j.witness(v,s,b)
        assert abs(a-w)<=max(1e-5,abs(a)*2e-6)
    return a>=j.TH
def sequential(v,order,full=False):
    hit=False;n=0
    for s in order:
        hit=check(v,s) or hit;n+=1
        if hit and not full:break
    return dict(hit=hit,checks=n,epochs=0,extra=0,drain=0,late=0)
def vectorized(v,order):
    n,w=v.shape;ss=np.asarray(order)
    ds=np.rint(np.arange(n)[None,:]*ss[:,None]/(n-1)).astype(int)
    idx=np.arange(w)[None,:]+ds[:,0,None]
    acc=np.zeros((len(ss),w),dtype=v.dtype)
    valid=np.ones((len(ss),w),bool)
    for i in range(n):
        idx=np.arange(w)[None,:]+ds[:,i,None]
        valid &= (idx>=0)&(idx<w)
        acc+=v[i,np.clip(idx,0,w-1)]
    vals=np.where(valid,acc/np.sqrt(n),-np.inf)
    hit=bool(np.any(vals>=j.TH))
    if hit:
        g,b=np.argwhere(vals>=j.TH)[0]
        assert check(v,int(ss[g]))
    return dict(hit=hit,checks=len(ss),epochs=0,extra=0,drain=0,late=0)
def gated(v,order,pool,policy):
    assert len(set(order))==65
    cap=S["cap"];quota=S["quota"];idx=0;done=set();epochs=0;extra=0;drain=0;late=0;hit=False;receipts=[];rate=None
    while idx<len(order) and not hit:
        epochs+=1;start=idx;limit=min(len(order),start+quota);inflight={};closed=False;last=time.perf_counter()
        while True:
            # Admission is synchronous and bounded independently of telemetry.
            lead=0
            if policy=="undina" and rate is not None:
                lead=min(quota-1,math.ceil(rate*S["latency_p95_seconds"])+S["margin"])
            elif policy=="sham":lead=1
            if idx>=limit or (idx>start and idx+lead>=limit):
                closed=True
            if closed and policy=="undina" and idx>=limit:late+=1
            if not closed and not hit:
                while idx<limit and len(inflight)<cap:
                    if policy=="fixed" and inflight:break
                    if idx>start and idx+lead>=limit:closed=True;break
                    s=order[idx];idx+=1
                    assert s not in done and s not in inflight.values()
                    inflight[pool.submit(check,v,s)]=s
                # Fixed batching fills a wave only when previous batch drained.
                if policy=="fixed" and len(inflight)==1:
                    while idx<limit and len(inflight)<cap:
                        s=order[idx];idx+=1;inflight[pool.submit(check,v,s)]=s
            if not inflight:break
            was_closed=closed or hit
            completed,_=cf.wait(inflight,return_when=cf.ALL_COMPLETED if policy=="fixed" else cf.FIRST_COMPLETED)
            now=time.perf_counter();instant=len(completed)/max(now-last,1e-9)
            rate=instant if rate is None else .8*rate+.2*instant;last=now
            new_hit=False
            for f in completed:
                s=inflight.pop(f);assert s not in done;done.add(s)
                yes=f.result();new_hit|=yes
                if was_closed:drain+=1
                receipts.append((epochs,s,yes))
            if new_hit and not hit:
                hit=True;extra+=len(inflight)+len(completed)-1;closed=True
            if idx>=limit:closed=True
        assert not inflight and len(done)==idx and idx>start
        # Two separate quiescence ledger checks, no pending futures.
        assert len({s for _,s,_ in receipts})==idx
        assert set(order[:idx])==done
    if not hit:assert len(done)==65
    return dict(hit=hit,checks=len(done),epochs=epochs,extra=extra,drain=drain,late=late)
def bench(data,markov,pool):
    mo=tuple(s for g in markov for s in j.GROUPS[g])
    methods={
        "full":lambda v:sequential(v,j.ASC,True),
        "center":lambda v:sequential(v,j.CENTER),
        "juxtapose":lambda v:sequential(v,ORDER),
        "fixed":lambda v:gated(v,ORDER,pool,"fixed"),
        "undina":lambda v:gated(v,ORDER,pool,"undina"),
        "sham":lambda v:gated(v,ORDER,pool,"sham"),
        "center_undina":lambda v:gated(v,j.CENTER,pool,"undina"),
        "markov":lambda v:sequential(v,mo),
        "markov_undina":lambda v:gated(v,mo,pool,"undina"),
        "vectorized":lambda v:vectorized(v,j.ASC)}
    combined={k:[] for k in methods};ref=None
    for rep in range(S["repeats"]):
        names=list(methods);rng=np.random.default_rng(31000+rep);rng.shuffle(names)
        for k in names:
            wall=time.perf_counter();cpu=time.process_time()
            res=[methods[k](j.normalized(x)) for x,_,_,_ in data]
            elapsed=time.perf_counter()-wall;cput=time.process_time()-cpu
            pred=[r["hit"] for r in res]
            if ref is None:ref=pred
            assert pred==ref,("decision mismatch",k)
            combined[k].append({"wall_seconds":elapsed,"cpu_seconds":cput,**{a:sum(r[a] for r in res) for a in ("checks","epochs","extra","drain","late")}})
    return {k:{"repeats":v,"wall_median":statistics.median(r["wall_seconds"] for r in v),"cpu_median":statistics.median(r["cpu_seconds"] for r in v),"checks_median":statistics.median(r["checks"] for r in v)} for k,v in combined.items()}
def main():
    lock=calibrate();t=time.perf_counter();m=json.loads((OUT/"manifest.json").read_text())
    assert m["protocol_sha256"]==hashbytes(P.read_bytes())
    assert len(m["selected"])==6
    selected=m["selected"][2:];assert len({q["target"] for q in selected})==4
    assert not set(q["target"] for q in selected)&set(hold.S["exclude"])
    rows={};load=0
    with cf.ThreadPoolExecutor(max_workers=S["cap"]) as pool:
        for q in selected:
            a=time.perf_counter();data=gather(q);load+=time.perf_counter()-a
            rows[q["target"]]=bench(data,lock["markov_order"],pool)
            emit("R18BC_TARGET",{"target":q["target"],"metrics":rows[q["target"]]})
            save("partial.json",rows)
    totals={k:{a:sum(v[k][a] for v in rows.values()) for a in ("wall_median","cpu_median","checks_median")} for k in next(iter(rows.values()))}
    summary={"targets":list(rows),"samples":1024,"decision_mismatches":0,"unknown":0,"totals":totals,"per_target":rows,"lock":lock,"load_seconds":load,"run_seconds":time.perf_counter()-t,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"scope":"Native in-process pilot; fixed settings; lead-free comparison and fixed lead sham, not matched-time sham. Fresh archive targets, synthetic injections. No production promotion.","production_adoption":False}
    save("result.json",summary);emit("R18BC_SUMMARY",{k:v for k,v in summary.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="scan":main()
    else:raise ValueError("phase")
