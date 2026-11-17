#!/usr/bin/env python3
"""Coarse vector-wave R18B2 with receipt-grounded context and conditional rank ablations."""
import hashlib,json,pathlib,statistics,sys,time,resource
import numpy as np
import juxtapose_r18a as j
import juxtapose_undina_r18bc as old
import r17c_frozen_holdout as hold
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r18b2_vector_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r18b2");OUT.mkdir(parents=True,exist_ok=True)
hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))}
hold.OUT=OUT;hold.P=P;eng.OUT=OUT;old.OUT=OUT;old.S={**old.S,"seed":S["seed"]}
R=tuple(S["rank"]);C=tuple(range(8));TH=j.TH
assert TH==S["threshold"]
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2)+"\n")
def emit(k,x):print(k,json.dumps(x),flush=True)
def flat(groups):return tuple(s for g in groups for s in j.GROUPS[g])
def normalize(X):return np.asarray([j.normalized(x) for x in X])
def feature(v):return float(np.max(np.median(v,axis=0)))
def vectors(V,shifts):
    n,w=V.shape[1:];ss=np.asarray(shifts)
    ds=np.rint(np.arange(n)[None,:]*ss[:,None]/(n-1)).astype(int)
    acc=np.zeros((len(V),len(ss),w),dtype=V.dtype);valid=np.ones((len(ss),w),bool)
    for i in range(n):
        idx=np.arange(w)[None,:]+ds[:,i,None]
        valid&=(idx>=0)&(idx<w)
        acc+=V[:,i,np.clip(idx,0,w-1)]
    scores=np.where(valid[None],acc/np.sqrt(n),-np.inf)
    assert np.isfinite(scores.max(axis=(1,2))).all()
    maxima=scores.max(axis=(1,2));pred=maxima>=TH
    # Recompute threshold-near cases with unchanged scorer; verify actual positive witness.
    for i in np.flatnonzero(pred | (np.abs(maxima-TH)<1e-4)):
        g,k=np.unravel_index(int(np.argmax(scores[i])),scores[i].shape)
        a,b=j.score(V[i],int(ss[g]))
        if abs(maxima[i]-TH)<1e-4:
            pred[i]=any(j.score(V[i],s)[0]>=TH for s in shifts)
        elif pred[i]:
            assert a>=TH and abs(a-j.witness(V[i],int(ss[g]),b))<=max(1e-5,abs(a)*2e-6)
    return pred
def vector_route(V,rank,width,policy=None):
    result=np.zeros(len(V),bool);checks=waves=0;width_counts={}
    for base in range(0,len(V),S["tile_batch"]):
        active=np.arange(base,min(len(V),base+S["tile_batch"]));pos=0
        while len(active) and pos<8:
            # Synchronous admission: no next wave until all current receipts finalized.
            remaining=8-pos
            g=min(remaining,width if policy is None else policy["large" if len(active)>=8 else "small"])
            shifts=flat(rank[pos:pos+g]);a=vectors(V[active],shifts)
            checks+=len(active)*len(shifts);waves+=1
            width_counts[str(g)]=width_counts.get(str(g),0)+1
            result[active[a]]=True
            active=active[~a];pos+=g
        if len(active):assert pos==8
    return result,checks,waves,width_counts
def conditional(cache):
    sizes=np.array([len(g) for g in j.GROUPS]);hits=cache<=sizes[None,:]
    survived=np.ones(len(cache),bool);remain=set(C);rank=[]
    while remain:
        g=min(remain,key=lambda g:(-(int(hits[survived,g].sum())+1)/(int(survived.sum())+2)/sizes[g],R.index(g)))
        rank.append(g);remain.remove(g);survived&=~hits[:,g]
    return tuple(rank)
def count_order(cache,rank):
    sizes=np.array([len(g) for g in j.GROUPS]);active=np.ones(len(cache),bool);cost=0
    for g in rank:
        cost+=int(np.minimum(cache[active,g],sizes[g]).sum())
        active&=cache[:,g]>sizes[g]
    return cost
def calibrate():
    timer=time.perf_counter();X,y,_,targets,_=j.build(j.S["source_artifacts"][0])
    V=normalize(X);cache=j.cache_cost(X);features=np.array([feature(v) for v in V]);boundary=float(np.median(features))
    buckets={};routes={}
    for name,mask in [("low",features<=boundary),("high",features>boundary)]:
        c=cache[mask];base=count_order(c,R)
        best,_=j.tournament(c)
        gain=(base-best[0])/base
        routes[name]=tuple(best[2]) if len(c)>=64 and gain>=.05 else R
        buckets[name]={"n":int(len(c)),"global_cost":base,"best_cost":best[0],"gain":gain,"rank":routes[name]}
    markov=conditional(cache)
    # Fixed and coarse telemetry-conditioned policies compared on identical R9 powers.
    widths=S["groups_per_wave"];policies=[{"large":a,"small":b} for a in widths for b in widths]
    tests={}
    ref=vector_route(V,R,8)[0]
    for p in policies:
        k=f'{p["large"]}:{p["small"]}';times=[];count=None
        for _ in range(2):
            t=time.perf_counter();r=vector_route(V,R,8,p);times.append(time.perf_counter()-t)
            assert np.array_equal(r[0],ref);count=r[1]
        tests[k]={"policy":p,"median_seconds":statistics.median(times),"checks":count}
    fixed=min((tests[f"{a}:{a}"] for a in widths),key=lambda d:(d["median_seconds"],d["checks"]))
    adaptive=min(tests.values(),key=lambda d:(d["median_seconds"],d["checks"]))
    lock={"protocol_sha256":hashlib.sha256(P.read_bytes()).hexdigest(),"training_targets":sorted(targets),"n":len(X),"routes":routes,"context_boundary":boundary,"context_buckets":buckets,"markov_rank":markov,"fixed":fixed,"adaptive":adaptive,"tuning":tests,"calibration_seconds":time.perf_counter()-timer}
    save("lock.json",lock);emit("R18B2_LOCK",lock);return lock
def routed(V,lock):
    feats=np.array([feature(v) for v in V]);pred=np.zeros(len(V),bool);checks=waves=0
    for name,ix in [("low",np.flatnonzero(feats<=lock["context_boundary"])),("high",np.flatnonzero(feats>lock["context_boundary"]))]:
        if len(ix):
            r=vector_route(V[ix],lock["routes"][name],lock["fixed"]["policy"]["large"]);pred[ix]=r[0];checks+=r[1];waves+=r[2]
    return pred,checks,waves,{}
def scalar(V):
    pred=[];checks=0
    for v in V:
        r=old.sequential(v,flat(R));pred.append(r["hit"]);checks+=r["checks"]
    return np.asarray(pred),checks,0,{}
def bench(data,lock):
    X=np.asarray([x for x,_,_,_ in data])
    modes={"scalar_j":scalar,"vector_full":lambda V:vector_route(V,C,8),"vector_center_2":lambda V:vector_route(V,C,2)}
    modes.update({f"vector_j_{w}":(lambda V,w=w:vector_route(V,R,w)) for w in S["groups_per_wave"]})
    modes["vector_j_fixed_calibrated"]=lambda V:vector_route(V,R,lock["fixed"]["policy"]["large"])
    modes["vector_undina"]=lambda V:vector_route(V,R,8,lock["adaptive"]["policy"])
    modes["vector_keymaster"]=lambda V:routed(V,lock)
    modes["vector_residual_markov"]=lambda V:vector_route(V,lock["markov_rank"],lock["fixed"]["policy"]["large"])
    ref=np.array([j.decision(x,j.ASC,True)[0] for x in X])
    allrows={k:[] for k in modes}
    for rep in range(S["repeats"]):
        keys=list(modes);np.random.default_rng(18200+rep).shuffle(keys)
        for k in keys:
            t=time.perf_counter();cpu=time.process_time();V=normalize(X);r=modes[k](V)
            dt=time.perf_counter()-t;dc=time.process_time()-cpu
            assert np.array_equal(r[0],ref),k
            allrows[k].append({"wall_seconds":dt,"cpu_seconds":dc,"checks":r[1],"waves":r[2],"width_counts":r[3]})
    return {k:{"wall_median":statistics.median(x["wall_seconds"] for x in a),"cpu_median":statistics.median(x["cpu_seconds"] for x in a),"checks_median":statistics.median(x["checks"] for x in a),"repeats":a} for k,a in allrows.items()}
def main():
    lock=calibrate();m=json.loads((OUT/"manifest.json").read_text())
    assert m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest() and len(m["selected"])==6
    rows={};load=0;t=time.perf_counter()
    for q in m["selected"][2:]:
        assert q["target"] not in hold.S["exclude"]
        a=time.perf_counter();data=old.gather(q);load+=time.perf_counter()-a
        rows[q["target"]]=bench(data,lock);save("partial.json",rows)
        emit("R18B2_TARGET",{"target":q["target"],"metrics":rows[q["target"]]})
    assert len(rows)==4
    total={k:{a:sum(v[k][a] for v in rows.values()) for a in ("wall_median","cpu_median","checks_median")} for k in next(iter(rows.values()))}
    best=min(("vector_full","vector_center_2","vector_j_fixed_calibrated"),key=lambda k:total[k]["wall_median"])
    summary={"targets":list(rows),"samples":1024,"decision_mismatches":0,"unknown":0,"totals":total,"per_target":rows,"calibration":lock,"best_baseline_descriptive":best,"load_seconds":load,"run_seconds":time.perf_counter()-t,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"scope":"synchronous vector adaptation, not original native asynchronous gate; calibration context and Markov rank only; synthetic injections in fresh archive OFF data","production_adoption":False}
    save("result.json",summary);emit("R18B2_SUMMARY",{k:v for k,v in summary.items() if k not in ("per_target","calibration")})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="scan":main()
    else:raise ValueError("phase")
