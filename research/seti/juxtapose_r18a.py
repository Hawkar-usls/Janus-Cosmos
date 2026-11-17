#!/usr/bin/env python3
"""Source-linked JUXTAPOSE: eight spectral pivot groups, exact existential decision."""
import hashlib,io,json,math,os,pathlib,re,zipfile
import numpy as np
import requests
import itertools,time,statistics
from r9_expanded_search import normalized
S=json.loads(pathlib.Path("research/seti/r15_r16_r17_protocol.json").read_text())["r17"]
OUT=pathlib.Path("runtime/juxtapose-r18a");OUT.mkdir(parents=True,exist_ok=True)
def emit(k,x):print(k,json.dumps(x),flush=True)
def download(a):
    u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip"
    buf=bytearray()
    with requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},stream=True,timeout=(15,60)) as r:
        r.raise_for_status()
        for b in r.iter_content(1024*1024):
            buf.extend(b);assert len(buf)<=160*1024*1024
    assert hashlib.sha256(buf).hexdigest()==a["digest"]
    z=zipfile.ZipFile(io.BytesIO(buf));assert sum(i.file_size for i in z.infolist())<=240*1024*1024
    names={}
    for i in z.infolist():
        b=pathlib.PurePosixPath(i.filename).name
        if b in ("r9-result.json","result.json") or re.fullmatch("[0-9a-f]{24}\\.npy",b):
            assert b not in names;names[b]=i.filename
    data=json.loads(z.read(names["r9-result.json" if a["split"]=="train_R9" else "result.json"]))
    return z,names,data
def build(a):
    z,names,data=download(a);X=[];Y=[];scores=[];rows=[];targets=set();library=[]
    for q in data["results"]:
        target=q["target"];targets.add(target)
        for o in q["observations"]:
            if not o["role"].startswith("OFF"):continue
            for t in o["tiles"]:
                x=np.load(io.BytesIO(z.read(names[t["artifact"]])),allow_pickle=False)
                assert list(x.shape)==t["shape"] and hashlib.sha256(x.tobytes()).hexdigest()==t["slice_sha256"]
                assert x.shape[0]==16
                rng=np.random.default_rng(int.from_bytes(hashlib.sha256((str(S["seed"])+"|"+a["split"]+"|"+o["name"]+"|"+str(t["center_mhz"])).encode()).digest()[:8],"big"))
                width=S["patch_channels"]
                # Four additional hard patches around the strongest stationary peaks.
                med=np.median(normalized(x),axis=0)
                order=np.argsort(-med,kind="stable"); centers=[]
                for k in order:
                    if width//2<=k<x.shape[1]-width//2 and all(abs(int(k)-c)>=width for c in centers):centers.append(int(k))
                    if len(centers)==4:break
                assert len(centers)==4
                for trial in range(S["trials_per_tile"]+4):
                    hard=trial>=S["trials_per_tile"]
                    start=centers[trial-S["trials_per_tile"]]-width//2 if hard else int(rng.integers(0,x.shape[1]-width+1))
                    raw=x[:,start:start+width].copy()
                    shift=int(rng.integers(-32,33));ds=np.rint(np.arange(16)*shift/15).astype(int);c=int(rng.integers(max(16,16-ds.min()),min(width-16,width-16-ds.max())))
                    amp=S["amplitude"][trial%len(S["amplitude"])]
                    scale=1.4826*np.median(np.abs(raw-np.median(raw,axis=1)[:,None]),axis=1)
                    for label in (0,1):
                        p=raw.copy()
                        if label:p[np.arange(16),c+ds]+=amp*scale
                        X.append(p);Y.append(label)
                        # Known injected trajectory is NOT used by this baseline; blind max over all shifts.
                        rows.append({"target":target,"source":o["name"],"tile_sha256":t["slice_sha256"],"center_mhz":t["center_mhz"],"trial":trial,"hard_control":bool(hard),"start":start,"label":label,"amplitude":amp if label else 0})
    for h in data.get("candidates",data.get("forward",[]))+data.get("reverse",[]):
        if h["status"]=="CONTROL_TRACK_PRESENT_RFI_CONCERN" or h.get("cadence")=="OFF_ON_OFF":
            library.append({"source_run":a["run"],"candidate":h,"label":"RFI_CONCERN_OR_REVERSE_CONTROL_NOT_CONFIRMED_EMITTER"})
    z.close()
    return np.asarray(X,np.float32),np.asarray(Y,np.int64),rows,targets,library


# JUXTAPOSE sources: original 250x250 tournament + cached exact receipts, ESSENCE v3.6 exact stop.
# Re-express pivot units for spectral detection; CNF elimination is not applicable here.
PROTO=pathlib.Path("research/seti/juxtapose_r18a_protocol.json")
P=json.loads(PROTO.read_text())
TH=float(P["threshold"])
ASC=tuple(range(-32,33))
CENTER=tuple(sorted(ASC,key=lambda s:(abs(s),s)))
GROUPS=tuple(tuple(int(v) for v in a) for a in np.array_split(CENTER,8))
assert len(GROUPS)==8 and sorted(s for g in GROUPS for s in g)==list(ASC)
def score(v,s):
    n,w=v.shape
    ds=np.rint(np.arange(n)*s/(n-1)).astype(int)
    lo=max(0,-int(ds.min()));hi=min(w,w-int(ds.max()))
    vals=sum(v[i,lo+k:hi+k] for i,k in enumerate(ds))/np.sqrt(n)
    j=int(np.argmax(vals))
    return float(vals[j]),int(lo+j)
def witness(v,s,b):
    ds=np.rint(np.arange(len(v))*s/(len(v)-1)).astype(int)
    return float(sum(float(v[i,b+k]) for i,k in enumerate(ds))/np.sqrt(len(v)))
def cache_cost(train):
    rows=[]
    for p in train:
        v=normalized(p)
        rows.append([next(( (i+1) for i,s in enumerate(g) if score(v,s)[0]>=TH),len(g)+1) for g in GROUPS])
    return np.array(rows,np.int16)
def tournament(cache):
    # All 40320 order replays on cached exact group receipts.
    # Every negative pays full cost; positives stop at first verified group hit.
    orders=list(itertools.permutations(range(8)));assert len(orders)==40320
    sizes=np.array([len(g) for g in GROUPS],np.int16)
    best=None;results=[]
    for off in range(0,len(orders),256):
        chunk=np.asarray(orders[off:off+256],np.int16)
        active=np.ones((len(chunk),len(cache)),bool)
        costs=np.zeros(active.shape,np.int16)
        for j in range(8):
            g=chunk[:,j];n=sizes[g,None];c=cache[:,g].T
            costs+=active*np.minimum(c,n)
            active&=c>n
        for perm,total,peak in zip(chunk,costs.sum(1),costs.max(1)):
            key=(int(total),int(peak),tuple(int(x) for x in perm))
            results.append({"order":list(key[2]),"total_checks":key[0],"peak_checks":key[1]})
            if best is None or key<best:best=key
    return best,results
def decision(p,order,full=False):
    v=normalized(p);hit=None;count=0
    for s in order:
        a,b=score(v,s);count+=1
        if a>=TH and hit is None:
            # independent scalar reconstruction; tolerate FP rounding, do not invent threshold hit.
            w=witness(v,s,b)
            assert abs(a-w)<=max(1e-5,abs(a)*2e-6)
            hit=(s,b,a)
            if not full:break
    return bool(hit),count,hit
def measure(X,y,order,full=False):
    runs=[];last=None
    for _ in range(3):
        start=time.perf_counter();res=[decision(p,order,full) for p in X];runs.append(time.perf_counter()-start)
        signature=[(r[0],r[1]) for r in res]
        if last is not None:assert signature==last
        last=signature
    pred=np.asarray([r[0] for r in res]);counts=np.asarray([r[1] for r in res])
    m={"n":len(y),"checks":int(counts.sum()),"mean_checks":float(counts.mean()),"positive_mean_checks":float(counts[pred].mean()) if pred.any() else None,"negative_mean_checks":float(counts[~pred].mean()) if (~pred).any() else None,"recall":float(pred[y==1].mean()),"control_alert_fraction":float(pred[y==0].mean()),"median_seconds":statistics.median(runs),"timing_repeats":runs}
    return m,pred
def main():
    t=time.perf_counter();X,y,meta,targets,_=build(S["source_artifacts"][0]);load_train=time.perf_counter()-t
    t=time.perf_counter();cache=cache_cost(X);cache_sec=time.perf_counter()-t
    t=time.perf_counter();champ,orders=tournament(cache);tourney_sec=time.perf_counter()-t
    ranking=champ[2];order=tuple(s for g in ranking for s in GROUPS[g])
    lock={"protocol_sha256":hashlib.sha256(PROTO.read_bytes()).hexdigest(),"groups":GROUPS,"champion":ranking,"ordered_shifts":order,"train_targets":sorted(targets),"train_samples":len(y),"orders_replayed":len(orders),"training_checks_upper_bound":len(y)*65,"tournament_seconds":tourney_sec,"cache_seconds":cache_sec,"train_load_and_patch_seconds":load_train}
    (OUT/"frozen_rank.json").write_text(json.dumps(lock,indent=2)+"\n");(OUT/"all_order_scores.json").write_text(json.dumps(orders)+"\n");emit("JUXTAPOSE_LOCK",lock)
    t=time.perf_counter();T,Y,tm,tt,_=build(S["source_artifacts"][1]);load_test=time.perf_counter()-t
    assert not targets&tt and not({r["tile_sha256"] for r in meta}&{r["tile_sha256"] for r in tm})
    results={};pred={}
    for name,o,full in [("full_exact",ASC,True),("ascending_early",ASC,False),("center_out_early",CENTER,False),("juxtapose",order,False)]:
        results[name],pred[name]=measure(T,Y,o,full)
        assert np.array_equal(pred[name],pred["full_exact"])
        emit("JUXTAPOSE_METHOD",{"method":name,**results[name]})
    for name,m in results.items():
        m["exact_decision_equal"]=True
        m["full_search_check_reduction"]=1-m["checks"]/results["full_exact"]["checks"]
        m["full_search_runtime_ratio"]=results["full_exact"]["median_seconds"]/m["median_seconds"]
    fast=min(("ascending_early","center_out_early"),key=lambda k:results[k]["median_seconds"])
    delta=(results[fast]["median_seconds"]-results["juxtapose"]["median_seconds"])/len(Y)
    cost=load_train+cache_sec+tourney_sec
    summary={"lock":lock,"test_targets":sorted(tt),"test_samples":len(Y),"methods":results,"best_simple_early_baseline":fast,"training_total_seconds":cost,"test_load_and_patch_seconds":load_test,"break_even_samples_vs_best_early":math.ceil(cost/delta) if delta>0 else None,"decision_mismatches":0,"scope":"Exact finite threshold existence on 128-channel 16-row patches, no persistence requirement; no guarantee of same strongest track or top-k if stopping early. Synthetic injections in real OFF data; R11 already viewed, development only. No independent source classification or production change.","production_adoption":False}
    (OUT/"result.json").write_text(json.dumps(summary,indent=2)+"\n")
    emit("JUXTAPOSE_SUMMARY",summary)
if __name__=="__main__":main()
