#!/usr/bin/env python3
"""Exact triangle join: cached common supports and a frozen cost selector."""
import hashlib,io,json,os,pathlib,statistics,time,zipfile
import numpy as np
import requests
import r19bcd_evidence_relations as old
P=pathlib.Path("research/seti/r20b_join_dispatch_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r20b");OUT.mkdir(parents=True,exist_ok=True)
FT=S["thresholds"]["frequency"];DT=S["thresholds"]["drift"]
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def emit(k,d):print(k,json.dumps(d),flush=True)
def load():
    a=S["source_artifact"];b=bytearray()
    with requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},stream=True,timeout=(15,60)) as r:
        r.raise_for_status()
        for c in r.iter_content(65536):b.extend(c);assert len(b)<8*1024**2
    assert hashlib.sha256(b).hexdigest()==a["digest"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        names=[n for n in z.namelist() if pathlib.PurePosixPath(n).name=="rfi_catalog.json"]
        assert len(names)==1 and z.getinfo(names[0]).file_size<8*1024**2
        return json.loads(z.read(names[0]))["records"]
def dense(M):
    x,y,z=M
    return np.argwhere(x[:,:,None]&y[:,None,:]&z[None,:,:])
def supports(M):
    x,y,z=M
    X=x&z.any(axis=1)[None,:];Y=y&z.any(axis=0)[None,:]
    bx=X.sum(axis=1,dtype=np.int64);cy=Y.sum(axis=1,dtype=np.int64)
    work=int(np.dot(bx,cy))
    return X,Y,work
def cached(M,sup=None):
    x,y,z=M;X,Y,work=supports(M) if sup is None else sup
    rows=[]
    for i in np.flatnonzero(X.any(axis=1)&Y.any(axis=1)):
        b=np.flatnonzero(X[i]);c=np.flatnonzero(Y[i])
        k,l=np.nonzero(z[np.ix_(b,c)])
        if len(k):rows.append(np.column_stack((np.full(len(k),i),b[k],c[l])))
    return (np.concatenate(rows,axis=0) if rows else np.empty((0,3),dtype=np.int64)),work
def route(M,policy):
    x,y,z=M;cube=x.shape[0]*x.shape[1]*y.shape[1]
    if cube<=policy["min_cube"]:return dense(M),"dense",cube
    sup=supports(M)
    if sup[2]<=cube*policy["max_ratio"]:
        a,w=cached(M,sup);return a,"cached",w
    return dense(M),"dense",cube
def run(A,method,policy=None):
    count=int(np.prod([len(a) for a in A]))
    if method=="original_semijoin":
        a,w=old.semi(A,FT,DT);return a,"original_semijoin",w
    M=old.mats(A,FT,DT)
    if method=="dense":return dense(M),"dense",count
    if method=="cached_semijoin":
        a,w=cached(M);return a,"cached",w
    return route(M,policy)
def datasets(seed,reps):
    rng=np.random.default_rng(seed)
    for n in S["sizes"]:
        for family in S["families"]:
            for rep in range(reps):
                # Unequal dimensions challenge actual catalog shape; dense case deliberately includes many witnesses.
                sizes=[n,max(2,n-rng.integers(0,max(1,n//4))), max(2,n-rng.integers(0,max(1,n//3)))]
                span={"sparse":100000,"medium":2000,"dense":150,"hostile":1000}[family]
                A=[np.column_stack((rng.uniform(0,span,k),rng.normal(0,.1,k))) for k in sizes]
                if family=="hostile":A[2][:,0]+=5000
                else:
                    for a in A:a[0]=[10,.001]
                yield {"n":n,"family":family,"rep":rep},A
def measure(A,method,policy,reference,repeats):
    ts=[];cp=[];chosen=None;work=None
    for _ in range(repeats):
        t=time.perf_counter();cpu=time.process_time();a,chosen,work=run(A,method,policy)
        ts.append(time.perf_counter()-t);cp.append(time.process_time()-cpu)
        assert np.array_equal(a,reference),method
    return {"wall":statistics.median(ts),"cpu":statistics.median(cp),"route":chosen,"expansion":work}
def calibrate():
    T=time.perf_counter()
    policies=[{"min_cube":n,"max_ratio":r} for n in S["router_grid"]["min_cube"] for r in S["router_grid"]["max_ratio"]]
    sums=[0.0]*len(policies);cases=0
    for meta,A in datasets(S["cal_seed"],S["cal_reps"]):
        ref=run(A,"dense")[0];cases+=1
        order=np.random.default_rng(S["cal_seed"]+cases).permutation(len(policies))
        for i in order:sums[i]+=measure(A,"adaptive",policies[i],ref,2)["wall"]
    best=min(range(len(policies)),key=lambda i:(sums[i],i))
    return {"policy":policies[best],"cases":cases,"tuning":[{"policy":p,"wall":s} for p,s in zip(policies,sums)],"calibration_seconds":time.perf_counter()-T}
def bench(meta,A,policy):
    ref=run(A,"dense")[0];cube=int(np.prod([len(a) for a in A]))
    # Independent brute check where manageable, including empty and boundary cases.
    if cube<=4096:
        r=[]
        for i,a in enumerate(A[0]):
            for j,b in enumerate(A[1]):
                for k,c in enumerate(A[2]):
                    if max(a[0],b[0],c[0])-min(a[0],b[0],c[0])<=FT and max(a[1],b[1],c[1])-min(a[1],b[1],c[1])<=DT:r.append((i,j,k))
        assert np.array_equal(np.asarray(r,dtype=np.int64).reshape(-1,3),ref)
    modes=list(S["modes"])
    # Randomize method order per case without result-dependent choices.
    seed=int.from_bytes(hashlib.sha256(json.dumps(meta,sort_keys=True).encode()).digest()[:4],"big")
    np.random.default_rng(seed).shuffle(modes)
    return {**meta,"cube":cube,"matches":len(ref),"modes":{k:measure(A,k,policy,ref,S["repeats"]) for k in modes}}
def real(records,policy):
    rows=[]
    for f in sorted(set(r["payload"]["window_mhz"] for r in records)):
        for cadence in ("ON","OFF"):
            A=[]
            for i in (1,2,3):
                a=[[r["payload"]["projected_hz"],r["payload"]["drift"]] for r in records if r["payload"]["window_mhz"]==f and r["payload"]["role"]==cadence+str(i)]
                A.append(np.asarray(a,dtype=float).reshape(-1,2))
            rows.append(bench({"family":"archive","center_mhz":f,"cadence":cadence},A,policy))
    return rows
def total(rows):
    return {"cases":len(rows),"cube":sum(r["cube"] for r in rows),"matches":sum(r["matches"] for r in rows),"modes":{k:{"wall":sum(r["modes"][k]["wall"] for r in rows),"cpu":sum(r["modes"][k]["cpu"] for r in rows),"expansion":sum(r["modes"][k]["expansion"] for r in rows),"routes":{p:sum(r["modes"][k]["route"]==p for r in rows) for p in ("dense","cached","original_semijoin")}} for k in S["modes"]}}
def main():
    lock=calibrate();save("lock.json",lock);emit("R20B_LOCK",lock)
    rows=[bench(m,A,lock["policy"]) for m,A in datasets(S["test_seed"],S["test_reps"])]
    t=time.perf_counter();rec=load();loadsec=time.perf_counter()-t
    replay=real(rec,lock["policy"])
    d={"status":"EXACT_MATCH","mismatches":0,"lock":lock,"holdout":total(rows),"by_family":{k:total([x for x in rows if x["family"]==k]) for k in S["families"]},"archive":total(replay),"archive_load_seconds":loadsec,"rows":rows,"archive_rows":replay,"scope":"Same declared tolerance predicates only; candidate catalog archive is a diagnostic replay, not new observation. Calibration cost separate. All cases include actual routing cost.","production_adoption":False}
    save("result.json",d);emit("R20B_SUMMARY",{k:v for k,v in d.items() if k not in ("rows","archive_rows","lock")})
if __name__=="__main__":main()
