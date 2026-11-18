#!/usr/bin/env python3
"""R35F: same-objective native JIT acceleration; no sensitivity/threshold change."""
import hashlib,json,pathlib,time,statistics,resource,io,zipfile,os
import numpy as np
from numba import njit
import requests
import r35b_grand_benchmark as ref
OUT=pathlib.Path("runtime/r35f");OUT.mkdir(parents=True,exist_ok=True)
ART={"id":14405617460,"sha256":"fa60d60f8c961d521fe9bfb0f915447b327a05f03c5756d5b743a965ce46a9b0"}
D=ref.D.astype(np.int64)
@njit(cache=False,fastmath=False)
def kernel(x,D,guard):
    B,n,w=x.shape
    best=np.full(B,-np.inf,np.float64);base=np.zeros(B,np.int64);which=np.zeros(B,np.int64)
    # One contiguous scratch vector, no full drift-plane allocation, exact order of additions.
    acc=np.empty(w,np.float64)
    for b in range(B):
        for k in range(len(D)):
            lo=max(0,-D[k].min())+guard;hi=min(w,w-D[k].max())-guard
            if lo>=hi:continue
            for j in range(lo,hi):acc[j]=0.
            for i in range(n):
                d=D[k,i]
                for j in range(lo,hi):acc[j]+=x[b,i,j+d]
            for j in range(lo,hi):
                s=acc[j]/4.
                if s>best[b]:
                    best[b]=s;base[b]=j;which[b]=k
    p=np.empty((B,n),np.int64)
    for b in range(B):
        for i in range(n):p[b,i]=base[b]+D[which[b],i]
    return best,p
def native(V,kind):
    x=np.ascontiguousarray(V,np.float64);guard=0
    if kind=="clip":x=np.clip(x,-3,3)
    if kind=="neighborhood":
        v=x;x=.5*v.copy();x[:,:,1:]+=.25*v[:,:,:-1];x[:,:,:-1]+=.25*v[:,:,1:];guard=1
    return kernel(x,D,guard)
def measure(V,kind):
    a,p=ref.template(V,kind);b,q=native(V,kind)
    err=float(np.max(abs(a-b)));m=int(np.sum(np.any(p!=q,axis=1)))
    assert np.array_equal(a,b) and m==0,("exactness",kind,err,m)
    ts=[];ns=[]
    for i in range(3):
        jobs=[(ref.template,ts),(native,ns)]
        if i%2:jobs.reverse()
        for f,l in jobs:
            t=time.perf_counter();x,y=f(V,kind);l.append(time.perf_counter()-t)
            assert np.array_equal(x,a) and np.array_equal(y,p)
    return {"examples":len(V),"score_exact":True,"path_mismatches":m,"max_abs_error":err,"numpy_median_seconds":statistics.median(ts),"jit_median_seconds":statistics.median(ns),"numpy_repeats":ts,"jit_repeats":ns}
def save(n,v):(OUT/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def main():
    rng=np.random.default_rng(3551);x=rng.normal(size=(4,16,128));x[1]=0.;x[2]=1.;x[3]=3
    t=time.perf_counter();native(x,"raw");warm=time.perf_counter()-t
    for kind in ("raw","clip","neighborhood"):
        a,p=ref.template(x,kind);v,q=native(x,kind);assert np.array_equal(a,v) and np.array_equal(p,q)
    u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{ART['id']}/zip"
    r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=60);r.raise_for_status();b=r.content
    assert len(b)<64*1024**2 and hashlib.sha256(b).hexdigest()==ART["sha256"]
    z=zipfile.ZipFile(io.BytesIO(b));assert len(z.namelist())==6 and sum(i.file_size for i in z.infolist())<64*1024**2
    rows=[]
    for f in sorted(k for k in z.namelist() if k.startswith("eval_")):
        with np.load(io.BytesIO(z.read(f)),allow_pickle=False) as d:
            V=d["V"].astype(np.float64)
            for kind in ("raw","clip","neighborhood"):
                row={"target":f,"kind":kind,**measure(V,kind)};rows.append(row);print("R35F_TARGET",json.dumps(row),flush=True)
                save("partial.json",rows)
    z.close();assert len(rows)==12
    a=sum(r["numpy_median_seconds"] for r in rows);b=sum(r["jit_median_seconds"] for r in rows)
    d={"status":"PASS_EXACT_OBJECTIVE_REPLAY" if b<a else "HOLD_NO_SPEED_GAIN","exam_examples":1024,"template_bank":len(D),"evaluations":sum(r["examples"] for r in rows),"score_mismatches":0,"path_mismatches":0,"numpy_seconds":a,"jit_seconds":b,"warm_speedup":a/b,"jit_compile_seconds":warm,"including_compile_speedup":a/(b+warm),"peak_rss_kib_process":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"per_target":rows,"scope":"Historical exact numerical engineering benchmark only, no independent scientific gain. Same CPU, one thread, 3 repeats; no production changes.","production_adoption":False}
    save("result.json",d);print("R35F_SUMMARY",json.dumps({k:v for k,v in d.items() if k!="per_target"}),flush=True)
if __name__=="__main__":main()
