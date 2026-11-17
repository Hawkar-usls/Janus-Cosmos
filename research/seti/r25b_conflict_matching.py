#!/usr/bin/env python3
"""R25B: maximum-cardinality/minimum-distance matching, factored by exact components."""
import json,pathlib,time,statistics,itertools
import numpy as np
from scipy.optimize import linear_sum_assignment
import r20b_join_dispatch as archive
P=pathlib.Path("research/seti/r25b_matching_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r25b");OUT.mkdir(parents=True,exist_ok=True)
FT=S["frequency_tolerance_hz"];DT=S["drift_tolerance_hz_s"]
def matrix(A,B):
    d=np.abs(A[:,None,:]-B[None,:,:]);ok=(d[:,:,0]<=FT)&(d[:,:,1]<=DT)
    cost=np.rint((d[:,:,0]/FT+d[:,:,1]/DT)*1e6).astype(np.int64)
    return ok,cost
def exact(ok,c):
    n,m=ok.shape
    if not n or not m:return []
    # one unmatched dummy per left node. Penalty exceeds every feasible total distance.
    M=(min(n,m)+1)*2000001
    C=np.full((n,m+n),M,dtype=np.int64);C[:,:m]=np.where(ok,c,2*M)
    a,b=linear_sum_assignment(C)
    return [(int(i),int(j)) for i,j in zip(a,b) if j<m and ok[i,j]]
def objective(q,c):return (len(q),sum(int(c[i,j]) for i,j in q))
def components(ok):
    n,m=ok.shape;parent=list(range(n+m))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in np.argwhere(ok):
        a,b=root(int(i)),root(n+int(j));parent[b]=a
    groups={}
    for i in range(n+m):groups.setdefault(root(i),[]).append(i)
    return [(np.array([i for i in g if i<n],int),np.array([i-n for i in g if i>=n],int)) for g in groups.values() if any(i<n for i in g) and any(i>=n for i in g)]
def factored(ok,c):
    groups=components(ok);r=[]
    for a,b in groups:r.extend((int(a[i]),int(b[j])) for i,j in exact(ok[np.ix_(a,b)],c[np.ix_(a,b)]))
    return sorted(r),len(groups)
def greedy(ok,c):
    q=[];a=set();b=set()
    for i,j in sorted(map(tuple,np.argwhere(ok)),key=lambda k:(c[k],k)):
        if i not in a and j not in b:q.append((int(i),int(j)));a.add(i);b.add(j)
    return q
def verify(q,ok):
    assert len({i for i,j in q})==len({j for i,j in q})==len(q) and all(ok[i,j] for i,j in q)
def brute(ok,c):
    n,m=ok.shape;best=(0,0)
    def walk(i,used,k,s):
        nonlocal best
        if i==n:
            if (k,-s)>(best[0],-best[1]):best=(k,s)
            return
        walk(i+1,used,k,s)
        for j in np.flatnonzero(ok[i]):
            if int(j) not in used:walk(i+1,used|{int(j)},k+1,s+int(c[i,j]))
    walk(0,set(),0,0);return best
def selftest():
    rng=np.random.default_rng(25020)
    for i in range(64):
        n=int(rng.integers(1,7));m=int(rng.integers(1,7));o=rng.random((n,m))<.5;c=rng.integers(0,10,(n,m))
        a=exact(o,c);b,_=factored(o,c);verify(a,o);verify(b,o);assert objective(a,c)==objective(b,c)==brute(o,c)
    o=np.array([[1,1],[1,0]],bool);c=np.array([[0,1],[1,0]])
    assert len(exact(o,c))==2 and len(greedy(o,c))==1
    return {"brute_force_cases":64,"greedy_counterexample":"PASS"}
def bench(A,B,truth=None):
    ok,c=matrix(A,B);out={};qs={}
    for k in ("global","factored","greedy"):
        ts=[]
        for _ in range(3):
            t=time.perf_counter();o,d=matrix(A,B)
            if k=="global":q=exact(o,d)
            elif k=="factored":q,_=factored(o,d)
            else:q=greedy(o,d)
            ts.append(time.perf_counter()-t)
        verify(q,ok);qs[k]=q;o=objective(q,c)
        out[k]={"pairs":o[0],"distance_cost":o[1],"wall":statistics.median(ts)}
        if truth is not None:
            tp=len(set(q)&truth);out[k].update(tp=tp,fp=len(q)-tp,fn=len(truth)-tp)
    assert objective(qs["global"],c)==objective(qs["factored"],c)
    # Alternate optimal assignment diagnostic; not an identity certificate.
    ambiguous=0;checked=0
    if max(ok.shape,default=0)<=32:
        for i,j in qs["global"]:
            o=ok.copy();o[i,j]=False;ambiguous+=int(objective(exact(o,c),c)==objective(qs["global"],c));checked+=1
    return {"shape":list(ok.shape),"compatible_edges":int(ok.sum()),"components":len(components(ok)),"alternative_optimal_selected_edges":ambiguous,"alternative_edges_checked":checked,"modes":out,"all_compatible_edges":np.argwhere(ok).tolist() if truth is None else None,"selected_pairs":qs if truth is None else None}
def synthetic(rng,n,fam):
    span={"sparse":n*1000,"crowded":n*18,"dropout":n*100,"clutter":n*35}[fam]
    A=np.column_stack((rng.uniform(0,span,n),rng.uniform(-.2,.2,n)))
    idx=np.arange(n)
    if fam=="dropout":idx=idx[rng.random(n)>.35]
    B=A[idx]+np.column_stack((rng.uniform(-18,18,len(idx)),rng.uniform(-.022,.022,len(idx))))
    ids=list(idx)
    if fam in ("dropout","clutter","crowded"):
        extra=max(1,n//2);B=np.concatenate((B,np.column_stack((rng.uniform(0,span,extra),rng.uniform(-.2,.2,extra)))));ids+=[-1]*extra
    perm=rng.permutation(len(B));B=B[perm];truth={(int(ids[j]),i) for i,j in enumerate(perm) if ids[j]>=0}
    return A,B,truth
def main():
    check=selftest();rng=np.random.default_rng(S["seed"]);rows=[]
    for n,fam,rep in itertools.product(S["sizes"],S["families"],range(S["replicates"])):
        A,B,T=synthetic(rng,n,fam);rows.append({"n":n,"family":fam,"rep":rep,**bench(A,B,T)})
    records=archive.load();r=[]
    for f in sorted({x["payload"]["window_mhz"] for x in records}):
        for cad in ("ON","OFF"):
            for i,j in ((1,2),(2,3),(1,3)):
                data=[]
                for k in (i,j):
                    rr=[x for x in records if x["payload"]["window_mhz"]==f and x["payload"]["role"]==cad+str(k)]
                    data.append(rr)
                A,B=[np.array([[x["payload"]["projected_hz"],x["payload"]["drift"]] for x in rr],float).reshape(-1,2) for rr in data]
                if not len(A) or not len(B):continue
                r.append({"frequency_mhz":f,"roles":[cad+str(i),cad+str(j)],"record_ids":[[x["record_id"] for x in rr] for rr in data],**bench(A,B)})
    def total(a,synthetic):
        keys=["pairs","distance_cost","wall"]+(["tp","fp","fn"] if synthetic else [])
        return {"cases":len(a),"edges":sum(x["compatible_edges"] for x in a),"alternative_edges_checked":sum(x["alternative_edges_checked"] for x in a),"alternative_optimal_selected_edges":sum(x["alternative_optimal_selected_edges"] for x in a),"modes":{k:{s:sum(x["modes"][k][s] for x in a) for s in keys} for k in ("global","factored","greedy")}}
    d={"status":"RESEARCH_ONLY","selftests":check,"synthetic":total(rows,True),"by_family":{k:total([x for x in rows if x["family"]==k],True) for k in S["families"]},"archive":total(r,False),"objective_mismatches":0,"rows":rows,"archive_rows":r,"scope":"Two-scan one-to-one assignment under fixed compatibility. True identity is not established by optimal objective. No triple/source inference; overlapping pairs across pairwise comparisons not independent. Raw edges retained.","production_adoption":False}
    (OUT/"result.json").write_text(json.dumps(d,indent=2)+"\n");print("R25B_SUMMARY",json.dumps({k:v for k,v in d.items() if k not in ("rows","archive_rows")}),flush=True)
if __name__=="__main__":main()
