#!/usr/bin/env python3
"""R27A original Taylor core vs native JANUS, common CPU preprocessing."""
import json,pathlib,hashlib,sys,types,time
import numpy as np
import requests
import r26c_reference_alignment as r
import r24_odontocrystal as od
import r17c_frozen_holdout as hold
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r27a_benchmark_protocol.json");S=json.loads(P.read_text());OUT=pathlib.Path("runtime/r27a");OUT.mkdir(parents=True,exist_ok=True)
r.S=S;r.OUT=OUT;od.S=S;hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))};hold.OUT=OUT;hold.P=P;eng.OUT=OUT
def emit(k,d):print(k,json.dumps(d),flush=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def load_tree():
    # Exact upstream sources; only their dependency namespaces are isolated.
    spec=S["comparators"]["turbo"];root="https://raw.githubusercontent.com/UCBerkeleySETI/turbo_seti/"+spec["ref"]+"/"
    for name,(path,sha) in zip(("turbo_seti.find_doppler.kernels._bitrev","_janus_external_taylor"),spec["files"]):
        a=requests.get(root+path,timeout=30);a.raise_for_status();b=a.content
        assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()==sha
        m=types.ModuleType(name);sys.modules[name]=m;exec(compile(b,path,"exec"),m.__dict__)
    return sys.modules["_janus_external_taylor"].flt,sys.modules["turbo_seti.find_doppler.kernels._bitrev"].bitrev
def init_tree():
    global flt,bitrev,paths
    t=time.perf_counter();flt,bitrev=load_tree();n=16;stride=288;probe=90
    paths=np.zeros((n,n),int)
    for i in range(n):
        a=np.zeros((n,stride),np.float64);a[i]=np.arange(stride)
        b=flt(a.ravel(),n).reshape(n,stride)
        paths[:,i]=np.rint([b[bitrev(k,4),probe]-probe for k in range(n)]).astype(int)
    assert np.array_equal(paths[:,0],np.zeros(n,int)) and np.array_equal(paths[:,-1],np.arange(n))
    assert np.max(abs(np.diff(paths,axis=1)))<=1
    return time.perf_counter()-t
def tree(V):
    B,n,w=V.shape;assert n==16;pad=64;stride=w+2*pad+2*n;t=np.arange(n);best=np.full(B,-np.inf);pp=np.zeros((B,n),int)
    for b in (-3,-2,-1,0,1,2):
        for i,x in enumerate(V):
            col=np.arange(stride)[None,:]-pad+b*t[:,None];ok=(col>=0)&(col<w)
            a=np.where(ok,x[t[:,None],np.clip(col,0,w-1)],0).astype(np.float64)
            z=flt(a.ravel(),n).reshape(n,stride)
            for k in range(n):
                d=paths[k]+b*t
                if abs(d[-1])>32:continue
                lo=max(0,-int(d.min()));hi=min(w,w-int(d.max()))
                zz=z[bitrev(k,4),pad+lo:pad+hi]/4
                if not len(zz):continue
                j=int(zz.argmax());sc=float(zz[j])
                if sc>best[i]:best[i]=sc;pp[i]=lo+j+d
    replay=V[np.arange(B)[:,None],t[None,:],pp].sum(1)/4
    assert np.allclose(replay,best,atol=1e-8,rtol=1e-10)
    return best,pp
def straight(V):
    B,n,w=V.shape;t=np.arange(n);best=np.full(B,-np.inf);p=np.zeros((B,n),int)
    for sh in range(-32,33):
        d=np.rint(t*sh/(n-1)).astype(int);lo=max(0,-int(d.min()));hi=min(w,w-int(d.max()))
        v=sum(V[:,i,lo+k:hi+k] for i,k in enumerate(d))/4;j=v.argmax(1);q=v[np.arange(B),j];m=q>best
        best[m]=q[m];p[m]=lo+j[m,None]+d
    return best,p
def arm(X,R,k):
    t=time.perf_counter();V=r.norm(X)
    if k=="linear_reference":s,p=straight(V)
    elif k=="turbo_tree":s,p=tree(V)
    else:
        if k=="janus_registered":V=r.register(V,r.norm(R))[0]
        s,p,_=od.sliced(V,0,None)
    return (s,p,time.perf_counter()-t)
BASE=["janus_exact","janus_registered","linear_reference","turbo_tree"]
def evaluate(rows,R):
    X=np.stack([x[0] for x in rows]);a={k:arm(X,R,k) for k in BASE}
    return a
def joint(a,L):
    s,p,t=a["janus_exact"];v,q,w=a["linear_reference"]
    z=(s-L["pool"]["janus_exact"])/L["scale"]["janus_exact"];u=(v-L["pool"]["linear_reference"])/L["scale"]["linear_reference"]
    return np.maximum(z,u),np.where((z>=u)[:,None],p,q),t+w
def cut(x):return r.cut(x)
def calibrate(cal):
    y=np.concatenate([np.array([x[1] for x in rows]) for rows,a in cal]);pool={};scale={}
    for k in BASE:
        s=np.concatenate([a[k][0] for rows,a in cal]);pool[k]=cut(s[y==0]);scale[k]=max(.001,float(np.std(s[y==0])))
    L={"pool":pool,"scale":scale};z=np.concatenate([joint(a,L)[0] for rows,a in cal]);pool["janus_joint"]=cut(z[y==0]);return L
def met(rows,a,L):
    a["janus_joint"]=joint(a,L);out={}
    for k,(s,p,t) in a.items():
        pred=s>=L["pool"][k];m=np.zeros_like(p,bool);out[k]={**od.measure(pred,p,m,rows),"seconds":t}
    return out
def selftest():
    rng=np.random.default_rng(2701);V=rng.normal(size=(8,16,128));s,p=tree(V);best=np.full(len(V),-np.inf)
    for b in (-3,-2,-1,0,1,2):
        for d in paths+b*np.arange(16):
            if abs(d[-1])>32:continue
            lo=max(0,-int(d.min()));hi=min(128,128-int(d.max()))
            z=sum(V[:,i,lo+k:hi+k] for i,k in enumerate(d))/4;best=np.maximum(best,z.max(1))
    assert np.allclose(best,s,atol=1e-10)
    return {"upstream_source_hashes":"PASS","tree_map_exhaustive_8":"PASS","tree_witness":"PASS"}
def main():
    jit=init_tree();checks=selftest();emit("R27A_SELFTEST",checks)
    m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest() and len(m["selected"])==6
    cal=[];seen=set();load=0
    for q in m["selected"][:2]:
        t=time.perf_counter();rows,R,_,h=r.gather(q);load+=time.perf_counter()-t;assert not (h&seen);seen|=h;cal.append((rows,evaluate(rows,R)))
    L=calibrate(cal);L["targets"]=[q["target"] for q in m["selected"][:2]];save("lock.json",L);emit("R27A_LOCK",L)
    out={}
    for q in m["selected"][2:]:
        t=time.perf_counter();rows,R,_,h=r.gather(q);load+=time.perf_counter()-t;assert not (h&seen);seen|=h
        out[q["target"]]=met(rows,evaluate(rows,R),L);save("partial.json",out);emit("R27A_TARGET",{"target":q["target"],"metrics":out[q["target"]]})
    keys=("detected","localized","positives","control_alerts","controls","seconds");total={k:{a:sum(x[k][a] for x in out.values()) for a in keys} for k in S["modes"]}
    for k in total:
        for f,vals in (("amplitudes",S["amplitudes"]),("gaps",S["gap_levels"]),("models",S["injection_models"])):
            total[k][f]={str(v):{a:sum(x[k][f][str(v)][a] for x in out.values()) for a in keys[:5]} for v in vals}
    d={"status":"FRESH_ARCHIVE_KERNEL_COMPARISON","targets":list(out),"selftests":checks,"lock":L,"modes":total,"per_target":out,"samples":1024,"external_fetch_and_jit_seconds":jit,"load_seconds_all_targets":load,"scope":S["scope"],"production_adoption":False};save("result.json",d);emit("R27A_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="selftest":emit("R27A_SELFTEST",{"init_seconds":init_tree(),"checks":selftest()})
    else:main()
