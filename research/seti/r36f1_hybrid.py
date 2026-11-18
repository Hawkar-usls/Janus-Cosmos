#!/usr/bin/env python3
"""R36F1: reproduce frozen R28 registered detector before proposing optimized implementation."""
import concurrent.futures,hashlib,json,pathlib,resource,time
import numpy as np
import r36_galaxy_swarm as g
import r28b_router as rb
P=pathlib.Path("research/seti/r36f1_hybrid_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r36f1");OUT.mkdir(parents=True,exist_ok=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def lin(V):return rb.a.straight(V)
def od(V):return rb.od.sliced(V,0,None)[:2]
def get(r,k):
    ids=r[1][:,k];st=r[2][:,k];assert np.all(ids>=0)
    return r[0][:,k],g.D[ids]+st[:,None]
def same(a,b):
    assert np.array_equal(a[0],b[0]) and np.array_equal(a[1],b[1]),"score/witness mismatch"
def cut(v):
    x=np.sort(v)[::-1];return float(np.nextafter(x[int(len(x)*.05)],np.inf))
def fuse(a,b,L):
    s,p=a;z,q=b
    u=(s-L["linear_threshold"])/L["linear_scale"];v=(z-L["odonto_threshold"])/L["odonto_scale"]
    m=v>u;return np.where(m,v,u),np.where(m[:,None],q,p)
def inputs():
    d=g.fetch();assert S["source"]["input_artifact"]==g.S["source"]["artifact"]
    return [(n,q,h,np.ascontiguousarray(q["V"],dtype=np.float64)) for n,q,h in d]
def main():
    t=time.perf_counter();data=inputs();fetch=time.perf_counter()-t
    x=data[0][3][:2];t=time.perf_counter();cpu=time.process_time()
    r=g.scan(x,g.D,g.R,g.L,False);same(lin(x),get(r,0))
    compile={"wall":time.perf_counter()-t,"cpu":time.process_time()-cpu}
    t=time.perf_counter();cal=[(lin(V),od(V),q["y"]) for n,q,h,V in data if n.startswith("cal_")];fit_time=time.perf_counter()-t
    y=np.concatenate([v[2] for v in cal]);L={}
    for k,i in (("linear",0),("odonto",1)):
        z=np.concatenate([v[i][0] for v in cal])[y==0];th=cut(z);L[k+"_threshold"]=th;L[k+"_scale"]=max(.1,th-float(np.quantile(z,.8)))
    th=cut(np.concatenate([fuse(a,b,L)[0] for a,b,y in cal])[y==0]);L["hybrid_threshold"]=th
    save("lock.json",L);lh=g.sha((OUT/"lock.json").read_bytes());rows=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for n,q,h,V in data:
            if not n.startswith("eval_"):continue
            def original():return fuse(lin(V),od(V),L)
            def compiled():return fuse(get(g.scan(V,g.D,g.R,g.L,False),0),od(V),L)
            def twin():
                a=pool.submit(g.scan,V,g.D,g.R,g.L,False);b=pool.submit(od,V);return fuse(get(a.result(),0),b.result(),L)
            def duplicate():
                a=g.separate(V);return fuse(get(a,0),od(V),L),get(a,1)
            def shared():
                a=pool.submit(g.scan,V,g.D,g.R,g.ALL,False);b=pool.submit(od,V);r=a.result();return fuse(get(r,0),b.result(),L),get(r,1)
            fs={"original":original,"compiled":compiled,"twin":twin,"duplicate_aux":duplicate,"shared_aux":shared}
            out={};times={k:[] for k in fs}
            for i in range(3):
                for k in list(fs)[i:]+list(fs)[:i]:
                    out[k],t=g.timed(fs[k]);times[k].append(t)
                for k in ("compiled","twin"):same(out["original"],out[k])
                for k in ("duplicate_aux","shared_aux"):same(out["original"],out[k][0])
                same(out["duplicate_aux"][1],out["shared_aux"][1])
            z,p=out["original"];pred=z>=th;act=q["active"];loc=((abs(p-q["truth"])<=2)&act).sum(1)>=.75*act.sum(1)
            assert np.array_equal(pred,q["registered_hybrid_alert"]),"archived alert mismatch"
            assert np.array_equal(loc,q["registered_hybrid_localized"]),"archived localization mismatch"
            hit=pred&loc&(q["y"]==1)
            row={"target":n,"array_sha256":h,"localized":int(hit.sum()),"weak":int((hit&(q["amp"]<=1.25)).sum()),"controls":int((pred&(q["y"]==0)).sum()),"positives":int((q["y"]==1).sum()),"negatives":int((q["y"]==0).sum()),"exact":True,"timing":{k:{j:float(np.median([v[j] for v in t])) for j in ("wall","cpu")} for k,t in times.items()},"repeats":times}
            rows.append(row);save("partial.json",rows);print("R36F1_TARGET",json.dumps(row),flush=True)
    assert len(rows)==4 and g.sha((OUT/"lock.json").read_bytes())==lh
    m={k:sum(r[k] for r in rows) for k in ("localized","weak","controls","positives","negatives")};assert m["positives"]==m["negatives"]==512
    assert m["localized"]==164 and m["weak"]==20 and m["controls"]==12
    timing={k:{j:sum(r["timing"][k][j] for r in rows) for j in ("wall","cpu")} for k in rows[0]["timing"]}
    gates={k:{"warm_wall":timing[k]["wall"]<timing["original"]["wall"],"warm_cpu":timing[k]["cpu"]<timing["original"]["cpu"],"cold_wall":timing[k]["wall"]+compile["wall"]<timing["original"]["wall"]} for k in ("compiled","twin")}
    r={"status":"PASS_EXACT_ARCHIVED_HYBRID_REPLAY","metrics":m,"timing":timing,"gates":gates,"compile":compile,"archive_fetch_seconds":fetch,"reference_fit_seconds":fit_time,"lock_sha256":lh,"exact_archive_alert_and_localization":True,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"scope":"Registered linear+ODONTO detector, historical arrays, no I/O normalization or HANNAH performance claim; warm CPU and wall separate; extra smooth bank advisory only","r36f":"HOLD_FRESH_NATIVE_AND_PROTECTED_REAL_EXAM","production_adoption":False}
    save("result.json",r);print("R36F1_SUMMARY",json.dumps(r),flush=True)
if __name__=="__main__":main()
