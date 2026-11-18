#!/usr/bin/env python3
"""Use unchanged native adapters on one authenticated common-input artifact."""
import hashlib,io,json,os,pathlib,resource,sys,time,zipfile
import numpy as np
import requests
MODE=sys.argv[1];assert MODE in ("bliss","fastdd","turbo")
OUT=pathlib.Path("runtime/r35g_native")/MODE;OUT.mkdir(parents=True,exist_ok=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def sha(b):return hashlib.sha256(b).hexdigest()
def restore():
    aid=os.environ["INPUT_ID"];digest=os.environ["INPUT_DIGEST"].removeprefix("sha256:");assert aid.isdigit() and len(digest)==64
    u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{aid}/zip";r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=120);r.raise_for_status()
    raw=r.content;assert len(raw)<128*1024**2 and sha(raw)==digest
    z=zipfile.ZipFile(io.BytesIO(raw));assert sum(i.file_size for i in z.infolist())<512*1024**2
    names={pathlib.PurePosixPath(i.filename).name:i.filename for i in z.infolist() if not i.is_dir()}
    assert len(names)==len([i for i in z.infolist() if not i.is_dir()])
    rec=json.loads(z.read(names["receipt.json"]));assert rec["status"]=="PASS_COMMON_REPLAY" and len(rec["records"])==10
    return z,names,rec,{"id":int(aid),"digest":digest}
def init():
    t=time.perf_counter()
    if MODE=="bliss":
        import r35a_bliss_native as m
        m.OUT=OUT
        fn=lambda x,tag:m.run(x,tag)
    elif MODE=="fastdd":
        import r35c_fastdd_native as m
        m.OUT=OUT;f=m.build();D,R,use=m.recover_paths(f)
        fn=lambda x,tag:m.score(f,x,D,R,use)
    else:
        import r35d_turbo_native as m
        m.OUT=OUT;k=m.Kernels(False,1,0)
        fn=lambda x,tag:m.run(x,tag,k)
    elapsed=time.perf_counter()-t
    rng=np.random.default_rng(35071);x=rng.normal(size=(3,16,128)).astype(np.float32);truth=[]
    for i,s in enumerate((-15,0,15)):
        p=64+np.rint(np.arange(16)*s/15).astype(int);truth.append(p);x[i,np.arange(16),p]+=30
    s,p,_=fn(x,"canary");assert np.all(s>5) and np.all((abs(p-np.stack(truth))<=2).sum(1)>=12)
    return fn,elapsed
def main():
    fn,compile_seconds=init();z,n,receipt,art=restore()
    def load(r):
        raw=z.read(n[r["name"]]);assert sha(raw)==r["sha256"]
        q=np.load(io.BytesIO(raw),allow_pickle=False)
        assert sha(q["V"].tobytes())==r["array_sha256"]
        return q
    cal=[]
    for r in receipt["records"]:
        if not r["name"].startswith("cal_"):continue
        with load(r) as q:
            s,p,t=fn(q["V"],pathlib.Path(r["name"]).stem)
            cal.append(s[q["y"]==0])
    assert len(cal)==2
    c=np.concatenate(cal);assert len(c)==256 and np.isfinite(c).all()
    sortedc=np.sort(c)[::-1];th=float(np.nextafter(sortedc[int(len(c)*.05)],np.inf));assert int((c>=th).sum())<=12
    save("lock.json",{"threshold":th,"calibration_control_alerts":int((c>=th).sum()),"calibration_controls":256,"source":art,"parent_lock":receipt["lock_sha256"]})
    rs=[]
    for r in receipt["records"]:
        if r["name"].startswith("cal_"):continue
        with load(r) as q:
            s,p,t=fn(q["V"],pathlib.Path(r["name"]).stem)
            y=q["y"];a=q["amp"];pred=s>=th
            loc=((abs(p-q["truth"])<=2)&q["active"]).sum(1)>=.75*q["active"].sum(1)
            hit=pred&loc&(y==1);base=q["hybrid_r28_localized"]&(y==1);plus=q["hybrid_plus_smooth_localized"]&(y==1)
            row={"target":r["name"],"localized":int(hit.sum()),"weak":int((hit&(a<=1.25)).sum()),"very_weak":int((hit&(a==.75)).sum()),"control_alerts":int((pred&(y==0)).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum()),"added_vs_r28":int((hit&~base).sum()),"lost_vs_r28":int((base&~hit).sum()),"added_vs_plus":int((hit&~plus).sum()),"lost_vs_plus":int((plus&~hit).sum()),"models":{str(k):int((hit&(q["model"]==k)).sum()) for k in np.unique(q["model"])}, "timing":t}
            rs.append(row);save("partial.json",rs)
            np.savez_compressed(OUT/(pathlib.Path(r["name"]).stem+"_scores.npz"),score=s,path=p,alert=pred,localized_positive=hit)
    def total(prefix):
        xx=[r for r in rs if r["target"].startswith(prefix)];assert len(xx)==4
        d={k:sum(r[k] for r in xx) for k in ("localized","weak","very_weak","control_alerts","positives","controls","added_vs_r28","lost_vs_r28","added_vs_plus","lost_vs_plus")}
        d["models"]={k:sum(v["models"].get(k,0) for v in xx) for k in sorted({i for v in xx for i in v["models"]})}
        d["timing"]={k:sum(v["timing"][k] for v in xx) for k in xx[0]["timing"]}
        return d
    evals=total("eval_");shared=total("shared_");assert evals["positives"]==evals["controls"]==512 and shared["positives"]==128
    d={"status":"NATIVE_COMMON_INPUT_COMPARISON","mode":MODE,"source":art,"threshold":th,"eval":evals,"shared":shared,"freshness_status":receipt["freshness_status"],"r32_overlap":receipt["r32_overlap"],"compile_setup_seconds":compile_seconds,"peak_rss_kib_process":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"scope":"Pinned native CPU adapters, registered float32 arrays with synthetic coordinate headers. turboSETI DAT path is approximate. Separate native timing and bridge; not GPU or blind origin test.","production_adoption":False}
    save("result.json",d);print("R35G_NATIVE_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
