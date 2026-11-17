#!/usr/bin/env python3
"""Prospective archive holdout of frozen R17B, no training or production adoption."""
import hashlib,io,json,math,os,pathlib,time,zipfile
import numpy as np
import requests
import torch
from torch import nn
import r11_control_calibration as r11
import r9_expanded_search as eng
import r4_chunk_replay as base
import r17b_hybrid as hybrid
P=pathlib.Path("research/seti/r15c_r16c_r17c_protocol.json")
S=json.loads(P.read_text())["r17c"]
OUT=pathlib.Path("runtime/r17c");OUT.mkdir(parents=True,exist_ok=True)
eng.OUT=OUT
eng.SPEC={**r11.S,"window_centers_mhz":S["windows_mhz"]}
def emit(k,x):print(k,json.dumps(x),flush=True)
def save(name,x):(OUT/name).write_text(json.dumps(x,indent=2)+"\n")
def sha(b):return hashlib.sha256(b).hexdigest()
def freeze():
    raw,url=base.fetch(r11.S["selection"]["index"],limit=32*base.MIB)
    (OUT/"target-index.json").write_bytes(raw)
    pool=set(r11.names(json.loads(raw)))-set(S["exclude"])
    order=sorted(pool,key=lambda t:sha((S["seed"]+"|"+t).encode()))
    m={"protocol_sha256":sha(P.read_bytes()),"index_sha256":sha(raw),"candidate_order":order[:S["max_queries"]],"selected":[],"attempts":[]}
    for target in m["candidate_order"]:
        try:
            fs,q=eng.query(target);fs.sort(key=lambda f:(f["day"],f["seconds"],f["url"]))
            cand=[]
            for i in range(len(fs)-3):
                a=fs[i:i+4]
                if [x["position"] for x in a]!=["S","R","S","R"] or len({x["day"] for x in a})!=1:continue
                if not all(0<y["seconds"]-x["seconds"]<=600 for x,y in zip(a,a[1:])):continue
                cand.append(a)
            cand.sort(key=lambda a:(a[0]["day"]<58575,a[0]["day"],a[0]["seconds"],a[0]["url"]))
            for a in cand[:5]:
                try:
                    hs=[eng.hdr(f) for f in a]
                    assert all(eng.valid(h,f,target) and h["shape"][0]==16 for h,f in zip(hs,a))
                    assert len({(h["shape"][0],round(h["ts"],7),round(h["df"],12)) for h in hs})==1
                    coords=[(h["src_raj"],h["src_dej"]) for h in hs]
                    assert coords[0]==coords[2] and coords[1]==coords[3] and coords[0]!=coords[1] and all(float(c[0])!=0 or float(c[1])!=0 for c in coords)
                    m["selected"].append({"target":target,"split":"calibration" if len(m["selected"])<2 else "evaluation","query":q,"files":[{**f,"header":h,"role":r} for f,h,r in zip(a,hs,("ON1","OFF1","ON2","OFF2"))]})
                    break
                except Exception as e:m["attempts"].append({"target":target,"day":a[0]["day"],"error":str(e)})
        except Exception as e:m["attempts"].append({"target":target,"error":str(e)})
        if len(m["selected"])==6:break
    save("manifest.json",m)
    emit("R17C_MANIFEST",{"selected":[{"target":q["target"],"split":q["split"],"files":[f["name"] for f in q["files"]]} for q in m["selected"]],"index_sha256":m["index_sha256"],"attempts":m["attempts"]})
def load_model():
    a=S["weights_artifact"];u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip"
    b=bytearray()
    with requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},stream=True,timeout=(15,60)) as r:
        r.raise_for_status()
        for c in r.iter_content(1024*1024):
            b.extend(c);assert len(b)<8*1024**2
    assert sha(b)==a["digest"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        assert sum(i.file_size for i in z.infolist())<8*1024**2
        paths=[i for i in z.namelist() if pathlib.PurePosixPath(i).name=="hybrid_weights.pt"]
        assert len(paths)==1
        w=z.read(paths[0])
    assert sha(w)==a["weights_sha256"]
    obj=torch.load(io.BytesIO(w),weights_only=True,map_location="cpu")
    model=nn.Sequential(nn.Linear(12,32),nn.ReLU(),nn.Linear(32,2))
    model.load_state_dict(obj["weights"],strict=True);model.eval()
    return model,np.asarray(obj["mean"],np.float32),np.asarray(obj["scale"],np.float32),obj
def classical(v):
    n,w=v.shape;best=-float("inf")
    for shift in range(-32,33):
        ds=np.rint(np.arange(n)*shift/(n-1)).astype(int)
        lo=max(0,-ds.min());hi=min(w,w-ds.max())
        val=sum(v[i,lo+k:hi+k] for i,k in enumerate(ds))/np.sqrt(n)
        best=max(best,float(val.max()))
    return best
def gather(q,model,mu,sd):
    rows=[]
    for f in q["files"]:
        if not f["role"].startswith("OFF"):continue
        for t in eng.get_tiles(f):
            x=np.load(OUT/t["artifact"],allow_pickle=False)
            assert sha(x.tobytes())==t["slice_sha256"] and list(x.shape)==t["shape"] and x.shape[0]==16
            rng=np.random.default_rng(int.from_bytes(hashlib.sha256((S["seed"]+"|"+f["name"]+"|"+str(t["center_mhz"])).encode()).digest()[:8],"big"))
            w=128;centers=[]
            for k in np.argsort(-np.median(eng.normalized(x),axis=0),kind="stable"):
                if w//2<=k<x.shape[1]-w//2 and all(abs(int(k)-c)>=w for c in centers):centers.append(int(k))
                if len(centers)==4:break
            assert len(centers)==4
            for i in range(16):
                hard=i>=12
                start=centers[i-12]-w//2 if hard else int(rng.integers(0,x.shape[1]-w+1))
                raw=x[:,start:start+w].copy()
                shift=int(rng.integers(-32,33));ds=np.rint(np.arange(16)*shift/15).astype(int)
                c=int(rng.integers(max(16,16-ds.min()),min(w-16,w-16-ds.max())))
                amp=(3,4,6,8)[i%4]
                scale=1.4826*np.median(np.abs(raw-np.median(raw,axis=1)[:,None]),axis=1)
                for label in (0,1):
                    p=raw.copy()
                    if label:p[np.arange(16),c+ds]+=amp*scale
                    a=time.perf_counter();v=eng.normalized(p);score=classical(v);classic_seconds=time.perf_counter()-a
                    a=time.perf_counter();fe=np.asarray(hybrid.track_features(eng.normalized(p)),np.float32)
                    with torch.no_grad():prob=float(model(torch.from_numpy(np.clip((fe-mu)/sd,-20,20))[None]).softmax(1)[0,1])
                    hybrid_seconds=time.perf_counter()-a
                    assert abs(score-float(fe[0]))<=max(.0001,abs(score)*1e-6)
                    rows.append({"target":q["target"],"split":q["split"],"source":f["name"],"tile_sha256":t["slice_sha256"],"center_mhz":t["center_mhz"],"start":start,"trial":i,"hard":hard,"label":label,"amplitude":amp if label else 0,"classical":score,"hybrid":prob,"classical_seconds":classic_seconds,"hybrid_seconds":hybrid_seconds})
    assert len(rows)==256
    save("target-"+q["target"]+".json",rows)
    emit("R17C_TARGET",{"target":q["target"],"split":q["split"],"samples":len(rows)})
    return rows
def threshold(rows,k,recall):
    positives=sorted((r[k] for r in rows if r["label"]==1),reverse=True)
    return float(positives[max(0,math.ceil(len(positives)*recall)-1)])
def met(rows,k,t):
    y=np.array([r["label"] for r in rows]);p=np.array([r[k]>=t for r in rows])
    return {"n":len(y),"positives":int(sum(y==1)),"negatives":int(sum(y==0)),"detected":int(sum(p[y==1])),"control_alerts":int(sum(p[y==0])),"recall":float(np.mean(p[y==1])),"control_alert_fraction":float(np.mean(p[y==0])),"accuracy":float(np.mean(p==y)),"seconds":sum(r[k+"_seconds"] for r in rows)}
def report(rows,thresholds):
    return {k:met(rows,k,thresholds[k]) for k in ("classical","hybrid")}
def scan():
    torch.set_num_threads(1);m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==sha(P.read_bytes())
    assert len(m["selected"])==6 and len({q["target"] for q in m["selected"]})==6
    model,mu,sd,obj=load_model()
    cal=[];test=[]
    for q in m["selected"][:2]:cal+=gather(q,model,mu,sd)
    assert len(cal)==512
    thresholds={k:threshold(cal,k,S["calibration_recall"]) for k in ("classical","hybrid")}
    lock={"weights_sha256":S["weights_artifact"]["weights_sha256"],"manifest_sha256":sha((OUT/"manifest.json").read_bytes()),"thresholds":thresholds,"calibration":report(cal,thresholds),"original_thresholds":{"classical":obj["classical_threshold"],"hybrid":obj["threshold"]},"status":"FROZEN_BEFORE_EVALUATION"}
    save("calibration_lock.json",lock);emit("R17C_LOCK",lock)
    for q in m["selected"][2:]:test+=gather(q,model,mu,sd)
    assert len(test)==1024 and not({r["tile_sha256"] for r in cal}&{r["tile_sha256"] for r in test})
    fixed=report(test,thresholds)
    # descriptive only: thresholds computed on evaluation positives, not a deployable result
    post={k:threshold(test,k,S["calibration_recall"]) for k in thresholds}
    summary={"calibration_targets":[q["target"] for q in m["selected"][:2]],"evaluation_targets":[q["target"] for q in m["selected"][2:]],"calibration_samples":len(cal),"evaluation_samples":len(test),"lock":lock,"prospective":fixed,"hard":report([r for r in test if r["hard"]],thresholds),"ordinary":report([r for r in test if not r["hard"]],thresholds),"by_target":{q["target"]:report([r for r in test if r["target"]==q["target"]],thresholds) for q in m["selected"][2:]},"descriptive_test_recall_matched":{"thresholds":post,"metrics":report(test,post),"not_prospective":True},"scope":"Frozen R17B weights; synthetic injections vs unmodified real OFF controls; correlated samples; archive holdout, not external blind test; calibration recall matched, test recall not guaranteed; timing includes normalization and scoring, excludes I/O; not source-origin classification","production_adoption":False}
    save("result.json",summary);emit("R17C_SUMMARY",summary)
if __name__=="__main__":
    import sys
    {"freeze":freeze,"scan":scan}[sys.argv[1]]()
