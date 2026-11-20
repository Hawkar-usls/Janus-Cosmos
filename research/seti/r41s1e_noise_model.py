#!/usr/bin/env python3
"""Archive-only conditional contrast noise diagnostics, not source inference."""
import hashlib, io, json, os, pathlib, time, zipfile
import numpy as np
import requests
from scipy.signal import savgol_filter
from scipy.ndimage import median_filter
from sklearn.metrics import roc_auc_score
O=pathlib.Path("runtime/r41s1e");O.mkdir(parents=True,exist_ok=True)
P=pathlib.Path("research/seti/r41s1e_noise_model_protocol.json");L=json.loads(P.read_text())
def sha(x):return hashlib.sha256(x).hexdigest()
def save(k,d): (O/k).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def rms(x):return float(np.sqrt(np.mean(x*x)))
def robust(x):return float(1.4826*np.median(abs(x-np.median(x))))
def corr(x,y):
    x=x-x.mean();y=y-y.mean();d=np.linalg.norm(x)*np.linalg.norm(y)
    return float(np.clip(x@y/d,-1,1)) if d>1e-12 else None
def ratio(a,b):return float(a/b) if b>1e-12 else None
def scalar(x):
    assert x.shape==(16,) and np.isfinite(x).all() and np.max(abs(x))<=4.000000001
    med=x-median_filter(x,size=3,mode="nearest");sg=x-savgol_filter(x,5,2,mode="interp")
    y=x-x.mean();e=abs(np.fft.rfft(y))**2;e[1:8]*=2
    return {"mean":float(x.mean()),"std":float(x.std()),"median_residual_rms":rms(med),"median_residual_mad":robust(med),"savgol_residual_rms":rms(sg),"hf_fraction":ratio(e[4:].sum(),e[1:].sum()),"haar_detail_fraction":ratio(((y[::2]-y[1::2])**2).sum()/2,(y*y).sum()),"residual_lag1":corr(sg[:-1],sg[1:]),"difference_rms":rms(np.diff(x)),"clipped_fraction":float((abs(x)>=4).mean())}
def pred(a,b):
    # Equal-variance Gaussian working approximation, scale comes from held-out
    # complements' OFF differences only. OFF can contain a real signal or RFI.
    scores={k:[] for k in ("zero","target_only","shared","separate")}
    scales=[]
    for j in (0,1):
        tr=np.arange(16)%2==j;te=~tr;x=a[tr];y=b[tr]
        sigma=max(.5,robust(np.diff(y))/np.sqrt(2));scales.append(float(sigma))
        shared=float(np.r_[x,y].mean())
        preds={"zero":(0.,0.),"target_only":(max(0.,float(x.mean())),0.),"shared":(max(0.,shared),)*2,"separate":(float(x.mean()),float(y.mean()))}
        for k,(m,n) in preds.items():
            residual=np.r_[a[te]-m,b[te]-n]
            scores[k].append(float(np.mean(-.5*np.log(2*np.pi*sigma*sigma)-.5*(residual/sigma)**2)))
    s={k:float(np.mean(v)) for k,v in scores.items()}
    return {"target_vs_zero":s["target_only"]-s["zero"],"shared_vs_zero":s["shared"]-s["zero"],"target_vs_shared":s["target_only"]-s["shared"],"separate_vs_shared":s["separate"]-s["shared"],"off_train_sigma":float(np.mean(scales))}
def diagnose(v,fixed):
    a,b,c,d=np.asarray(v,float);assert v.shape==(4,16)
    fs={}
    for k,x in zip(("primary_on","primary_off","repeat_on_crossfit","repeat_off_crossfit"),(a,b,c,d)):
        fs.update({k+"_"+j:z for j,z in scalar(x).items()})
    fs.update({"on_crossfit_disagreement":rms(a-c),"off_crossfit_disagreement":rms(b-d),"on_crossfit_corr":corr(a,c),"off_crossfit_corr":corr(b,d),"on_fixed_block_disagreement":rms(np.asarray(fixed[0])-a.reshape(2,8).mean(1)),"off_fixed_block_disagreement":rms(np.asarray(fixed[1])-b.reshape(2,8).mean(1))})
    fs.update({"primary_"+k:x for k,x in pred(a,b).items()})
    assert all(x is None or np.isfinite(x) for x in fs.values())
    return fs
def med(rows):
    return {k:float(np.median(v)) if (v:=[r["features"][k] for r in rows if r["features"][k] is not None]) else None for k in rows[0]["features"]}
def syn():
    rng=np.random.default_rng(413301);out=[];t=np.linspace(0,1,16)
    for i in range(64):
        x=rng.normal(size=(4,16));ar=x.copy()
        for j in range(1,16):ar[:,j]=.75*ar[:,j-1]+np.sqrt(1-.75**2)*x[:,j]
        smooth=2+0.35*np.sin(2*np.pi*t);jagged=2+1.1*(-1.)**np.arange(16)
        for k,v in (("noise_iid",x),("noise_ar1",ar),("smooth_target",x+np.array([smooth,0*smooth,smooth,0*smooth])),("jagged_target",x+np.array([jagged,0*jagged,jagged,0*jagged])),("weak_target",x+np.array([smooth,0*smooth,smooth,0*smooth])*.35),("shared_narrow",x+smooth),("transient_target",x+np.array([smooth*(t<.5),0*smooth,0*smooth,0*smooth]))):
            v=np.clip(v,-4,4);f=diagnose(v,v[2:].reshape(2,2,8).mean(2));out.append({"pair":i,"class":k,"features":f})
    return {"rows":out,"medians":{k:med([r for r in out if r["class"]==k]) for k in sorted({r["class"] for r in out})}}
def archive():
    x=L["parent"];h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"};root="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions"
    def get(u):
        r=requests.get(u,headers=h,timeout=120);r.raise_for_status();return r
    run=get(f"{root}/runs/{x['run']}").json();assert run["head_sha"]==x["head_sha"] and run["conclusion"]=="success"
    md=get(f"{root}/artifacts/{x['artifact']}").json();assert md["workflow_run"]["id"]==x["run"] and md["digest"]=="sha256:"+x["archive_sha256"]
    raw=get(f"{root}/artifacts/{x['artifact']}/zip").content;assert sha(raw)==x["archive_sha256"] and len(raw)<128*1024**2
    z=zipfile.ZipFile(io.BytesIO(raw));names=z.namelist();assert len(names)==len(set(names)) and sum(i.file_size for i in z.infolist())<512*1024**2
    assert all(not pathlib.PurePosixPath(n).is_absolute() and ".." not in pathlib.PurePosixPath(n).parts for n in names)
    root=pathlib.PurePosixPath(next(n for n in names if pathlib.PurePosixPath(n).name=="result.json")).parent
    def rd(k):return z.read(str(root/k))
    v=rd("result.json");assert sha(v)==x["result_sha256"];d=json.loads(v);assert d["qualified"]==891 and d["decision_changes"]==0 and len(d["trace_sha256"])==22
    for k,h in d["trace_sha256"].items():assert sha(rd(k))==h
    out=[]
    for i in range(1,12):
        g=json.loads(rd(f"group-{i:02d}.json"));assert g["parity"] and not g["decision_changes"]
        count=0
        for kind in ("eval","shared"):
            rr=[v for v in g["records"] if v["kind"]==kind];tr=json.loads(rd(f"trace-{i:02d}-{kind}.json"));assert tr["decision_changes"]==0 and [v["row"] for v in rr]==[v["row"] for v in tr["records"]]
            for v,w in zip(rr,tr["records"]):
                assert v["supported"]==w["advanced"] and not w["decision_changed"] and v["lane"]==w["lane"]
                if not w["advanced"]:continue
                e=w["states"][-1]["snapshot"]["evidence"];a=e["differentiation"];assert a["features"]==v["features"]
                detail=a["details"];cur=np.array([detail[k]["primary_morphology"]["trace"] for k in ("on","off")]+[detail[k]["template_shift"]["crossfit_morphology"]["trace"] for k in ("on","off")],float)
                fixed=[e["advisory_repeat"]["D2"]["on_off"]["candidate"][k] for k in ("on","off")]
                f=diagnose(cur,fixed)
                out.append({"group":i,"fit_viewed":g["fit_viewed"],"kind":kind,"row":v["row"],"lane":v["lane"],"source_tile_sha256":w["source_tile_sha256"],"parent_state_sha256":w["states"][-1]["state_sha256"],"features":f,"label_evaluation_only":v["label_evaluation_only"]});count+=1
        assert count==g["qualified"];print("R41S1E_GROUP",i,count,flush=True)
    assert len(out)==891 and sum(not x["fit_viewed"] for x in out)==409
    return out
def table(rows):
    def auc(rs,k):
        x=[(r["features"][k],1-r["label_evaluation_only"]) for r in rs if r["features"][k] is not None]
        return float(roc_auc_score([a[1] for a in x],[a[0] for a in x])) if len({a[1] for a in x})==2 else None
    return {"qualified":len(rows),"control":med([r for r in rows if not r["label_evaluation_only"]]),"positive":med([r for r in rows if r["label_evaluation_only"]]),"auc":{k:{"pooled":auc(rows,k),"group_median":float(np.median(v)) if (v:=[z for g in sorted({r["group"] for r in rows}) if (z:=auc([r for r in rows if r["group"]==g],k)) is not None]) else None,"groups":len(v)} for k in rows[0]["features"]}}
def main():
    t=time.perf_counter();y=syn();save("synthetic.json",y);print("R41S1E_SYNTHETIC",len(y["rows"]),flush=True)
    x=archive();save("records.json",x)
    res={"status":"PASS_ARCHIVED_NOISE_DIAGNOSTICS__HOLD_SCIENTIFIC_SCREEN","source":L["parent"],"qualified":len(x),"features":len(x[0]["features"]),"synthetic_count":len(y["rows"]),"decision_changes":0,"baseline":{"localized":318,"weak":69,"controls":53},"tables":{"fit":table([v for v in x if v["fit_viewed"]]),"audit":table([v for v in x if not v["fit_viewed"]]),"audit_parent":table([v for v in x if not v["fit_viewed"] and v["lane"]=="parent"]),"audit_rescue":table([v for v in x if not v["fit_viewed"] and v["lane"]=="rescue"])},"wall_s":time.perf_counter()-t,"records_sha256":sha((O/"records.json").read_bytes()),"synthetic_sha256":sha((O/"synthetic.json").read_bytes()),"external_validation":False,"production_adoption":False,"claim":"descriptive clipped contrast traces; no calibrated likelihood, no raw noise model, no physical source inference"}
    save("result.json",res);print("R41S1E_SUMMARY",json.dumps({k:v for k,v in res.items() if k!="tables"}),flush=True)
if __name__=="__main__":main()
