#!/usr/bin/env python3
"""R17B: hybrid features and hard control patches; development evaluation only."""
import hashlib,io,json,math,os,pathlib,re,zipfile
import numpy as np
import requests
import torch
from torch import nn
from r9_expanded_search import normalized
S=json.loads(pathlib.Path("research/seti/r15_r16_r17_protocol.json").read_text())["r17"]
OUT=pathlib.Path("runtime/r17b");OUT.mkdir(parents=True,exist_ok=True)
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
                        v=normalized(p);features=track_features(v);X.append(features);Y.append(label)
                        # Known injected trajectory is NOT used by this baseline; blind max over all shifts.
                        scores.append(features[0]);rows.append({"target":target,"source":o["name"],"tile_sha256":t["slice_sha256"],"center_mhz":t["center_mhz"],"trial":trial,"hard_control":bool(hard),"start":start,"label":label,"amplitude":amp if label else 0})
    for h in data.get("candidates",data.get("forward",[]))+data.get("reverse",[]):
        if h["status"]=="CONTROL_TRACK_PRESENT_RFI_CONCERN" or h.get("cadence")=="OFF_ON_OFF":
            library.append({"source_run":a["run"],"candidate":h,"label":"RFI_CONCERN_OR_REVERSE_CONTROL_NOT_CONFIRMED_EMITTER"})
    z.close()
    return np.asarray(X,np.float32),np.asarray(Y,np.int64),np.array(scores),rows,targets,library

def track_features(v):
    n,w=v.shape; candidates=[]
    for shift in range(-32,33):
        ds=np.rint(np.arange(n)*shift/(n-1)).astype(int)
        lo=max(0,-ds.min());hi=min(w,w-ds.max())
        vals=sum(v[i,lo+k:hi+k] for i,k in enumerate(ds))/np.sqrt(n)
        for j in np.argsort(vals)[-3:]:
            b=int(lo+j);q=v[np.arange(n),b+ds]
            candidates.append((float(vals[j]),b,shift,float(np.mean(q>3)),float(np.std(q)),float(np.median(q))))
    candidates.sort(reverse=True)
    best,b,shift,persist,std,med=candidates[0]
    second=next((q[0] for q in candidates if abs(q[1]-b)>5 and abs(q[1]+q[2]-b-shift)>5),0.)
    stationary=max(q[0] for q in candidates if q[2]==0)
    median=np.median(v,axis=0);occup=np.mean(v>3,axis=0)
    return [best,persist,abs(shift),std,med,second,best-second,stationary,float(np.max(median)),float(np.mean(occup>=.5)),float(np.mean(v>3)),float(np.max(occup))]
def metrics(y,p):
    y=np.asarray(y);p=np.asarray(p)
    return {"n":len(y),"accuracy":float(np.mean(y==p)),"injected_recall":float(np.mean(p[y==1]==1)),"unmodified_alert_fraction":float(np.mean(p[y==0]==1))}
def main():
    torch.set_num_threads(2);torch.manual_seed(S["seed"]);np.random.seed(S["seed"])
    tr,ty,ts,trm,tt,lib=build(S["source_artifacts"][0])
    mu=tr.mean(0);sd=tr.std(0);sd=np.maximum(sd,1e-5)
    tx=torch.from_numpy(np.clip((tr-mu)/sd,-20,20).astype(np.float32));labels=torch.from_numpy(ty)
    model=nn.Sequential(nn.Linear(len(mu),32),nn.ReLU(),nn.Linear(32,2))
    opt=torch.optim.AdamW(model.parameters(),lr=.003);loss=nn.CrossEntropyLoss()
    for epoch in range(40):
        for idx in torch.randperm(len(tx)).split(32):
            opt.zero_grad();l=loss(model(tx[idx]),labels[idx]);l.backward();opt.step()
    threshold=max(np.linspace(4,60,225),key=lambda t:(np.mean((ts>=t)==ty),-t))
    model.eval()
    with torch.no_grad(): train_p=model(tx).softmax(1)[:,1].numpy()
    hybrid_thr=max(np.linspace(.05,.95,91),key=lambda t:(np.mean((train_p>=t)==ty),-t))
    torch.save({"weights":model.state_dict(),"mean":mu.tolist(),"scale":sd.tolist(),"threshold":float(hybrid_thr),"classical_threshold":float(threshold)},OUT/"hybrid_weights.pt")
    emit("R17B_FROZEN",{"weights_sha256":hashlib.sha256((OUT/"hybrid_weights.pt").read_bytes()).hexdigest(),"threshold":float(hybrid_thr),"baseline_threshold":float(threshold)})
    te,ey,es,tem,et,lib2=build(S["source_artifacts"][1])
    assert not(tt&et) and len(tt)==len(et)==4
    assert not({x["tile_sha256"] for x in trm}&{x["tile_sha256"] for x in tem})
    with torch.no_grad(): probs=model(torch.from_numpy(np.clip((te-mu)/sd,-20,20).astype(np.float32))).softmax(1)[:,1].numpy()
    pred=probs>=hybrid_thr;baseline=es>=threshold
    def group(mask):
        return {"hybrid":metrics(ey[mask],pred[mask]),"baseline":metrics(ey[mask],baseline[mask])}
    report={"train_targets":sorted(tt),"development_targets":sorted(et),"train_samples":len(ty),"development_samples":len(ey),"model_parameters":sum(p.numel() for p in model.parameters()),"overall":group(np.ones(len(ey),bool)),"ordinary":group(np.array([not r["hard_control"] for r in tem])),"hard":group(np.array([r["hard_control"] for r in tem])),"by_target":{t:group(np.array([r["target"]==t for r in tem])) for t in sorted(et)},"weights_sha256":hashlib.sha256((OUT/"hybrid_weights.pt").read_bytes()).hexdigest(),"scope":"Injection detector with hard OFF backgrounds; no independent RFI-versus-ETI labels. R11 previously evaluated in R17, so development only.","production_adoption":False}
    for r,p,b in zip(tem,probs,baseline):r.update(hybrid_probability=float(p),baseline_alert=bool(b))
    (OUT/"result.json").write_text(json.dumps(report,indent=2)+"\n")
    (OUT/"sample_manifest.json").write_text(json.dumps({"train":trm,"development":tem},indent=2)+"\n")
    (OUT/"rfi_library.json").write_text(json.dumps(lib+lib2,indent=2)+"\n")
    emit("R17B_SUMMARY",report)
if __name__=="__main__":main()
