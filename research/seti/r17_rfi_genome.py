#!/usr/bin/env python3
"""R17: train/test separated real OFF backgrounds; injection recognition, not origin classifier."""
import hashlib,io,json,math,os,pathlib,re,zipfile
import numpy as np
import requests
import torch
from torch import nn
from r9_expanded_search import normalized
S=json.loads(pathlib.Path("research/seti/r15_r16_r17_protocol.json").read_text())["r17"]
OUT=pathlib.Path("runtime/r17");OUT.mkdir(parents=True,exist_ok=True)
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
                for trial in range(S["trials_per_tile"]):
                    start=int(rng.integers(0,x.shape[1]-width+1));raw=x[:,start:start+width].copy()
                    shift=int(rng.integers(-32,33));ds=np.rint(np.arange(16)*shift/15).astype(int);c=int(rng.integers(max(16,16-ds.min()),min(width-16,width-16-ds.max())))
                    amp=S["amplitude"][trial%len(S["amplitude"])]
                    scale=1.4826*np.median(np.abs(raw-np.median(raw,axis=1)[:,None]),axis=1)
                    for label in (0,1):
                        p=raw.copy()
                        if label:p[np.arange(16),c+ds]+=amp*scale
                        v=normalized(p);X.append(np.clip(v,-8,16));Y.append(label)
                        # Known injected trajectory is NOT used by this baseline; blind max over all shifts.
                        best=-1e10
                        for s in range(-32,33):
                            offs=np.rint(np.arange(16)*s/15).astype(int);lo=max(0,-offs.min());hi=min(width,width-offs.max())
                            val=sum(v[i,lo+k:hi+k] for i,k in enumerate(offs))/4
                            best=max(best,float(val.max()))
                        scores.append(best);rows.append({"target":target,"source":o["name"],"tile_sha256":t["slice_sha256"],"center_mhz":t["center_mhz"],"trial":trial,"start":start,"label":label,"amplitude":amp if label else 0})
    for h in data.get("candidates",data.get("forward",[]))+data.get("reverse",[]):
        if h["status"]=="CONTROL_TRACK_PRESENT_RFI_CONCERN" or h.get("cadence")=="OFF_ON_OFF":
            library.append({"source_run":a["run"],"candidate":h,"label":"RFI_CONCERN_OR_REVERSE_CONTROL_NOT_CONFIRMED_EMITTER"})
    z.close()
    return np.asarray(X,np.float32),np.asarray(Y,np.int64),np.array(scores),rows,targets,library
class Tiny(nn.Module):
    def __init__(self):
        super().__init__();self.embed=nn.Linear(128,64);self.pos=nn.Parameter(torch.zeros(1,16,64))
        self.enc=nn.TransformerEncoder(nn.TransformerEncoderLayer(64,4,128,dropout=0.,batch_first=True),2)
        self.head=nn.Linear(64,2)
    def forward(self,x):return self.head(self.enc(self.embed(x)+self.pos).mean(1))
def metrics(y,p):
    y=np.asarray(y);p=np.asarray(p);return {"n":len(y),"accuracy":float(np.mean(y==p)),"injected_recall":float(np.mean(p[y==1]==1)),"unmodified_alert_fraction":float(np.mean(p[y==0]==1))}
def main():
    torch.set_num_threads(2);torch.manual_seed(S["seed"]);np.random.seed(S["seed"])
    tr,ty,ts,trm,tt,lib=build(S["source_artifacts"][0])
    te,ey,es,tem,et,lib2=build(S["source_artifacts"][1])
    assert not (tt&et) and len(tt)==len(et)==4
    assert not ({x["tile_sha256"] for x in trm}&{x["tile_sha256"] for x in tem})
    model=Tiny();opt=torch.optim.AdamW(model.parameters(),lr=.001);loss=nn.CrossEntropyLoss()
    tx=torch.from_numpy(tr);labels=torch.from_numpy(ty)
    for epoch in range(12):
        model.train()
        for idx in torch.randperm(len(tx)).split(32):
            opt.zero_grad();l=loss(model(tx[idx]),labels[idx]);l.backward();opt.step()
        emit("R17_EPOCH",{"epoch":epoch+1,"last_train_loss":float(l.detach())})
    # Freeze weights before evaluation; baseline threshold selected only on training.
    torch.save(model.state_dict(),OUT/"experimental_weights.pt")
    threshold=max(np.linspace(4,40,145),key=lambda t:(np.mean((ts>=t)==ty),-t))
    model.eval()
    with torch.no_grad(): pred=model(torch.from_numpy(te)).argmax(1).numpy()
    report={"train_targets":sorted(tt),"heldout_targets":sorted(et),"train_samples":len(ty),"test_samples":len(ey),"model_parameters":sum(p.numel() for p in model.parameters()),"transformer":metrics(ey,pred),"baseline":metrics(ey,es>=threshold),"baseline_threshold":float(threshold),"library_entries":len(lib+lib2),"weight_sha256":hashlib.sha256((OUT/"experimental_weights.pt").read_bytes()).hexdigest(),"scope":"synthetic injection recognition in disjoint archived OFF backgrounds, not true RFI-vs-ETI classification; paired samples correlated; no external blind custodian","production_adoption":False}
    (OUT/"result.json").write_text(json.dumps(report,indent=2)+"\n")
    (OUT/"rfi_library.json").write_text(json.dumps(lib+lib2,indent=2)+"\n")
    (OUT/"sample_manifest.json").write_text(json.dumps({"train":trm,"heldout":tem},indent=2)+"\n")
    emit("R17_SUMMARY",report)
if __name__=="__main__":main()
