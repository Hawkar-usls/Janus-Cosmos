#!/usr/bin/env python3
"""Freeze Joint Min score; certify replay only, leave fresh exam HOLD."""
import hashlib,io,json,os,pathlib,zipfile
import numpy as np
import requests
P=pathlib.Path("research/seti/r41s1g_frozen_noise_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r41s1g");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+"\n")
def download(spec):
    h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    root="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions"
    r=requests.get(f"{root}/artifacts/{spec['artifact']}",headers=h,timeout=30);r.raise_for_status();m=r.json()
    assert m["workflow_run"]["id"]==spec["run"] and m["digest"]=="sha256:"+spec["archive_sha256"]
    r=requests.get(f"{root}/artifacts/{spec['artifact']}/zip",headers=h,timeout=120);r.raise_for_status()
    assert len(r.content)<128*1024**2 and sha(r.content)==spec["archive_sha256"]
    z=zipfile.ZipFile(io.BytesIO(r.content));assert len(z.namelist())==len(set(z.namelist())) and sum(x.file_size for x in z.infolist())<512*1024**2
    return z
def score(f,lock):
    if f is None:return None
    keys=lock["features"];x=np.asarray([np.nan if f[k] is None else f[k] for k in keys],float)
    x=np.where(np.isfinite(x),x,np.asarray(lock["imputation"]))[~np.asarray(lock["dropped_all_missing"],bool)]
    return float(((x-np.asarray(lock["scaler_mean"]))/np.asarray(lock["scaler_scale"]))@np.asarray(lock["coefficients"])+lock["intercept"])
def main():
    z=download(S["source"]);rb=z.read("result.json");assert sha(rb)==S["source"]["result_sha256"]
    d=json.loads(rb);m=d["models"]["raw_noise"]["rescue"]
    lock={k:m[k] for k in ("features","imputation","dropped_all_missing","scaler_mean","scaler_scale","coefficients","intercept","train_positive_threshold")}
    assert len(lock["features"])==40
    save("frozen_noise_model.json",lock)
    rows={}
    for i in range(1,12):
        g=json.loads(z.read(f"group-{i:02d}.json"))
        if g["fit_viewed"]:continue
        for v in g["records"]:
            if v["lane"]=="rescue":rows[(i,v["kind"],v["row"])]=v
    replay=[];err=0.
    for p in m["predictions"]:
        v=rows[(p["group"],p["kind"],p["row"])];x=score(v["features"],lock)
        err=max(err,abs(x-p["score"]));assert abs(x-p["score"])<=1e-9
        keep=x>=lock["train_positive_threshold"];assert keep==p["would_keep"]
        replay.append({"group":p["group"],"kind":p["kind"],"row":p["row"],"score":x,"would_keep":keep,"label_evaluation_only":v["label_evaluation_only"],"localized":v["localized_evaluation_only"],"weak":v["weak"]})
    assert len(replay)==105
    lost=[p for p in replay if not p["would_keep"]]
    assert sum(p["label_evaluation_only"]==0 and p["kind"]=="eval" for p in lost)==2 and sum(p["localized"] for p in lost)==0
    save("historical_replay.json",replay)
    prior=download(S["prior_qualification"]);rec=[]
    for name in prior.namelist():
        if name.endswith(".json"):
            v=json.loads(prior.read(name));rec.append({"name":name,"sha256":sha(prior.read(name)),"data":v})
    assert rec
    save("prior_qualification.json",rec)
    result={"status":"PASS_FROZEN_SCORE_REPLAY__HOLD_NEW_INDEPENDENT_GROUPS","model_sha256":sha((O/"frozen_noise_model.json").read_bytes()),"historical_replayed":len(replay),"max_score_error":err,"historical_control_removed":2,"historical_localized_lost":0,"fresh_groups_measured":0,"external_validation":False,"production_adoption":False}
    save("result.json",result);print("R41S1G_SUMMARY",json.dumps(result),flush=True)
if __name__=="__main__":main()
