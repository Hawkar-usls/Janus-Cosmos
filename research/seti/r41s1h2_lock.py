"""Freeze H1 analytic score and development-only operating points before any real evaluation."""
import hashlib,io,json,os,pathlib,zipfile
from collections import Counter
import numpy as np
import requests
P=pathlib.Path("research/seti/r41s1h2_operating_point_protocol.json")
S=json.loads(P.read_text());Q=S["source"];T=S["threshold_policy"]
O=pathlib.Path("runtime/r41s1h2");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def write(name,d):
    b=(json.dumps(d,indent=2,sort_keys=True,allow_nan=False)+"\n").encode();(O/name).write_bytes(b);return sha(b)
def restore():
    api="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions"
    h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    def get(u):
        r=requests.get(u,headers=h,timeout=90);r.raise_for_status();return r
    run=get(f"{api}/runs/{Q['run']}").json();assert run["status"]=="completed" and run["conclusion"]=="success"
    a=get(f"{api}/artifacts/{Q['artifact']}").json()
    assert a["workflow_run"]["id"]==Q["run"] and a["digest"]=="sha256:"+Q["archive_sha256"] and not a["expired"]
    b=get(f"{api}/artifacts/{Q['artifact']}/zip").content
    assert len(b)<4*1024**2 and sha(b)==Q["archive_sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        files=[x for x in z.infolist() if not x.is_dir()]
        assert sum(x.file_size for x in files)<16*1024**2
        names={pathlib.PurePosixPath(x.filename).name:x.filename for x in files}
        assert len(names)==len(files)
        result=json.loads(z.read(names["result.json"]));raw=z.read(names["rows.json"]);rows=json.loads(raw)
    assert result["protocol_sha256"]==Q["protocol_sha256"] and result["code_sha256"]==Q["code_sha256"]
    assert result["heldout_examples"]==1152 and result["development_examples"]==768
    assert result["historical_tuning"] is False and result["decision_changes"]==0
    assert sha(pathlib.Path("research/seti/r41s1h1_physical_profile.py").read_bytes())==Q["code_sha256"]
    assert sha(pathlib.Path("research/seti/r41s1h1_physical_protocol.json").read_bytes())==Q["protocol_sha256"]
    return rows,sha(raw)
def verify(rows):
    counts=Counter((r["split"],r["seed"],r["class"]) for r in rows)
    expected={(split,seed,k):48 for split,seeds in (("development",T["development_seeds"]),("heldout",T["heldout_seeds"])) for seed in seeds for k in T["positive_classes"]+T["negative_classes"]}
    assert counts==expected and len(rows)==1920
    assert len({(r["seed"],r["class"],r["index"]) for r in rows})==1920
    for r in rows:
        assert r["positive"]==(r["class"] in T["positive_classes"])
        assert r["split"] in ("development","heldout") and 0<=r["index"]<48
        assert np.isfinite([r["physical"],r["fixed"]]).all()
def metrics(rows,key,cut):
    pos=[r for r in rows if r["positive"]];neg=[r for r in rows if not r["positive"]]
    tp=sum(r[key]>=cut for r in pos);fp=sum(r[key]>=cut for r in neg)
    return {"positives":len(pos),"negatives":len(neg),"tp":tp,"fn":len(pos)-tp,"fp":fp,"tn":len(neg)-fp,"tpr":tp/len(pos),"fpr":fp/len(neg),"classes":{k:{"n":len(rr),"kept":sum(r[key]>=cut for r in rr)} for k in sorted({r["class"] for r in rows}) if (rr:=[r for r in rows if r["class"]==k])}}
def main():
    rows,rawsha=restore();verify(rows)
    dev=[r for r in rows if r["split"]=="development"]
    cuts={k:float(min(r[k] for r in dev if r["positive"])) for k in ("physical","fixed")}
    assert all(metrics(dev,k,v)["fn"]==0 for k,v in cuts.items())
    lock={"schema":"janus.cosmos.r41s1h2.lock.v1","status":"LOCKED_SYNTHETIC_OPERATING_POINT","source":Q,"rows_sha256":rawsha,"h2_protocol_sha256":sha(P.read_bytes()),"thresholds":cuts,"rule":"keep score >= threshold; unsupported keep; parent only; no threshold changes","calibration":"development seeds only","native_validation":False,"historical_tuning":False,"production_adoption":False}
    locksha=write("lock.json",lock)
    # Heldout is accessed only after the numeric lock has been serialized.
    held=[r for r in rows if r["split"]=="heldout"]
    d={"status":"LOCKED_SYNTHETIC_OPERATING_POINT__HOLD_NATIVE_VALIDATION","lock_sha256":locksha,"thresholds":cuts,"development":{k:metrics(dev,k,v) for k,v in cuts.items()},"heldout":{k:metrics(held,k,v) for k,v in cuts.items()},"heldout_per_seed":{str(seed):{k:metrics([r for r in held if r["seed"]==seed],k,v) for k,v in cuts.items()} for seed in T["heldout_seeds"]},"independent_real_groups":0,"native_test_run":False,"production_adoption":False,"external_validation":False}
    assert sha((O/"lock.json").read_bytes())==locksha
    write("result.json",d);print("R41S1H2_LOCK",json.dumps({k:v for k,v in d.items() if k!="heldout_per_seed"}),flush=True)
if __name__=="__main__":main()
