"""Authenticated score-only historical audit; no thresholds selected or deployed."""
import hashlib,io,json,os,pathlib,zipfile
import numpy as np
import requests
P=pathlib.Path("research/seti/r42a1_score_audit_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r42a1");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(k,d):(O/k).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def archive(q):
    u="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/"+str(q["id"])
    h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    r=requests.get(u,headers=h,timeout=90);r.raise_for_status();a=r.json()
    assert a["workflow_run"]["id"]==S["run"] and a["digest"]=="sha256:"+q["digest"] and not a["expired"]
    r=requests.get(u+"/zip",headers=h,timeout=120);r.raise_for_status();b=r.content
    assert len(b)<128*1024**2 and sha(b)==q["digest"]
    z=zipfile.ZipFile(io.BytesIO(b));f=[i for i in z.infolist() if not i.is_dir()]
    assert sum(i.file_size for i in f)<512*1024**2
    names={pathlib.PurePosixPath(i.filename).name:i for i in f};assert len(names)==len(f)
    return lambda n:z.read(names[n])
def counts(alert,hit,y,w,shared):
    return {"localized":int(hit.sum()),"weak":int((hit&w).sum()),"controls":int((alert&~y).sum()),"ordinary":int((hit&~shared).sum()),"shared":int((hit&shared).sum())}
def main():
    c=archive(S["common"]);r=json.loads(c("receipt.json"));assert len(r["records"])==22 and sum(x["rows"] for x in r["records"])==3168
    rows=[]
    for x in r["records"]:
        b=c(x["name"]);assert sha(b)==x["sha256"]
        with np.load(io.BytesIO(b),allow_pickle=False) as a:
            assert sha(a["V"].tobytes())==x["array_sha256"]
            if not x["fit"]:rows.append((x,{k:a[k] for k in ("y","truth","active","weak","minimal_alert","minimal_hit","noise_alert","noise_hit")}))
    assert len(rows)==10
    cat=lambda k:np.concatenate([x[1][k] for x in rows])
    y=cat("y").astype(bool);w=cat("weak").astype(bool);shared=np.concatenate([np.full(len(a["y"]),x["kind"]=="shared",bool) for x,a in rows])
    assert y.sum()==800 and (~y).sum()==640
    out={"status":"HISTORICAL_SCORE_AUDIT_ONLY","positives":800,"negatives":640,"reference":{},"native":{},"production_adoption":False,"external_validation":False}
    for k in ("minimal","noise"):
        d=counts(cat(k+"_alert"),cat(k+"_hit"),y,w,shared);assert all(d[v]==r["audit"][k][v] for v in d);out["reference"][k]=d
    out["reference"]["curve_status"]="HOLD_CONTINUOUS_FULL_PIPELINE_SCORES_NOT_EXPORTED"
    for k,src in S["competitors"].items():
        rd=archive(src);res=json.loads(rd("result.json"));lock=json.loads(rd("lock.json"))
        assert res["input_id"]==lock["input_id"]==str(S["common"]["id"]) and res["input_digest"]==lock["input_digest"]==S["common"]["digest"]
        scores=[];locs=[];alerts={a:[] for a in res["arms"]};hits={a:[] for a in res["arms"]}
        for x,a in rows:
            with np.load(io.BytesIO(rd(x["name"][:-4]+"-scores.npz")),allow_pickle=False) as z:
                s=z["score"];p=z["path"];assert np.isfinite(s).all() and p.shape==a["truth"].shape
                loc=((abs(p-a["truth"])<=2)&a["active"]).sum(1)>=.75*a["active"].sum(1)
                scores.append(s);locs.append(loc)
                for arm,th in lock["thresholds"].items():
                    al=s>=th;hit=al&loc&(a["y"]==1);assert np.array_equal(al,z[arm+"_alert"]) and np.array_equal(hit,z[arm+"_hit"])
                    alerts[arm].append(al);hits[arm].append(hit)
        s=np.concatenate(scores);loc=np.concatenate(locs)
        locked={}
        for arm in alerts:
            d=counts(np.concatenate(alerts[arm]),np.concatenate(hits[arm]),y,w,shared);assert all(d[v]==res["arms"][arm][v] for v in d);locked[arm]=d
        points=[{"controls":0,"localized":0,"weak":0,"ordinary":0,"shared":0,"fpr":0.,"tpr":0.,"fnr":1.}]
        for t in np.unique(s)[::-1]:
            a=s>=t;d=counts(a,a&loc&y,y,w,shared)
            points.append({**d,"fpr":d["controls"]/640,"tpr":d["localized"]/800,"fnr":1-d["localized"]/800})
        table=[{"cap":cap,**max((d for d in points if d["fpr"]<=cap+1e-12),key=lambda d:(d["localized"],-d["controls"]))} for cap in S["fpr_caps"]]
        out["native"][k]={"locked":locked,"table":table,"curve":points}
        print("R42A1_METHOD",k,json.dumps({"locked":locked,"table":table}),flush=True)
    save("result.json",out)
    print("R42A1_SUMMARY",json.dumps({"status":out["status"],"reference":out["reference"],"methods":list(out["native"])}),flush=True)
if __name__=="__main__":main()
