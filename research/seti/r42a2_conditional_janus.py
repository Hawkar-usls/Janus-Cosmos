"""Historical conditional NOISE margins on fixed Minimal alerts; not full detector ROC."""
import pathlib,json,hashlib,io,zipfile,os
import numpy as np
import requests
P=pathlib.Path("research/seti/r42a2_conditional_janus_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r42a2");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def archive(q,anchor):
    h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    u="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/"+str(q["id"])
    def get(u):
        r=requests.get(u,headers=h,timeout=120);r.raise_for_status();return r
    a=get(u).json();assert a["workflow_run"]["id"]==q["run"] and a["digest"]=="sha256:"+q["digest"]
    b=get(u+"/zip").content;assert len(b)<256*1024**2 and sha(b)==q["digest"]
    z=zipfile.ZipFile(io.BytesIO(b));files=[x for x in z.infolist() if not x.is_dir()]
    assert sum(x.file_size for x in files)<1024*1024**2 and len({x.filename for x in files})==len(files)
    m=[x for x in files if pathlib.PurePosixPath(x.filename).name==anchor];assert len(m)==1
    root=pathlib.PurePosixPath(m[0].filename).parent
    return lambda n:z.read(str(root/n))
def main():
    rd=archive(S["noise"],"result.json");b=rd("result.json");assert sha(b)==S["noise"]["result_sha256"]
    f=json.loads(b);models=f["models"]["raw_noise"];pred={}
    for lane in ("parent","rescue"):
        m=models[lane]
        for r in m["predictions"]:
            k=(r["group"],r["kind"],r["row"]);assert k not in pred
            margin=r["score"]-m["train_positive_threshold"];assert np.isfinite(margin) and (margin>=0)==r["would_keep"]
            pred[k]=(lane,margin)
    c=archive(S["common"],"receipt.json");rec=json.loads(c("receipt.json"));assert len(rec["records"])==22 and sum(x["rows"] for x in rec["records"])==3168
    rows=[];used=set()
    for x in rec["records"]:
        raw=c("inputs/"+x["name"]);assert sha(raw)==x["sha256"]
        with np.load(io.BytesIO(raw),allow_pickle=False) as a:
            assert sha(a["V"].tobytes())==x["array_sha256"]
            if x["fit"]:continue
            g=json.loads(rd(f'group-{x["group"]:02d}.json'));assert not g["fit_viewed"]
            old={r["row"]:r for r in g["records"] if r["kind"]==x["kind"]}
            assert set(old)==set(np.flatnonzero(a["minimal_alert"]))
            for i in range(len(a["y"])):
                k=(x["group"],x["kind"],i);r=old.get(i);v=pred.get(k)
                if v is not None:
                    assert r and r["lane"]==v[0];used.add(k)
                if r:assert (r["label_evaluation_only"],r["weak"],r["localized_evaluation_only"])==(int(a["y"][i]),bool(a["weak"][i]),bool(a["minimal_hit"][i]))
                margin=None if v is None else v[1]
                keep=bool(a["minimal_alert"][i] and (margin is None or margin>=0))
                assert keep==bool(a["noise_alert"][i]) and (keep and bool(a["minimal_hit"][i]))==bool(a["noise_hit"][i])
                rows.append({"group":x["group"],"kind":x["kind"],"row":i,"y":bool(a["y"][i]),"weak":bool(a["weak"][i]),"base":bool(a["minimal_alert"][i]),"hit":bool(a["minimal_hit"][i]),"frozen":keep,"lane":r["lane"] if r else None,"margin":margin})
    assert used==set(pred) and len(rows)==1440 and sum(r["y"] for r in rows)==800
    def count(keep):
        rr=[r for r,k in zip(rows,keep) if k]
        hit=sum(r["hit"] for r in rr);weak=sum(r["hit"] and r["weak"] for r in rr)
        return {"localized":hit,"weak":weak,"controls":sum(not r["y"] for r in rr),"parent_controls":sum(not r["y"] and r["lane"]=="parent" for r in rr),"rescue_controls":sum(not r["y"] and r["lane"]=="rescue" for r in rr),"lost":318-hit,"weak_lost":69-weak,"shared":sum(r["hit"] and r["kind"]=="shared" for r in rr)}
    ref=count([r["frozen"] for r in rows]);assert (ref["localized"],ref["weak"],ref["controls"],ref["parent_controls"],ref["rescue_controls"])==(318,69,51,38,13)
    out={"status":"HISTORICAL_CONDITIONAL_CURVE_ONLY","reference":ref,"arms":{},"production_adoption":False,"external_validation":False,"full_janus_curve":False}
    for arm in S["arms"]:
        def eligible(r):
            return r["margin"] is not None and (arm=="both_lanes" or arm.startswith(r["lane"]+"_only"))
        vals=sorted({r["margin"] for r in rows if eligible(r)}|{0.})
        # Include strictest endpoint; all non-alerts stay rejected.
        vals.append(float(np.nextafter(vals[-1],np.inf)))
        points=[]
        for t in vals:
            d=count([r["margin"]>=t if eligible(r) else r["frozen"] for r in rows])
            if t==0:assert d==ref
            points.append(d|{"fpr":d["controls"]/640,"tpr":d["localized"]/800})
        table=[]
        for cap in S["fpr_caps"]:
            sub=[d for d in points if d["fpr"]<=cap+1e-12];table.append({"cap":cap,"best":max(sub,key=lambda d:(d["localized"],d["weak"],-d["controls"])) if sub else None})
        safe=min((d for d in points if d["lost"]==d["weak_lost"]==0),key=lambda d:d["controls"])
        detail={"supported":sum(eligible(r) for r in rows),"floor":min(d["controls"] for d in points),"zero_loss":safe,"table":table}
        out["arms"][arm]=detail|{"points":points}
        print("R42A2_ARM",arm,json.dumps(detail),flush=True)
    save("margins.json",rows);save("result.json",out)
    print("R42A2_SUMMARY",json.dumps({k:v for k,v in out.items() if k!="arms"}),flush=True)
if __name__=="__main__":main()
