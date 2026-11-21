"""R42A5: fixed-route historical fusion, rigorous archive/row parity, fit-only calibration."""
import hashlib,io,json,os,pathlib,zipfile
import requests
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

P=pathlib.Path("research/seti/r42a5_combined_score_protocol.json");Q=json.loads(P.read_text())
O=pathlib.Path("runtime/r42a5");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def archive(q):
    url="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/"+str(q["artifact"])
    h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    def get(u):
        r=requests.get(u,headers=h,timeout=180);r.raise_for_status();return r
    m=get(url).json();assert m["workflow_run"]["id"]==q["run"] and m["digest"]=="sha256:"+q["digest"] and not m["expired"]
    b=get(url+"/zip").content;assert len(b)<256*1024**2 and sha(b)==q["digest"]
    z=zipfile.ZipFile(io.BytesIO(b));files=[x for x in z.infolist() if not x.is_dir()]
    assert sum(x.file_size for x in files)<1024*1024**2
    names={pathlib.PurePosixPath(x.filename).name:x for x in files};assert len(names)==len(files)
    return lambda n:z.read(names[n])
def checked(read,n,digest):
    raw=read(n);assert sha(raw)==digest,n;return raw
def npz(raw):
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:return {k:z[k] for k in z.files}
def equal(a,b):assert np.array_equal(a,b)
def squish(x):return x/(1+np.abs(x))
def load():
    read={k:archive(v) for k,v in Q["sources"].items()}
    mr,nr,cr,orr=[read[k] for k in ("minimal","noise","common","models")]
    orig=json.loads(checked(orr,"result.json",Q["sources"]["models"]["result_sha256"]))
    models=orig["models"]["raw_noise"]
    assert orig["baseline"]["controls"]==53 and models["hypothetical"]["controls"]==51
    ns=json.loads(nr("result.json"));assert ns["rows"]==3168 and ns["score_supported"]==3143
    assert ns["original_qualified"]==891
    cs=json.loads(cr("receipt.json"));assert len(cs["records"])==22
    common={x["name"]:x for x in cs["records"]}
    rows=[]
    for g in range(1,12):
        num=f"{g:02d}";mg=json.loads(mr("group-"+num+".json"));ng=json.loads(nr("group-"+num+".json"))
        fit=g in Q["split"]["fit"];assert mg["group"]==ng["group"]==g and mg["fit"]==ng["fit"]==fit
        assert mg["rows"]==ng["rows"]==288
        for kind in ("eval","shared"):
            tag=num+"-"+kind
            a=next(x for x in mg["records"] if x["kind"]==kind)
            b=next(x for x in ng["records"] if x["kind"]==kind)
            c=common[tag+".npz"];assert c["group"]==g and c["fit"]==fit and c["array_sha256"]==a["registered_sha256"]
            x=npz(checked(mr,tag+"-scores.npz",a["score_file_sha256"]))
            n=npz(checked(nr,tag+"-noise-all.npz",b["file_sha256"]))
            z=npz(checked(cr,tag+".npz",c["sha256"]))
            y=z["y"].astype(bool);w=z["weak"].astype(bool)
            route=x["parent_route"].astype(bool);valid=n["supported"].astype(bool);eligible=n["frozen_eligible"].astype(bool)
            m=x["score"].astype(float);raw=n["score"].astype(float)
            assert len(m)==len(y)==a["rows"]==b["rows"]==c["rows"]
            equal(x["source_tile_sha"],n["tile_sha"]);equal(z["tile_sha"],n["tile_sha"])
            equal(route,n["parent_route"]);equal(x["path"],np.where(route[:,None],n["parent_path"],n["rescue_path"]))
            equal(m>=0,x["alert"]);equal(x["alert"],z["minimal_alert"]);equal(n["original_minimal_alert"],x["alert"])
            equal(n["original_minimal_hit"],z["minimal_hit"])
            loc=((abs(x["path"]-z["truth"])<=2)&z["active"]).sum(1)>=.75*z["active"].sum(1)
            equal((m>=0)&loc&y,z["minimal_hit"])
            assert np.isfinite(m).all() and np.array_equal(np.isfinite(raw),valid) and np.all(~eligible|valid)
            cuts=np.where(route,float(models["parent"]["train_positive_threshold"]),float(models["rescue"]["train_positive_threshold"]))
            nm=np.zeros(len(m));nm[valid]=squish(raw[valid]-cuts[valid])
            ref=np.where(eligible,np.minimum(m,nm),m)
            score=np.where(valid&(eligible|~x["alert"]),np.minimum(m,nm),m)
            equal(score>=0,ref>=0)
            if not fit:
                equal(ref>=0,z["noise_alert"]);equal((ref>=0)&loc&y,z["noise_hit"])
            rows.append(dict(group=np.full(len(m),g),row=np.arange(len(m)),shared=np.full(len(m),kind=="shared"),
                m=m,n=nm,valid=valid,eligible=eligible,route=route,baseline=m>=0,reference=ref>=0,
                primary=score,y=y,weak=w,loc=loc))
    d={k:np.concatenate([x[k] for x in rows]) for k in rows[0]}
    assert len(d["m"])==3168 and d["valid"].sum()==3143 and d["eligible"].sum()==891
    return d
def metrics(a,d,mask):
    a=a&mask;hit=a&d["y"]&d["loc"];ref=d["reference"]&mask;rh=ref&d["y"]&d["loc"];false=a&~d["y"];rf=ref&~d["y"]
    return {"localized":int(hit.sum()),"weak":int((hit&d["weak"]).sum()),"controls":int(false.sum()),
      "ordinary":int((hit&~d["shared"]).sum()),"shared":int((hit&d["shared"]).sum()),
      "reference_hits_lost":int((rh&~hit).sum()),"reference_weak_lost":int((rh&~hit&d["weak"]).sum()),
      "reference_hits_gained":int((hit&~rh).sum()),"reference_false_removed":int((rf&~false).sum()),
      "reference_false_added":int((false&~rf).sum())}
def threshold(s,d,fit):
    null=np.sort(s[fit&~d["y"]])[::-1];assert len(null)==768
    t=float(np.nextafter(null[42],np.inf));assert np.isfinite(t) and ((s>=t)&fit&~d["y"]).sum()<=42
    return t
def ids(mask,d):
    return [{"group":int(d["group"][i]),"kind":"shared" if d["shared"][i] else "eval","row":int(d["row"][i])} for i in np.flatnonzero(mask)]
def main():
    d=load();fit=d["group"]<=6;audit=~fit
    assert fit.sum()==1728 and audit.sum()==1440 and (audit&~d["y"]).sum()==640 and (audit&d["y"]).sum()==800
    base=metrics(d["baseline"],d,audit);ref=metrics(d["reference"],d,audit)
    for name,v in (("minimal",base),("noise",ref)):
        assert all(v[k]==x for k,x in Q["reference"][name].items()),name
    assert ((d["baseline"]&d["loc"]&d["y"]&audit)&~d["reference"]).sum()==0
    # Primary needs no new training. Secondary is a separate predeclared fixed-hyperparameter head.
    X=np.column_stack([d["m"],d["n"],~d["valid"],d["route"],d["m"]*d["n"]]).astype(float)
    target=d["y"]&d["loc"]
    sc=StandardScaler().fit(X[fit]);lr=LogisticRegression(C=.1,class_weight="balanced",random_state=42,max_iter=2000).fit(sc.transform(X[fit]),target[fit])
    scores={"primary_min":d["primary"],"secondary_logistic":lr.decision_function(sc.transform(X))}
    assert all(np.isfinite(s).all() for s in scores.values())
    locks={k:threshold(s,d,fit) for k,s in scores.items()}
    md={"protocol_sha256":sha(P.read_bytes()),"code_sha256":sha(pathlib.Path(__file__).read_bytes()),
        "fit_groups":Q["split"]["fit"],"audit_groups":Q["split"]["audit"],"fit_null_budget":42,"thresholds":locks,
        "secondary_model":{"features":["M","N","missing","parent","M*N"],"mean":sc.mean_.tolist(),"scale":sc.scale_.tolist(),"coef":lr.coef_[0].tolist(),"intercept":float(lr.intercept_[0])},
        "fit_metrics":{k:metrics(s>=locks[k],d,fit) for k,s in scores.items()}}
    save("fit_lock.json",md);print("R42A5_FIT_LOCK",json.dumps(md),flush=True)
    results={};curves={}
    for k,s in scores.items():
        alert=s>=locks[k];v=metrics(alert,d,audit)
        rh=d["reference"]&d["loc"]&d["y"]&audit;hit=alert&d["loc"]&d["y"]&audit
        result={"fit":md["fit_metrics"][k],"audit":v,"threshold":locks[k],
            "historical_goal_met":v["controls"]<=35 and v["reference_hits_lost"]==0 and v["reference_weak_lost"]==0 and v["localized"]>=318 and v["weak"]>=69,
            "lost_ids":ids(rh&~hit,d),"gained_ids":ids(hit&~rh,d),
            "per_group":[{"group":g,**metrics(alert,d,d["group"]==g)} for g in range(7,12)]}
        results[k]=result
        # Descriptive curve only; not a second threshold selection or an adopted decision.
        ts=np.unique(s[audit])[::-1];points=[]
        for t in ts:
            z=metrics(s>=t,d,audit);points.append({"threshold":float(t),**z})
        curves[k]=points
        print("R42A5_ARM",json.dumps({"arm":k,**{j:x for j,x in result.items() if j not in ("lost_ids","gained_ids")}}),flush=True)
    np.savez_compressed(O/"all_scores.npz",**d,**scores)
    save("descriptive_curves.json",curves)
    out={"status":"PASS_HISTORICAL_SCORE_CONTRACT__HOLD_INDEPENDENT_VALIDATION","rows":3168,"noise_supported":3143,"missing":25,"reference_minimal":base,"reference_noise":ref,
       "arms":results,"fit_lock_sha256":sha((O/"fit_lock.json").read_bytes()),"full_rerouting":False,"audit_previously_viewed":True,
       "production_adoption":False,"external_validation":False,"independent_pointing_hold":True}
    save("result.json",out)
    print("R42A5_SUMMARY",json.dumps({k:v for k,v in out.items() if k!="arms"}),flush=True)
if __name__=="__main__":main()
