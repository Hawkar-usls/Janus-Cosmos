"""Historical descriptive diagnostics, never select a deployed threshold."""
import hashlib,io,json,os,pathlib,zipfile,collections
import numpy as np
import requests
from sklearn.metrics import roc_curve,roc_auc_score
import r41s1h2_lock as h2
P=pathlib.Path("research/seti/r41s1h2a_diagnostics_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r41s1h2a");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def write(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def archive(q,anchor='result.json'):
    api="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions"
    h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    def get(u):
        r=requests.get(u,headers=h,timeout=90);r.raise_for_status();return r
    a=get(f"{api}/artifacts/{q['id']}").json();assert a["workflow_run"]["id"]==q["run"] and a["digest"]=="sha256:"+q["digest"] and not a["expired"]
    b=get(f"{api}/artifacts/{q['id']}/zip").content;assert len(b)<64*1024**2 and sha(b)==q["digest"]
    z=zipfile.ZipFile(io.BytesIO(b));files=[x for x in z.infolist() if not x.is_dir()]
    assert sum(x.file_size for x in files)<256*1024**2
    assert len({x.filename for x in files})==len(files)
    a=[x for x in files if pathlib.PurePosixPath(x.filename).name==anchor]
    assert len(a)==1,(anchor,len(a))
    root=pathlib.PurePosixPath(a[0].filename).parent
    def rd(name):
        exact=str(root/name)
        m=[x for x in files if x.filename==exact]
        assert len(m)==1,(exact,len(m))
        return z.read(m[0])
    return rd
def curves():
    rows,rawsha=h2.restore();h2.verify(rows)
    rd=archive(S["source_h2"]);b=rd("lock.json");assert sha(b)==S["source_h2"]["lock_sha256"];lock=json.loads(b)
    assert lock["rows_sha256"]==rawsha
    held=[r for r in rows if r["split"]=="heldout"];q=S["roc"]["empirical_fpr_caps"]
    out={}
    for group in ["pooled"]+list(map(str,h2.T["heldout_seeds"])):
        rr=held if group=="pooled" else [r for r in held if str(r["seed"])==group]
        y=np.array([r["positive"] for r in rr]);g={}
        for k in ("physical","fixed"):
            sc=np.array([r[k] for r in rr]);fpr,tpr,_=roc_curve(y,sc,drop_intermediate=False)
            # No thresholds stored; descriptive measured envelope only.
            points=[{"fpr":float(f),"tpr":float(t),"fnr":float(1-t)} for f,t in zip(fpr,tpr)]
            table=[]
            for cap in q:
                ok=np.where(fpr<=cap+1e-12)[0];j=ok[np.argmax(tpr[ok])]
                table.append({"fpr_cap":cap,"actual_fpr":float(fpr[j]),"tpr":float(tpr[j]),"tp":int(round(tpr[j]*y.sum())),"fp":int(round(fpr[j]*(~y).sum()))})
            g[k]={"auc":float(roc_auc_score(y,sc)),"positives":int(y.sum()),"negatives":int((~y).sum()),"locked_metrics":h2.metrics(rr,k,lock["thresholds"][k]),"table":table,"curve":points}
        out[group]=g
    write("roc_det.json",out)
    print("H2A_ROC",json.dumps({k:{"auc":v["auc"],"table":v["table"]} for k,v in out["pooled"].items()}),flush=True)
    return out
def morph():
    rd=archive(S["source_f"]);d=rd("result.json");assert sha(d)==S["source_f"]["result_sha256"];f=json.loads(d)
    assert f["baseline"]["controls"]==53 and f["models"]["raw_noise"]["hypothetical"]["controls"]==51
    h=archive(S["source_h"],"h_result.json");res=json.loads(h("h_result.json"));assert res["groups"]==11 and res["decision_changes"]==0
    pred={(p["group"],p["kind"],p["row"]):p["would_keep"] for lane in ("parent","rescue") for p in f["models"]["raw_noise"][lane]["predictions"]}
    controls=[];tot=collections.Counter()
    names=("broad","zero","rival")
    def profile(v,role):
        if not v:return "unsupported"
        a={k:v.get("h_"+role+"_"+k+"_minus_narrow_mean") for k in names}
        if any(x is None or not np.isfinite(x) for x in a.values()):return "unsupported"
        a["narrow"]=0.;lo=min(a.values());best=[k for k,x in a.items() if x==lo]
        return best[0] if len(best)==1 else "ambiguous"
    for i in range(1,12):
        ff=json.loads(rd(f"group-{i:02d}.json"));hh=json.loads(h(f"group-{i:02d}.json"))
        assert ff["group"]==hh["group"]==i and ff["fit_viewed"]==hh["fit_viewed"] and ff["decision_changes"]==hh["decision_changes"]==0
        assert len(ff["records"])==len(hh["records"])
        for a,b in zip(ff["records"],hh["records"]):
            assert all(a[k]==b[k] for k in ("row","kind","lane","supported","weak","label_evaluation_only","localized_evaluation_only"))
            if ff["fit_viewed"] or a["kind"]!="eval" or a["label_evaluation_only"]!=0 or not pred.get((i,a["kind"],a["row"]),True):continue
            tot[a["lane"]]+=1
            if a["lane"]!="parent":continue
            v=b["features"];roles={k:profile(v,k) for k in ("primary_on","repeat_on","primary_off","repeat_off")}
            shared=bool(v and all(v.get("h_"+k+"_zero_minus_narrow_mean") is not None and v["h_"+k+"_zero_minus_narrow_mean"]>0 for k in roles))
            controls.append({"group":i,"row":a["row"],"supported":a["supported"],"profiles":roles,"on_repeat_agree":roles["primary_on"]==roles["repeat_on"],"shared_like":shared,"features":{k:x for k,x in (v or {}).items() if k.startswith("h_")}})
    assert tot=={"parent":38,"rescue":13} and len(controls)==38,dict(tot)
    r={"status":"HISTORICAL_PROFILE_DESCRIPTIVE_ONLY","controls":len(controls),"primary_on":dict(collections.Counter(x["profiles"]["primary_on"] for x in controls)),"repeat_on":dict(collections.Counter(x["profiles"]["repeat_on"] for x in controls)),"on_repeat_agreement":sum(x["on_repeat_agree"] and x["profiles"]["primary_on"]!="unsupported" for x in controls),"shared_like":sum(x["shared_like"] for x in controls),"supported":sum(x["profiles"]["primary_on"]!="unsupported" for x in controls),"rows":controls,"causal_labels":False}
    write("historical_morphology.json",r)
    print("H2A_MORPH",json.dumps({k:v for k,v in r.items() if k!="rows"}),flush=True)
    return r
def main():
    r=curves();m=morph()
    write("result.json",{"status":"DESCRIPTIVE_ROC_AND_HISTORICAL_MORPHOLOGY__HOLD_NATIVE_VALIDATION","protocol_sha256":sha(P.read_bytes()),"code_sha256":sha(pathlib.Path(__file__).read_bytes()),"locked_h2_unchanged":True,"historical_parent":m["controls"],"production_adoption":False,"external_validation":False})
if __name__=="__main__":main()
