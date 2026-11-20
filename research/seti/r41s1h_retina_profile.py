#!/usr/bin/env python3
"""Conditional held-out physical profile test for native registered pixel layers."""
import sys
ARGS=sys.argv[1:]
import pathlib,json
import numpy as np
from sklearn.metrics import roc_auc_score
import r41s1f_unclipped_joint as f
P=pathlib.Path("research/seti/r41s1h_retina_profile_protocol.json")
O=f.O
SIDE=f.SIDE
X=np.arange(-24,25,dtype=float)
def fit(y,t):
    v=t-t.mean();a=max(0.,float((y-y.mean())@v/(v@v))) if v@v>0 else 0.
    return y.mean()+a*v
def profile(cube,p):
    assert cube.shape==(4,16,208) and np.all((p>=24)&(p<=183))
    strips=cube[:,np.arange(16)[:,None],p[:,None]+np.arange(-24,25)[None,:]].astype(float)
    narrow=np.exp(-X**2/2);broad=np.exp(-X**2/50);zero=np.zeros(49)
    rivals=[np.exp(-(X-v)**2/2) for v in list(range(-12,-5))+list(range(6,13))]
    out={};details={}
    for k,strip in zip(f.retina.FACES,strips):
        vals=[]
        for j in range(4):
            te=np.arange(16)//4==j;tr=~te
            shoulder=strip[tr][:,SIDE];scale=f.mad(shoulder-np.median(shoulder,axis=1,keepdims=True))
            if not np.isfinite(scale) or scale<=0:return None
            data=(strip-np.median(strip[:,SIDE],axis=1,keepdims=True))/scale
            y=data[tr].mean(axis=0)
            pn,pb,pz=[fit(y,t) for t in (narrow,broad,zero)]
            pr=min((fit(y,t) for t in rivals),key=lambda q:np.mean((data[tr]-q)**2))
            losses=[float(np.mean((data[te]-q)**2)) for q in (pn,pb,pz,pr)]
            vals.append([losses[1]-losses[0],losses[2]-losses[0],losses[3]-losses[0]])
        a=np.asarray(vals)
        for i,name in enumerate(("broad_minus_narrow","zero_minus_narrow","rival_minus_narrow")):
            out[k+"_"+name+"_mean"]=float(a[:,i].mean());out[k+"_"+name+"_min"]=float(a[:,i].min())
        details[k]=vals
    for name in ("broad_minus_narrow","zero_minus_narrow","rival_minus_narrow"):
        out["on_off_"+name]=out["primary_on_"+name+"_mean"]-out["primary_off_"+name+"_mean"]
    assert np.isfinite(list(out.values())).all()
    return {"features":out,"folds":details}
def synthetic():
    rng=np.random.default_rng(413601);p=90+np.arange(16)//2;rr=np.arange(16);rows=[]
    for i in range(64):
        z=rng.normal(size=(4,16,208))
        def add(amp,sigma,roles):
            v=z.copy()
            for k in roles:v[k]+=amp*np.exp(-(np.arange(208)[None,:]-p[:,None])**2/(2*sigma*sigma))
            return v
        narrow=add(4,1,[0]);weak=add(1,1,[0]);broad=add(4,5,[0]);shared=add(4,1,[0,1,2,3]);scr=narrow.copy()
        for j,d in enumerate((8,-9,12,-11)):scr[0,j*4:j*4+4]=np.roll(narrow[0,j*4:j*4+4],d,axis=1)
        for k,v in (("noise",z),("narrow",narrow),("weak",weak),("broad",broad),("shared",shared),("scrambled",scr)):
            a=profile(v,p);assert a is not None;rows.append({"pair":i,"class":k,"features":a["features"]})
    med={k:{j:float(np.median([r["features"][j] for r in rows if r["class"]==k])) for j in rows[0]["features"]} for k in sorted({r["class"] for r in rows})}
    k="primary_on_broad_minus_narrow_mean";checks={"narrow_gt_broad":med["narrow"][k]>med["broad"][k],"narrow_gt_scramble":med["narrow"]["primary_on_zero_minus_narrow_mean"]>med["scrambled"]["primary_on_zero_minus_narrow_mean"],"shared_off_gt_narrow_off":med["shared"]["primary_off_zero_minus_narrow_mean"]>med["narrow"]["primary_off_zero_minus_narrow_mean"]}
    f.save("h_synthetic.json",{"rows":rows,"medians":med,"direction_checks":checks});print("R41S1H_SYNTHETIC",json.dumps(checks),flush=True)
def prepare():
    f.prepare();synthetic()
    f.save("h_protocol_receipt.json",{"sha256":f.s.sha(P.read_bytes()),"code_sha256":f.s.sha(pathlib.Path(__file__).read_bytes())})
def group(i):
    receipt=f.read("h_protocol_receipt.json");assert receipt["sha256"]==f.s.sha(P.read_bytes()) and receipt["code_sha256"]==f.s.sha(pathlib.Path(__file__).read_bytes())
    original=f.raw_diagnose
    def wrapped(cube,p):
        raw=original(cube,p)
        if raw is not None:
            d=profile(cube,p);raw["profile_h"]=d
            if d is not None:raw["features"].update({"h_"+k:v for k,v in d["features"].items()})
        return raw
    f.raw_diagnose=wrapped
    f.group(i)
def summary():
    gr=[f.read(f"group-{i:02d}.json") for i in range(1,12)]
    assert all(g["parity"] and g["decision_changes"]==0 for g in gr)
    rows=[{**v,"group":g["group"],"fit_viewed":g["fit_viewed"]} for g in gr for v in g["records"] if v["lane"]=="parent" and v["features"] is not None and any(k.startswith("h_") for k in v["features"])]
    assert rows
    keys=sorted(k for k in rows[0]["features"] if k.startswith("h_"))
    def auc(rs,k):
        yy=[1-v["label_evaluation_only"] for v in rs]
        return float(roc_auc_score(yy,[v["features"][k] for v in rs])) if len(set(yy))==2 else None
    tables={}
    for name,rs in (("fit",[v for v in rows if v["fit_viewed"]]),("audit",[v for v in rows if not v["fit_viewed"]])):
        tables[name]={"supported":len(rs),"auc_control":{k:{"pooled":auc(rs,k),"group_median":float(np.median(a)) if (a:=[x for i in sorted({v["group"] for v in rs}) if (x:=auc([v for v in rs if v["group"]==i],k)) is not None]) else None} for k in keys}}
    r={"status":"PASS_HELDOUT_PROFILE_DIAGNOSTIC__HOLD_FRESH_VALIDATION","groups":11,"parent_supported":len(rows),"features":len(keys),"synthetic_checks":f.read("h_synthetic.json")["direction_checks"],"decision_changes":0,"tables":tables,"production_adoption":False,"external_validation":False,"baseline":{"localized":318,"weak":69,"controls":53}}
    f.save("h_result.json",r);print("R41S1H_SUMMARY",json.dumps(r),flush=True)
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="group":group(int(ARGS[1]))
    elif ARGS[0]=="summary":summary()
    else:raise ValueError(ARGS)
