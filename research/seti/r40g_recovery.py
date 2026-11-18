#!/usr/bin/env python3
"""R40G: a second, rival-trajectory witness for rejected R40F alerts."""
import json,pathlib,time
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier
import r40f_screen as f
c=f.c;t=f.t
P=pathlib.Path("research/seti/r40g_recovery_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40g");O.mkdir(parents=True,exist_ok=True)
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def screen(q,L):
    F,_=f.features(q)
    v=np.clip(F,np.array(L["lo"]),np.array(L["hi"]))
    return ((v-np.array(L["mean"]))/np.array(L["scale"]))@np.array(L["coef"])+L["intercept"]
def features(q):
    V=q["V"].astype(float);n,T,W=V.shape;p=q["r40c_fused_path"].astype(int)
    med=np.median(V,axis=2,keepdims=True);mad=1.4826*np.median(abs(V-med),axis=2,keepdims=True)
    X=np.clip((V-med)/np.maximum(.5,mad),-4,4)
    offsets=np.array([d for d in range(-32,33) if abs(d)>3])
    paths=p[:,None,:]+offsets[None,:,None]
    ok=((paths>=0)&(paths<W)).all(2);assert (ok.sum(1)>=1).all()
    trace=X[np.arange(n)[:,None],np.arange(T)[None,:],p]
    rival=X[np.arange(n)[:,None,None],np.arange(T)[None,None,:],np.clip(paths,0,W-1)]
    comb=lambda x: x.mean(-1)
    a=comb(trace);v=comb(rival);v=np.where(ok,v,np.nan)
    center=np.nanmedian(v,axis=1);disp=np.nanmedian(abs(v-center[:,None]),axis=1)
    z=(a-center)/np.maximum(.1,1.4826*disp)
    blocks=trace.reshape(n,4,4).mean(2)
    block_r=rival.reshape(n,len(offsets),4,4).mean(3)
    bc=np.nanmedian(np.where(ok[:,:,None],block_r,np.nan),axis=1)
    adj=blocks-bc
    odds=trace[:,1::2].mean(1);evens=trace[:,::2].mean(1)
    base,_=f.features(q)
    vals=[a,center,disp,z,a-np.nanmax(v,axis=1),np.sum((v<=a[:,None])&ok,axis=1)/ok.sum(1),np.min(adj,axis=1),np.median(adj,axis=1),np.max(adj,axis=1),np.std(adj,axis=1),odds,evens,np.minimum(odds,evens),abs(odds-evens),np.median(trace,axis=1),np.mean(trace>0,axis=1),np.max(abs(trace),axis=1),np.max(abs(trace),axis=1)/np.maximum(.01,np.sum(abs(trace),axis=1)),np.mean(abs(V)>8,axis=(1,2)),np.mean(mad,axis=(1,2))]
    out=np.concatenate([np.stack(vals,axis=1),base],axis=1)
    assert np.isfinite(out).all()
    return out
def compare(q,old,probs,thresholds):
    alert=q["r40c_fused_alert"].astype(bool);hit=q["r40c_fused_localized"].astype(bool);y=q["y"];weak=np.isin(q["amp"],t.S["weak_amplitudes"]);null=y==0
    keeps={"r40c_fused":alert,"r40f_screen":alert&old}
    for k,p in probs.items():keeps[k]=alert&(old|(p>=thresholds[k]))
    out={}
    for k,a in keeps.items():
        h=a&hit;out[k]={"localized":int(h.sum()),"weak":int((h&weak).sum()),"control_alerts":int((a&null).sum()),"lost":int((hit&~a).sum()),"weak_lost":int((hit&~a&weak).sum()),"recovered":int((h&~old).sum()),"restored_controls":int((a&~old&null).sum()),"positives":int((y==1).sum()),"controls":int(null.sum())}
    return out,keeps
def main():
    st=time.perf_counter()
    z,n=c.artifact(S["source"],16);raw=z.read(n["lock.json"]);assert c.sha(raw)==S["source"]["lock_sha256"];sf=json.loads(raw)
    mz,mn=c.artifact(f.S["model"],8);mr=mz.read(mn["model.json"]);L=json.loads(mr);assert c.sha(mr)==f.S["model"]["model_sha256"]==sf["parent_model_sha256"]
    hz,hn,rec,_=t.parent.restore();train=[];cal=[];roots=set()
    for r in rec["records"]:
        x=t.parent.load(hz,hn,r);roots.update(str(v) for v in x["tile_sha"]);q=f.historical(x,L)
        row=(q,features(q),screen(q,sf)>=sf["threshold"])
        (cal if r["name"].startswith("cal_") else train).append(row)
    assert len(train)==8 and len(cal)==2
    X=np.concatenate([v[q["r40c_fused_alert"]] for q,v,o in train]);y=np.concatenate([q["r40c_fused_localized"][q["r40c_fused_alert"]] for q,v,o in train])
    assert y.any() and (~y).any()
    sc=StandardScaler().fit(X)
    lr=LogisticRegression(**S["models"]["logistic"]).fit(sc.transform(X),y)
    tree=ExtraTreesClassifier(**S["models"]["tree"]).fit(X,y)
    assert lr.n_iter_.max()<S["models"]["logistic"]["max_iter"]
    def preds(X):return {"r40g_logistic":lr.decision_function(sc.transform(X)),"r40g_tree":tree.predict_proba(X)[:,1]}
    out={k:[] for k in S["arms"][2:]};fallback={k:[] for k in out}
    for q,v,o in train:
        h=q["r40c_fused_localized"].astype(bool)
        for k,p in preds(v).items():
            out[k].extend(p[h&~o].tolist());fallback[k].extend(p[h].tolist())
    th={k:float(min(v or fallback[k])) for k,v in out.items()}
    lock={"protocol_sha256":c.sha(P.read_bytes()),"parent_screen_sha256":c.sha(raw),"train_alerts":len(y),"train_localized":int(y.sum()),"thresholds":th,"training_rejected_positives":{k:len(v) for k,v in out.items()},"logistic":{"mean":sc.mean_.tolist(),"scale":sc.scale_.tolist(),"coef":lr.coef_.tolist(),"intercept":lr.intercept_.tolist()},"tree_seed_and_params":S["models"]["tree"]}
    # Serialize trained tree as numeric arrays; no unsafe pickle.
    lock["trees"]=[{"children_left":e.tree_.children_left.tolist(),"children_right":e.tree_.children_right.tolist(),"feature":e.tree_.feature.tolist(),"threshold":e.tree_.threshold.tolist(),"value":e.tree_.value.tolist()} for e in tree.estimators_]
    for q,v,o in cal:
        m,_=compare(q,o,preds(v),th);assert m["r40f_screen"]["lost"]==0
        assert all(m[k]["lost"]==0 for k in th)
    save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes());emit("R40G_LOCK",{k:v for k,v in lock.items() if k not in ("trees","logistic")})
    rows=[]
    for target,a in f.S["artifacts"].items():
        zz,nn=c.artifact(a,128);res=json.loads(zz.read(nn["result.json"]));assert res["target"]==target and res["model_sha256"]==f.S["model"]["model_sha256"] and not roots.intersection(res["source_tile_hashes"])
        for kind in ("eval","shared"):
            q=f.load(zz,nn,kind+".npz");old=screen(q,sf)>=sf["threshold"]
            prior=f.load(z,n,target+"_"+kind+".npz")
            assert np.array_equal(old&q["r40c_fused_alert"],prior["r40f_screen_priority_alert"])
            pr=preds(features(q));m,keeps=compare(q,old,pr,th)
            # Preserve original flags, candidate path, and all scores; no rejection.
            np.savez_compressed(O/(target+"_"+kind+".npz"),original_alert=q["r40c_fused_alert"],original_path=q["r40c_fused_path"],screen_score=screen(q,sf),**pr,**{k+"_priority_alert":v for k,v in keeps.items()})
            row={"target":target,"kind":kind,"metrics":m};rows.append(row);emit("R40G_GROUP",row)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def total(kind):
        rr=[v["metrics"] for v in rows if v["kind"]==kind];assert len(rr)==5
        return {k:{m:sum(x[k][m] for x in rr) for m in rr[0][k]} for k in S["arms"]}
    e=total("eval");sh=total("shared")
    assert [e["r40f_screen"][k] for k in ("localized","weak","control_alerts")]==[187,26,23]
    gates={}
    for k in th:
        v=e[k];w=sh[k];gates[k]={"recovered":v["recovered"],"remaining_lost":v["lost"],"control_alerts":v["control_alerts"],"weak_lost":v["weak_lost"],"shared_lost":w["lost"],"gate":v["lost"]==0 and w["lost"]==0 and v["weak_lost"]==0 and v["control_alerts"]<=41}
    d={"status":"DEVELOPMENT_REPLAY_ONLY","eval":e,"shared":sh,"gates":gates,"rows":rows,"lock_sha256":lh,"wall_seconds":time.perf_counter()-st,"production_adoption":False,"discovery_relaunch_allowed":False}
    save("result.json",d);emit("R40G_SUMMARY",{k:v for k,v in d.items() if k!="rows"})
if __name__=="__main__":main()
