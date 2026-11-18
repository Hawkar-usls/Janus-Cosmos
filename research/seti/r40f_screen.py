#!/usr/bin/env python3
"""R40F: protected priority screen, fixed historical fit and calibration, viewed transfer replay."""
import io,json,pathlib,sys,time
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
# r40cde uses argv[1] for output root; no main runs on import
sys.argv[1:]=["screen"]
import r40cde as c
import r40b_transfer as t
P=pathlib.Path("research/seti/r40f_screen_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40f");O.mkdir(parents=True,exist_ok=True)
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def signedlog(x):return np.sign(x)*np.log1p(abs(x))
def features(q):
    V=np.asarray(q["V"],float);n,T,W=V.shape;assert (T,W)==(16,128)
    p=np.asarray(q["r40c_fused_path"],int);b=np.asarray(q["r28_frozen_path"],int)
    assert p.shape==b.shape==(n,T) and np.all((p>=0)&(p<W)) and np.all((b>=0)&(b<W))
    ix=np.arange(n)[:,None];ti=np.arange(T)[None,:]
    agreement=(abs(p-b)<=2).sum(1)
    vals=[signedlog(q["r40c_fused_score"]),signedlog(q["r28_frozen_score"]),agreement/16, q["r28_frozen_alert"].astype(float)]
    profile=np.median(V,axis=1)
    vals.extend([signedlog(np.max(profile,axis=1)),np.log1p(np.std(profile,axis=1)),(abs(V)>8).mean((1,2)),np.log1p(np.max(abs(V),axis=(1,2)))])
    for path in (p,b):
        trace=V[ix,ti,path];delta=np.diff(path,axis=1)
        sides=path[:,:,None]+np.array([-5,-4,-3,3,4,5])
        ok=(sides>=0)&(sides<W)
        fl=V[ix[:,:,None],ti[:,:,None],np.clip(sides,0,W-1)]
        bg=np.nanmedian(np.where(ok,fl,np.nan),axis=2)
        d=trace-bg
        block=trace.reshape(n,4,4).mean(2)
        vals.extend([signedlog(trace.mean(1)),np.log1p(trace.std(1)),signedlog(np.median(trace,axis=1)),signedlog(np.quantile(trace,.25,axis=1)),signedlog(np.quantile(trace,.75,axis=1)),signedlog(trace.max(1)),(trace>0).mean(1),signedlog(d.mean(1)),signedlog(np.median(d,axis=1)),signedlog(block.min(1)),signedlog(block.max(1)),np.log1p(block.std(1)),np.abs(delta).mean(1),np.abs(np.diff(delta,axis=1)).mean(1)])
    F=np.stack(vals,axis=1);assert np.isfinite(F).all()
    return F,agreement
def historical(q,L):
    F,p=t.feature(q);pr=c.base_predictions(q,F,p,L)
    vals={"r28_frozen":(q["hybrid_r28_score"],q["hybrid_r28_paths"]),"r40c_fused":c.fused(q,*pr["r40c_fused"],L["margins"]["r40c_fused"])}
    d={k:q[k] for k in ("V","y","amp","truth","active","tile_sha")}
    for k,(s,pp) in vals.items():
        _,h,a=t.metrics(q,s,pp,L["thresholds"][k])
        d.update({k+"_score":s,k+"_path":pp,k+"_alert":a,k+"_localized":h})
    return d
def load(z,n,filename):
    with np.load(io.BytesIO(z.read(n[filename])),allow_pickle=False) as x:return {k:x[k].copy() for k in x.files}
def compare(q,stat,th):
    F,agree=features(q);alert=q["r40c_fused_alert"].astype(bool);hit=q["r40c_fused_localized"].astype(bool)
    keep=stat>=th;prot=q["r28_frozen_alert"].astype(bool)&(agree>=12)
    y=q["y"];weak=np.isin(q["amp"],t.S["weak_amplitudes"]);null=y==0
    arms={"r28_frozen":(q["r28_frozen_alert"],q["r28_frozen_localized"]),"r40c_fused":(alert,hit),"r40f_screen":(alert&keep,hit&keep),"r40f_protected":(alert&(keep|prot),hit&(keep|prot))}
    out={}
    for k,(a,h) in arms.items():
        out[k]={"localized":int(h.sum()),"weak":int((h&weak).sum()),"control_alerts":int((a&null).sum()),"lost_vs_r40c":int((hit&~h).sum()),"weak_lost_vs_r40c":int((hit&~h&weak).sum()),"controls_removed_vs_r40c":int((alert&~a&null).sum()),"controls_added_vs_r40c":int((a&~alert&null).sum()),"positives":int((y==1).sum()),"controls":int(null.sum())}
    return out,arms,agree
def main():
    st=time.perf_counter();z,n=c.artifact(S["model"],8);raw=z.read(n["model.json"]);assert c.sha(raw)==S["model"]["model_sha256"];L=json.loads(raw)
    hz,hn,rec,old=t.parent.restore()
    train=[];cal=[];roots=set()
    for r in rec["records"]:
        q=t.parent.load(hz,hn,r);roots.update(str(v) for v in q["tile_sha"])
        d=historical(q,L);F,_=features(d)
        if r["name"].startswith("cal_"):cal.append((d,F))
        else:train.append((d,F))
    assert len(train)==8 and len(cal)==2
    X=np.concatenate([F[q["r40c_fused_alert"]] for q,F in train]);y=np.concatenate([q["r40c_fused_localized"][q["r40c_fused_alert"]] for q,F in train]);assert y.any() and (~y).any()
    lo,hi=np.quantile(X,[.01,.99],axis=0);X=np.clip(X,lo,hi);sc=StandardScaler().fit(X);opts=S["model_options"]
    model=LogisticRegression(**opts).fit(sc.transform(X),y)
    assert model.n_iter_.max()<opts["max_iter"]
    def predict(F):return model.decision_function(sc.transform(np.clip(F,lo,hi)))
    v=np.concatenate([predict(F)[q["r40c_fused_localized"]] for q,F in cal])
    th=float(v.min()) if len(v) else -1e10
    lock={"parent_model_sha256":c.sha(raw),"protocol_sha256":c.sha(P.read_bytes()),"train_alerts":len(y),"train_localized":int(y.sum()),"cal_localized":len(v),"threshold":th,"lo":lo.tolist(),"hi":hi.tolist(),"mean":sc.mean_.tolist(),"scale":sc.scale_.tolist(),"coef":model.coef_[0].tolist(),"intercept":float(model.intercept_[0])}
    for q,F in cal:
        m,_,_=compare(q,predict(F),th);assert m["r40f_screen"]["lost_vs_r40c"]==0 and m["r40f_protected"]["lost_vs_r40c"]==0
    save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes());emit("R40F_LOCK",{"sha256":lh,"train_alerts":len(y),"train_localized":int(y.sum()),"cal_localized":len(v),"threshold":th})
    rows=[]
    for target,a in S["artifacts"].items():
        z,n=c.artifact(a,128);res=json.loads(z.read(n["result.json"]));assert res["target"]==target and res["model_sha256"]==S["model"]["model_sha256"] and not roots.intersection(res["source_tile_hashes"])
        for kind in ("eval","shared"):
            q=load(z,n,kind+".npz");F,_=features(q);score=predict(F);m,arms,agree=compare(q,score,th)
            for k in ("r28_frozen","r40c_fused"):
                for name in ("localized","weak","control_alerts"):assert m[k][name]==res["results"][kind][k][name]
            # Compare actual graph scores, paths and decisions, not just counts.
            graph={}
            for k in ("r40d_static","r40d_adaptive"):
                graph[k]={"score_equal":bool(np.array_equal(q[k+"_score"],q["r40c_fused_score"])),"path_equal":bool(np.array_equal(q[k+"_path"],q["r40c_fused_path"])),"alert_equal":bool(np.array_equal(q[k+"_alert"],q["r40c_fused_alert"])),"max_score_delta":float(np.max(abs(q[k+"_score"]-q["r40c_fused_score"])))}
            row={"target":target,"kind":kind,"metrics":m,"graph":graph};rows.append(row)
            np.savez_compressed(O/(target+"_"+kind+".npz"),screen_score=score,agreement=agree,original_alert=q["r40c_fused_alert"],original_path=q["r40c_fused_path"],original_score=q["r40c_fused_score"],**{k+"_priority_alert":v[0] for k,v in arms.items()})
            emit("R40F_GROUP",row)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def total(kind):
        rr=[v["metrics"] for v in rows if v["kind"]==kind];assert len(rr)==5
        return {k:{m:sum(x[k][m] for x in rr) for m in rr[0][k]} for k in S["arms"]}
    e=total("eval");sh=total("shared")
    assert e["r40c_fused"]["localized"]==198 and e["r40c_fused"]["weak"]==27 and e["r40c_fused"]["control_alerts"]==60 and sh["r40c_fused"]["localized"]==37
    gates={}
    for k in S["arms"][2:]:
        gates[k]={"reduction":e["r40c_fused"]["control_alerts"]-e[k]["control_alerts"],"lost":e[k]["lost_vs_r40c"],"weak_lost":e[k]["weak_lost_vs_r40c"],"shared_lost":sh[k]["lost_vs_r40c"],"gate":e[k]["control_alerts"]<60 and e[k]["lost_vs_r40c"]==0 and sh[k]["lost_vs_r40c"]==0 and e[k]["weak"]>=27}
    d={"status":"DEVELOPMENT_REPLAY_ONLY","eval":e,"shared":sh,"gates":gates,"all_graph_alerts_equal":all(g["alert_equal"] for row in rows for g in row["graph"].values()),"all_graph_paths_equal":all(g["path_equal"] for row in rows for g in row["graph"].values()),"all_graph_scores_equal":all(g["score_equal"] for row in rows for g in row["graph"].values()),"rows":rows,"wall_seconds":time.perf_counter()-st,"lock_sha256":lh,"production_adoption":False,"discovery_relaunch_allowed":False}
    save("result.json",d);emit("R40F_SUMMARY",{k:v for k,v in d.items() if k!="rows"})
if __name__=="__main__":main()
