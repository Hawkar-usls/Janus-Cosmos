#!/usr/bin/env python3
"""R40B: learned path evidence on historical source-separated folds; no production promotion."""
import hashlib,json,pathlib,time
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import r40a_sensitivity as parent
P=pathlib.Path("research/seti/r40b_transfer_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r40b");OUT.mkdir(parents=True,exist_ok=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def feature(q):
    V=q["V"].astype(float);n,t,w=V.shape;assert (t,w)==(16,128)
    keys=S["paths"];paths=np.stack([q[k+"_paths"] for k in keys],axis=1).astype(int)
    assert np.all((paths>=0)&(paths<w))
    ix=np.arange(n)[:,None];ti=np.arange(t)[None,:];out=[]
    for j,k in enumerate(keys):
        p=paths[:,j];trace=V[ix,ti,p];d=np.diff(p,axis=1)
        vote=(abs(paths-p[:,None,:])<=2).sum(2)>=12
        vote[:,j]=False
        side=p[:,:,None]+np.array([-5,-4,-3,3,4,5])
        ok=(side>=0)&(side<w)
        v=V[ix[:,:,None],ti[:,:,None],np.clip(side,0,w-1)]
        bg=np.nanmedian(np.where(ok,v,np.nan),axis=2)
        dif=trace-bg
        block=trace.reshape(n,4,4).mean(2)
        vals=[q[k+"_score"],vote.sum(1),trace.mean(1),trace.std(1),np.median(trace,axis=1),np.quantile(trace,.25,axis=1),np.quantile(trace,.75,axis=1),trace.max(1), (trace>0).mean(1),dif.mean(1),np.median(dif,axis=1),block.min(1),block.max(1),block.std(1),np.mean(abs(d),axis=1),np.mean(abs(np.diff(d,axis=1)),axis=1),abs(p[:,-1]-p[:,0]),np.max(np.median(V,axis=1),axis=1),np.std(np.median(V,axis=1),axis=1)]
        # Include arm identity; the learner may distinguish fixed path generators.
        vals += [np.full(n,float(a==j)) for a in range(len(keys))]
        out.append(np.stack(vals,axis=1))
    F=np.stack(out,axis=1);assert np.isfinite(F).all()
    return F,paths
def targets(q,p):
    active=q["active"].astype(bool);truth=q["truth"]
    return (((abs(p-truth[:,None,:])<=2)&active[:,None,:]).sum(2)>=.75*active.sum(1)[:,None])&(q["y"][:,None]==1)
def predict(F,p,scale,model):
    n,k,d=F.shape;a=model.decision_function(scale.transform(F.reshape(-1,d))).reshape(n,k);j=a.argmax(1)
    return a[np.arange(n),j],p[np.arange(n),j]
def fuse(s,p,b,q,L):
    a=(b-L["base_center"])/L["base_scale"];c=(s-L["learned_center"])/L["learned_scale"];choose=c>a
    return np.maximum(a,c),np.where(choose[:,None],p,q["hybrid_r28_paths"])
def metrics(q,s,p,t):
    y=q["y"];act=q["active"].astype(bool);truth=q["truth"]
    pred=s>=t;hit=pred&(((abs(p-truth)<=2)&act).sum(1)>=.75*act.sum(1))&(y==1)
    weak=np.isin(q["amp"],S["weak_amplitudes"])
    return {"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"control_alerts":int((pred&(y==0)).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum())},hit,pred
def group(r):return r["name"].split("_",1)[1].removesuffix(".npz")
def main():
    start=time.perf_counter();z,n,rec,old=parent.restore()
    assert S["source"]==parent.S["source"]
    entries={}
    for r in rec["records"]:
        q=parent.load(z,n,r);F,p=feature(q)
        prefix=r["name"].split("_",1)[0];g=group(r)
        assert (prefix,g) not in entries
        entries[prefix,g]=(q,F,p)
    groups=sorted(g for k,g in entries if k=="eval");cal=sorted(g for k,g in entries if k=="cal")
    assert len(groups)==4 and len(cal)==2 and not set(groups)&set(cal)
    assert all(("shared",g) in entries for g in groups)
    # Source hashes must be separated at the group level, including shared-reference stress rows.
    roots={}
    for (k,g),(q,F,p) in entries.items():
        roots.setdefault(g,set()).update(str(v) for v in q["tile_sha"])
    assert all(roots[a].isdisjoint(roots[b]) for i,a in enumerate(roots) for b in list(roots)[i+1:])
    base_threshold=float(old["thresholds"]["hybrid_r28"])
    rows=[]
    for held in groups:
        train=[entries[k,g] for g in groups if g!=held for k in ("eval","shared")]
        X=np.concatenate([F.reshape(-1,F.shape[-1]) for q,F,p in train]);y=np.concatenate([targets(q,p).reshape(-1) for q,F,p in train])
        assert y.any() and (~y).any()
        scale=StandardScaler().fit(X);opts=S["model"];model=LogisticRegression(C=opts["C"],class_weight=opts["class_weight"],max_iter=opts["max_iter"],random_state=opts["random_state"]).fit(scale.transform(X),y)
        assert model.n_iter_.max()<opts["max_iter"]
        ca=[]
        for g in cal:
            q,F,p=entries["cal",g];s,pp=predict(F,p,scale,model);ca.append((q,s,pp))
        bs=np.concatenate([q["hybrid_r28_score"][q["y"]==0] for q,s,p in ca]);cs=np.concatenate([s[q["y"]==0] for q,s,p in ca])
        assert len(bs)==len(cs)==256
        bc=parent.cut(bs);ct=parent.cut(cs);assert abs(bc-base_threshold)<1e-6
        L={"base_center":bc,"base_scale":max(.1,bc-float(np.quantile(bs,.8))),"learned_center":ct,"learned_scale":max(.1,ct-float(np.quantile(cs,.8)))}
        fs=np.concatenate([fuse(s,p,q["hybrid_r28_score"],q,L)[0][q["y"]==0] for q,s,p in ca]);ft=parent.cut(fs)
        lock={"held":held,"train":[g for g in groups if g!=held],"cal":cal,"thresholds":{"r28_frozen":base_threshold,"learned_path":ct,"r28_plus_learned":ft},"scaler_mean":scale.mean_.tolist(),"scaler_scale":scale.scale_.tolist(),"coef":model.coef_.tolist(),"intercept":model.intercept_.tolist(),"margin":L,"cal_alerts":{"r28":int((bs>=base_threshold).sum()),"learned":int((cs>=ct).sum()),"fused":int((fs>=ft).sum())},"train_path_labels":int(y.sum()),"train_path_count":len(y)}
        assert max(lock["cal_alerts"].values())<=12
        name="lock_"+held+".json";save(name,lock);lh=sha((OUT/name).read_bytes())
        for kind in ("eval","shared"):
            q,F,p=entries[kind,held];st=time.perf_counter();s,pp=predict(F,p,scale,model);f,fp=fuse(s,pp,q["hybrid_r28_score"],q,L)
            vals={"r28_frozen":(q["hybrid_r28_score"],q["hybrid_r28_paths"]),"learned_path":(s,pp),"r28_plus_learned":(f,fp)}
            per={};hits={};alerts={};arr={"y":q["y"],"truth":q["truth"],"active":q["active"],"amp":q["amp"]}
            for k,(ss,path) in vals.items():
                per[k],hits[k],alerts[k]=metrics(q,ss,path,lock["thresholds"][k]);arr.update({k+"_score":ss,k+"_path":path,k+"_alert":alerts[k],k+"_localized":hits[k]})
            assert np.array_equal(hits["r28_frozen"],q["hybrid_r28_localized"]) and np.array_equal(alerts["r28_frozen"],q["hybrid_r28_alert"])
            b=hits["r28_frozen"];a=alerts["r28_frozen"];null=q["y"]==0;weak=np.isin(q["amp"],S["weak_amplitudes"])
            for k in per:
                h=hits[k];v=alerts[k];per[k].update({"added":int((h&~b).sum()),"lost":int((b&~h).sum()),"weak_added":int((h&~b&weak).sum()),"weak_lost":int((b&~h&weak).sum()),"controls_added":int((v&~a&null).sum()),"controls_removed":int((a&~v&null).sum())})
            np.savez_compressed(OUT/(kind+"_"+held+".npz"),**arr)
            assert sha((OUT/name).read_bytes())==lh
            row={"held":held,"kind":kind,"metrics":per,"lock_sha256":lh,"prediction_seconds":time.perf_counter()-st};rows.append(row);save("partial.json",rows);emit("R40B_GROUP",row)
    def total(kind):
        rr=[v["metrics"] for v in rows if v["kind"]==kind];assert len(rr)==4
        return {k:{m:sum(x[k][m] for x in rr) for m in rr[0][k]} for k in S["arms"]}
    e=total("eval");sh=total("shared");b=e["r28_frozen"];bs=sh["r28_frozen"]
    assert b["positives"]==b["controls"]==512 and bs["positives"]==128
    gates={}
    for k in S["arms"][1:]:
        v=e[k];w=sh[k]
        gates[k]={"total_gain":v["localized"]-b["localized"],"weak_gain":v["weak"]-b["weak"],"control_reduction":b["control_alerts"]-v["control_alerts"],"shared_gain":w["localized"]-bs["localized"],"dual_gate":v["localized"]>=b["localized"] and v["weak"]>b["weak"] and v["control_alerts"]<b["control_alerts"] and w["localized"]>=bs["localized"]}
    d={"status":"HISTORICAL_GROUP_HELD_OUT_DEVELOPMENT_ONLY","eval":e,"shared":sh,"gates":gates,"rows":rows,"wall_seconds":time.perf_counter()-start,"root_groups":{g:len(v) for g,v in roots.items()},"production_adoption":False,"discovery_relaunch_allowed":False}
    save("result.json",d);emit("R40B_SUMMARY",{k:v for k,v in d.items() if k!="rows"})
if __name__=="__main__":main()
