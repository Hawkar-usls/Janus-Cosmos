#!/usr/bin/env python3
"""R40CDE: freeze one historical model, new-frequency injection transfer, and bounded graph/null audits."""
import hashlib,io,json,os,pathlib,sys,time,zipfile
import numpy as np
import requests
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import r40b_transfer as t
P=pathlib.Path("research/seti/r40cde_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r40cde")/sys.argv[1];OUT.mkdir(parents=True,exist_ok=True)
def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def artifact(a,cap=128):
    u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip"
    r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=120)
    r.raise_for_status();raw=r.content;assert len(raw)<cap*1024**2 and sha(raw)==a["sha256"]
    z=zipfile.ZipFile(io.BytesIO(raw));entries=[i for i in z.infolist() if not i.is_dir()]
    assert sum(i.file_size for i in entries)<512*1024**2
    n={pathlib.PurePosixPath(i.filename).name:i.filename for i in entries};assert len(n)==len(entries)
    return z,n
def graph(logits,paths,iters):
    # Boundary potential 1 at source, 0 at sink. Path nodes connected through soft geometric agreement.
    G=S["graph"];n,k=logits.shape;assert k==8
    prior=np.exp(np.clip(logits-logits.max(1,keepdims=True),-12,0))+.01
    sim=np.mean(abs(paths[:,:,None,:]-paths[:,None,:,:])<=2,axis=3)
    A=G["similarity_floor"]+.5*sim
    A[:,np.arange(k),np.arange(k)]=0
    D=prior.copy();E=np.ones_like(D)
    for j in range(iters+1):
        M=-A.copy();M[:,np.arange(k),np.arange(k)]=D+E+A.sum(2)
        v=np.linalg.solve(M,D[...,None])[...,0]
        Q=D*(1-v);R=E*v;F=A*(v[:,:,None]-v[:,None,:])
        assert np.max(abs(Q-R-F.sum(2)))<1e-8 and np.all(v>=-1e-8) and np.all(v<=1+1e-8)
        if j==iters:break
        def update(d,f):return np.maximum(G["conductance_floor"],d+G["dt"]*(G["growth"]*abs(f)-G["decay"]*d))
        D=update(D,Q);E=update(E,R);A=update(A,F);A[:,np.arange(k),np.arange(k)]=0
    score=logits+G["score_log_weight"]*np.log(np.maximum(Q,1e-12)/np.maximum(Q.max(1,keepdims=True),1e-12))
    i=score.argmax(1);return score[np.arange(n),i],paths[np.arange(n),i]
def base_predictions(q,F,p,lock):
    m=np.asarray(lock["mean"]);d=np.asarray(lock["scale"]);c=np.asarray(lock["coef"]);a=((F-m)/d)@c+lock["intercept"]
    assert np.isfinite(a).all()
    j=a.argmax(1);ix=np.arange(len(j))
    return {"r40c_fused":(a[ix,j],p[ix,j]),"r40d_static":graph(a,p,0),"r40d_adaptive":graph(a,p,S["graph"]["iterations"])}
def fused(q,s,p,L):
    return t.fuse(s,p,q["hybrid_r28_score"],q,L)
def freeze():
    z,n,rec,old=t.parent.restore()
    assert S["source"]==t.S["source"]
    ms=json.loads(z.read(n["manifest.json"]))
    assert not {v["target"] for v in ms["selected"]}&set(S["targets"])
    assert not {int(v["day"]) for v in ms["selected"]}&set(S["days"])
    entries={}
    for r in rec["records"]:
        q=t.parent.load(z,n,r);F,p=t.feature(q);g=t.group(r);k=r["name"].split("_",1)[0];entries[k,g]=(q,F,p)
    groups=sorted(g for k,g in entries if k=="eval");cal=sorted(g for k,g in entries if k=="cal");assert len(groups)==4 and len(cal)==2
    roots={}
    for (k,g),(q,F,p) in entries.items():roots.setdefault(g,set()).update(str(v) for v in q["tile_sha"])
    assert all(roots[a].isdisjoint(roots[b]) for i,a in enumerate(roots) for b in list(roots)[i+1:])
    tr=[entries[k,g] for g in groups for k in ("eval","shared")]
    X=np.concatenate([F.reshape(-1,F.shape[-1]) for q,F,p in tr]);y=np.concatenate([t.targets(q,p).reshape(-1) for q,F,p in tr])
    opts=t.S["model"];sc=StandardScaler().fit(X);mdl=LogisticRegression(C=opts["C"],class_weight=opts["class_weight"],max_iter=opts["max_iter"],random_state=opts["random_state"]).fit(sc.transform(X),y)
    assert mdl.n_iter_.max()<opts["max_iter"]
    L={"mean":sc.mean_.tolist(),"scale":sc.scale_.tolist(),"coef":mdl.coef_[0].tolist(),"intercept":float(mdl.intercept_[0]),"train_groups":groups,"cal_groups":cal,"train_rows":len(y),"source":S["source"],"protocol_sha256":sha(P.read_bytes()),"original":old,"source_roots":sorted(set().union(*roots.values()))}
    tests=[]
    for g in cal:
        q,F,p=entries["cal",g];pred=base_predictions(q,F,p,L)
        a=mdl.decision_function(sc.transform(F.reshape(-1,F.shape[-1]))).reshape(F.shape[:2])
        assert np.allclose(a.max(1),pred["r40c_fused"][0],atol=1e-10)
        tests.append((g,q,pred))
    bs=np.concatenate([q["hybrid_r28_score"][q["y"]==0] for g,q,p in tests])
    bt=t.parent.cut(bs);assert abs(bt-old["thresholds"]["hybrid_r28"])<1e-6
    L["thresholds"]={"r28_frozen":float(old["thresholds"]["hybrid_r28"])};L["margins"]={};L["calibration"]={};L["cal_group_maxima"]={}
    for k in S["arms"][1:]:
        cs=np.concatenate([p[k][0][q["y"]==0] for g,q,p in tests]);ct=t.parent.cut(cs)
        m={"base_center":bt,"base_scale":max(.1,bt-float(np.quantile(bs,.8))),"learned_center":ct,"learned_scale":max(.1,ct-float(np.quantile(cs,.8)))}
        fs=np.concatenate([fused(q,*p[k],m)[0][q["y"]==0] for g,q,p in tests]);th=t.parent.cut(fs)
        L["thresholds"][k]=th;L["margins"][k]=m;L["calibration"][k]=int((fs>=th).sum())
        L["cal_group_maxima"][k]=[float(fused(q,*p[k],m)[0][q["y"]==0].max()) for g,q,p in tests]
    L["cal_group_maxima"]["r28_frozen"]=[float(q["hybrid_r28_score"][q["y"]==0].max()) for g,q,p in tests]
    assert max(L["calibration"].values())<=12
    save("model.json",L);modelsha=sha((OUT/"model.json").read_bytes())
    emit("R40C_FREEZE",{"model_sha256":modelsha,"train_groups":groups,"cal_groups":cal,"train_rows":len(y),"thresholds":L["thresholds"],"calibration":L["calibration"],"graph_flow_conservation":"PASS","single_model":True})
    # R40E: group maxima empirical ranks, not valid calibrated probabilities absent exchangeability.
    ncal=len(cal);planned=S["r40e"]["planned_campaign_groups"];alpha=S["r40e"]["alpha"]
    # exact disjoint-uniform intervals and identical events both obey union bound
    u=np.arange(10000);events=np.stack([(u>=i*100)&(u<(i+1)*100) for i in range(planned)])
    correlated=np.repeat(events[:1],planned,axis=0)
    assert events.any(0).sum()==500 and correlated.any(0).sum()==100
    e={"status":"HOLD_INSUFFICIENT_EXCHANGEABLE_NULL_GROUPS","alpha":alpha,"planned_groups":planned,"per_group_alpha":alpha/planned,"calibration_groups":ncal,"min_empirical_p":1/(ncal+1),"minimum_groups_for_rank_resolution":int(np.ceil(planned/alpha))-1,"dependency_canary":{"disjoint_union":float(events.any(0).mean()),"identical_union":float(correlated.any(0).mean()),"sum_bound":float(events.mean(1).sum()),"PASS":True},"formal_FWER_proved_for_archive":False,"reason":"Null labels and exchangeability of real observations not established. Do not multiply votes or overlapping patches. Resolution alone is not sufficient for validity.","production_adoption":False}
    save("campaign_gate.json",e);emit("R40E_SUMMARY",e)
def evaluate():
    target=sys.argv[2];assert target in S["targets"]
    af={"id":int(os.environ["MODEL_ID"]),"sha256":os.environ["MODEL_DIGEST"].removeprefix("sha256:")}
    z,n=artifact(af,8);raw=z.read(n["model.json"]);L=json.loads(raw);lh=sha(raw)
    assert L["protocol_sha256"]==sha(P.read_bytes()) and L["source"]==S["source"]
    (OUT/"model.json").write_bytes(raw)
    # Import existing exact R35G feature-generating pipeline; no new detector configuration.
    import r35b_grand_benchmark as b
    b.S["seed"]=S["seed"];b.ref.S=b.S;b.ref.OUT=OUT;b.sh.eng.OUT=OUT
    b.sh.eng.SPEC={**b.sh.eng.SPEC,"window_centers_mhz":S["windows_mhz"],"window_half_hz":S["half_hz"],"limits":{**b.sh.eng.SPEC["limits"],"spectrum_budget_per_file_mib":S["max_mib_per_scan"]}}
    assert b.sh.eng.SPEC["window_centers_mhz"]==S["windows_mhz"]
    mz,mn=artifact(S["manifest"],32);m=json.loads(mz.read(mn["manifest.json"]))
    assert m["new_obs_power_read"] is False
    group=[v for v in m["selected"]["Parkes"] if v["target"]==target];assert len(group)==1
    q=group[0];assert int(q["day"])==S["days"][S["targets"].index(target)]
    assert [v["role"] for v in q["files"]]==["ON1","OFF1","ON2","OFF2"]
    assert target not in L["train_groups"]+L["cal_groups"]
    save("manifest.json",q)
    b.taylor.init_tree()
    rows,R,ch,hs=b.ref.gather(q);assert not set(L["source_roots"])&hs and len(rows)==256 and len(ch)==32
    def encode(rows,R):
        d=b.panel(rows,R);arms=b.base_arms(d,L["original"]["component_fit"])
        v={"V":d["Y"].astype(np.float32),"y":b.labels(d),"truth":np.stack([r[7] for r in rows]),"active":np.stack([r[8] for r in rows]),"amp":np.asarray([r[4] for r in rows]),"tile_sha":np.asarray([r[6] for r in rows])}
        for k in t.S["paths"]:
            v[k+"_score"],v[k+"_paths"]=arms[k][:2]
        return v
    results={}
    for kind,qq in (("eval",encode(rows,R)),("shared",encode([v[0] for v in ch],np.stack([v[1] for v in ch])))):
        F,p=t.feature(qq);pr=base_predictions(qq,F,p,L)
        vals={"r28_frozen":(qq["hybrid_r28_score"],qq["hybrid_r28_paths"])}
        vals.update({k:fused(qq,*pp,L["margins"][k]) for k,pp in pr.items()})
        per={};hits={};alerts={};arr={k:v for k,v in qq.items() if k in ("y","truth","active","amp","tile_sha","V")}
        for k,(s,p) in vals.items():
            per[k],hits[k],alerts[k]=t.metrics(qq,s,p,L["thresholds"][k]);arr.update({k+"_score":s,k+"_path":p,k+"_alert":alerts[k],k+"_localized":hits[k]})
            if kind=="eval":
                mx=float(s[qq["y"]==0].max());cal=L["cal_group_maxima"][k];per[k]["group_max_control_score"]=mx;per[k]["descriptive_group_tail_rank"]=(1+sum(v>=mx for v in cal))/(1+len(cal))
        base=hits["r28_frozen"];ba=alerts["r28_frozen"];weak=np.isin(qq["amp"],t.S["weak_amplitudes"]);null=qq["y"]==0
        for k in per:
            h=hits[k];a=alerts[k];per[k].update({"added":int((h&~base).sum()),"lost":int((base&~h).sum()),"weak_added":int((h&~base&weak).sum()),"weak_lost":int((base&~h&weak).sum()),"controls_added":int((a&~ba&null).sum()),"controls_removed":int((ba&~a&null).sum())})
        np.savez_compressed(OUT/(kind+".npz"),**arr);results[kind]=per
    assert sha((OUT/"model.json").read_bytes())==lh
    d={"status":"NEW_FREQUENCY_INJECTIONS_REVIEW_ONLY","target":target,"model_sha256":lh,"model_artifact":af,"new_windows_mhz":S["windows_mhz"],"source_tile_hashes":sorted(hs),"results":results,"scope":S["limits"],"production_adoption":False}
    save("result.json",d);emit("R40CDE_TARGET",d)
if __name__=="__main__":
    if sys.argv[1]=="freeze":freeze()
    elif sys.argv[1]=="evaluate":evaluate()
    else:raise ValueError("invalid mode")
