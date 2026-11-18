#!/usr/bin/env python3
"""R40J guarded stationary secondary track; no replacement of parent decisions."""
import sys
sys.argv[1:]=["secondary"]
import json,pathlib,time
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import r40h_frozen as h
c=h.c;t=h.t;f=h.f
P=pathlib.Path("research/seti/r40j_secondary_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40j");O.mkdir(parents=True,exist_ok=True)
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def emit(k,v):print(k,json.dumps(v,allow_nan=False),flush=True)
def features(V):
    V=np.asarray(V,float);N,T,W=V.shape;assert (T,W)==(16,128)
    med=np.median(V,axis=2,keepdims=True);mad=1.4826*np.median(abs(V-med),axis=2,keepdims=True)
    X=np.clip((V-med)/np.maximum(.5,mad),-4,4)
    centers=np.arange(5,123);p=X[:,:,centers];side=X[:,:,centers[:,None]+np.array([-5,-4,-3,3,4,5])]
    local=p-np.median(side,axis=3);local=np.clip(local,-4,4)
    block=local.reshape(N,4,4,len(centers)).mean(2)
    prev=np.median(X,axis=1)[:,centers]
    F=np.stack([p.mean(1),p.std(1),np.median(p,axis=1),np.quantile(p,.25,axis=1),np.quantile(p,.75,axis=1),(p>0).mean(1),local.mean(1),np.median(local,axis=1),np.sort(block,axis=1)[:,1:3].mean(1),block.min(1),block.max(1),block.std(1),local[:,::2].mean(1),local[:,1::2].mean(1),prev,np.max(abs(p),axis=1)/np.maximum(.01,np.sum(abs(p),axis=1))],axis=2)
    fixed=.5*local.mean(1)+.5*np.sort(block,axis=1)[:,1:3].mean(1)
    assert np.isfinite(F).all() and np.isfinite(fixed).all()
    return F,fixed,centers
def candidate(V,scale,mdl):
    a,b,centers=features(V);n,k,d=a.shape
    learned=mdl.decision_function(scale.transform(a.reshape(-1,d))).reshape(n,k)
    out={}
    for name,s in (("fixed",b),("learned",learned)):
        j=s.argmax(1);out[name]=(s[np.arange(n),j],np.repeat(centers[j,None],16,axis=1))
    return out
def base(q,L):
    x=f.historical(q,L["r40c"]);_,_,_,a=h.replay(x,L)
    return a,x["r40c_fused_path"],a&x["r40c_fused_localized"]
def correct(q,p):
    act=q["active"].astype(bool);return (((abs(p-q["truth"])<=2)&act).sum(1)>=.75*act.sum(1))&(q["y"]==1)
def met(q,a,p,hit,cands,th):
    weak=np.isin(q["amp"],t.S["weak_amplitudes"]);null=q["y"]==0
    out={};details={}
    def item(aa,hh,second):
        return {"localized":int(hh.sum()),"weak":int((hh&weak).sum()),"control_alerts":int((aa&null).sum()),"lost_parent":int((hit&~hh).sum()),"added_localized":int((hh&~hit).sum()),"added_weak":int((hh&~hit&weak).sum()),"additional_candidate_paths":int(second.sum()),"additional_control_paths":int((second&null).sum()),"positives":int((~null).sum()),"controls":int(null.sum())}
    out["r40g"]=item(a,hit,np.zeros(len(a),bool))
    for k,(s,pp) in cands.items():
        # On rejected rows only: exactly one output path per alerted tile.
        rescue=(~a)&(s>=th[k]["primary"])
        aa=a|rescue;hh=hit|(rescue&correct(q,pp))
        out[k+"_primary"]=item(aa,hh,rescue);details[k+"_primary"]=(aa,hh,rescue,pp,s)
        # Two-path review lane never replaces primary; account for each added candidate.
        different=(abs(pp-p)>2).sum(1)>=5
        extra=(s>=th[k]["secondary"])&((~a)|different)
        aa=a|extra;hh=hit|(extra&correct(q,pp))
        out[k+"_secondary"]=item(aa,hh,extra);details[k+"_secondary"]=(aa,hh,extra,pp,s)
    return out,details
def main():
    st=time.perf_counter();L=h.readlocks();z,n,rec,_=t.parent.restore()
    train=[];cal=[]
    for rr in rec["records"]:
        q=t.parent.load(z,n,rr)
        if rr["name"].startswith("cal_"):cal.append(q)
        else:train.append(q)
    assert len(train)==8 and len(cal)==2
    X=[];Y=[]
    for q in train:
        F,_,centers=features(q["V"]);act=q["active"].astype(bool)
        yes=(((abs(centers[None,:,None]-q["truth"][:,None,:])<=2)&act[:,None,:]).sum(2)>=.75*act.sum(1)[:,None])&(q["y"][:,None]==1)
        X.append(F.reshape(-1,F.shape[-1]));Y.append(yes.ravel())
    X=np.concatenate(X);Y=np.concatenate(Y);assert Y.any() and (~Y).any()
    scale=StandardScaler().fit(X);mdl=LogisticRegression(C=.1,class_weight="balanced",max_iter=1000,random_state=401).fit(scale.transform(X),Y)
    assert mdl.n_iter_.max()<1000
    nulls={k:{"primary":[],"secondary":[]} for k in ("fixed","learned")};base_controls=0
    for q in cal:
        a,p,hit=base(q,L);base_controls+=int((a&(q["y"]==0)).sum())
        for k,(s,pp) in candidate(q["V"],scale,mdl).items():
            nulls[k]["primary"].extend(s[(q["y"]==0)&~a].tolist())
            nulls[k]["secondary"].extend(s[q["y"]==0].tolist())
    assert base_controls==10
    th={k:{j:float(np.nextafter(max(v),np.inf)) for j,v in x.items()} for k,x in nulls.items()}
    lock={"protocol_sha256":c.sha(P.read_bytes()),"thresholds":th,"base_calibration_alerts":base_controls,"train_columns":len(Y),"positive_columns":int(Y.sum()),"mean":scale.mean_.tolist(),"scale":scale.scale_.tolist(),"coef":mdl.coef_.tolist(),"intercept":mdl.intercept_.tolist()};save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes());emit("R40J_LOCK",{k:v for k,v in lock.items() if k not in ("mean","scale","coef","intercept")})
    for q in cal:
        a,p,hit=base(q,L);m,_=met(q,a,p,hit,candidate(q["V"],scale,mdl),th)
        assert all(x["additional_control_paths"]==0 for x in m.values())
    rows=[]
    def run(q,a,p,hit,target,kind,source):
        ca=candidate(q["V"],scale,mdl);m,detail=met(q,a,p,hit,ca,th)
        assert all(v["lost_parent"]==0 for v in m.values())
        row={"source":source,"target":target,"kind":kind,"metrics":m};rows.append(row);emit("R40J_GROUP",row)
        np.savez_compressed(O/(source+"_"+target+"_"+kind+".npz"),parent_alert=a,parent_path=p,parent_hit=hit,**{k+"_"+j:v for k,vals in detail.items() for j,v in zip(("alert","localized","extra","path","score"),vals)})
    for target,art in f.S["artifacts"].items():
        zz,nn=c.artifact(art,128)
        for kind in ("eval","shared"):
            q=f.load(zz,nn,kind+".npz");_,_,_,a=h.replay(q,L);p=q["r40c_fused_path"];hit=a&q["r40c_fused_localized"]
            run(q,a,p,hit,target,kind,"r40cde")
    zz,nn=c.artifact(S["source"],192);rec=json.loads(zz.read(nn["receipt.json"]))
    for rr in rec["records"]:
        raw=zz.read(nn[rr["name"]]);assert c.sha(raw)==rr["sha256"]
        q=f.load(zz,nn,rr["name"]);assert c.sha(q["V"].tobytes())==rr["array_sha256"]
        kind,target=rr["name"].removesuffix(".npz").split("_",1)
        run(q,q["r40g_alert"],q["r40g_path"],q["r40g_hit"],target,kind,"r40h")
    def tot(source,kind):
        rr=[v["metrics"] for v in rows if v["source"]==source and v["kind"]==kind];assert len(rr)==5
        return {k:{j:sum(x[k][j] for x in rr) for j in rr[0][k]} for k in rr[0]}
    totals={source:{kind:tot(source,kind) for kind in ("eval","shared")} for source in ("r40cde","r40h")}
    e=totals["r40h"]["eval"];sh=totals["r40h"]["shared"];assert [e["r40g"][k] for k in ("localized","weak","control_alerts")]==[213,41,31] and sh["r40g"]["weak"]==18
    gates={k:{"shared_weak":sh[k]["weak"],"controls":e[k]["control_alerts"],"gate":sh[k]["weak"]>=21 and e[k]["control_alerts"]<=35 and e[k]["lost_parent"]==0 and sh[k]["lost_parent"]==0} for k in e if k!="r40g"}
    assert c.sha((O/"lock.json").read_bytes())==lh
    d={"status":"DEVELOPMENT_REPLAY_ONLY","totals":totals,"gates":gates,"lock_sha256":lh,"rows":rows,"seconds":time.perf_counter()-st,"production_adoption":False}
    save("result.json",d);emit("R40J_SUMMARY",{k:v for k,v in d.items() if k!="rows"})
if __name__=="__main__":main()
