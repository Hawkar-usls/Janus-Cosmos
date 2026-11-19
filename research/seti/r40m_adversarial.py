#!/usr/bin/env python3
"""R40M: explicit rival hypotheses on frozen compact R40L evidence."""
import json,hashlib,io,os,pathlib,time,zipfile
import numpy as np
import requests
from sklearn.ensemble import ExtraTreesClassifier
P=pathlib.Path("research/seti/r40m_adversarial_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40m");O.mkdir(parents=True,exist_ok=True)
def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def artifact():
    a=S["source"];r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=120)
    r.raise_for_status();raw=r.content;assert len(raw)<8*1024**2 and sha(raw)==a["sha256"]
    z=zipfile.ZipFile(io.BytesIO(raw));entries=[i for i in z.infolist() if not i.is_dir()];assert sum(i.file_size for i in entries)<32*1024**2
    n={pathlib.PurePosixPath(i.filename).name:i.filename for i in entries};assert len(n)==len(entries)
    return z,n
def load(z,n,name):
    with np.load(io.BytesIO(z.read(n[name])),allow_pickle=False) as q:d={k:q[k].copy() for k in q.files}
    F=d["F"];assert F.shape==(len(d["y"]),30) and np.isfinite(F).all()
    assert np.array_equal(F[:,24],d["s"]) and np.array_equal(F[:,25],d["support"])
    return d
def rivals(F):
    # Stored features: ON indices 0..11, OFF 12..23, joint evidence 24..29.
    # Larger values are evidence against the proposed rivals, not proof of identity.
    concentration=np.maximum(F[:,9]-(F[:,7]+F[:,8])/2,F[:,21]-(F[:,19]+F[:,20])/2)
    return np.stack([np.min(F[:,[10,11,22,23]],axis=1),np.minimum(F[:,7],F[:,19]),-F[:,27],F[:,28],-concentration],axis=1)
def mask(d):return (~d["parent_alert"].astype(bool))&(d["y"]==0)
def fixed(d,L):
    a=d["parent_alert"].astype(bool);base=d["parent_hit"].astype(bool);res=(~a)&(d["s"]>=L["parent_joint_min_threshold"])
    assert np.array_equal(res,d["joint_rescue"])
    assert np.array_equal(base|(res&d["correct"]),d["joint_min_hit"])
    F=d["F"];ls=((F-np.asarray(L["mean"]))/np.asarray(L["scale"]))@np.asarray(L["coef"])[0]+L["intercept"][0]
    assert np.allclose(ls,d["learned_score"],atol=1e-10,rtol=1e-10)
    top=d["support"]>=L["thresholds"]["topa05"];learn=ls>=L["thresholds"]["learned05"]
    assert np.array_equal(a|(res&top),d["topa05_alert"]) and np.array_equal(a|(res&learn),d["learned05_alert"])
    return a,base,res,top,learn
def main():
    st=time.perf_counter();z,n=artifact();raw=z.read(n["lock.json"]);assert sha(raw)==S["source"]["lock_sha256"];L=json.loads(raw)
    train=[];cal=[];ev=[]
    for name in sorted(n):
        if not name.endswith(".npz"):continue
        if name.startswith("historical_cal_"):cal.append((name,load(z,n,name)))
        elif name.startswith(("historical_eval_","historical_shared_")):train.append((name,load(z,n,name)))
        elif name.startswith(("eval_","shared_")):ev.append(name)
    assert len(train)==8 and len(cal)==2 and len(ev)==10
    for name,d in train+cal:fixed(d,L)
    null=np.concatenate([rivals(d["F"])[mask(d)] for _,d in train]);assert len(null)>0
    center=np.median(null,axis=0);scale=np.maximum(.01,np.quantile(null,.75,axis=0)-np.quantile(null,.25,axis=0))
    def gs(d):return np.min((rivals(d["F"])-center)/scale,axis=1)
    def feats(d):return np.column_stack([d["F"],rivals(d["F"])])
    X=np.concatenate([feats(d)[~d["parent_alert"]] for _,d in train]);Y=np.concatenate([d["correct"][~d["parent_alert"]] for _,d in train])
    assert Y.any() and (~Y).any()
    mdl=ExtraTreesClassifier(n_estimators=128,max_depth=4,min_samples_leaf=8,class_weight="balanced",random_state=401,n_jobs=1).fit(X,Y)
    def ts(d):return mdl.predict_proba(feats(d))[:,1]
    def pos(d):return d["joint_rescue"]&d["correct"]
    gp=np.concatenate([gs(d)[pos(d)] for _,d in train]);tp=np.concatenate([ts(d)[pos(d)] for _,d in train]);gc=np.concatenate([gs(d)[mask(d)] for _,d in cal])
    assert len(gp)==len(tp)>0 and len(gc)==246
    th={"ghost_preserve":float(np.nextafter(gp.min(),-np.inf)),"tree_preserve":float(np.nextafter(tp.min(),-np.inf)),"ghost_null025":float(np.nextafter(np.sort(gc)[::-1][int(.025*len(gc))],np.inf))}
    trees=[{"children_left":e.tree_.children_left.tolist(),"children_right":e.tree_.children_right.tolist(),"feature":e.tree_.feature.tolist(),"threshold":e.tree_.threshold.tolist(),"value":e.tree_.value.tolist()} for e in mdl.estimators_]
    lock={"protocol_sha256":sha(P.read_bytes()),"source_lock_sha256":sha(raw),"thresholds":th,"center":center.tolist(),"scale":scale.tolist(),"train_rows":len(Y),"train_positive":int(Y.sum()),"train_joint_localized_rescues":len(gp),"calibration_parent_rejected_nulls":len(gc),"trees":trees}
    save("lock.json",lock);lh=sha((O/"lock.json").read_bytes());emit("R40M_LOCK",{k:v for k,v in lock.items() if k!="trees"})
    def met(d):
        a,base,res,top,learn=fixed(d,L);g=gs(d);tree=ts(d);p=g>=th["ghost_preserve"];n=g>=th["ghost_null025"];tr=tree>=th["tree_preserve"]
        guards={"r40g":np.zeros(len(a),bool),"joint_min":np.ones(len(a),bool),"topa05":top,"learned05":learn,"topa_learned":top&learn,"ghost_preserve":p,"learned_ghost_preserve":learn&p,"topa_learned_ghost_preserve":top&learn&p,"ghost_null025":n,"learned_ghost_null025":learn&n,"ghost_tree_preserve":tr}
        y=d["y"];weak=d["weak"];old=d["joint_min_hit"];out={};arr={"ghost_score":g,"tree_score":tree,**d}
        for k,v in guards.items():
            rescue=res&v;aa=a|rescue;hh=base|(rescue&d["correct"])
            out[k]={"localized":int(hh.sum()),"weak":int((hh&weak).sum()),"control_alerts":int((aa&(y==0)).sum()),"lost_parent":int((base&~hh).sum()),"lost_joint":int((old&~hh).sum()),"lost_joint_weak":int((old&~hh&weak).sum()),"additional_control_alerts":int((rescue&(y==0)).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum())}
            assert out[k]["lost_parent"]==0
            arr[k+"_alert"]=aa;arr[k+"_hit"]=hh
        return out,arr
    for name,d in cal:
        m,_=met(d);assert all(v["additional_control_alerts"]==0 for v in m.values())
    rows=[]
    for name in ev:
        d=load(z,n,name);m,arr=met(d);row={"record":name,"metrics":m};rows.append(row);emit("R40M_GROUP",row);np.savez_compressed(O/name,**arr)
    def total(kind):
        rr=[v["metrics"] for v in rows if v["record"].startswith(kind+"_")];assert len(rr)==5
        return {k:{j:sum(x[k][j] for x in rr) for j in rr[0][k]} for k in rr[0]}
    e=total("eval");sh=total("shared")
    assert [e["joint_min"][k] for k in ("localized","weak","control_alerts")]==[215,43,43] and [sh["joint_min"][k] for k in ("localized","weak")]==[111,36]
    assert e["topa05"]["control_alerts"]==39 and e["learned05"]["control_alerts"]==40 and sh["learned05"]["localized"]==111
    gates={k:{"strict":e[k]["localized"]>=215 and e[k]["weak"]>=43 and e[k]["control_alerts"]<=35 and sh[k]["localized"]>=111 and sh[k]["weak"]>=36 and e[k]["lost_joint"]==sh[k]["lost_joint"]==0,"minimum":e[k]["localized"]>=214 and e[k]["weak"]>=41 and e[k]["control_alerts"]<=35 and sh[k]["weak"]>=21 and e[k]["lost_parent"]==sh[k]["lost_parent"]==0} for k in e}
    assert sha((O/"lock.json").read_bytes())==lh
    d={"status":"DEVELOPMENT_REPLAY_ONLY","eval":e,"shared":sh,"gates":gates,"lock_sha256":lh,"rows":rows,"seconds":time.perf_counter()-st,"production_adoption":False}
    save("result.json",d);emit("R40M_SUMMARY",{k:v for k,v in d.items() if k!="rows"})
if __name__=="__main__":main()
