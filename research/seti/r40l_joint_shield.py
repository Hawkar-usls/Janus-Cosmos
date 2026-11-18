#!/usr/bin/env python3
"""R40L: protected TOPA/KUSTO and learned joint witness. Research replay."""
import sys
sys.argv[1:]=["joint_shield"]
import json,pathlib,time
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
import r40k_separate as k
h=k.h;c=k.c;t=k.t;f=k.f;b=k.b
P=pathlib.Path("research/seti/r40l_joint_shield_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40l");O.mkdir(parents=True,exist_ok=True);k.O=O
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def cut(v,rate):
    z=np.sort(np.asarray(v,float))[::-1];assert len(z)>0 and np.isfinite(z).all()
    return float(np.nextafter(z[int(len(z)*rate)],np.inf))
def evidence(rows,R,q):
    ca=k.regen(rows,R,q);s,p=ca["joint_min"];X=b.ref.norm(np.stack([v[0] for v in rows]));D=b.ref.norm(R)
    i=np.arange(len(X))[:,None];j=np.arange(16)[None,:]
    def stats(V):
        v=V[i,j,p];side=V[i[:,:,None],j[:,:,None],p[:,:,None]+np.array([-5,-4,-3,3,4,5])]
        u=np.clip(v-np.median(side,axis=2),-4,4)
        block=np.sort(u.reshape(-1,4,4).mean(2),axis=1)
        return u,[u.mean(1),u.std(1),np.median(u,axis=1),np.quantile(u,.25,axis=1),np.quantile(u,.75,axis=1),(u>0).mean(1),block[:,0],block[:,1],block[:,2],block[:,3],u[:,::2].mean(1),u[:,1::2].mean(1)]
    x,xf=stats(X);d,df=stats(D)
    support=np.minimum(xf[7],df[7])
    regime=np.quantile(np.maximum(np.median(D,axis=1),0),.95,axis=1)
    # Independent selected-path temporal pieces; no truth or injection strength in features.
    F=np.stack(xf+df+[s,support,regime,abs(x.mean(1)-d.mean(1)),np.minimum(x,d).mean(1),np.abs(x-d).mean(1)],axis=1)
    assert np.isfinite(F).all()
    return s,p,support,regime,F
def parent(q,L):
    hh=f.historical(q,L["r40c"]);_,_,_,a=h.replay(hh,L)
    return a.astype(bool),hh["r40c_fused_path"],a&hh["r40c_fused_localized"]
def correct(q,p):
    return (((abs(p-q["truth"])<=2)&q["active"].astype(bool)).sum(1)>=.75*q["active"].sum(1))&(q["y"]==1)
def encode(rows,R,q,a,p,hit):
    s,path,support,regime,F=evidence(rows,R,q)
    return {"s":s,"p":path,"support":support,"regime":regime,"F":F,"parent_alert":a,"parent_hit":hit,"parent_path":p,"y":q["y"],"weak":np.isin(q["amp"],t.S["weak_amplitudes"]),"correct":correct(q,path)}
def selftest():
    assert cut([4,4,3,2],.25)>4
    a=np.array([True,False,False]);v=np.array([False,True,True]);g=np.array([True,False,True])
    assert np.array_equal(a|(v&g),[True,False,True])
def main():
    st=time.perf_counter();selftest();L=h.readlocks()
    kz,kn=c.artifact(S["source"]["r40k"],32);raw=kz.read(kn["lock.json"])
    assert c.sha(raw)=="469f0172868ace2303e2a46b766a21c9c23e97f3cad733a4f2a8ee4bc1fe0cf3"
    parentlock=json.loads(raw);threshold=parentlock["thresholds"]["joint_min"]
    z,n,rec,_=t.parent.restore();m=json.loads(z.read(n["manifest.json"]))
    train=[];cal=[]
    records={v["name"]:v for v in rec["records"]}
    for group in m["selected"]:
        split=group["split"];name=group["target"]
        if split=="component_fit":continue
        rows,R,ch,hs=k.gather(group,"historical_"+name)
        pairs=[("cal" if split=="null_calibration" else "eval",rows,R)]
        if split=="evaluation":pairs.append(("shared",[v[0] for v in ch],np.stack([v[1] for v in ch])))
        else:assert sorted(hs)==next(v["raw_tile_hashes"] for v in parentlock["calibration"] if v["target"]==name)
        for kind,rr,ref in pairs:
            q=t.parent.load(z,n,records[kind+"_"+name+".npz"]);a,p,hit=parent(q,L)
            d=encode(rr,ref,q,a,p,hit)
            (cal if kind=="cal" else train).append((name,kind,d))
            emit("R40L_HISTORICAL",{"target":name,"kind":kind,"rows":len(a),"registered_exact":True})
    assert len(train)==8 and len(cal)==2
    X=np.concatenate([d["F"][~d["parent_alert"]] for _,_,d in train]);Y=np.concatenate([d["correct"][~d["parent_alert"]] for _,_,d in train])
    assert Y.any() and (~Y).any()
    sc=StandardScaler().fit(X);mdl=LogisticRegression(C=.1,class_weight="balanced",max_iter=1000,random_state=401).fit(sc.transform(X),Y);assert mdl.n_iter_.max()<1000
    refcut=float(np.median(np.concatenate([d["regime"][d["y"]==0] for _,_,d in train])))
    def score(d):return mdl.decision_function(sc.transform(d["F"]))
    def cat(key):return np.concatenate([d[key][(d["y"]==0)&~d["parent_alert"]] for _,_,d in cal])
    cs=cat("support");cr=cat("regime");cl=np.concatenate([score(d)[(d["y"]==0)&~d["parent_alert"]] for _,_,d in cal])
    targets=np.concatenate([np.repeat(name,int(((d["y"]==0)&~d["parent_alert"]).sum())) for name,_,d in cal])
    nc=sum(int((d["y"]==0).sum()) for _,_,d in cal);npalert=sum(int(((d["y"]==0)&d["parent_alert"]).sum()) for _,_,d in cal)
    assert nc==256 and npalert==10 and len(cs)==246 and (cat("s")<threshold).all()
    thresholds={"topa05":cut(cs,.05),"topa10":cut(cs,.10),"learned05":cut(cl,.05),"learned10":cut(cl,.10)}
    regime=cr>refcut;stats={str(i):{"controls":int((regime==i).sum()),"targets":len(set(targets[regime==i]))} for i in (0,1)}
    fallback=any(v["controls"]<128 or v["targets"]<3 for v in stats.values())
    kt=[thresholds["topa10"]]*2 if fallback else [cut(cs[regime==i],.10) for i in (0,1)]
    lock={"protocol_sha256":c.sha(P.read_bytes()),"parent_lock_sha256":c.sha(raw),"parent_joint_min_threshold":threshold,"thresholds":thresholds,"kusto_thresholds":kt,"kusto_fallback":fallback,"strata":stats,"reference_cut":refcut,"calibration_controls":nc,"parent_control_alerts":npalert,"train_rows":len(Y),"train_positive":int(Y.sum()),"mean":sc.mean_.tolist(),"scale":sc.scale_.tolist(),"coef":mdl.coef_.tolist(),"intercept":mdl.intercept_.tolist()}
    save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes());emit("R40L_LOCK",{k:v for k,v in lock.items() if k not in ("mean","scale","coef","intercept")})
    def metrics(d):
        a=d["parent_alert"];hit=d["parent_hit"];y=d["y"];weak=d["weak"];res=(~a)&(d["s"]>=threshold);sl=score(d);sup=d["support"];reg=d["regime"]>refcut
        guards={"r40g":np.zeros(len(a),bool),"joint_min":np.ones(len(a),bool),"topa05":sup>=thresholds["topa05"],"topa10":sup>=thresholds["topa10"],"kusto_topa10":sup>=np.asarray(kt)[reg.astype(int)],"learned05":sl>=thresholds["learned05"],"learned10":sl>=thresholds["learned10"]}
        baseline=hit|(res&d["correct"]);out={};arr={"joint_rescue":res,"learned_score":sl,**d}
        for key,guard in guards.items():
            rescue=res&guard;aa=a|rescue;hh=hit|(rescue&d["correct"])
            out[key]={"localized":int(hh.sum()),"weak":int((hh&weak).sum()),"control_alerts":int((aa&(y==0)).sum()),"lost_parent":int((hit&~hh).sum()),"lost_joint":int((baseline&~hh).sum()),"lost_joint_weak":int((baseline&~hh&weak).sum()),"additional_control_alerts":int((rescue&(y==0)).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum())}
            assert out[key]["lost_parent"]==0
            arr[key+"_alert"]=aa;arr[key+"_hit"]=hh
        return out,arr
    for name,kind,d in cal:
        m,_=metrics(d);assert all(v["additional_control_alerts"]==0 for v in m.values())
    # Retain compact features and scores for all historical groups, not spectra.
    for name,kind,d in train+cal:
        m,arr=metrics(d);np.savez_compressed(O/("historical_"+kind+"_"+name+".npz"),**arr)
    zz,nn=c.artifact(S["source"]["r40h"],192);rec=json.loads(zz.read(nn["receipt.json"]));m=json.loads(zz.read(nn["manifest.json"]))
    assert len(m["selected"])==5 and rec["manifest_sha256"]==c.sha(zz.read(nn["manifest.json"]))
    b.S["seed"]=h.S["seed"];b.ref.S=b.S;b.sh.eng.SPEC={**b.sh.eng.SPEC,"window_centers_mhz":h.S["windows_mhz"],"window_half_hz":h.S["half_hz"],"limits":{**b.sh.eng.SPEC["limits"],"spectrum_budget_per_file_mib":h.S["max_mib_per_scan"]}}
    results=[]
    for group in m["selected"]:
        name=group["target"];rows,R,ch,hs=k.gather(group,name);expected=next(v for v in rec["rows"] if v["target"]==name);assert sorted(hs)==expected["raw_tile_hashes"]
        for kind,rr,ref in (("eval",rows,R),("shared",[v[0] for v in ch],np.stack([v[1] for v in ch]))):
            key=kind+"_"+name+".npz";ar=next(v for v in rec["records"] if v["name"]==key);assert c.sha(zz.read(nn[key]))==ar["sha256"];q=f.load(zz,nn,key);assert c.sha(q["V"].tobytes())==ar["array_sha256"]
            d=encode(rr,ref,q,q["r40g_alert"].astype(bool),q["r40g_path"],q["r40g_hit"].astype(bool));out,arr=metrics(d)
            old=f.load(kz,kn,key);assert np.array_equal(d["s"],old["joint_min_score"]) and np.array_equal(d["p"],old["joint_min_path"]) and np.array_equal(arr["joint_rescue"],old["joint_min_rescue"]) and np.array_equal(arr["joint_min_hit"],old["joint_min_localized"])
            row={"target":name,"kind":kind,"registered_exact":True,"joint_min_exact":True,"metrics":out};results.append(row);save("partial.json",results);emit("R40L_GROUP",row);np.savez_compressed(O/key,**arr)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def tot(kind):
        rr=[r["metrics"] for r in results if r["kind"]==kind];assert len(rr)==5
        return {k:{j:sum(x[k][j] for x in rr) for j in rr[0][k]} for k in rr[0]}
    e=tot("eval");sh=tot("shared");assert [e["joint_min"][k] for k in ("localized","weak","control_alerts")]==[215,43,43] and sh["joint_min"]["localized"]==111 and sh["joint_min"]["weak"]==36
    gates={k:{"minimum":e[k]["localized"]>=214 and e[k]["weak"]>=41 and sh[k]["weak"]>=21 and e[k]["control_alerts"]<=35 and e[k]["lost_parent"]==sh[k]["lost_parent"]==0,"joint_preservation":sh[k]["weak"]>=36 and sh[k]["localized"]>=111} for k in e if k not in ("r40g","joint_min")}
    d={"status":"DEVELOPMENT_REPLAY_ONLY","eval":e,"shared":sh,"gates":gates,"rows":results,"lock_sha256":lh,"seconds":time.perf_counter()-st,"production_adoption":False}
    save("result.json",d);emit("R40L_SUMMARY",{k:v for k,v in d.items() if k!="rows"})
if __name__=="__main__":main()
