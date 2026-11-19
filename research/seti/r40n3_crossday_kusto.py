#!/usr/bin/env python3
"""R40N3: authenticated eight-day CROSSDAY bank and conditional KUSTO replay."""
import sys
sys.argv[1:]=["controls"]
import json,pathlib,time,io
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier
import r40n1_crossday as x
n=x.n;c=x.c;b=x.b;h=x.h;l=x.l;f=x.f
P=pathlib.Path("research/seti/r40n3_crossday_kusto_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40n3");O.mkdir(parents=True,exist_ok=True);n.O=O;n.k.O=O;n.l.O=O
makebank=x.makebank;cross=x.cross
def save(name,d):(O/name).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def archive(k):return c.artifact(S["sources"][k],64)
def main():
    start=time.perf_counter();x.selftest()
    az,an=archive("n1");bz,bn=archive("n2");rz,rn=archive("replay")
    a=n.read(az,an,"controls.json");am=n.read(az,an,"manifest.json")
    bb=n.read(bz,bn,"result.json");bm=n.read(bz,bn,"manifest.json")
    assert c.sha(az.read(an["manifest.json"]))==a["manifest_sha256"] and c.sha(bz.read(bn["manifest.json"]))==bb["manifest_sha256"]
    assert a["groups"]==4 and a["controls"]==512 and bb["total_groups"]==8 and bb["total_controls"]==1024 and set(bb["instruments"])=={"Parkes"}
    assert c.sha(rz.read(rn["lock.json"]))==S["sources"]["replay"]["lock_sha256"];rl=n.read(rz,rn,"lock.json")
    lz,ln=n.archive("r40l");raw=lz.read(ln["lock.json"]);assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];old=json.loads(raw)
    hm,exam,fr,sel=n.manifests();kz,kn=n.archive("r40k");kr=n.read(kz,kn,"result.json");kl=n.read(kz,kn,"lock.json")
    prior_days={int(g["day"]) for g in hm["selected"]+exam["selected"]}
    groups=am["selected"]+bm["selected"];assert len(groups)==8 and len({g["target"] for g in groups})==8 and len({int(g["day"]) for g in groups})==8 and not prior_days&{int(g["day"]) for g in groups}
    seen=set(fr["excluded_raw_hashes"]);seen.update(v for r in kr["rows"] for v in r["raw_tile_hashes"]);seen.update(v for r in kl["calibration"] for v in r["raw_tile_hashes"])
    for name in rn:
        if name.endswith(".npz"):
            with np.load(io.BytesIO(rz.read(rn[name])),allow_pickle=False) as q:seen.update(q["source_tile_sha"].tolist())
    initial=set(seen);bank=[];controls=[];records=[]
    cfg=json.loads(pathlib.Path("research/seti/r40n2_instrument_controls_protocol.json").read_text())["selection"]
    n.set_config(cfg["seed"]+"|Parkes",cfg["parkes_windows_mhz"],cfg["half_hz"],cfg["max_mib_per_scan"])
    for i,g in enumerate(groups):
        name=g["target"];day=int(g["day"])
        rec=next(v for v in (a["records"] if i<4 else bb["records"]) if v["target"]==name);hs=set(rec["raw_tile_hashes"])
        assert len(hs)==16 and not hs&seen;seen.update(hs)
        if i<4:
            d=f.load(az,an,"control_"+name+".npz");raw=az.read(an["bank_"+name+".npz"]);(O/("bank_"+name+".npz")).write_bytes(raw)
            with np.load(io.BytesIO(raw),allow_pickle=False) as q:
                grid=q["index"];assert grid.shape==(16,4)
                items=[(int(v[0]),int(v[1]),float(v[2]),float(v[3]),q["tile_"+str(j)].copy()) for j,v in enumerate(grid)]
        else:
            assert g["telescope"]=="Parkes";d=f.load(bz,bn,"control_Parkes_"+name+".npz")
            rr,R,ch,hs2,meta,cm,pre=n.gather(g,name);assert hs2==hs
            ix=np.array([v[1]==0 for v in rr]);assert ix.sum()==128
            rows=[v for v,q in zip(rr,ix) if q];mm=[v for v,q in zip(meta,ix) if q]
            fresh=n.encode(rows,R[ix],mm,pre)
            for key in ("s","p","support","regime","F","raw_features","repeat_features","frequency_mhz","source_tile_sha","y"):assert np.array_equal(fresh[key],d[key]),"N2 parity "+key
            items=makebank(pre)
            np.savez_compressed(O/("bank_"+name+".npz"),index=np.asarray([[slot,wi,f0,df] for slot,wi,f0,df,z in items]),**{"tile_"+str(j):v[4] for j,v in enumerate(items)})
        assert len(items)==16 and {(v[0],v[1]) for v in items}=={(slot,wi) for slot in range(4) for wi in range(4)}
        assert all(np.isfinite(v[4]).all() and v[3]!=0 for v in items)
        assert len(d["y"])==128 and not d["y"].any() and set(d["source_tile_sha"].tolist())<=hs
        assert int(d["parent_alert"].sum())==rec["parent_alerts"] and int(d["joint_rescue"].sum())==rec.get("extra_joint_alerts",rec.get("joint_extra_alerts"))
        assert np.array_equal(d["joint_rescue"],(~d["parent_alert"])&(d["s"]>=old["parent_joint_min_threshold"]))
        bank.append((day,items));controls.append((name,day,d));row={"target":name,"day":day,"controls":128,"raw_tile_hashes":sorted(hs),"bank_source":"N1_AUTHENTICATED" if i<4 else "N2_RAW_REGENERATED_EXACT"};records.append(row);emit("R40N3_BANK",row);save("partial.json",records)
    assert len(seen-initial)==128 and len(bank)==8
    bank.sort(key=lambda v:v[0]);controls.sort(key=lambda v:v[1])
    for name,day,d in controls:
        d["crossday_features"]=cross(d["frequency_mhz"],day,bank);assert (d["crossday_features"][:,0]<=7).all()
        np.savez_compressed(O/("control_"+name+".npz"),**d)
    rg=np.concatenate([d["regime"]>old["reference_cut"] for _,_,d in controls]);pa=np.concatenate([d["parent_alert"] for _,_,d in controls]);ta=np.repeat([day for _,day,_ in controls],128)
    stats={str(i):{"controls":int(((rg==i)&~pa).sum()),"days":len(set(ta[(rg==i)&~pa].tolist()))} for i in (0,1)}
    assert len(rg)==1024 and int(pa.sum())==50 and stats==bb["instruments"]["Parkes"]["strata"]=={"0":{"controls":487,"days":8},"1":{"controls":487,"days":8}}
    qualified=True
    def raw_score(d):
        v=rl["models"]["raw_logistic"];X=d["raw_features"];return ((X-np.asarray(v["mean"]))/np.asarray(v["scale"]))@np.asarray(v["coef"])[0]+v["intercept"][0]
    fixed={}
    for key,z in (("joint",np.concatenate([d["s"] for _,_,d in controls])),("raw_logistic",np.concatenate([raw_score(d) for _,_,d in controls]))):
        pooled=l.cut(z[~pa],.025);conditional=[l.cut(z[(rg==i)&~pa],.025) for i in (0,1)];v=bb["instruments"]["Parkes"]["thresholds_development_only"][key]
        assert np.allclose([pooled]+conditional,[v["pooled"]]+v["conditional"],rtol=0,atol=1e-10)
        fixed[key]={"pooled":v["pooled"],"conditional":v["conditional"]}
    save("fixed_n2_thresholds.json",fixed)
    manifest={"n1_manifest_sha256":a["manifest_sha256"],"n2_manifest_sha256":bb["manifest_sha256"],"records":records};save("manifest.json",manifest);mh=c.sha((O/"manifest.json").read_bytes())
    receipt={"status":"PASS_KUSTO_COUNTS_AND_BANK_PARITY","manifest_sha256":mh,"controls":1024,"groups":8,"strata":stats,"parent_alerts":50,"fallback":False,"fixed_threshold_parity":True,"calibration_bank_days":7,"train_test_bank_days":8,"formal_fap":False,"physical_null_verified":False,"production_adoption":False}
    save("controls.json",receipt);emit("R40N3_CONTROL_SUMMARY",receipt)
    daymap={g["target"]:int(g["day"]) for g in hm["selected"]+exam["selected"]}
    def load(name):
        d=f.load(rz,rn,name);target=name.removesuffix(".npz").split("_")[-1];assert daymap[target] not in {day for day,_ in bank};d["crossday_features"]=cross(d["frequency_mhz"],daymap[target],bank);return d
    train=[load(name) for name in sorted(rn) if name.startswith(("historical_eval_","historical_shared_")) and name.endswith(".npz")];assert len(train)==8
    models={};thresholds={};serial={}
    def use(d):return ~d["parent_alert"]
    def null(d):return use(d)&(d["y"]==0)
    def pos(d):return d["joint_rescue"]&d["correct"]
    y=np.concatenate([d["correct"][use(d)] for d in train]);assert y.any() and (~y).any()
    for feature in ("raw","crossday"):
        def fv(d,feature=feature):return d["raw_features"] if feature=="raw" else np.column_stack([d["raw_features"],d["crossday_features"][:,1:]])
        X=np.concatenate([fv(d)[use(d)] for d in train]);sc=StandardScaler().fit(X)
        for kind in ("logistic","tree"):
            model=LogisticRegression(C=.1,class_weight="balanced",max_iter=1000,random_state=401).fit(sc.transform(X),y) if kind=="logistic" else ExtraTreesClassifier(n_estimators=128,max_depth=4,min_samples_leaf=8,class_weight="balanced",random_state=401,n_jobs=1).fit(X,y)
            def score(d,model=model,sc=sc,fv=fv,kind=kind):return model.decision_function(sc.transform(fv(d))) if kind=="logistic" else model.predict_proba(fv(d))[:,1]
            name=feature+"_"+kind;models[name]=score
            gp=np.concatenate([score(d)[pos(d)] for d in train]);gc=np.concatenate([score(d)[null(d)] for _,_,d in controls]);reg=rg[~pa];assert len(gp)>0 and len(gc)==len(reg)
            pooled=l.cut(gc,.025);conditional=[l.cut(gc[reg==i],.025) for i in (0,1)] if qualified else [pooled,pooled]
            thresholds[name]={"preserve":float(np.nextafter(gp.min(),-np.inf)),"pooled":pooled,"kusto":conditional}
            serial[name]={"mean":sc.mean_.tolist(),"scale":sc.scale_.tolist(),"coef":model.coef_.tolist(),"intercept":model.intercept_.tolist()} if kind=="logistic" else {"trees":[{"children_left":e.tree_.children_left.tolist(),"children_right":e.tree_.children_right.tolist(),"feature":e.tree_.feature.tolist(),"threshold":e.tree_.threshold.tolist(),"value":e.tree_.value.tolist()} for e in model.estimators_]}
    lock={"protocol_sha256":c.sha(P.read_bytes()),"manifest_sha256":mh,"source_lock_sha256":S["sources"]["replay"]["lock_sha256"],"fixed_n2_thresholds":fixed,"thresholds":thresholds,"models":serial,"kusto_qualified":qualified,"bank_days":[day for day,_ in bank],"bank_sha256":{p.name:c.sha(p.read_bytes()) for p in sorted(O.glob("bank_*.npz"))},"train_rows":len(y)}
    save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes());emit("R40N3_LOCK",{k:v for k,v in lock.items() if k!="models"})
    def measure(d):
        a=d["parent_alert"];base=d["parent_hit"];res=d["joint_rescue"];y=d["y"];weak=d["weak"];correct=d["correct"];oldhit=base|(res&correct);reg=(d["regime"]>old["reference_cut"]).astype(int)
        guards={"joint_min":np.ones(len(y),bool),"r40g":np.zeros(len(y),bool)}
        if "learned05_alert" in d:guards["learned05"]=d["learned05_alert"]&~a
        for key,v in (("joint",d["s"]),("raw_logistic",raw_score(d))):
            th=fixed[key];guards["n2_"+key+"_pooled"]=v>=th["pooled"];guards["n2_"+key+"_kusto"]=v>=np.asarray(th["conditional"])[reg]
            if "learned05" in guards:guards["n2_"+key+"_kusto_learned"]=guards["n2_"+key+"_kusto"]&guards["learned05"]
        for name,score in models.items():
            v=score(d);th=thresholds[name]
            for cut in ("preserve","pooled","kusto"):guards[name+"_"+cut]=v>=(np.asarray(th[cut])[reg] if cut=="kusto" else th[cut])
        out={}
        for name,g in guards.items():
            rescue=res&g;aa=a|rescue;hh=base|(rescue&correct)
            out[name]={"localized":int(hh.sum()),"weak":int((hh&weak).sum()),"control_alerts":int((aa&(y==0)).sum()),"lost_parent":int((base&~hh).sum()),"lost_joint":int((oldhit&~hh).sum()),"lost_joint_weak":int((oldhit&~hh&weak).sum()),"false_paths_on_positive":int((aa&(y==1)&~hh).sum())}
            d[name+"_n3_alert"]=aa;d[name+"_n3_hit"]=hh
        return out
    calibration=[]
    for name,day,d in controls:
        row={"target":name,"metrics":measure(d),"covered":int((d["crossday_features"][:,0]>0).sum()),"min_covered_days":int(d["crossday_features"][:,0].min()),"max_covered_days":int(d["crossday_features"][:,0].max())};calibration.append(row);emit("R40N3_CALIBRATION",row)
    save("calibration.json",calibration)
    results=[]
    for name in sorted(rn):
        if not name.startswith(("eval_","shared_")) or not name.endswith(".npz"):continue
        d=load(name);row={"record":name,"metrics":measure(d),"crossday_covered":int((d["crossday_features"][:,0]>0).sum()),"min_covered_days":int(d["crossday_features"][:,0].min()),"max_covered_days":int(d["crossday_features"][:,0].max())};results.append(row);save("replay_partial.json",results);emit("R40N3_GROUP",row);np.savez_compressed(O/name,**d)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def total(kind):
        rr=[r["metrics"] for r in results if r["record"].startswith(kind+"_")];assert len(rr)==5
        return {a:{k:sum(v[a][k] for v in rr) for k in rr[0][a]} for a in rr[0]}
    e=total("eval");sh=total("shared");assert [e["joint_min"][k] for k in ("localized","weak","control_alerts")]==[215,43,43] and [sh["joint_min"][k] for k in ("localized","weak")]==[111,36]
    gates={a:qualified and e[a]["localized"]>=215 and e[a]["weak"]>=43 and e[a]["control_alerts"]<=35 and sh[a]["localized"]>=111 and sh[a]["weak"]>=36 and e[a]["lost_joint"]==sh[a]["lost_joint"]==0 for a in e}
    result={"status":"DEVELOPMENT_REPLAY_ONLY","controls":receipt,"calibration":calibration,"eval":e,"shared":sh,"gates":gates,"lock_sha256":lh,"rows":results,"seconds":time.perf_counter()-start,"production_adoption":False};save("result.json",result);emit("R40N3_SUMMARY",{k:v for k,v in result.items() if k not in ("rows","controls","calibration")})
if __name__=="__main__":main()
