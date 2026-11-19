#!/usr/bin/env python3
"""R40N5: matched-frequency recurrence bank and nondegenerate protected replay."""
import sys
sys.argv[1:]=["controls"]
import json,pathlib,time,io
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import roc_auc_score
import r40n1_crossday as x
n=x.n;c=x.c;b=x.b;l=x.l;f=x.f
P=pathlib.Path("research/seti/r40n5_matched_frequency_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40n5");O.mkdir(parents=True,exist_ok=True)
def save(name,d):(O/name).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def archive(k):return c.artifact(S["sources"][k],64)
# Each day yields matched receiver evidence, not one independent sample per channel.
# Native channel offsets are justified by fixed Parkes resolution, no frequency remapping.
RADII=(2,8,32)
def neighborhood(freq,bank,omit=None):
    out=np.zeros((len(freq),1+7*len(RADII)),float)
    eligible=[(day,items) for day,items in bank if day!=omit]
    for i,fq in enumerate(freq):
        vals=[]
        for day,items in eligible:
            aligned={}
            for slot,wi,f0,df,z in items:
                ix=int(round((fq-f0)*1e6/df))-5
                if 32<=ix<len(z)-32:
                    aligned[slot,wi]=(z[ix-32:ix+33].astype(float)[::(1 if df>0 else -1)])
            pairs=[(aligned[slot,wi],aligned[slot+1,wi]) for slot in (0,2) for wi in range(4) if (slot,wi) in aligned and (slot+1,wi) in aligned]
            if not pairs:continue
            vals.append(pairs)
        if not vals:continue
        out[i,0]=len(vals)/len(eligible)
        for j,r in enumerate(RADII):
            per=[]
            for daypairs in vals:
                # Highest simultaneous local contrast across ON/OFF on this day.
                # Compare signed physical offsets (all grids use same |df|).
                pieces=[]
                for a,d in daypairs:
                    aa=a[32-r:33+r];dd=d[32-r:33+r];joint=np.minimum(aa,dd);k=int(np.argmax(joint))
                    pieces.append((float(joint[k]),float(dd.max()),float((joint>=.75).mean()),float((dd>=.75).mean()),k-r))
                per.append(max(pieces,key=lambda v:v[0]))
            v=np.asarray(per);offset=v[:,4];k=1+7*j
            out[i,k:k+7]=[np.median(v[:,0]),v[:,0].max(),np.mean(v[:,1]),np.mean(v[:,2]),np.mean(v[:,3]),np.mean(v[:,0]>=.75),np.mean(abs(offset-np.median(offset))) / r]
    assert np.isfinite(out).all()
    return out
def placebo(freq,bank,omit=None):
    # Fixed sign derived from physical frequency quantization, no labels.
    df=abs(bank[0][1][0][3]);sign=np.where(np.floor(freq*1e6/df).astype(np.int64)%2==0,1,-1)
    return neighborhood(freq+sign*256*df/1e6,bank,omit)
def seven_day(feat,freq,bank,day):
    if day in {d for d,_ in bank}:return feat(freq,bank,day)
    return np.mean([feat(freq,bank,d) for d,_ in bank],axis=0)
def central(freq,bank,omit=None):
    return x.cross(freq,omit,bank)[:,1:]
def selftest():
    z=np.zeros(100);z[54]=2.;items=[(i,0,3000.,2.8,z) for i in range(4)];bank=[(1,items),(2,items)]
    fq=np.array([3000.+(50+5)*2.8/1e6]);v=neighborhood(fq,bank,1);assert v[0,0]==1 and v[0,1]==0 and v[0,8]==2
    assert not neighborhood(np.array([3100.]),bank).any()
    assert not neighborhood(fq,[(1,items)],1).any()
    neg=[(i,0,3000.,-2.8,z) for i in range(4)];w=neighborhood(np.array([3000.-55*2.8/1e6]),[(1,neg),(2,neg)],1);assert np.array_equal(v,w)
    emit("R40N5_SELFTEST",{"offset_recurrence":"PASS","signed_native_grid":"PASS","day_exclusion":"PASS","offband_abstention":"PASS"})
def model_fit(X,y,kind):
    if kind=="logistic":
        sc=StandardScaler().fit(X);model=LogisticRegression(C=.1,class_weight="balanced",max_iter=1000,random_state=401).fit(sc.transform(X),y)
        assert model.n_iter_.max()<1000
        def score(v):return model.decision_function(sc.transform(v))
        serial={"mean":sc.mean_.tolist(),"scale":sc.scale_.tolist(),"coef":model.coef_.tolist(),"intercept":model.intercept_.tolist()}
    else:
        model=ExtraTreesClassifier(n_estimators=128,max_depth=4,min_samples_leaf=8,class_weight="balanced",random_state=401,n_jobs=1).fit(X,y)
        def score(v):return model.predict_proba(v)[:,1]
        serial={"trees":[{"children_left":e.tree_.children_left.tolist(),"children_right":e.tree_.children_right.tolist(),"feature":e.tree_.feature.tolist(),"threshold":e.tree_.threshold.tolist(),"value":e.tree_.value.tolist()} for e in model.estimators_]}
    return score,serial
def col_audit(X,y):
    out=[]
    for i,v in enumerate(X.T):
        distinct=len(np.unique(v));out.append({"column":i,"std":float(v.std()),"min":float(v.min()),"max":float(v.max()),"unique":distinct,"auc":float(roc_auc_score(y,v)) if distinct>1 else .5})
    return out
def matched_placebo(freq,bank,omit=None):
    df=abs(bank[0][1][0][3]);sign=np.where(np.floor(freq*1e6/df).astype(np.int64)%2==0,1,-1)
    v=neighborhood(freq+sign*256*df/1e6,bank,omit)
    missing=v[:,0]<1
    if missing.any():v[missing]=neighborhood(freq[missing]-sign[missing]*256*df/1e6,bank,omit)
    return v
def coverage(fq,centers):
    q=np.asarray(fq);dist=abs(q[:,None]-np.asarray(centers)[None,:]);idx=dist.argmin(1)
    assert np.isfinite(q).all() and np.all(dist.min(1)<=.0501),"Candidate frequency outside specified windows"
    return {"minimum_mhz":float(q.min()),"maximum_mhz":float(q.max()),"counts":{str(v):int((idx==j).sum()) for j,v in enumerate(centers)}}
def main():
    start=time.perf_counter();selftest();z,nn=archive("n3");rz,rn=archive("replay")
    assert c.sha(z.read(nn["lock.json"]))==S["sources"]["n3"]["lock_sha256"] and c.sha(rz.read(rn["lock.json"]))==S["sources"]["replay"]["lock_sha256"]
    lock=n.read(z,nn,"lock.json");manifest=n.read(z,nn,"manifest.json");receipt=n.read(z,nn,"controls.json")
    assert c.sha(z.read(nn["manifest.json"]))==lock["manifest_sha256"] and receipt["controls"]==1024 and receipt["groups"]==8
    bank=[];days={}
    for rec in manifest["records"]:
        name=rec["target"];day=int(rec["day"]);fn="bank_"+name+".npz";raw=z.read(nn[fn]);assert c.sha(raw)==lock["bank_sha256"][fn]
        with np.load(io.BytesIO(raw),allow_pickle=False) as q:
            a=q["index"];assert a.shape==(16,4);items=[(int(v[0]),int(v[1]),float(v[2]),float(v[3]),q["tile_"+str(j)].copy()) for j,v in enumerate(a)]
        assert len({round(abs(v[3]),10) for v in items})==1
        bank.append((day,items));days[name]=day
    bank.sort(key=lambda v:v[0]);assert len(bank)==8 and len(set(days.values()))==8
    assert len({round(abs(it[3]),10) for _,items in bank for it in items})==1
    hm,exam,fr,sel=n.manifests();daymap={g["target"]:int(g["day"]) for g in hm["selected"]+exam["selected"]}
    lz,ln=n.archive("r40l");raw=lz.read(ln["lock.json"]);assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];old=json.loads(raw)
    # Authenticated selection manifests, no new target selection.
    az,an=archive("n1");bz,bn=archive("n2")
    ma=n.read(az,an,"manifest.json");mb=n.read(bz,bn,"manifest.json")
    aa=n.read(az,an,"controls.json");bb=n.read(bz,bn,"result.json")
    assert c.sha(az.read(an["manifest.json"]))==aa["manifest_sha256"]==manifest["n1_manifest_sha256"]
    assert c.sha(bz.read(bn["manifest.json"]))==bb["manifest_sha256"]==manifest["n2_manifest_sha256"]
    groups0=ma["selected"]+mb["selected"];assert {g["target"]:int(g["day"]) for g in groups0}==days
    assert not set(days.values())&set(daymap.values())
    W=S["windows"];historical=W["historical_mhz"];evaluation=W["evaluation_mhz"]
    seen=set(fr["excluded_raw_hashes"])
    for v in manifest["records"]:seen.update(v["raw_tile_hashes"])
    kz,kn=n.archive("r40k");kr=n.read(kz,kn,"result.json");kl=n.read(kz,kn,"lock.json")
    for v in kr["rows"]+kl["calibration"]:seen.update(v["raw_tile_hashes"])
    for name in rn:
        if name.endswith(".npz"):
            with np.load(io.BytesIO(rz.read(rn[name])),allow_pickle=False) as q:seen.update(q["source_tile_sha"].tolist())
    newbank=[];newrecords=[];bank_hashes=dict(lock["bank_sha256"])
    n.set_config("R40N5-RAW-REFERENCE-ONLY",historical,W["historical_read_half_hz"],W["max_mib_per_scan"])
    for g in sorted(groups0,key=lambda q:int(q["day"])):
        name=g["target"];day=int(g["day"]);w=O/"raw"/name;w.mkdir(parents=True,exist_ok=True)
        n.ref.OUT=w;b.sh.eng.OUT=w
        assert [v["role"] for v in g["files"]]==["ON1","OFF1","ON2","OFF2"]
        pre=[n.ref.tiles(v) for v in g["files"]]
        assert len(pre)==4 and all(len(v)==4 for v in pre)
        assert all([a["center_mhz"] for a,V in v]==historical for v in pre)
        hs={a["slice_sha256"] for v in pre for a,V in v};assert len(hs)==16 and not hs&seen;seen.update(hs)
        items=x.makebank(pre);assert len(items)==16 and all(np.isfinite(v[4]).all() for v in items)
        assert len({round(abs(v[3]),10) for v in items})==1 and round(abs(items[0][3]),10)==round(abs(bank[0][1][0][3]),10)
        path=O/("historical_bank_"+name+".npz")
        np.savez_compressed(path,index=np.asarray([[slot,wi,f0,df] for slot,wi,f0,df,z0 in items]),**{"tile_"+str(j):v[4] for j,v in enumerate(items)})
        bank_hashes[path.name]=c.sha(path.read_bytes());newbank.append((day,items))
        rec={"target":name,"day":day,"raw_tile_hashes":sorted(hs),"windows_mhz":historical,"bank_sha256":bank_hashes[path.name]};newrecords.append(rec);save("bank_manifest.json",newrecords);emit("R40N5_BANK",rec)
        for p0 in w.glob("*.npy"):p0.unlink()
    assert len(newbank)==8 and len({v for rec in newrecords for v in rec["raw_tile_hashes"]})==128
    emit("R40N5_BANK_SUMMARY",{"days":8,"historical_raw_tiles":128,"retained_raw_tiles":128,"windows_mhz":historical+evaluation})
    def augment(d,day,band):
        fq=d["frequency_mhz"];centers=historical if band=="historical" else evaluation
        bk=newbank if band=="historical" else bank
        d["central"]=seven_day(central,fq,bk,day)
        d["neighborhood"]=seven_day(neighborhood,fq,bk,day);d["placebo"]=seven_day(matched_placebo,fq,bk,day)
        assert np.all(d["central"][:,0]==1) and np.all(d["neighborhood"][:,0]==1) and np.all(d["placebo"][:,0]==1),"Incomplete matched-frequency support; do not fit/impute"
        return d
    controls=[]
    for name,day in sorted(days.items()):
        d=f.load(z,nn,"control_"+name+".npz");assert np.array_equal(d["crossday_features"],x.cross(d["frequency_mhz"],day,bank))
        assert len(d["y"])==128 and not d["y"].any();coverage(d["frequency_mhz"],evaluation);controls.append((name,augment(d,day,"evaluation")))
    train=[]
    for name in sorted(rn):
        if not name.startswith(("historical_eval_","historical_shared_")) or not name.endswith(".npz"):continue
        target=name.removesuffix(".npz").split("_")[-1];d=f.load(rz,rn,name);assert daymap[target] not in set(days.values())
        d["central_original"]=x.cross(d["frequency_mhz"],daymap[target],bank);assert not d["central_original"][:,0].any();cv=coverage(d["frequency_mhz"],historical);train.append((target,augment(d,daymap[target],"historical")));emit("R40N5_TRAIN_COVERAGE",{"record":name,"frequency":cv,"old_covered":0,"new_covered":len(d["frequency_mhz"])})
    assert len(train)==8 and len(set(v for v,_ in train))==4
    def use(d):return ~d["parent_alert"]
    def null(d):return use(d)&(d["y"]==0)
    def pos(d):return d["joint_rescue"]&d["correct"]
    reg=np.concatenate([(d["regime"]>old["reference_cut"])[null(d)] for _,d in controls]);assert (reg==0).sum()==(reg==1).sum()==487
    y=np.concatenate([d["correct"][use(d)] for _,d in train]);groups=np.concatenate([[name]*int(use(d).sum()) for name,d in train]);assert len(y)==882
    audit={}
    for field in ("central_original","central","neighborhood","placebo"):
        X=np.concatenate([d[field][use(d)] for _,d in train]);audit[field]=col_audit(X,y)
    audit["previous_crossday_coefficients"]=lock["models"]["crossday_logistic"]["coef"][0][-7:]
    audit["previous_crossday_scales"]=lock["models"]["crossday_logistic"]["scale"][-7:]
    save("audit.json",audit);emit("R40N5_AUDIT",audit)
    variable={field:sum(v["std"]>1e-12 for v in audit[field][1:]) for field in ("central","neighborhood","placebo")}
    save("variation_gate.json",variable);emit("R40N5_VARIATION",variable)
    assert all(v>=1 for v in variable.values()),"Degenerate recurrence representation; HOLD before fit"
    models={};serial={};cuts={};auc={}
    def fv(d,feature):
        if feature=="raw":return d["raw_features"]
        return np.column_stack([d["raw_features"],d[feature]])
    for feature in ("raw","central","neighborhood","placebo"):
        X=np.concatenate([fv(d,feature)[use(d)] for _,d in train])
        for kind in ("logistic","tree"):
            name=feature+"_"+kind;scores=np.zeros(len(y));per={}
            for group in sorted(set(groups)):
                test=groups==group;sc,_=model_fit(X[~test],y[~test],kind);scores[test]=sc(X[test])
                per[group]=float(roc_auc_score(y[test],scores[test])) if y[test].any() and (~y[test]).any() else None
            auc[name]={"pooled_oof_auc":float(roc_auc_score(y,scores)),"by_group":per}
            sc,s=model_fit(X,y,kind)
            def score(d,sc=sc,feature=feature):return sc(fv(d,feature))
            models[name]=score;serial[name]=s
            gp=np.concatenate([score(d)[pos(d)] for _,d in train]);gc=np.concatenate([score(d)[null(d)] for _,d in controls]);assert len(gc)==974 and len(gp)>0
            # Prespecified calibration ranks, not a search on the seen evaluation set.
            cuts[name]={"preserve":float(np.nextafter(gp.min(),-np.inf)),"pooled025":l.cut(gc,.025),"kusto025":[l.cut(gc[reg==i],.025) for i in (0,1)],"pooled005":l.cut(gc,.005),"kusto005":[l.cut(gc[reg==i],.005) for i in (0,1)]}
    save("oof.json",auc);emit("R40N5_OOF",auc)
    lk={"protocol_sha256":c.sha(P.read_bytes()),"source_lock_sha256":S["sources"]["n3"]["lock_sha256"],"bank_sha256":bank_hashes,"new_bank_manifest_sha256":c.sha((O/"bank_manifest.json").read_bytes()),"bank_days":sorted(days.values()),"models":serial,"thresholds":cuts,"historical_oof_auc":auc,"reference_cut":old["reference_cut"],"train_rows":len(y),"calibration_controls":1024,"calibration_parent_rejected":974,"calibration_strata":[487,487],"bank_support":"seven days, own day excluded, external average of 8 seven-day subsets"}
    save("lock.json",lk);lh=c.sha((O/"lock.json").read_bytes());emit("R40N5_LOCK",{k:v for k,v in lk.items() if k!="models"})
    def measure(d,baseline=False):
        a=d["parent_alert"];base=d["parent_hit"];res=d["joint_rescue"];y=d["y"];weak=d["weak"];correct=d["correct"];oldhit=base|(res&correct);reg=(d["regime"]>old["reference_cut"]).astype(int)
        guards={"joint_min":np.ones(len(y),bool),"r40g":np.zeros(len(y),bool)}
        if baseline:
            guards["learned05"]=d["learned05_alert"]&~a
            guards["n3_tree_kusto"]=d["crossday_tree_kusto_n3_alert"]&~a
        for name,score in models.items():
            v=score(d)
            for cut,th in cuts[name].items():guards[name+"_"+cut]=v>=(np.asarray(th)[reg] if isinstance(th,list) else th)
        out={}
        for name,g in guards.items():
            rescue=res&g;aa=a|rescue;hh=base|(rescue&correct)
            out[name]={"localized":int(hh.sum()),"weak":int((hh&weak).sum()),"control_alerts":int((aa&(y==0)).sum()),"lost_parent":int((base&~hh).sum()),"lost_joint":int((oldhit&~hh).sum()),"lost_joint_weak":int((oldhit&~hh&weak).sum()),"false_paths_on_positive":int((aa&(y==1)&~hh).sum())}
            d[name+"_n5_alert"]=aa;d[name+"_n5_hit"]=hh
        return out
    cal=[]
    for name,d in controls:
        cal.append({"target":name,"metrics":measure(d),"neighborhood_coverage":float(d["neighborhood"][:,0].mean()),"placebo_coverage":float(d["placebo"][:,0].mean())})
    save("calibration.json",cal)
    # Viewed development panel; never use it to fit or choose cuts.
    results=[]
    for name in sorted(nn):
        if not name.startswith(("eval_","shared_")) or not name.endswith(".npz"):continue
        target=name.removesuffix(".npz").split("_")[-1];d=f.load(z,nn,name);day=daymap[target];assert day not in set(days.values())
        assert np.array_equal(d["crossday_features"],x.cross(d["frequency_mhz"],day,bank))
        cv=coverage(d["frequency_mhz"],evaluation);d=augment(d,day,"evaluation");row={"record":name,"metrics":measure(d,True),"neighborhood_coverage":float(d["neighborhood"][:,0].mean()),"placebo_coverage":float(d["placebo"][:,0].mean())};results.append(row);save("partial.json",results);np.savez_compressed(O/name,**d)
        assert c.sha((O/"lock.json").read_bytes())==lh
        emit("R40N5_GROUP",{"record":name,"neighborhood_coverage":row["neighborhood_coverage"],"placebo_coverage":row["placebo_coverage"]})
    def total(kind):
        rr=[v["metrics"] for v in results if v["record"].startswith(kind+"_")];assert len(rr)==5
        return {a:{j:sum(v[a][j] for v in rr) for j in rr[0][a]} for a in rr[0]}
    e=total("eval");sh=total("shared")
    assert [e["joint_min"][k] for k in ("localized","weak","control_alerts")]==[215,43,43] and [sh["joint_min"][k] for k in ("localized","weak")]==[111,36]
    assert [e["n3_tree_kusto"][k] for k in ("localized","weak","control_alerts")]==[215,43,40] and [sh["n3_tree_kusto"][k] for k in ("localized","weak")]==[111,36]
    gates={a:e[a]["localized"]>=215 and e[a]["weak"]>=43 and e[a]["control_alerts"]<=35 and sh[a]["localized"]>=111 and sh[a]["weak"]>=36 and e[a]["lost_joint"]==sh[a]["lost_joint"]==0 for a in e}
    result={"status":"DEVELOPMENT_REPLAY_ONLY","nonconstant_recurrence_columns":variable,"bank_windows":historical+evaluation,"eval":e,"shared":sh,"gates":gates,"oof":auc,"rows":results,"lock_sha256":lh,"seconds":time.perf_counter()-start,"production_adoption":False};save("result.json",result);emit("R40N5_SUMMARY",{k:v for k,v in result.items() if k!="rows"})
if __name__=="__main__":main()
