#!/usr/bin/env python3
"""R40N1 frozen expanded metadata census, cross-day witness and KUSTO research replay."""
import sys
sys.argv[1:]=["controls"]
import json,pathlib,time,hashlib
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier
import r40n_raw_rivals as n
c=n.c;b=n.b;h=n.h;l=n.l;f=n.f;t=n.t
P=pathlib.Path("research/seti/r40n1_crossday_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r40n1");O.mkdir(parents=True,exist_ok=True);n.O=O;n.k.O=O;n.l.O=O
def save(name,d):(O/name).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def coord(g):
    hd=g["files"][0]["header"];return (float(hd["src_raj"]),float(hd["src_dej"]))
def select():
    hm,m,fr,sel=n.manifests();rz,rn=n.archive("r37_manifest");rm=n.read(rz,rn,"manifest.json")
    prior=hm["selected"]+m["selected"]+rm["selected"]["Parkes"]+rm["selected"].get("GBT",[])
    cde=json.loads(pathlib.Path("research/seti/r40cde_protocol.json").read_text())
    names=set(sel["excluded_names"])|{g["target"] for g in prior}|set(cde["targets"])|set(b.sh.hold.S["exclude"])
    days=set(sel["excluded_days"])|{int(g["day"]) for g in prior}|set(cde["days"])|set(b.sh.C["exclude_days"])
    coords={coord(g) for g in prior if g.get("files") and "header" in g["files"][0]}
    urls={v["url"] for g in prior for v in g.get("files",[]) if "url" in v}
    filenames={v["name"] for g in prior for v in g.get("files",[]) if "name" in v}
    zz,nn=c.artifact(S["sources"]["partial"],32);old=n.read(zz,nn,"manifest.json")
    assert [g["target"] for g in old["selected"]]==S["selection"]["retain"]
    eng=b.sh.eng;eng.OUT=O;cfg=n.S["new_controls"];n.set_config(cfg["seed"],cfg["windows_mhz"],cfg["half_hz"],cfg["max_mib_per_scan"])
    def validate(target,fs,q):
        day=int(fs[0]["day"])
        assert target not in names and day not in days
        assert [v["position"] for v in fs]==["S","R","S","R"] and len({v["day"] for v in fs})==1
        assert all(0<v["seconds"]-u["seconds"]<=600 for u,v in zip(fs,fs[1:]))
        assert len({v["url"] for v in fs})==len({v["name"] for v in fs})==4 and not any(v["url"] in urls or v["name"] in filenames for v in fs)
        hd=[eng.hdr(v) for v in fs]
        assert all(eng.valid(x,v,target) and x["shape"][0]==16 for x,v in zip(hd,fs))
        assert len({(x["shape"][0],round(x["ts"],7),round(x["df"],12)) for x in hd})==1
        co=[(float(x["src_raj"]),float(x["src_dej"])) for x in hd]
        assert co[0]==co[2] and co[1]==co[3] and co[0]!=co[1] and all(x!=(0.,0.) and np.isfinite(x).all() for x in co) and co[0] not in coords
        return {"target":target,"day":day,"split":"null_calibration","query":q,"files":[{**v,"header":x,"role":role} for v,x,role in zip(fs,hd,("ON1","OFF1","ON2","OFF2"))]}
    raw,url=b.sh.hold.base.fetch(b.sh.hold.r11.S["selection"]["index"],limit=32*b.sh.hold.base.MIB)
    (O/"target-index.json").write_bytes(raw);pool=set(b.sh.hold.r11.names(json.loads(raw)))-names
    order=sorted(pool,key=lambda x:(x in set(old["candidate_order"]),hashlib.sha256((S["selection"]["seed"]+"|"+x).encode()).hexdigest()))
    manifest={"protocol_sha256":c.sha(P.read_bytes()),"index_sha256":c.sha(raw),"index_url":url,"pool_count":len(order),"candidate_order":order[:S["selection"]["max_queries"]],"selected":[],"attempts":[],"prior_partial_sha256":c.sha(zz.read(nn["manifest.json"]))}
    def add(g):
        manifest["selected"].append(g);names.add(g["target"]);days.add(g["day"]);coords.add(coord(g))
        urls.update(v["url"] for v in g["files"]);filenames.update(v["name"] for v in g["files"]);save("manifest.json",manifest)
        emit("R40N1_SELECTED",{"target":g["target"],"day":g["day"],"count":len(manifest["selected"])})
    for g in old["selected"]:
        new=validate(g["target"],g["files"],g["query"])
        assert all(v["header"]==w["header"] for v,w in zip(g["files"],new["files"]));add(new)
    for target in manifest["candidate_order"]:
        if len(manifest["selected"])==S["selection"]["target_groups"]:break
        if target in names:continue
        attempt={"target":target,"cadences":0,"errors":[]};manifest["attempts"].append(attempt)
        try:
            fs,q=eng.query(target);fs.sort(key=lambda x:(x["day"],x["seconds"],x["url"]))
            groups=[fs[i:i+4] for i in range(len(fs)-3) if [v["position"] for v in fs[i:i+4]]==["S","R","S","R"] and len({v["day"] for v in fs[i:i+4]})==1 and fs[i]["day"] not in days and all(0<v["seconds"]-u["seconds"]<=600 for u,v in zip(fs[i:i+4],fs[i+1:i+4]))]
            groups.sort(key=lambda x:(x[0]["day"]<58575,x[0]["day"],x[0]["seconds"],x[0]["url"]))
            attempt["cadences"]=len(groups)
            for fs in groups[:S["selection"]["max_header_candidates_per_target"]]:
                try:add(validate(target,fs,q));break
                except Exception as e:attempt["errors"].append({"day":fs[0]["day"],"error":type(e).__name__+": "+str(e)[:200]})
        except Exception as e:attempt["errors"].append({"error":type(e).__name__+": "+str(e)[:200]})
        save("manifest.json",manifest)
    assert len({g["day"] for g in manifest["selected"]})==len({coord(g) for g in manifest["selected"]})==len(manifest["selected"])
    save("manifest.json",manifest);emit("R40N1_MANIFEST",{"selected":[[g["target"],g["day"]] for g in manifest["selected"]],"pool":len(order),"queries":len(manifest["attempts"]),"sha256":c.sha((O/"manifest.json").read_bytes())})
    return manifest,fr,hm,m
def makebank(pre):
    out=[]
    for slot,tiles in enumerate(pre):
        for wi,(a,V) in enumerate(tiles):
            Z=b.ref.norm(V.astype(float));ix=np.arange(5,Z.shape[1]-5);u=Z[:,ix]-np.median(Z[:,ix[:,None]+np.array([-5,-4,-3,3,4,5])],axis=2)
            z=np.clip(u,-4,4).mean(0)
            out.append((slot,wi,float(a["channel0_mhz"]),float(a["df_hz"]),z.astype(np.float32)))
    return out
def cross(freq,day,bank):
    out=np.zeros((len(freq),8),float)
    eligible=[(d,items) for d,items in bank if d!=day]
    for i,fq in enumerate(freq):
        vals=[]
        for d,items in eligible:
            values={};res=[]
            for slot,wi,f0,df,z in items:
                x=(fq-f0)*1e6/df;ix=int(round(x));err=abs(x-ix)*abs(df)
                if 5<=ix<5+len(z) and err<=abs(df)*.5000001:
                    values[slot,wi]=(float(z[ix-5]),err)
            ps=[]
            for slot in (0,2):
                for wi in range(4):
                    if (slot,wi) in values and (slot+1,wi) in values:
                        a,e=values[slot,wi];v,ee=values[slot+1,wi];ps.append((min(a,v),v,max(e,ee)))
            if ps:vals.append((max(v[0] for v in ps),max(v[1] for v in ps),max(v[2] for v in ps)))
        if vals:
            v=np.asarray(vals);out[i]=[len(v),len(v)/max(1,len(eligible)),v[:,0].mean(),v[:,0].max(),v[:,1].mean(),v[:,1].max(),(v[:,0]>=.75).mean(),v[:,2].max()]
    assert np.isfinite(out).all();return out
def selftest():
    z=np.zeros(100);z[20]=3.;items=[(i,0,3000.,-2.8,z) for i in range(4)]
    a=cross(np.array([3000.-25*2.8/1e6]),1,[(1,items),(2,items)])
    assert a[0,0]==1 and a[0,2]==3 and a[0,7]<1e-6
    assert not cross(np.array([3000.]),1,[(1,items)]).any()
    emit("R40N1_SELFTEST",{"cross_day_exclusion":"PASS","signed_grid":"PASS"})
def main():
    start=time.perf_counter();selftest();m,fr,hm,exam=select();mh=c.sha((O/"manifest.json").read_bytes())
    cfg=n.S["new_controls"];n.set_config(cfg["seed"],cfg["windows_mhz"],cfg["half_hz"],cfg["max_mib_per_scan"])
    L=h.readlocks();b.taylor.init_tree();lz,ln=n.archive("r40l");raw=lz.read(ln["lock.json"]);assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];old=json.loads(raw)
    kz,kn=n.archive("r40k");kr=n.read(kz,kn,"result.json");kl=n.read(kz,kn,"lock.json")
    rz,rn=c.artifact(S["sources"]["replay"],64);assert c.sha(rz.read(rn["lock.json"]))==S["sources"]["replay"]["lock_sha256"]
    seen=set(fr["excluded_raw_hashes"]);seen.update(v for r in kr["rows"] for v in r["raw_tile_hashes"]);seen.update(v for r in kl["calibration"] for v in r["raw_tile_hashes"])
    for name in rn:
        if name.endswith(".npz"):
            q=f.load(rz,rn,name);seen.update(q["source_tile_sha"].tolist())
    bank=[];controls=[];records=[]
    for g in m["selected"]:
        name=g["target"];rr,R,ch,hs,meta,cm,pre=n.gather(g,name);assert len(hs)==16 and not seen&hs;seen.update(hs)
        ix=np.array([v[1]==0 for v in rr]);assert ix.sum()==128
        rows=[v for v,x in zip(rr,ix) if x];mm=[v for v,x in zip(meta,ix) if x];R=R[ix];d=n.encode(rows,R,mm,pre)
        pan=b.panel(rows,R);arms=b.base_arms(pan,L["r40c"]["original"]["component_fit"])
        q={"V":pan["Y"].astype(np.float32),"y":np.zeros(128,int),"truth":np.stack([v[7] for v in rows]),"active":np.stack([v[8] for v in rows]),"amp":np.array([v[4] for v in rows]),"tile_sha":np.array([v[6] for v in rows])}
        for a in t.S["paths"]:q[a+"_score"],q[a+"_paths"]=arms[a][:2]
        hh=f.historical(q,L["r40c"]);_,_,_,a=h.replay(hh,L);a=a.astype(bool)
        d.update(parent_alert=a,parent_hit=np.zeros(128,bool),joint_rescue=(~a)&(d["s"]>=old["parent_joint_min_threshold"]),regime_high=d["regime"]>old["reference_cut"])
        items=makebank(pre);bank.append((int(g["day"]),items))
        # Export independent day reference contrast grids, not synthetic injection samples.
        np.savez_compressed(O/("bank_"+name+".npz"),index=np.asarray([[slot,wi,f0,df] for slot,wi,f0,df,z in items]),**{"tile_"+str(i):v[4] for i,v in enumerate(items)})
        controls.append((name,int(g["day"]),d));row={"target":name,"day":g["day"],"controls":128,"parent_alerts":int(a.sum()),"extra_joint_alerts":int(d["joint_rescue"].sum()),"raw_tile_hashes":sorted(hs)};records.append(row);save("partial.json",records);emit("R40N1_CONTROL",row)
        assert c.sha((O/"manifest.json").read_bytes())==mh
    for name,day,d in controls:
        d["crossday_features"]=cross(d["frequency_mhz"],day,bank);np.savez_compressed(O/("control_"+name+".npz"),**d)
    assert controls
    rg=np.concatenate([d["regime_high"] for _,_,d in controls]);pa=np.concatenate([d["parent_alert"] for _,_,d in controls]);ta=np.repeat([day for _,day,_ in controls],128)
    stats={str(i):{"controls":int(((rg==i)&~pa).sum()),"days":len(set(ta[(rg==i)&~pa].tolist()))} for i in (0,1)}
    full=len(controls)==8 and len(rg)==1024;qualified=full and all(v["controls"]>=128 and v["days"]>=3 for v in stats.values())
    receipt={"status":"PASS_KUSTO_COUNTS_ONLY" if qualified else "HOLD_KUSTO_PARTIAL_OR_SPARSE","manifest_sha256":mh,"controls":len(rg),"groups":len(controls),"strata":stats,"fallback":not qualified,"physical_null_verified":False,"alias_completeness_proven":False,"formal_fap":False,"records":records,"production_adoption":False}
    save("controls.json",receipt);emit("R40N1_CONTROL_SUMMARY",{k:v for k,v in receipt.items() if k!="records"})
    if not full:
        save("result.json",receipt);return
    daymap={g["target"]:int(g["day"]) for g in hm["selected"]+exam["selected"]}
    def load(name):
        d=f.load(rz,rn,name);target=name.removesuffix(".npz").split("_")[-1];d["crossday_features"]=cross(d["frequency_mhz"],daymap[target],bank);return d
    train=[load(name) for name in sorted(rn) if name.startswith(("historical_eval_","historical_shared_")) and name.endswith(".npz")];assert len(train)==8
    models={};thresholds={};serial={}
    def use(d):return ~d["parent_alert"]
    def null(d):return use(d)&(d["y"]==0)
    def pos(d):return d["joint_rescue"]&d["correct"]
    y=np.concatenate([d["correct"][use(d)] for d in train]);assert y.any() and (~y).any()
    for feature in ("raw","crossday"):
        def fv(d,feature=feature):return d["raw_features"] if feature=="raw" else np.column_stack([d["raw_features"],d["crossday_features"]])
        X=np.concatenate([fv(d)[use(d)] for d in train]);sc=StandardScaler().fit(X)
        for kind in ("logistic","tree"):
            model=LogisticRegression(C=.1,class_weight="balanced",max_iter=1000,random_state=401).fit(sc.transform(X),y) if kind=="logistic" else ExtraTreesClassifier(n_estimators=128,max_depth=4,min_samples_leaf=8,class_weight="balanced",random_state=401,n_jobs=1).fit(X,y)
            def score(d,model=model,sc=sc,fv=fv,kind=kind):return model.decision_function(sc.transform(fv(d))) if kind=="logistic" else model.predict_proba(fv(d))[:,1]
            name=feature+"_"+kind;models[name]=score
            gp=np.concatenate([score(d)[pos(d)] for d in train]);gc=np.concatenate([score(d)[null(d)] for _,_,d in controls]);reg=rg[~pa];assert len(gp)>0 and len(gc)==len(reg)
            pooled=l.cut(gc,.025);conditional=[l.cut(gc[reg==i],.025) for i in (0,1)] if qualified else [pooled,pooled]
            thresholds[name]={"preserve":float(np.nextafter(gp.min(),-np.inf)),"pooled":pooled,"kusto":conditional}
            serial[name]={"mean":sc.mean_.tolist(),"scale":sc.scale_.tolist(),"coef":model.coef_.tolist(),"intercept":model.intercept_.tolist()} if kind=="logistic" else {"trees":[{"children_left":e.tree_.children_left.tolist(),"children_right":e.tree_.children_right.tolist(),"feature":e.tree_.feature.tolist(),"threshold":e.tree_.threshold.tolist(),"value":e.tree_.value.tolist()} for e in model.estimators_]}
    lock={"protocol_sha256":c.sha(P.read_bytes()),"manifest_sha256":mh,"source_lock_sha256":S["sources"]["replay"]["lock_sha256"],"thresholds":thresholds,"models":serial,"kusto_qualified":qualified,"bank_days":[day for day,_ in bank],"bank_sha256":{p.name:c.sha(p.read_bytes()) for p in sorted(O.glob("bank_*.npz"))},"train_rows":len(y)}
    save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes());emit("R40N1_LOCK",{k:v for k,v in lock.items() if k!="models"})
    def measure(d):
        a=d["parent_alert"];base=d["parent_hit"];res=d["joint_rescue"];y=d["y"];weak=d["weak"];correct=d["correct"];oldhit=base|(res&correct);reg=(d["regime"]>old["reference_cut"]).astype(int)
        guards={"joint_min":np.ones(len(y),bool),"r40g":np.zeros(len(y),bool)}
        if "learned05_alert" in d:guards["learned05"]=d["learned05_alert"]&~a
        for name,score in models.items():
            v=score(d);th=thresholds[name]
            for cut in ("preserve","pooled","kusto"):guards[name+"_"+cut]=v>=(np.asarray(th[cut])[reg] if cut=="kusto" else th[cut])
        out={}
        for name,g in guards.items():
            rescue=res&g;aa=a|rescue;hh=base|(rescue&correct)
            out[name]={"localized":int(hh.sum()),"weak":int((hh&weak).sum()),"control_alerts":int((aa&(y==0)).sum()),"lost_parent":int((base&~hh).sum()),"lost_joint":int((oldhit&~hh).sum()),"lost_joint_weak":int((oldhit&~hh&weak).sum()),"false_paths_on_positive":int((aa&(y==1)&~hh).sum())}
            d[name+"_n1_alert"]=aa;d[name+"_n1_hit"]=hh
        return out
    for name,day,d in controls:emit("R40N1_CALIBRATION",{"target":name,"metrics":measure(d)})
    results=[]
    for name in sorted(rn):
        if not name.startswith(("eval_","shared_")) or not name.endswith(".npz"):continue
        d=load(name);row={"record":name,"metrics":measure(d),"crossday_covered":int((d["crossday_features"][:,0]>0).sum())};results.append(row);save("replay_partial.json",results);emit("R40N1_GROUP",row);np.savez_compressed(O/name,**d)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def total(kind):
        rr=[r["metrics"] for r in results if r["record"].startswith(kind+"_")];assert len(rr)==5
        return {a:{k:sum(v[a][k] for v in rr) for k in rr[0][a]} for a in rr[0]}
    e=total("eval");sh=total("shared");assert [e["joint_min"][k] for k in ("localized","weak","control_alerts")]==[215,43,43] and [sh["joint_min"][k] for k in ("localized","weak")]==[111,36]
    gates={a:qualified and e[a]["localized"]>=215 and e[a]["weak"]>=43 and e[a]["control_alerts"]<=35 and sh[a]["localized"]>=111 and sh[a]["weak"]>=36 and e[a]["lost_joint"]==sh[a]["lost_joint"]==0 for a in e}
    result={"status":"DEVELOPMENT_REPLAY_ONLY","controls":receipt,"eval":e,"shared":sh,"gates":gates,"lock_sha256":lh,"rows":results,"seconds":time.perf_counter()-start,"production_adoption":False};save("result.json",result);emit("R40N1_SUMMARY",{k:v for k,v in result.items() if k not in ("rows","controls")})
if __name__=="__main__":main()
