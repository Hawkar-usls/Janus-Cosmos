#!/usr/bin/env python3
"""R40N6: disjoint same-band training and frozen new-day holdout."""
import sys
sys.argv[1:]=["controls"]
import json,pathlib,time,hashlib,re,io
import numpy as np
from sklearn.metrics import roc_auc_score
import r40n5a_verified_frequency as a
import r40n2_instrument_controls as q
n=a.n;c=a.c;b=a.b;l=a.l;f=a.f;h=n.h;t=n.t
qual=q.qual;scout=q.scout;eng=q.eng;index_query=q.index_query
x=a.x
model_fit=a.model_fit;col_audit=a.col_audit;seven_day=a.seven_day;central=a.central;neighborhood=a.neighborhood;matched_placebo=a.matched_placebo;coverage=a.coverage
P=pathlib.Path("research/seti/r40n6_same_band_protocol.json");S=json.loads(P.read_text());C=S["selection"]
O=pathlib.Path("runtime/r40n6");O.mkdir(parents=True,exist_ok=True);n.O=O;n.k.O=O;n.l.O=O;eng.OUT=O
def save(name,d):(O/name).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def archive(k):return c.artifact(S["sources"][k],64)
def select(old):
    names,days,ns,urls,coords=qual.history()
    hm,m,fr,sel=n.manifests();zz,nn=n.archive("r37_manifest");rm=n.read(zz,nn,"manifest.json")
    cde=json.loads(pathlib.Path("research/seti/r40cde_protocol.json").read_text())
    names.update(sel["excluded_names"]);names.update(cde["targets"]);days.update(sel["excluded_days"]);days.update(cde["days"])
    prior=hm["selected"]+m["selected"]+rm["selected"]["Parkes"]+rm["selected"].get("GBT",[])+old["selected"]
    # Identical or ambiguous historical declinations use R37B conservative sky exclusion.
    for g in prior:
        names.add(g["target"]);days.add(int(g["day"]))
        for v in g["files"]:ns.add(v["name"]);urls.add(v["url"])
        hd=g["files"][0]["header"];mt=g["files"][0].get("metadata",{})
        coords.append({"target":g["target"],**qual.history_pointing(hd,mt)})
    raw,u=scout.body(scout.S["index"]);strings=set(scout.walk(json.loads(raw)))
    (O/"target-index.json").write_bytes(raw)
    aliases={x[:-2].upper() for x in strings if re.fullmatch(r"[A-Za-z0-9_-]{2,48}_S",x,re.I)}
    gb={x[:-4].upper() for x in strings if re.fullmatch(r"[A-Za-z0-9_-]{2,48}_OFF",x,re.I)}
    published=set(scout.S["gbt_targets"])
    order=lambda pool:sorted(pool-names,key=lambda x:(x not in published,hashlib.sha256((C["seed"]+"|"+x).encode()).hexdigest()))
    queues={"Parkes":order(aliases)[:C["max_targets"]]}
    # Override only metadata windows, not cadence or pointing checks.
    qual.S["parkes"]={**qual.S["parkes"],"windows_mhz":C["windows_mhz"],"half_hz":C["half_hz"]}
    manifest={"protocol_sha256":c.sha(P.read_bytes()),"target_index_sha256":c.sha(raw),"old_manifest_sha256":c.sha(json.dumps(old,sort_keys=True).encode()),"candidate_counts":{k:len(v) for k,v in queues.items()},"candidate_order":queues,"selected":[],"attempts":[]}
    save("manifest.json",manifest)
    for tel in ("Parkes",):
        for target in queues[tel]:
            if len(manifest["selected"])>=C["groups"]:break
            if target in names:continue
            a={"telescope":tel,"target":target,"errors":[]};manifest["attempts"].append(a)
            try:
                rr=[index_query(tel,target,q) for q in ((target,target+"_OFF") if tel=="GBT" else (target,))]
                fs=qual.entries(rr,tel,target,days,ns,urls);cand=qual.cadence(fs);a["files"]=len(fs);a["cadences"]=len(cand);a["index_sha256"]=[v["index_sha256"] for v in rr]
                for group in cand[:C["header_attempts"]]:
                    try:
                        hs,pp=qual.qualify(tel,target,group,coords)
                        day=int(group[0]["day"]);assert day not in days and target not in names
                        g={"telescope":tel,"target":target,"day":day,"pointing":pp[0],"split":"unassigned","files":[{**v,"header":hd,"role":role} for v,hd,role in zip(group,hs,("ON1","OFF1","ON2","OFF2"))]}
                        manifest["selected"].append(g);names.add(target);days.add(day);coords.append({"target":target,"pointing":pp[0],"status":"FULL_HEADER_POINTING"});ns.update(v["name"] for v in group);urls.update(v["url"] for v in group)
                        emit("R40N6_SELECTED",{"telescope":tel,"target":target,"day":day,"count":len(manifest["selected"])});break
                    except Exception as e:a["errors"].append({"day":group[0]["day"],"error":type(e).__name__+": "+str(e)[:200]})
            except Exception as e:a["errors"].append({"error":type(e).__name__+": "+str(e)[:200]})
            save("manifest.json",manifest)
    assert len({g["day"] for g in manifest["selected"]})==len({g["target"] for g in manifest["selected"]})==len(manifest["selected"])
    if len(manifest["selected"])>=C["min_groups"]:
        for i,g in enumerate(manifest["selected"]):g["split"]="fit" if i<len(manifest["selected"])-C["holdout_groups"] else "holdout"
    save("manifest.json",manifest);emit("R40N6_MANIFEST",{"targets":[[g["telescope"],g["target"],g["day"]] for g in manifest["selected"]],"queries":len(manifest["attempts"]),"candidate_counts":manifest["candidate_counts"],"sha256":c.sha((O/"manifest.json").read_bytes())})
    return manifest,fr

def main():
    start=time.perf_counter();a.selftest();n.selftest()
    z,nn=archive("n3");rz,rn=archive("replay");uz,un=archive("n5a");az,an=archive("n1");bz,bn=archive("n2")
    for key,zz,ns in (("n3",z,nn),("replay",rz,rn),("n5a",uz,un)):assert c.sha(zz.read(ns["lock.json"]))==S["sources"][key]["lock_sha256"]
    lock=n.read(z,nn,"lock.json");man=n.read(z,nn,"manifest.json")
    assert c.sha(z.read(nn["manifest.json"]))==lock["manifest_sha256"]
    ma=n.read(az,an,"manifest.json");mb=n.read(bz,bn,"manifest.json")
    assert c.sha(az.read(an["manifest.json"]))==man["n1_manifest_sha256"] and c.sha(bz.read(bn["manifest.json"]))==man["n2_manifest_sha256"]
    oldgroups=ma["selected"]+mb["selected"];assert len(oldgroups)==8
    bank=[];days={};bank_hashes=dict(lock["bank_sha256"])
    for g in oldgroups:
        name=g["target"];day=int(g["day"]);fn="bank_"+name+".npz";raw=z.read(nn[fn]);assert c.sha(raw)==bank_hashes[fn]
        with np.load(io.BytesIO(raw),allow_pickle=False) as v:
            ids=v["index"];items=[(int(p[0]),int(p[1]),float(p[2]),float(p[3]),v["tile_"+str(i)].copy()) for i,p in enumerate(ids)]
        bank.append((day,items));days[name]=day
    assert len(set(days.values()))==8
    bank.sort(key=lambda v:v[0])
    manifest,fr=select({"selected":oldgroups});mh=c.sha((O/"manifest.json").read_bytes())
    groups0=manifest["selected"];assert not {int(g["day"]) for g in groups0}&set(days.values())
    if len(groups0)<C["min_groups"]:
        result={"status":"HOLD_INSUFFICIENT_DISJOINT_METADATA_GROUPS","groups":len(groups0),"required":C["min_groups"],"queries":len(manifest["attempts"]),"manifest_sha256":mh,"production_adoption":False}
        save("result.json",result);emit("R40N6_SUMMARY",result);return
    L=h.readlocks();b.taylor.init_tree()
    lz,ln=n.archive("r40l");raw=lz.read(ln["lock.json"]);assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];old=json.loads(raw)
    kz,kn=n.archive("r40k");kr=n.read(kz,kn,"result.json");kl=n.read(kz,kn,"lock.json")
    seen=set(fr["excluded_raw_hashes"])
    for v in man["records"]+kr["rows"]+kl["calibration"]:seen.update(v["raw_tile_hashes"])
    for v in n.read(uz,un,"bank_manifest.json")+n.read(uz,un,"unused_bank_manifest.json")["records"]:seen.update(v["raw_tile_hashes"])
    for name in rn:
        if name.endswith(".npz"):
            with np.load(io.BytesIO(rz.read(rn[name])),allow_pickle=False) as v:seen.update(v["source_tile_sha"].tolist())
    def augment(d,day):
        fq=d["frequency_mhz"];coverage(fq,C["windows_mhz"])
        d["central"]=seven_day(central,fq,bank,day)
        d["neighborhood"]=seven_day(neighborhood,fq,bank,day);d["placebo"]=seven_day(matched_placebo,fq,bank,day)
        assert all(np.all(d[k][:,0]==1) for k in ("central","neighborhood","placebo")),"Incomplete recurrence coverage"
        return d
    controls=[]
    for name,day in sorted(days.items()):
        d=f.load(z,nn,"control_"+name+".npz");assert np.array_equal(d["crossday_features"],x.cross(d["frequency_mhz"],day,bank))
        assert len(d["y"])==128 and not d["y"].any();controls.append((name,augment(d,day)))
    n.set_config(C["seed"],C["windows_mhz"],C["half_hz"],C["max_mib_per_scan"])
    def acquire(g):
        name=g["target"];day=int(g["day"])
        assert c.sha((O/"manifest.json").read_bytes())==mh
        rr,R,ch,hs,meta,cm,pre=n.gather(g,g["split"]+"_"+name);assert len(hs)==16 and not hs&seen;seen.update(hs)
        out=[]
        for kind,rows,refs,mm in (("eval",rr,R,meta),("shared",[v[0] for v in ch],np.stack([v[1] for v in ch]),cm)):
            d=n.encode(rows,refs,mm,pre);pan=b.panel(rows,refs);arms=b.base_arms(pan,L["r40c"]["original"]["component_fit"])
            inp={"V":pan["Y"].astype(np.float32),"y":np.array([v[1] for v in rows]),"truth":np.stack([v[7] for v in rows]),"active":np.stack([v[8] for v in rows]),"amp":np.array([v[4] for v in rows]),"tile_sha":np.array([v[6] for v in rows])}
            for path in t.S["paths"]:inp[path+"_score"],inp[path+"_paths"]=arms[path][:2]
            alert,path,hit=l.parent(inp,L);assert np.array_equal(inp["y"],d["y"])
            joint=(~alert)&(d["s"]>=old["parent_joint_min_threshold"])
            v=old["models"] if "models" in old else None
            d.update(parent_alert=alert,parent_hit=hit,parent_path=path,joint_rescue=joint,truth=inp["truth"],active=inp["active"],amp=inp["amp"])
            augment(d,day);out.append((kind,d))
        rec={"target":name,"day":day,"split":g["split"],"raw_tile_hashes":sorted(hs),"ordinary":256,"shared":32}
        emit("R40N6_ACQUIRED",rec)
        return rec,out
    train=[];records=[]
    for g in groups0:
        if g["split"]!="fit":continue
        rec,ds=acquire(g);records.append(rec)
        for kind,d in ds:
            train.append((g["target"],d));np.savez_compressed(O/("fit_"+kind+"_"+g["target"]+".npz"),**d)
        save("fit_records.json",records)
    assert len(train)==2*(len(groups0)-2)
    def use(d):return ~d["parent_alert"]
    def null(d):return use(d)&(d["y"]==0)
    def pos(d):return d["joint_rescue"]&d["correct"]
    reg=np.concatenate([(d["regime"]>old["reference_cut"])[null(d)] for _,d in controls]);assert (reg==0).sum()==(reg==1).sum()==487
    y=np.concatenate([d["correct"][use(d)] for _,d in train]);groups=np.concatenate([[name]*int(use(d).sum()) for name,d in train]);assert y.any() and (~y).any() and len(set(groups))>=4
    audit={}
    for field in ("central","neighborhood","placebo"):
        X=np.concatenate([d[field][use(d)] for _,d in train]);audit[field]=col_audit(X,y)
    save("audit.json",audit);emit("R40N6_AUDIT",audit)
    variable={field:sum(v["std"]>1e-12 for v in audit[field][1:]) for field in ("central","neighborhood","placebo")}
    save("variation_gate.json",variable);emit("R40N6_VARIATION",variable)
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
    save("oof.json",auc);emit("R40N6_OOF",auc)
    lk={"protocol_sha256":c.sha(P.read_bytes()),"source_lock_sha256":S["sources"]["n3"]["lock_sha256"],"bank_sha256":bank_hashes,"selection_manifest_sha256":mh,"fit_records_sha256":c.sha((O/"fit_records.json").read_bytes()),"bank_days":sorted(days.values()),"models":serial,"thresholds":cuts,"same_band_oof_auc":auc,"reference_cut":old["reference_cut"],"train_rows":len(y),"calibration_controls":1024,"calibration_parent_rejected":974,"calibration_strata":[487,487],"bank_support":"seven days, own day excluded, external average of 8 seven-day subsets"}
    save("lock.json",lk);lh=c.sha((O/"lock.json").read_bytes());emit("R40N6_LOCK",{k:v for k,v in lk.items() if k!="models"})
    def measure(d,baseline=False):
        a=d["parent_alert"];base=d["parent_hit"];res=d["joint_rescue"];y=d["y"];weak=d["weak"];correct=d["correct"];oldhit=base|(res&correct);reg=(d["regime"]>old["reference_cut"]).astype(int)
        guards={"joint_min":np.ones(len(y),bool),"r40g":np.zeros(len(y),bool)}
        if baseline:
            guards["learned05"]=d["learned05_alert"]&~a
        for name,score in models.items():
            v=score(d)
            for cut,th in cuts[name].items():guards[name+"_"+cut]=v>=(np.asarray(th)[reg] if isinstance(th,list) else th)
        out={}
        for name,g in guards.items():
            rescue=res&g;aa=a|rescue;hh=base|(rescue&correct)
            out[name]={"localized":int(hh.sum()),"weak":int((hh&weak).sum()),"control_alerts":int((aa&(y==0)).sum()),"lost_parent":int((base&~hh).sum()),"lost_joint":int((oldhit&~hh).sum()),"lost_joint_weak":int((oldhit&~hh&weak).sum()),"false_paths_on_positive":int((aa&(y==1)&~hh).sum())}
            d[name+"_n6_alert"]=aa;d[name+"_n6_hit"]=hh
        return out

    cal=[]
    for name,d in controls:cal.append({"target":name,"metrics":measure(d)})
    save("calibration.json",cal)
    results=[];holdout_records=[]
    for g in groups0:
        if g["split"]!="holdout":continue
        assert c.sha((O/"lock.json").read_bytes())==lh
        rec,ds=acquire(g);holdout_records.append(rec);save("holdout_records.json",holdout_records)
        for kind,d in ds:
            # R40L learned score is reproduced from its frozen logistic parameters.
            # Use the exact frozen R40L model/threshold keys, see source contract below.
            row={"record":kind+"_"+g["target"],"metrics":measure(d,True)}
            results.append(row);np.savez_compressed(O/("holdout_"+kind+"_"+g["target"]+".npz"),**d)
            save("partial.json",results);emit("R40N6_GROUP",row)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def total(kind):
        rr=[r["metrics"] for r in results if r["record"].startswith(kind+"_")];assert len(rr)==2
        return {a:{k:sum(v[a][k] for v in rr) for k in rr[0][a]} for a in rr[0]}
    e=total("eval");sh=total("shared");base=e["learned05"]["control_alerts"]
    gates={a:e[a]["lost_joint"]==sh[a]["lost_joint"]==0 and e[a]["lost_parent"]==sh[a]["lost_parent"]==0 and e[a]["control_alerts"]*640<=35*256 and e[a]["control_alerts"]*8<=7*base for a in e}
    result={"status":"NEW_DAY_DEVELOPMENT_HOLDOUT_ONLY","manifest_sha256":mh,"fit_groups":len(groups0)-2,"holdout_groups":2,"fit_records":records,"holdout_records":holdout_records,"eval":e,"shared":sh,"gates":gates,"oof":auc,"rows":results,"lock_sha256":lh,"seconds":time.perf_counter()-start,"production_adoption":False};save("result.json",result);emit("R40N6_SUMMARY",{k:v for k,v in result.items() if k not in ("rows","fit_records","holdout_records")})
if __name__=="__main__":main()
