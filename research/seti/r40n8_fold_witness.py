#!/usr/bin/env python3
"""R40N8: continuity and disjoint temporal-fold alternative witnesses."""
import sys
sys.argv[1:]=["controls"]
import json,pathlib,time,hashlib,re,io
import numpy as np
from sklearn.metrics import roc_auc_score
import r40n7_parent_trajectory as w
v=w.v;n=w.n;c=w.c;b=w.b;l=w.l;f=w.f;h=w.h;t=w.t
qual=w.qual;scout=w.scout;eng=w.eng;index_query=w.index_query
P=pathlib.Path("research/seti/r40n8_fold_witness_protocol.json");S=json.loads(P.read_text());C=S["selection"]
O=pathlib.Path("runtime/r40n8");O.mkdir(parents=True,exist_ok=True);n.O=O;n.k.O=O;n.l.O=O;eng.OUT=O
def save(name,d):(O/name).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def archive(k):return c.artifact(S["sources"][k],128)
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
                        emit("R40N8_SELECTED",{"telescope":tel,"target":target,"day":day,"count":len(manifest["selected"])});break
                    except Exception as e:a["errors"].append({"day":group[0]["day"],"error":type(e).__name__+": "+str(e)[:200]})
            except Exception as e:a["errors"].append({"error":type(e).__name__+": "+str(e)[:200]})
            save("manifest.json",manifest)
    assert len({g["day"] for g in manifest["selected"]})==len({g["target"] for g in manifest["selected"]})==len(manifest["selected"])
    if len(manifest["selected"])>=C["min_groups"]:
        for g in manifest["selected"]:g["split"]="holdout"
    save("manifest.json",manifest);emit("R40N8_MANIFEST",{"targets":[[g["telescope"],g["target"],g["day"]] for g in manifest["selected"]],"queries":len(manifest["attempts"]),"candidate_counts":manifest["candidate_counts"],"sha256":c.sha((O/"manifest.json").read_bytes())})
    return manifest,fr



FOLDS=[(np.arange(8),np.arange(8,16)),(np.arange(8,16),np.arange(8)),(np.arange(0,16,2),np.arange(1,16,2)),(np.arange(1,16,2),np.arange(0,16,2))]
def runstats(z,valid):
    vals=[]
    def longest(mask):
        a=np.zeros(len(mask),int);b=a.copy()
        for j in range(16):a=np.where(mask[:,j],a+1,0);b=np.maximum(b,a)
        return b
    for threshold in (.5,1.,2.):
        hit=(z>=threshold)&valid;gap=(~hit)&valid
        adj=valid[:,:-1]&valid[:,1:]
        isolated=hit&~np.pad(hit[:,:-1],((0,0),(1,0)))&~np.pad(hit[:,1:],((0,0),(0,1)))
        vals.extend([longest(hit)/16,longest(gap)/16,((hit[:,1:]!=hit[:,:-1])&adj).sum(1)/15,isolated.sum(1)/16,hit.sum(1)/16,hit.sum(1)/np.maximum(1,valid.sum(1))])
    return vals
def extra(X,D,p):
    N=len(X);ii=np.arange(N)[:,None];rr=np.arange(16)[None,:]
    valid=(p>=5)&(p<123);pp=np.clip(p,5,122);ok=(valid.sum(1)>=12)
    def mean(z,mask):return (z*mask).sum(1)/np.maximum(1,mask.sum(1))
    U=n.normalized_contrast(X);V=n.normalized_contrast(D);run=[];fold=[];paths=[];scores=[]
    for Z in (U,V):
        sel=Z[ii,rr,pp];run+=runstats(sel,valid)
        # Native ±2 residuals; strongest local bin relative to selected path.
        # Restrict all candidate neighbours to legal contrast support.
        nei=pp[:,:,None]+np.arange(-2,3);nv=(nei>=5)&(nei<123)
        z=Z[ii[:,:,None],rr[:,:,None],np.clip(nei,5,122)];z=np.where(nv,z,-np.inf)
        best=z.max(2);offset=z.argmax(2)-2;run += [mean(abs(offset),valid),mean(offset==0,valid),mean(best-sel,valid),mean(abs(np.diff(offset,axis=1)),valid[:,:-1]&valid[:,1:])]
        fps=[];fsc=[]
        for tr,te in FOLDS:
            vt=valid[:,tr];ve=valid[:,te];ok &= (vt.sum(1)>=4)&(ve.sum(1)>=4)
            best=np.full(N,-np.inf);bestp=np.zeros((N,16),int)
            for delta in b.D:
                lo=5-int(delta.min());hi=123-int(delta.max())
                if hi<=lo:continue
                starts=np.arange(lo,hi);path=starts[:,None]+delta
                # No test-row evidence enters the choice of alternative.
                distinct=(abs(path[None,:,tr]-p[:,None,tr])>2).sum(2)>=3
                sc=sum(Z[:,j,lo+delta[j]:hi+delta[j]] for j in tr)/len(tr)
                sc=np.where(distinct,sc,-np.inf);arg=sc.argmax(1);vv=sc[np.arange(N),arg];im=vv>best
                best[im]=vv[im];bestp[im]=path[arg[im]]
            assert np.isfinite(best).all()
            a=Z[ii,rr,bestp];check=a[:,tr].mean(1);held=a[:,te].mean(1)
            assert np.allclose(check,best,atol=1e-12,rtol=1e-12) and ((abs(bestp[:,tr]-p[:,tr])>2).sum(1)>=3).all()
            ct=mean(sel[:,tr],vt);ce=mean(sel[:,te],ve)
            fold.extend([check,held,ct-check,ce-held,check-held,ct-ce])
            fps.append(bestp);fsc.append(np.column_stack([check,held,ct,ce]))
        for j,k in ((0,1),(2,3)):
            fold += [np.mean(abs(fps[j]-fps[k]),axis=1),np.mean(abs(fps[j]-fps[k])<=2,axis=1)]
        paths.append(np.stack(fps,axis=1));scores.append(np.stack(fsc,axis=1))
    R=np.column_stack(run);F=np.column_stack(fold);P0=np.stack(paths,axis=1);Q=np.stack(scores,axis=1)
    assert R.shape==(N,44) and F.shape==(N,56) and np.isfinite(R).all() and np.isfinite(F).all()
    return R,F,ok,P0,Q
def selftest():
    rng=np.random.default_rng(808);X=rng.normal(0,.05,(2,16,128));D=X.copy();p=np.full((2,16),40);p[1,:8]=0
    X[0,:,40]+=3;D[0,:8,90]+=5;D[0,8:,100]+=5
    r,f0,ok,paths,sc=extra(X,D,p)
    assert ok.tolist()==[True,False] and r[0,0]==1 and r[0,1]==0
    D2=D.copy();D2[0,8:,100]+=100;D2[0,8:,90]+=100
    _,_,_,paths2,_=extra(X,D2,p)
    assert np.array_equal(paths[0,1,0],paths2[0,1,0]),"test rows changed train-selected rival"
    emit("R40N8_SELFTEST",{"run_gap_topology":"PASS","heldout_rival_selection_invariance":"PASS","edge_abstention":"PASS"})
def frozen_score(X,L,kind):
    if kind=="tree":return h.tree(X,L)
    return ((X-np.asarray(L["mean"]))/np.asarray(L["scale"]))@np.asarray(L["coef"])[0]+L["intercept"][0]
def main():
    start=time.perf_counter();w.selftest();selftest();n.selftest()
    zz,nn=archive("n6");z7,n7=archive("n7")
    for key,za,na in (("n6",zz,nn),("n7",z7,n7)):
        assert c.sha(za.read(na["lock.json"]))==S["sources"][key]["lock_sha256"] and c.sha(za.read(na["manifest.json"]))==S["sources"][key]["manifest_sha256"]
    nm=n.read(zz,nn,"manifest.json");m7=n.read(z7,n7,"manifest.json");lock7=n.read(z7,n7,"lock.json")
    fitrecs=n.read(zz,nn,"fit_records.json");holdrecs=n.read(zz,nn,"holdout_records.json");freshrecs=n.read(z7,n7,"fresh_records.json")
    priorrecs={r["target"]:r for r in fitrecs+holdrecs+freshrecs};assert len(priorrecs)==12
    az,an=archive("n1");bz,bn=archive("n2");uz,un=archive("n5a")
    assert c.sha(uz.read(un["lock.json"]))==S["sources"]["n5a"]["lock_sha256"]
    prior=n.read(az,an,"manifest.json")["selected"]+n.read(bz,bn,"manifest.json")["selected"]+nm["selected"]+m7["selected"]
    manifest,fr=select({"selected":prior});mh=c.sha((O/"manifest.json").read_bytes())
    fresh=manifest["selected"];assert not {int(g["day"]) for g in fresh}&{int(g["day"]) for g in prior}
    L=h.readlocks();b.taylor.init_tree()
    lz,ln=n.archive("r40l");raw=lz.read(ln["lock.json"]);assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];old=json.loads(raw)
    kz,kn=n.archive("r40k");kr=n.read(kz,kn,"result.json");kl=n.read(kz,kn,"lock.json")
    seen=set(fr["excluded_raw_hashes"])
    for row in kr["rows"]+kl["calibration"]+fitrecs+holdrecs+freshrecs+n.read(uz,un,"bank_manifest.json")+n.read(uz,un,"unused_bank_manifest.json")["records"]:seen.update(row["raw_tile_hashes"])
    for z0,n0 in ((az,an),(bz,bn)):
        for key in n0:
            if key.endswith(".npz"):
                with np.load(io.BytesIO(z0.read(n0[key])),allow_pickle=False) as vv:
                    if "source_tile_sha" in vv:seen.update(vv["source_tile_sha"].tolist())
    rrz,rrn=c.artifact(v.S["sources"]["replay"],64)
    for key in rrn:
        if key.endswith(".npz"):
            with np.load(io.BytesIO(rrz.read(rrn[key])),allow_pickle=False) as vv:seen.update(vv["source_tile_sha"].tolist())
    train=[];audit=[];dev=[];records=[];all_old=set()
    def oldkeep(d,k,cut):
        X=np.column_stack([d["temporal"],d["rival"]]);sc=frozen_score(X,lock7["models"][k],k.split("_")[-1]) if len(X) else np.empty(0)
        keep=np.ones(len(d["y"]),bool);keep[d["idx"]]=(~d["supported"])|(sc>=lock7["cuts"][k][cut]);return keep,sc
    def acquire(g,source):
        name=g["target"];day=int(g["day"]);split=g["split"]
        assert c.sha((O/"manifest.json").read_bytes())==mh
        n.set_config(v.C["seed"] if source=="n6" else w.C["seed"] if source=="n7" else C["seed"],C["windows_mhz"],C["half_hz"],C["max_mib_per_scan"])
        rr,R,ch,hs,meta,cm,pre=n.gather(g,source+"_"+name);assert len(hs)==16
        if source!="fresh":
            assert sorted(hs)==priorrecs[name]["raw_tile_hashes"] and not hs&all_old;all_old.update(hs)
        else:assert not hs&seen;seen.update(hs)
        ds=[]
        for kind,rows,refs,mm in (("eval",rr,R,meta),("shared",[q[0] for q in ch],np.stack([q[1] for q in ch]),cm)):
            q=n.encode(rows,refs,mm,pre);pan=b.panel(rows,refs);arms=b.base_arms(pan,L["r40c"]["original"]["component_fit"])
            inp={"V":pan["Y"].astype(np.float32),"y":np.array([r[1] for r in rows]),"truth":np.stack([r[7] for r in rows]),"active":np.stack([r[8] for r in rows]),"amp":np.array([r[4] for r in rows]),"tile_sha":np.array([r[6] for r in rows])}
            for path in t.S["paths"]:inp[path+"_score"],inp[path+"_paths"]=arms[path][:2]
            alert,path,hit=l.parent(inp,L);alert=alert.astype(bool);hit=hit.astype(bool)
            joint=(~alert)&(q["s"]>=old["parent_joint_min_threshold"])
            learned=((q["F"]-np.asarray(old["mean"]))/np.asarray(old["scale"]))@np.asarray(old["coef"])[0]+old["intercept"][0]
            learn=alert|(joint&(learned>=old["thresholds"]["learned05"]))
            d={"parent_alert":alert,"parent_hit":hit,"parent_path":path,"joint_rescue":joint,"learned05_alert":learn,"correct":q["correct"],"y":q["y"],"weak":q["weak"],"truth":inp["truth"],"active":inp["active"],"amp":inp["amp"],"source_tile_sha":inp["tile_sha"],"p":q["p"]}
            idx=np.flatnonzero(alert);pp=path[idx].astype(int)
            if len(idx):
                X=n.ref.norm(np.stack([rows[i][0] for i in idx]));D=n.ref.norm(refs[idx])
                st,rv,supported,a,off,rpaths,best=w.witness(X,D,pp);runs,folds,ok,fp,fs=extra(X,D,pp)
            else:
                st=np.empty((0,36));rv=np.empty((0,6));supported=np.empty(0,bool);a=np.empty((0,16));off=a.copy();rpaths=np.empty((0,2,16),int);best=np.empty((0,2));runs=np.empty((0,44));folds=np.empty((0,56));ok=np.empty(0,bool);fp=np.empty((0,2,4,16),int);fs=np.empty((0,2,4,4))
            d.update(idx=idx,temporal=st,rival=rv,supported=supported,selected_on=a,selected_off=off,rival_paths=rpaths,rival_scores=best)
            if source!="fresh":
                fn=("development_"+split+"_" if source=="n6" else "holdout_")+kind+"_"+name+".npz"
                before=f.load(z7,n7,fn)
                for k,value in d.items():assert np.array_equal(value,before[k]),"N7 parity "+name+" "+k
            d.update(runs=runs,folds=folds,fold_supported=ok,fold_paths=fp,fold_scores=fs)
            kt,ts=oldkeep(d,"rivals_tree","preserve");kl0,ls=oldkeep(d,"rivals_logistic","q01")
            d.update(n7_tree_keep=kt,n7_logistic_keep=kl0)
            if source=="n7":
                for k,keep in (("rivals_tree_preserve_joint",kt),("rivals_logistic_q01_joint",kl0)):
                    assert np.array_equal(alert&keep|joint,before[k+"_alert"]) and np.array_equal(hit&keep|(joint&d["correct"]),before[k+"_hit"]),"N7 locked mask parity"
                for j,i in enumerate(idx):
                    tags=[]
                    if not kt[i] and d["y"][i]==0:tags.append("tree_removed_control")
                    if not kl0[i] and hit[i]:tags.append("logistic_lost_localization")
                    if not tags:continue
                    _,_,f0,df=mm[i]
                    audit.append({"tags":tags,"target":name,"kind":kind,"day":day,"row":int(i),"source_sha256":str(inp["tile_sha"][i]),"parent_path":pp[j].tolist(),"frequency_mhz":(f0+pp[j]*df/1e6).tolist(),"truth_path":inp["truth"][i].tolist() if d["y"][i] else None,"active":inp["active"][i].tolist() if d["y"][i] else None,"amp":float(inp["amp"][i]) if d["y"][i] else None,"on_contrast":a[j].tolist(),"off_contrast":off[j].tolist(),"rival_paths":rpaths[j].tolist(),"rival_scores":best[j].tolist(),"old_tree_score":float(ts[j]),"old_logistic_score":float(ls[j]),"runs":runs[j].tolist(),"folds":folds[j].tolist(),"fold_paths":fp[j].tolist(),"fold_scores":fs[j].tolist(),"supported":bool(ok[j])})
            ds.append((kind,d))
        rec={"target":name,"day":day,"split":split,"source":source,"raw_tile_hashes":sorted(hs)}
        emit("R40N8_ACQUIRED",rec);return rec,ds
    for source,gg in (("n6",nm["selected"]),("n7",m7["selected"])):
        for g in gg:
            rec,ds=acquire(g,source);records.append(rec)
            for kind,d in ds:
                dev.append((g["target"],source,g["split"],kind,d))
                if source=="n6" and g["split"]=="fit":train.append((g["target"],d))
                np.savez_compressed(O/("development_"+source+"_"+kind+"_"+g["target"]+".npz"),**d)
            save("development_records.json",records)
    assert len(all_old)==192
    removed=[r for r in audit if "tree_removed_control" in r["tags"]];lost=[r for r in audit if "logistic_lost_localization" in r["tags"]]
    assert len(removed)==2 and {r["target"] for r in removed}=={"WISE0350-5658","GJ1005"} and len(lost)==1 and lost[0]["target"]=="WISE0350-5658"
    save("paired_failure_audit.json",audit);emit("R40N8_AUDIT",{"removed_controls":len(removed),"lost_localizations":len(lost),"cases":audit,"sha256":c.sha((O/"paired_failure_audit.json").read_bytes())})
    def fv(d,feature):
        base=np.column_stack([d["temporal"],d["rival"]])
        return np.column_stack([base]+([d["runs"]] if feature in ("runs","combined") else [])+([d["folds"]] if feature in ("folds","combined") else []))
    # Use same support definition on all matched new arms, report exact frozen N7 separately.
    def use(d):return d["supported"]&d["fold_supported"]
    y=np.concatenate([d["parent_hit"][d["idx"]][use(d)] for _,d in train]);groups=np.concatenate([[name]*int(use(d).sum()) for name,d in train])
    assert y.any() and (~y).any() and len(set(groups))==6
    models={};serial={};cuts={};oof={}
    for feature in ("base","runs","folds","combined"):
        X=np.concatenate([fv(d,feature)[use(d)] for _,d in train])
        for kind in ("logistic","tree"):
            name=feature+"_"+kind;ss=np.zeros(len(y))
            for group in sorted(set(groups)):
                te=groups==group;sc,_=v.model_fit(X[~te],y[~te],kind);ss[te]=sc(X[te])
            oof[name]=float(roc_auc_score(y,ss))
            sc,s=v.model_fit(X,y,kind);models[name]=(sc,feature);serial[name]=s;pos=sc(X[y])
            cuts[name]={"preserve":float(np.nextafter(pos.min(),-np.inf)),"q01":float(np.quantile(pos,.01)),"q05":float(np.quantile(pos,.05))}
    lk={"protocol_sha256":c.sha(P.read_bytes()),"source_n7_lock":S["sources"]["n7"]["lock_sha256"],"manifest_sha256":mh,"development_records_sha256":c.sha((O/"development_records.json").read_bytes()),"audit_sha256":c.sha((O/"paired_failure_audit.json").read_bytes()),"train_rows":len(y),"train_correct":int(y.sum()),"feature_shapes":{"base":42,"runs":44,"folds":56},"models":serial,"cuts":cuts,"oof":oof}
    save("lock.json",lk);lh=c.sha((O/"lock.json").read_bytes());emit("R40N8_LOCK",{k:v0 for k,v0 in lk.items() if k!="models"})
    def measure(d):
        aa=d["parent_alert"];hh=d["parent_hit"];res=d["joint_rescue"];lr=d["learned05_alert"]&~aa
        y=d["y"];weak=d["weak"];co=d["correct"];jm=hh|(res&co)
        arms={"r40g":(aa,hh),"joint_min":(aa|res,jm),"learned05":(aa|lr,hh|(lr&co))}
        for name,keep in (("n7_tree_preserve",d["n7_tree_keep"]),("n7_logistic_q01",d["n7_logistic_keep"])):
            arms[name]=(aa&keep|res,hh&keep|(res&co))
        for name,(sc,feature) in models.items():
            vals=sc(fv(d,feature)) if len(d["idx"]) else np.empty(0)
            for cut,th in cuts[name].items():
                keep=np.ones(len(y),bool);keep[d["idx"]]=(~use(d))|(vals>=th)
                arms[name+"_"+cut]=(aa&keep|res,hh&keep|(res&co))
        out={}
        for name,(al,hit) in arms.items():
            out[name]={"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"control_alerts":int((al&(y==0)).sum()),"lost_parent":int((hh&~hit).sum()),"lost_joint":int((jm&~hit).sum()),"lost_joint_weak":int((jm&~hit&weak).sum()),"false_paths_on_positive":int((al&(y==1)&~hit).sum())}
            d[name+"_hit"]=hit;d[name+"_alert"]=al
        return out
    view=[];case_out=[]
    for name,source,split,kind,d in dev:
        row={"target":name,"source":source,"split":split,"kind":kind,"metrics":measure(d)};view.append(row)
        for case in audit:
            if source=="n7" and name==case["target"] and kind==case["kind"]:
                i=case["row"];case_out.append({"target":name,"kind":kind,"row":i,"tags":case["tags"],"decisions":{key:bool(d[key+"_alert"][i]) for key in row["metrics"]}})
    save("development_metrics.json",view);save("case_decisions.json",case_out);emit("R40N8_CASE_DECISIONS",case_out)
    if len(fresh)<C["min_groups"]:
        result={"status":"HOLD_INSUFFICIENT_FRESH_GROUPS","fresh_groups":len(fresh),"audit_count":len(audit),"lock_sha256":lh,"oof":oof};save("result.json",result);emit("R40N8_SUMMARY",result);return
    out=[];freshrec=[]
    for g in fresh:
        assert c.sha((O/"lock.json").read_bytes())==lh
        rec,ds=acquire(g,"fresh");freshrec.append(rec)
        for kind,d in ds:
            row={"target":g["target"],"kind":kind,"metrics":measure(d)}
            out.append(row);np.savez_compressed(O/("holdout_"+kind+"_"+g["target"]+".npz"),**d);emit("R40N8_GROUP",row)
        save("fresh_records.json",freshrec);save("partial.json",out)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def total(kind):
        rows=[r["metrics"] for r in out if r["kind"]==kind];assert len(rows)==len(fresh)
        return {a:{k:sum(v[a][k] for v in rows) for k in rows[0][a]} for a in rows[0]}
    e=total("eval");sh=total("shared");nc=128*len(fresh);base=e["learned05"]["control_alerts"]
    gates={k:e[k]["lost_parent"]==e[k]["lost_joint"]==sh[k]["lost_parent"]==sh[k]["lost_joint"]==0 and e[k]["control_alerts"]*640<=35*nc and e[k]["control_alerts"]*8<=7*base for k in e}
    result={"status":"FRESH_DEVELOPMENT_HOLDOUT_ONLY","fresh_groups":len(fresh),"controls":nc,"audit_count":len(audit),"lock_sha256":lh,"oof":oof,"eval":e,"shared":sh,"gates":gates,"rows":out,"seconds":time.perf_counter()-start,"production_adoption":False};save("result.json",result);emit("R40N8_SUMMARY",{k:v for k,v in result.items() if k!="rows"})
if __name__=="__main__":main()
