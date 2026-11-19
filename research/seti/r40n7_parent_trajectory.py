#!/usr/bin/env python3
"""R40N7: protect parent detections with independent raw trajectory witnesses."""
import sys
sys.argv[1:]=["controls"]
import json,pathlib,time,hashlib,re,io
import numpy as np
from sklearn.metrics import roc_auc_score
import r40n6_same_band as v
n=v.n;c=v.c;b=v.b;l=v.l;f=v.f;h=v.h;t=v.t
qual=v.qual;scout=v.scout;eng=v.eng;index_query=v.index_query
P=pathlib.Path("research/seti/r40n7_parent_trajectory_protocol.json");S=json.loads(P.read_text());C=S["selection"]
O=pathlib.Path("runtime/r40n7");O.mkdir(parents=True,exist_ok=True);n.O=O;n.k.O=O;n.l.O=O;eng.OUT=O
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
                        emit("R40N7_SELECTED",{"telescope":tel,"target":target,"day":day,"count":len(manifest["selected"])});break
                    except Exception as e:a["errors"].append({"day":group[0]["day"],"error":type(e).__name__+": "+str(e)[:200]})
            except Exception as e:a["errors"].append({"error":type(e).__name__+": "+str(e)[:200]})
            save("manifest.json",manifest)
    assert len({g["day"] for g in manifest["selected"]})==len({g["target"] for g in manifest["selected"]})==len(manifest["selected"])
    if len(manifest["selected"])>=C["min_groups"]:
        for g in manifest["selected"]:g["split"]="holdout"
    save("manifest.json",manifest);emit("R40N7_MANIFEST",{"targets":[[g["telescope"],g["target"],g["day"]] for g in manifest["selected"]],"queries":len(manifest["attempts"]),"candidate_counts":manifest["candidate_counts"],"sha256":c.sha((O/"manifest.json").read_bytes())})
    return manifest,fr


def witness(X,D,p):
    # Only measured, legal bins; no wrapping or artificial edge values as evidence.
    assert X.shape==D.shape and X.shape[1:]==(16,128) and p.shape==(len(X),16)
    rows=np.arange(16)[None,:];ii=np.arange(len(X))[:,None]
    valid=(p>=5)&(p<123);count=valid.sum(1);supported=count>=12
    pp=np.clip(p,5,122);U=n.normalized_contrast(X);V=n.normalized_contrast(D)
    a=U[ii,rows,pp];d=V[ii,rows,pp]
    a=np.where(valid,a,0.);d=np.where(valid,d,0.)
    def mean(z,mask):
        return (z*mask).sum(1)/np.maximum(1,mask.sum(1))
    def stats(z):
        mu=mean(z,valid)
        block=np.stack([mean(z[:,j:j+4],valid[:,j:j+4]) for j in range(0,16,4)],1)
        return [mu,np.sqrt(mean((z-mu[:,None])**2,valid)),mean(z>0,valid),mean(z>1,valid),block.min(1),np.sort(block,axis=1)[:,1],block.max(1),mean(z[:,::2],valid[:,::2]),mean(z[:,1::2],valid[:,1::2])]
    # Recalculate raw ON and OFF separately on candidate and alternatives.
    best=np.full((len(X),2),-np.inf);rpaths=np.zeros((len(X),2,16),int);ix=np.arange(len(X))
    for delta in b.D:
        low=5-int(delta.min());high=123-int(delta.max())
        if high<=low:continue
        starts=np.arange(low,high);paths=starts[:,None]+delta
        distinct=(abs(paths[None,:,:]-p[:,None,:])>2).sum(2)>=5
        for k,Z in enumerate((U,V)):
            sc=sum(Z[:,j,low+dd:high+dd] for j,dd in enumerate(delta))/16
            sc=np.where(distinct,sc,-np.inf);arg=sc.argmax(1);vv=sc[ix,arg];m=vv>best[:,k]
            best[m,k]=vv[m];rpaths[m,k]=paths[arg[m]]
    assert np.isfinite(best).all()
    for k,Z in enumerate((U,V)):
        path=rpaths[:,k];assert ((abs(path-p)>2).sum(1)>=5).all()
        assert np.allclose(Z[ii,rows,path].mean(1),best[:,k],atol=1e-12,rtol=1e-12)
    st=np.column_stack(stats(a)+stats(d)+[count/16,mean(a-d,valid),mean(np.minimum(a,d),valid),mean(abs(a-d),valid),mean(a>0,valid)-mean(d>0,valid),np.mean(abs(np.diff(p,axis=1)),axis=1),np.ptp(p,axis=1),np.mean(abs(np.diff(p,n=2,axis=1)),axis=1)])
    # shape uses legal measured neighbour bins, not normalized contrast edges.
    st=np.column_stack([st,n.shape(X,p),n.shape(D,p)])
    rv=np.column_stack([best,mean(a,valid)-best[:,0],mean(d,valid)-best[:,1],best[:,0]-best[:,1],np.mean(abs(rpaths[:,0]-rpaths[:,1]),axis=1)])
    assert np.isfinite(st).all() and np.isfinite(rv).all()
    return st,rv,supported,a,d,rpaths,best
def selftest():
    rng=np.random.default_rng(701);X=rng.normal(0,.1,(2,16,128));D=X.copy();p=np.tile(np.arange(40,56),(2,1));X[0,np.arange(16),p[0]]+=5
    D[0,:,90]+=5; p[1,:8]=0
    st,rv,ok,a,d,rp,bs=witness(X,D,p)
    assert ok.tolist()==[True,False] and a[0].mean()>d[0].mean()+1 and np.count_nonzero(rp[0,1]==90)>=12
    emit("R40N7_SELFTEST",{"independent_off_rival":"PASS","edge_abstention":"PASS"})
def main():
    start=time.perf_counter();selftest();n.selftest()
    zz,nn=archive("n6");assert c.sha(zz.read(nn["lock.json"]))==S["sources"]["n6"]["lock_sha256"] and c.sha(zz.read(nn["manifest.json"]))==S["sources"]["n6"]["manifest_sha256"]
    nm=n.read(zz,nn,"manifest.json");fitrecs=n.read(zz,nn,"fit_records.json");holdrecs=n.read(zz,nn,"holdout_records.json");priorrecs={r["target"]:r for r in fitrecs+holdrecs}
    assert len(nm["selected"])==len(priorrecs)==8
    az,an=archive("n1");bz,bn=archive("n2");uz,un=archive("n5a")
    assert c.sha(uz.read(un["lock.json"]))==S["sources"]["n5a"]["lock_sha256"]
    prior=n.read(az,an,"manifest.json")["selected"]+n.read(bz,bn,"manifest.json")["selected"]+nm["selected"]
    manifest,fr=select({"selected":prior});mh=c.sha((O/"manifest.json").read_bytes())
    fresh=manifest["selected"];assert not {int(g["day"]) for g in fresh}&{int(g["day"]) for g in prior}
    L=h.readlocks();b.taylor.init_tree()
    lz,ln=n.archive("r40l");raw=lz.read(ln["lock.json"]);assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];old=json.loads(raw)
    kz,kn=n.archive("r40k");kr=n.read(kz,kn,"result.json");kl=n.read(kz,kn,"lock.json")
    seen=set(fr["excluded_raw_hashes"])
    for row in kr["rows"]+kl["calibration"]+fitrecs+holdrecs+n.read(uz,un,"bank_manifest.json")+n.read(uz,un,"unused_bank_manifest.json")["records"]:seen.update(row["raw_tile_hashes"])
    # Include all original calibration and historical replay raw hashes.
    for z0,n0 in ((az,an),(bz,bn)):
        for key in n0:
            if key.endswith(".npz"):
                with np.load(io.BytesIO(z0.read(n0[key])),allow_pickle=False) as vv:
                    if "source_tile_sha" in vv:seen.update(vv["source_tile_sha"].tolist())
    # R40N replay provenance, separately authenticated by original protocol.
    rrz,rrn=c.artifact(v.S["sources"]["replay"],64)
    for key in rrn:
        if key.endswith(".npz"):
            with np.load(io.BytesIO(rrz.read(rrn[key])),allow_pickle=False) as vv:seen.update(vv["source_tile_sha"].tolist())
    train=[];audit=[];dev=[];records=[];all_old=set()
    def acquire(g,previous=False):
        name=g["target"];day=int(g["day"]);split=g["split"]
        assert c.sha((O/"manifest.json").read_bytes())==mh
        n.set_config(v.C["seed"] if previous else C["seed"],C["windows_mhz"],C["half_hz"],C["max_mib_per_scan"])
        rr,R,ch,hs,meta,cm,pre=n.gather(g,("audit_" if previous else "new_")+name);assert len(hs)==16
        if previous:
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
            if previous:
                before=f.load(zz,nn,split+"_"+kind+"_"+name+".npz")
                for k,value in d.items():assert np.array_equal(value,before[k]),"N6 parity "+name+" "+k
            idx=np.flatnonzero(alert);pp=path[idx].astype(int)
            if len(idx):
                X=n.ref.norm(np.stack([rows[i][0] for i in idx]));D=n.ref.norm(refs[idx])
                st,rv,supported,a,off,rpaths,best=witness(X,D,pp)
            else:
                st=np.empty((0,36));rv=np.empty((0,6));supported=np.empty(0,bool);a=np.empty((0,16));off=a.copy();rpaths=np.empty((0,2,16),int);best=np.empty((0,2))
            # For parent-alert subset only; no decision for non-alert rows.
            d.update(idx=idx,temporal=st,rival=rv,supported=supported,selected_on=a,selected_off=off,rival_paths=rpaths,rival_scores=best)
            if previous and split=="holdout" and kind=="eval":
                for j,i in enumerate(idx):
                    if d["y"][i]!=0:continue
                    _,_,f0,df=mm[i];audit.append({"target":name,"day":day,"row":int(i),"source_sha256":str(inp["tile_sha"][i]),"parent_path":pp[j].tolist(),"frequency_mhz":(f0+pp[j]*df/1e6).tolist(),"on_contrast":a[j].tolist(),"off_contrast":off[j].tolist(),"rival_on_path":rpaths[j,0].tolist(),"rival_off_path":rpaths[j,1].tolist(),"rival_scores":best[j].tolist(),"temporal":st[j].tolist(),"rival":rv[j].tolist(),"supported":bool(supported[j])})
            ds.append((kind,d))
        rec={"target":name,"day":day,"split":split,"previous":previous,"raw_tile_hashes":sorted(hs)}
        emit("R40N7_ACQUIRED",rec);return rec,ds
    for g in nm["selected"]:
        rec,ds=acquire(g,True);records.append(rec)
        for kind,d in ds:
            dev.append((g["target"],g["split"],kind,d))
            if g["split"]=="fit":train.append((g["target"],d))
            np.savez_compressed(O/("development_"+g["split"]+"_"+kind+"_"+g["target"]+".npz"),**d)
        save("development_records.json",records)
    assert len(all_old)==128 and len(audit)==15
    save("parent_false_alert_audit.json",audit)
    emit("R40N7_AUDIT",{"controls":15,"by_target":{name:sum(x["target"]==name for x in audit) for name in sorted(set(x["target"] for x in audit))},"supported":sum(x["supported"] for x in audit),"sha256":c.sha((O/"parent_false_alert_audit.json").read_bytes())})
    def use(d):return d["supported"]
    def fv(d,feature):return d["temporal"] if feature=="temporal" else np.column_stack([d["temporal"],d["rival"]])
    y=np.concatenate([d["parent_hit"][d["idx"]][use(d)] for _,d in train]);groups=np.concatenate([[name]*int(use(d).sum()) for name,d in train])
    assert y.any() and (~y).any() and len(set(groups))==6
    models={};serial={};cuts={};oof={}
    for feature in ("temporal","rivals"):
        X=np.concatenate([fv(d,feature)[use(d)] for _,d in train])
        for kind in ("logistic","tree"):
            name=feature+"_"+kind;ss=np.zeros(len(y))
            for group in sorted(set(groups)):
                te=groups==group;sc,_=v.model_fit(X[~te],y[~te],kind);ss[te]=sc(X[te])
            oof[name]=float(roc_auc_score(y,ss))
            sc,s=v.model_fit(X,y,kind);models[name]=(sc,feature);serial[name]=s;pos=sc(X[y])
            cuts[name]={"preserve":float(np.nextafter(pos.min(),-np.inf)),"q01":float(np.quantile(pos,.01)),"q05":float(np.quantile(pos,.05))}
    lk={"protocol_sha256":c.sha(P.read_bytes()),"source_n6_lock":S["sources"]["n6"]["lock_sha256"],"manifest_sha256":mh,"development_records_sha256":c.sha((O/"development_records.json").read_bytes()),"audit_sha256":c.sha((O/"parent_false_alert_audit.json").read_bytes()),"train_rows":len(y),"train_correct":int(y.sum()),"models":serial,"cuts":cuts,"oof":oof}
    save("lock.json",lk);lh=c.sha((O/"lock.json").read_bytes());emit("R40N7_LOCK",{k:v0 for k,v0 in lk.items() if k!="models"})
    def measure(d):
        aa=d["parent_alert"];hh=d["parent_hit"];res=d["joint_rescue"];lr=d["learned05_alert"]&~aa
        y=d["y"];weak=d["weak"];co=d["correct"];jm=hh|(res&co)
        arms={"r40g":(aa,hh),"joint_min":(aa|res,jm),"learned05":(aa|lr,hh|(lr&co))}
        for name,(sc,feature) in models.items():
            vals=sc(fv(d,feature)) if len(d["idx"]) else np.empty(0)
            for cut,th in cuts[name].items():
                keep=np.ones(len(y),bool);keep[d["idx"]]=(~d["supported"])|(vals>=th)
                for rescue,rr in (("joint",res),("learned",lr)):
                    arms[name+"_"+cut+"_"+rescue]=(aa&keep|rr,hh&keep|(rr&co))
        out={}
        for name,(al,hit) in arms.items():
            out[name]={"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"control_alerts":int((al&(y==0)).sum()),"lost_parent":int((hh&~hit).sum()),"lost_joint":int((jm&~hit).sum()),"lost_joint_weak":int((jm&~hit&weak).sum()),"false_paths_on_positive":int((al&(y==1)&~hit).sum())}
            d[name+"_hit"]=hit;d[name+"_alert"]=al
        return out
    view=[]
    for name,split,kind,d in dev:
        row={"target":name,"split":split,"kind":kind,"metrics":measure(d)};view.append(row)
    save("development_metrics.json",view)
    if len(fresh)<C["min_groups"]:
        result={"status":"HOLD_INSUFFICIENT_FRESH_GROUPS","fresh_groups":len(fresh),"audit_count":len(audit),"lock_sha256":lh,"oof":oof}
        save("result.json",result);emit("R40N7_SUMMARY",result);return
    out=[];freshrec=[]
    for g in fresh:
        assert c.sha((O/"lock.json").read_bytes())==lh
        rec,ds=acquire(g);freshrec.append(rec)
        for kind,d in ds:
            row={"target":g["target"],"kind":kind,"metrics":measure(d)}
            out.append(row);np.savez_compressed(O/("holdout_"+kind+"_"+g["target"]+".npz"),**d);emit("R40N7_GROUP",row)
        save("fresh_records.json",freshrec);save("partial.json",out)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def total(kind):
        rows=[r["metrics"] for r in out if r["kind"]==kind];assert len(rows)==len(fresh)
        return {a:{k:sum(v[a][k] for v in rows) for k in rows[0][a]} for a in rows[0]}
    e=total("eval");sh=total("shared");nc=128*len(fresh);base=e["learned05"]["control_alerts"]
    gates={k:e[k]["lost_parent"]==e[k]["lost_joint"]==sh[k]["lost_parent"]==sh[k]["lost_joint"]==0 and e[k]["control_alerts"]*640<=35*nc and (e[k]["control_alerts"]*8<=7*base) for k in e}
    result={"status":"FRESH_DEVELOPMENT_HOLDOUT_ONLY","fresh_groups":len(fresh),"controls":nc,"audit_count":len(audit),"lock_sha256":lh,"oof":oof,"eval":e,"shared":sh,"gates":gates,"rows":out,"seconds":time.perf_counter()-start,"production_adoption":False}
    save("result.json",result);emit("R40N7_SUMMARY",{k:v for k,v in result.items() if k!="rows"})
if __name__=="__main__":main()
