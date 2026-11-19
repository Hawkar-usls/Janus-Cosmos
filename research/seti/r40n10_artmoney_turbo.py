#!/usr/bin/env python3
"""R40N10 exact scoring accelerator: staged checkpoints, no live repo mutation."""
import sys
ARGS=sys.argv[1:]
sys.argv[1:]=["controls"]
import json,pathlib,time,hashlib,re,io,shutil
import numpy as np
from sklearn.metrics import roc_auc_score
import r40n9_edge_invariance as old9
import r40n9_halo as halo
import r40n10_turbo as turbo
w=old9.w;v=w.v;n=w.n;c=w.c;b=w.b;l=w.l;f=w.f;h=w.h;t=w.t
qual=w.qual;scout=w.scout;eng=w.eng;index_query=w.index_query
P=pathlib.Path("research/seti/r40n10_artmoney_turbo_protocol.json");S=json.loads(P.read_text());C=S["selection"]
O=pathlib.Path("runtime/r40n10");O.mkdir(parents=True,exist_ok=True);CACHE=O/"cache";CACHE.mkdir(exist_ok=True)
PROG=O/"progress";PROG.mkdir(exist_ok=True)
n.O=O;n.k.O=O;n.l.O=O;eng.OUT=O
def write(path,d):path.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def save(name,d):write(O/name,d)
def load(name):return json.loads((O/name).read_text())
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def archive(k):return c.artifact(S["sources"][k],128)
def read_npz(path):
    with np.load(path,allow_pickle=False) as d:return {k:d[k] for k in d.files}
def checkpoint(name,d):
    write(PROG/(name+".json"),d);emit("R40N10_PROGRESS",d)
def checked_archive(key):
    z,nn=archive(key)
    if "lock_sha256" in S["sources"][key]:assert c.sha(z.read(nn["lock.json"]))==S["sources"][key]["lock_sha256"]
    if "manifest_sha256" in S["sources"][key]:assert c.sha(z.read(nn["manifest.json"]))==S["sources"][key]["manifest_sha256"]
    return z,nn
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
                        emit("R40N10_SELECTED",{"telescope":tel,"target":target,"day":day,"count":len(manifest["selected"])});break
                    except Exception as e:a["errors"].append({"day":group[0]["day"],"error":type(e).__name__+": "+str(e)[:200]})
            except Exception as e:a["errors"].append({"error":type(e).__name__+": "+str(e)[:200]})
            save("manifest.json",manifest)
    assert len({g["day"] for g in manifest["selected"]})==len({g["target"] for g in manifest["selected"]})==len(manifest["selected"])
    if len(manifest["selected"])>=C["min_groups"]:
        for g in manifest["selected"]:g["split"]="holdout"
    save("manifest.json",manifest);emit("R40N10_MANIFEST",{"targets":[[g["telescope"],g["target"],g["day"]] for g in manifest["selected"]],"queries":len(manifest["attempts"]),"candidate_counts":manifest["candidate_counts"],"sha256":c.sha((O/"manifest.json").read_bytes())})
    return manifest,fr




def prepare():
    start=time.perf_counter()
    paths=halo.candidates(b.D);bench=turbo.selftest(n.ref.norm,paths)
    assert bench["oracle_calls"]==28
    save("synthetic_parity.json",bench);emit("R40N10_SYNTHETIC",bench)
    # Authenticated N9 provides immutable earlier parent baseline arrays.
    arc={k:checked_archive(k) for k in ("n6","n7","n8","n9","n1","n2","n5a")}
    manifests={k:n.read(*arc[k],"manifest.json") for k in ("n6","n7","n8","n9","n1","n2")}
    z9,n9=arc["n9"]
    for k in ("n7","n8","n9"):
        z,nn=arc[k];(CACHE/(k+"_lock.json")).write_bytes(z.read(nn["lock.json"]))
    dev=[]
    for g in manifests["n6"]["selected"]:
        if g["split"]=="fit":dev.append({"source":"n6","group":g,"fit":True})
    dev.extend({"source":"n7","group":g,"fit":False} for g in manifests["n7"]["selected"] if g["target"]=="WISE0350-5658")
    dev.extend({"source":"n9","group":g,"fit":False} for g in manifests["n9"]["selected"])
    assert len(dev)==11 and sum(x["fit"] for x in dev)==6
    recs=n.read(z9,n9,"development_records.json")+n.read(z9,n9,"fresh_records.json")
    priorrecs={x["target"]:x for x in recs};assert len(priorrecs)==20
    for x in dev:
        name=x["group"]["target"];source=x["source"];x["raw_tile_hashes"]=priorrecs[name]["raw_tile_hashes"]
        for kind in ("eval","shared"):
            fn=("holdout_" if source=="n9" else "development_"+source+"_")+kind+"_"+name+".npz"
            (CACHE/fn).write_bytes(z9.read(n9[fn]))
    prior=sum([manifests[k]["selected"] for k in ("n1","n2","n6","n7","n8","n9")],[])
    manifest,fr=select({"selected":prior});mh=c.sha((O/"manifest.json").read_bytes())
    fresh=manifest["selected"]
    assert not {int(g["day"]) for g in fresh}&{int(g["day"]) for g in prior}
    # Freeze provenance exclusions before new powers.
    h.O=CACHE;L=h.readlocks();write(CACHE/"parent_locks.json",L)
    lz,ln=n.archive("r40l");raw=lz.read(ln["lock.json"]);assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];(CACHE/"joint_lock.json").write_bytes(raw)
    kz,kn=n.archive("r40k");kr=n.read(kz,kn,"result.json");kl=n.read(kz,kn,"lock.json")
    uz,un=arc["n5a"]
    seen=set(fr["excluded_raw_hashes"])
    for row in kr["rows"]+kl["calibration"]+recs+n.read(uz,un,"bank_manifest.json")+n.read(uz,un,"unused_bank_manifest.json")["records"]:seen.update(row["raw_tile_hashes"])
    for key in ("n1","n2"):
        z0,n0=arc[key]
        for fn in n0:
            if fn.endswith(".npz"):
                with np.load(io.BytesIO(z0.read(n0[fn])),allow_pickle=False) as vv:
                    if "source_tile_sha" in vv:seen.update(vv["source_tile_sha"].tolist())
    rrz,rrn=c.artifact(v.S["sources"]["replay"],64)
    for fn in rrn:
        if fn.endswith(".npz"):
            with np.load(io.BytesIO(rrz.read(rrn[fn])),allow_pickle=False) as vv:seen.update(vv["source_tile_sha"].tolist())
    plan={"manifest_sha256":mh,"protocol_sha256":c.sha(P.read_bytes()),"development":dev,"fresh":fresh,"excluded_raw_hashes":sorted(seen),"fresh_allowed":len(fresh)>=C["min_groups"]}
    save("plan.json",plan)
    checkpoint("00-prepare",{"stage":"PREPARE","completed_groups":0,"total_groups":len(dev)+(len(fresh) if plan["fresh_allowed"] else 0),"development_groups":len(dev),"fresh_groups":len(fresh),"manifest_sha256":mh,"seconds":time.perf_counter()-start,"synthetic_parity":bench})
def verified_plan():
    plan=load("plan.json");assert plan["protocol_sha256"]==c.sha(P.read_bytes()) and plan["manifest_sha256"]==c.sha((O/"manifest.json").read_bytes())
    return plan
def lock_cache(key):
    path=CACHE/(key+"_lock.json");assert c.sha(path.read_bytes())==S["sources"][key]["lock_sha256"];return json.loads(path.read_text())
def score(X,model,kind):
    return old9.frozen_score(X,model,kind) if len(X) else np.empty(0)
def field(d,key):
    base=np.column_stack([d["temporal"],d["rival"]]);ha=np.column_stack([d["halo"],d["artmoney"]])
    if key=="base":return base
    if key=="base_halo":return np.column_stack([base,d["halo"]])
    if key=="base_halo_art":return np.column_stack([base,ha])
    if key=="halo_art":return ha
    if key=="repeat":return d["repeat"]
    if key=="halo_art_repeat":return np.column_stack([ha,d["repeat"]])
    if key=="base_repeat":return np.column_stack([base,d["repeat"]])
    raise ValueError(key)
def fixed_mask(d,feature,kind,cut,L):
    vals=score(field(d,feature),L["models"][feature+"_"+kind],kind)
    keep=np.ones(len(d["y"]),bool);keep[d["idx"]]=(~d["halo_supported"])|(vals>=L["cuts"][feature+"_"+kind][cut])
    return keep
def arms9(d,L):
    return {feature+"_"+kind+"_"+cut:fixed_mask(d,feature,kind,cut,L) for feature in ("base","base_halo","base_halo_art","halo_art") for kind in ("logistic","tree") for cut in ("preserve","q01","q05")}
def parent_fresh(rows,refs,mm,pre,L,old,L7,L8):
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
        st,rv,supported,a,off,rpaths,best=w.witness(X,D,pp);runs,folds,ok,fp,fs=old9.extra(X,D,pp)
    else:
        st=np.empty((0,36));rv=np.empty((0,6));supported=np.empty(0,bool);a=np.empty((0,16));off=a.copy();rpaths=np.empty((0,2,16),int);best=np.empty((0,2));runs=np.empty((0,44));folds=np.empty((0,56));ok=np.empty(0,bool);fp=np.empty((0,2,4,16),int);fs=np.empty((0,2,4,4))
    d.update(idx=idx,temporal=st,rival=rv,supported=supported,selected_on=a,selected_off=off,rival_paths=rpaths,rival_scores=best,runs=runs,folds=folds,fold_supported=ok,fold_paths=fp,fold_scores=fs)
    base=np.column_stack([st,rv])
    for name,key,kind,cut in (("n7_tree_keep","rivals_tree","tree","preserve"),("n7_logistic_keep","rivals_logistic","logistic","q01")):
        sc=score(base,L7["models"][key],kind);keep=np.ones(len(alert),bool);keep[idx]=(~supported)|(sc>=L7["cuts"][key][cut]);d[name]=keep
    sc8=score(np.column_stack([base,runs,folds]),L8["models"]["combined_tree"],"tree")
    kk=np.ones(len(alert),bool);kk[idx]=(~(supported&ok))|(sc8>=L8["cuts"]["combined_tree"]["q01"])
    sr=score(np.column_stack([base,runs]),L8["models"]["runs_tree"],"tree")
    kr0=np.ones(len(alert),bool);kr0[idx]=(~(supported&ok))|(sr>=L8["cuts"]["runs_tree"]["q05"])
    d.update(n8_combined_keep=kk,n8_strict_runs_keep=kr0)
    return d
def verify_rows(d,rows):
    vals={"y":np.array([r[1] for r in rows]),"truth":np.stack([r[7] for r in rows]),"active":np.stack([r[8] for r in rows]),"amp":np.array([r[4] for r in rows]),"source_tile_sha":np.array([r[6] for r in rows])}
    for k,a in vals.items():assert np.array_equal(d[k],a),"Prior injection/row parity "+k
    assert np.array_equal(d["idx"],np.flatnonzero(d["parent_alert"]))
def measure(d,new,L9):
    aa=d["parent_alert"];hh=d["parent_hit"];res=d["joint_rescue"];lr=d["learned05_alert"]&~aa;y=d["y"];weak=d["weak"];co=d["correct"];jm=hh|(res&co)
    arms={"r40g":(aa,hh),"joint_min":(aa|res,jm),"learned05":(aa|lr,hh|(lr&co))}
    for name,keep in (("n7_tree_preserve",d["n7_tree_keep"]),("n8_combined_q01",d["n8_combined_keep"])):
        arms[name]=(aa&keep|res,hh&keep|(res&co))
    masks=arms9(d,L9)
    for name,keep in masks.items():arms["n9_"+name]=(aa&keep|res,hh&keep|(res&co))
    use=d["halo_supported"]&d["repeat_supported"]
    stress=d["repeat_stress_supported"]&use
    stress_d={**d,"repeat":d["repeat_stress"]}
    stress_m={}
    for feature in ("halo_art","repeat","halo_art_repeat","base_repeat"):
        for kind in ("logistic","tree"):
            name=feature+"_"+kind
            vals=score(field(d,feature),new["models"][name],kind);sv=score(field(stress_d,feature),new["models"][name],kind)
            for cut,th in new["cuts"][name].items():
                key=name+"_"+cut;keep=np.ones(len(y),bool);keep[d["idx"]]=(~use)|(vals>=th)
                arms[key]=(aa&keep|res,hh&keep|(res&co))
                stress_m[key]=int((hh[d["idx"]]&stress&(sv<th)).sum())
    out={}
    for name,(al,hit) in arms.items():
        out[name]={"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"control_alerts":int((al&(y==0)).sum()),"lost_parent":int((hh&~hit).sum()),"lost_joint":int((jm&~hit).sum()),"lost_joint_weak":int((jm&~hit&weak).sum()),"false_paths_on_positive":int((al&(y==1)&~hit).sum()),"repeat_stress_cases":int((hh[d["idx"]]&stress).sum()),"repeat_stress_lost":stress_m.get(name,int(((hh&~hit)[d["idx"]]&stress).sum()))}
        d[name+"_hit"]=hit;d[name+"_alert"]=al
    return out
def group(phase,i):
    start=time.perf_counter();plan=verified_plan();dev=phase=="dev"
    if dev:
        x=plan["development"][i];g=x["group"];source=x["source"]
    else:
        if not plan["fresh_allowed"] or i>=len(plan["fresh"]):
            checkpoint("fresh-"+str(i+1).zfill(2),{"stage":"FRESH_SKIPPED","slot":i+1,"reason":"not qualified"});return
        g=plan["fresh"][i];source="fresh"
        freeze=load("freeze_receipt.json");assert c.sha((O/"lock.json").read_bytes())==freeze["lock_sha256"] and freeze["plan_sha256"]==c.sha((O/"plan.json").read_bytes())
    name=g["target"];seed={"n6":v.C["seed"],"n7":w.C["seed"],"n9":old9.C["seed"],"fresh":C["seed"]}[source]
    n.set_config(seed,C["windows_mhz"],C["half_hz"],C["max_mib_per_scan"])
    paths=halo.candidates(b.D)
    # Warm both numba kernels before monitored real parity timing.
    turbo.bank_scores(np.zeros((2,16,176)),paths-32);turbo.legal_mask(paths,np.full((2,16),100))
    accelerator=turbo.Turbo(paths,6);halo.context_scores=accelerator
    rr,R,ch,hs,meta,cm,pre=n.gather(g,source+"_"+name);assert len(hs)==16
    if dev:assert sorted(hs)==x["raw_tile_hashes"]
    else:
        seen=set(plan["excluded_raw_hashes"])
        for fn in PROG.glob("fresh-*.json"):
            p0=json.loads(fn.read_text());seen.update(p0.get("raw_tile_hashes",[]))
        assert not hs&seen
    L9=lock_cache("n9")
    if not dev:
        L=json.loads((CACHE/"parent_locks.json").read_text());raw=(CACHE/"joint_lock.json").read_bytes();assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];old=json.loads(raw);L7=lock_cache("n7");L8=lock_cache("n8")
        b.taylor.init_tree();new=load("lock.json")
    cases=[];metrics=[];pair_count=0
    for kind,rows,refs,mm in (("eval",rr,R,meta),("shared",[q[0] for q in ch],np.stack([q[1] for q in ch]),cm)):
        if dev:
            fn=("holdout_" if source=="n9" else "development_"+source+"_")+kind+"_"+name+".npz"
            d=read_npz(CACHE/fn);verify_rows(d,rows)
        else:d=parent_fresh(rows,refs,mm,pre,L,old,L7,L8)
        idx=d["idx"];pp=d["parent_path"][idx].astype(int)
        ha=halo.extract(rows,refs,mm,pre,idx,pp,n.ref,c.sha,paths)
        if dev:
            for key,value in ha.items():assert np.array_equal(d[key],value),"R40N9 halo parity: "+name+" "+key
            pair_count+=len(idx)
        d.update(ha);d.update(turbo.repeat_witness(rows,mm,pre,idx,pp,n.ref,kind))
        if dev and source=="n9":
            masks=arms9(d,L9);alert=d["parent_alert"];hit=d["parent_hit"];res=d["joint_rescue"];co=d["correct"]
            for key,keep in masks.items():
                assert np.array_equal(alert&keep|res,d[key+"_alert"]) and np.array_equal(hit&keep|(res&co),d[key+"_hit"]),"N9 exact frozen decision parity "+key
            for j,k in enumerate(idx):
                tags=[]
                if hit[k] and not d["n8_combined_keep"][k]:tags.append("n8_lost")
                if hit[k] and not masks["base_halo_art_tree_preserve"][k]:tags.append("n9_art_tree_lost")
                if hit[k] and not masks["halo_art_tree_q05"][k]:tags.append("n9_halo_art_strict_lost")
                if d["y"][k]==0 and not masks["halo_art_tree_preserve"][k]:tags.append("n9_control_removed")
                if tags:cases.append({"target":name,"kind":kind,"row":int(k),"tags":tags,"repeat_supported":bool(d["repeat_supported"][j]),"repeat":d["repeat"][j].tolist(),"parent_path":pp[j].tolist(),"halo_traces":d["halo_traces"][j].tolist(),"source_tile_sha256":str(d["source_tile_sha"][k])})
        if dev and name=="WISE0350-5658" and kind=="eval":
            j=int(np.flatnonzero(idx==129)[0]);cases.append({"target":name,"kind":kind,"row":129,"tags":["wise_edge"],"repeat_supported":bool(d["repeat_supported"][j]),"repeat":d["repeat"][j].tolist(),"halo_traces":d["halo_traces"][j].tolist(),"halo_exact":True,"repeat_source_sha":d["repeat_source_sha"][j].tolist()})
        row={"target":name,"kind":kind,"parent_alerts":len(idx),"halo_supported":int(d["halo_supported"].sum()),"repeat_supported":int(d["repeat_supported"].sum()),"both_supported":int((d["repeat_supported"]&d["halo_supported"]).sum()),"repeat_stress_parent_hits":int((d["repeat_stress_supported"]&d["halo_supported"]&d["parent_hit"][idx]).sum())}
        if not dev:row["metrics"]=measure(d,new,L9)
        metrics.append(row);np.savez_compressed(O/(phase+"_"+str(i+1).zfill(2)+"_"+kind+".npz"),**d)
    # Do not silently accept a group without the full reference check.
    assert accelerator.checked==6 or accelerator.calls==0
    if not dev:assert c.sha((O/"lock.json").read_bytes())==freeze["lock_sha256"]
    record={"stage":phase.upper(),"slot":i+1,"target":name,"source":source,"day":int(g["day"]),"fit":bool(dev and x["fit"]),"completed_groups":i+1 if dev else len(plan["development"])+i+1,"total_groups":len(plan["development"])+(len(plan["fresh"]) if plan["fresh_allowed"] else 0),"raw_tile_hashes":sorted(hs),"exact_replayed_parent_rows":pair_count,"turbo":accelerator.report(),"cases":cases,"rows":metrics,"seconds":time.perf_counter()-start}
    checkpoint(phase+"-"+str(i+1).zfill(2),record)
def fit():
    start=time.perf_counter();plan=verified_plan();train=[];records=[]
    for i,x in enumerate(plan["development"]):
        rec=json.loads((PROG/("dev-"+str(i+1).zfill(2)+".json")).read_text());assert rec["target"]==x["group"]["target"];records.append(rec)
        if x["fit"]:
            for kind in ("eval","shared"):train.append((rec["target"],read_npz(O/("dev_"+str(i+1).zfill(2)+"_"+kind+".npz"))))
    hashes=[v0 for rec in records for v0 in rec["raw_tile_hashes"]];assert len(hashes)==len(set(hashes))==176
    def use(d):return d["halo_supported"]&d["repeat_supported"]
    y=np.concatenate([d["parent_hit"][d["idx"]][use(d)] for _,d in train]);groups=np.concatenate([[name]*int(use(d).sum()) for name,d in train])
    assert len(set(groups))==6 and y.any() and (~y).any()
    models={};cuts={};oof={}
    for feature in ("halo_art","repeat","halo_art_repeat","base_repeat"):
        X=np.concatenate([field(d,feature)[use(d)] for _,d in train])
        for kind in ("logistic","tree"):
            name=feature+"_"+kind;ss=np.zeros(len(y))
            for group0 in sorted(set(groups)):
                te=groups==group0;sc,_=v.model_fit(X[~te],y[~te],kind);ss[te]=sc(X[te])
            oof[name]=float(roc_auc_score(y,ss))
            sc,serial=v.model_fit(X,y,kind);models[name]=serial;pos=sc(X[y])
            cuts[name]={"preserve":float(np.nextafter(pos.min(),-np.inf)),"q01":float(np.quantile(pos,.01)),"q05":float(np.quantile(pos,.05))}
    lock={"source_n9_lock":S["sources"]["n9"]["lock_sha256"],"protocol_sha256":plan["protocol_sha256"],"manifest_sha256":plan["manifest_sha256"],"plan_sha256":c.sha((O/"plan.json").read_bytes()),"development_checkpoint_sha256":[c.sha((PROG/("dev-"+str(i+1).zfill(2)+".json")).read_bytes()) for i in range(len(records))],"train_rows":len(y),"train_correct":int(y.sum()),"oof":oof,"models":models,"cuts":cuts}
    save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes());save("freeze_receipt.json",{"lock_sha256":lh,"plan_sha256":lock["plan_sha256"]})
    checkpoint("12-freeze",{"stage":"FROZEN","completed_groups":11,"total_groups":11+(len(plan["fresh"]) if plan["fresh_allowed"] else 0),"lock_sha256":lh,"train_rows":len(y),"train_correct":int(y.sum()),"oof":oof,"seconds":time.perf_counter()-start})
def summarize():
    plan=verified_plan();freeze=load("freeze_receipt.json");assert c.sha((O/"lock.json").read_bytes())==freeze["lock_sha256"]
    dev=[json.loads((PROG/("dev-"+str(i+1).zfill(2)+".json")).read_text()) for i in range(len(plan["development"]))]
    fresh=[json.loads((PROG/("fresh-"+str(i+1).zfill(2)+".json")).read_text()) for i in range(len(plan["fresh"]))] if plan["fresh_allowed"] else []
    rows=[r for g in fresh for r in g["rows"]];allg=dev+fresh;fast=sum(g["turbo"]["paired_turbo_seconds"] for g in allg);ref=sum(g["turbo"]["paired_reference_seconds"] for g in allg)
    perf={"reference_calls":sum(g["turbo"]["oracle_calls"] for g in allg),"candidate_scores":sum(g["turbo"]["candidate_scores"] for g in allg),"paired_turbo_seconds":fast,"paired_reference_seconds":ref,"paired_speedup":ref/fast,"total_scoring_seconds":sum(g["turbo"]["turbo_seconds"] for g in allg),"total_group_seconds":sum(g["seconds"] for g in allg),"exact_replayed_parent_rows":sum(g["exact_replayed_parent_rows"] for g in dev)}
    out={"status":"FRESH_DEVELOPMENT_HOLDOUT_ONLY" if fresh else "HOLD_INSUFFICIENT_FRESH_GROUPS","fresh_groups":len(fresh),"controls":128*len(fresh),"lock_sha256":freeze["lock_sha256"],"manifest_sha256":plan["manifest_sha256"],"oof":load("lock.json")["oof"],"performance":perf,"synthetic_parity":load("synthetic_parity.json"),"audit_cases":[x for g in dev for x in g["cases"]],"rows":rows,"production_adoption":False}
    if fresh:
        def total(kind):
            r=[q["metrics"] for q in rows if q["kind"]==kind];assert len(r)==len(fresh)
            return {a:{k:sum(x[a][k] for x in r) for k in r[0][a]} for a in r[0]}
        e=total("eval");sh=total("shared");nc=128*len(fresh);base=e["learned05"]["control_alerts"]
        gates={k:e[k]["lost_parent"]==e[k]["lost_joint"]==sh[k]["lost_parent"]==sh[k]["lost_joint"]==0 and e[k]["control_alerts"]*640<=35*nc and e[k]["control_alerts"]*8<=7*base and e[k]["repeat_stress_lost"]==sh[k]["repeat_stress_lost"]==0 and e[k]["repeat_stress_cases"]+sh[k]["repeat_stress_cases"]>0 for k in e}
        out.update(eval=e,shared=sh,gates=gates)
    save("result.json",out);emit("R40N10_SUMMARY",{k:v for k,v in out.items() if k not in ("rows","audit_cases")})
    checkpoint("99-summary",{"stage":"COMPLETE","completed_groups":len(allg),"total_groups":len(allg),"status":out["status"],"performance":perf,"passed_arms":[k for k,v in out.get("gates",{}).items() if v]})
if __name__=="__main__":
    assert ARGS
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0] in ("dev","fresh"):group(ARGS[0],int(ARGS[1]))
    elif ARGS[0]=="fit":fit()
    elif ARGS[0]=="summarize":summarize()
    else:raise ValueError(ARGS)
