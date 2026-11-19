#!/usr/bin/env python3
"""R41 common archived developmental panel; feature ablation and cost replay."""
import sys
ARGS=sys.argv[1:]
import collections,itertools,json,pathlib,time,resource,io
import numpy as np
from sklearn.metrics import roc_auc_score
import r40n10_artmoney_turbo as a
n=a.n;c=a.c;b=a.b;h=a.h;w=a.w;old9=a.old9;halo=a.halo;turbo=a.turbo
P=pathlib.Path("research/seti/r41_consolidation_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r41");O.mkdir(parents=True,exist_ok=True);CACHE=O/"cache";CACHE.mkdir(exist_ok=True);PROG=O/"progress";PROG.mkdir(exist_ok=True)
a.O=O;a.CACHE=CACHE;a.PROG=PROG;n.O=O;n.k.O=O;n.l.O=O;a.eng.OUT=O
BLOCKS=tuple(S["blocks"]);KINDS=("eval","shared");CUTS=("preserve","q01","q05")
def save(k,d):a.write(O/k,d)
def load(k):return json.loads((O/k).read_text())
def checkpoint(k,d):
    a.write(PROG/(k+".json"),d);print("R41_PROGRESS",json.dumps(d,allow_nan=False),flush=True)
def prepare():
    z,nn=c.artifact(S["source"],128)
    for key in ("lock","manifest"):
        assert c.sha(z.read(nn[key+".json"]))==S["source"][key+"_sha256"]
    plan=json.loads(z.read(nn["plan.json"]));receipt=json.loads(z.read(nn["freeze_receipt.json"]))
    assert receipt["plan_sha256"]==c.sha(z.read(nn["plan.json"])) and receipt["lock_sha256"]==S["source"]["lock_sha256"]
    gs=plan["development"];assert len(gs)==11 and sum(x["fit"] for x in gs)==6
    assert len({x["group"]["day"] for x in gs})==11
    hashes=[k for x in gs for k in x["raw_tile_hashes"]];assert len(hashes)==len(set(hashes))==176
    for i,g in enumerate(gs):
        for kind in KINDS:
            fn="dev_"+str(i+1).zfill(2)+"_"+kind+".npz";raw=z.read(nn[fn]);(CACHE/fn).write_bytes(raw)
            d=a.read_npz(CACHE/fn);assert len(d["y"])==(256 if kind=="eval" else 32)
            assert np.array_equal(d["idx"],np.flatnonzero(d["parent_alert"]))
            for block in BLOCKS:assert len(d[block])==len(d["idx"]) and np.isfinite(d[block]).all()
    (CACHE/"n10_lock.json").write_bytes(z.read(nn["lock.json"]))
    for key in ("n7","n8","n9"):
        zz,nn0=a.checked_archive(key);(CACHE/(key+"_lock.json")).write_bytes(zz.read(nn0["lock.json"]))
    h.O=CACHE;L=h.readlocks();a.write(CACHE/"parent_locks.json",L)
    zz,nn0=n.archive("r40l");raw=zz.read(nn0["lock.json"]);assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];(CACHE/"joint_lock.json").write_bytes(raw)
    save("plan.json",{"protocol_sha256":c.sha(P.read_bytes()),"n10_plan_sha256":receipt["plan_sha256"],"development":gs,"source_hashes":sorted(hashes),"source_lock_sha256":S["source"]["lock_sha256"]})
    checkpoint("00-prepare",{"stage":"PREPARE","groups":11,"fit_groups":6,"audit_groups":5,"rows":3168,"source_tiles":176,"protocol_sha256":c.sha(P.read_bytes()),"external_validation":False})
def plan():
    p=load("plan.json");assert p["protocol_sha256"]==c.sha(P.read_bytes());return p
def archive_row(i,kind):return a.read_npz(CACHE/("dev_"+str(i+1).zfill(2)+"_"+kind+".npz"))
def source_lock():
    raw=(CACHE/"n10_lock.json").read_bytes();assert c.sha(raw)==S["source"]["lock_sha256"];return json.loads(raw)
def profile(i):
    x=plan()["development"][i];g=x["group"];name=g["target"];seed={"n6":a.v.C["seed"],"n7":w.C["seed"],"n9":old9.C["seed"]}[x["source"]]
    n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"])
    t0=time.perf_counter();p0=time.process_time();cold0=time.perf_counter()
    paths=halo.candidates(b.D);turbo.bank_scores(np.zeros((2,16,176)),paths-32);turbo.legal_mask(paths,np.full((2,16),100));b.taylor.init_tree()
    cold=time.perf_counter()-cold0
    t=time.perf_counter();rr,R,ch,hs,meta,cm,pre=n.gather(g,"r41_"+name);acq=time.perf_counter()-t
    assert sorted(hs)==x["raw_tile_hashes"]
    accel=turbo.Turbo(paths,6);halo.context_scores=accel
    L=json.loads((CACHE/"parent_locks.json").read_text());raw=(CACHE/"joint_lock.json").read_bytes();assert c.sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];old=json.loads(raw)
    L7=a.lock_cache("n7");L8=a.lock_cache("n8")
    records=[]
    for kind,rows,refs,mm in (("eval",rr,R,meta),("shared",[q[0] for q in ch],np.stack([q[1] for q in ch]),cm)):
        d=archive_row(i,kind);a.verify_rows(d,rows);cost={};capture={};hooks=[]
        def hook(mod,key,name,retain=False):
            original=getattr(mod,key)
            def wrapper(*args,**kwargs):
                t=time.perf_counter();ret=original(*args,**kwargs);cost[name]=cost.get(name,0.)+time.perf_counter()-t
                if retain:capture[name]=ret
                return ret
            hooks.append((mod,key,original));setattr(mod,key,wrapper)
        hook(n,"encode","encode_joint")
        hook(b,"panel","panel_all_bank")
        hook(b,"base_arms","bank_fusion",True)
        hook(a.l,"parent","parent_scoring")
        hook(a.l.f,"historical","historical_nested",True)
        hook(h,"replay","screen_replay_nested",True)
        hook(w,"witness","temporal_rival")
        hook(old9,"extra","runs_folds")
        t=time.perf_counter()
        try:replay=a.parent_fresh(rows,refs,mm,pre,L,old,L7,L8)
        finally:
            for mod,key,original in reversed(hooks):setattr(mod,key,original)
        cost["parent_fresh_total"]=time.perf_counter()-t
        for key in ("parent_alert","parent_hit","parent_path","joint_rescue","learned05_alert","correct","idx","temporal","rival","supported","runs","folds","fold_supported"):
            assert np.array_equal(replay[key],d[key]),"R41 raw replay parity "+key
        t=time.perf_counter();ref_before=accel.reference_s
        ha=halo.extract(rows,refs,mm,pre,d["idx"],d["parent_path"][d["idx"]].astype(int),n.ref,c.sha,paths)
        cost["halo_art_including_reference"]=time.perf_counter()-t;cost["halo_reference_extra"]=accel.reference_s-ref_before
        for key in ("halo","artmoney","halo_supported","intersection_counts","halo_source_sha"):assert np.array_equal(ha[key],d[key]),"R41 halo parity "+key
        t=time.perf_counter();rv=turbo.repeat_witness(rows,mm,pre,d["idx"],d["parent_path"][d["idx"]].astype(int),n.ref,kind)
        cost["repeat"]=time.perf_counter()-t
        for key in rv:assert np.array_equal(rv[key],d[key]),"R41 repeat parity "+key
        q=capture["historical_nested"];_,_,f_alert,g_alert=capture["screen_replay_nested"]
        assert np.array_equal(g_alert,d["parent_alert"])
        out={}
        for key in ("r28_frozen","r40c_fused"):
            out[key+"_alert"]=q[key+"_alert"];out[key+"_hit"]=q[key+"_localized"]
        out["r40f_priority_alert"]=f_alert;out["r40f_priority_hit"]=f_alert&q["r40c_fused_localized"]
        bank={}
        for key,(score,path,seconds) in capture["bank_fusion"].items():
            out["bank_"+key+"_score"]=score
            out["bank_"+key+"_localized"]=a.l.correct({"truth":d["truth"],"active":d["active"],"y":d["y"]},path)
            bank[key]=float(seconds)
        np.savez_compressed(O/("legacy_"+str(i+1).zfill(2)+"_"+kind+".npz"),**out)
        costs={**cost,"parent_required":sum(cost[k] for k in ("encode_joint","panel_all_bank","bank_fusion","parent_scoring")),"shared_witness_normalization_and_other":max(0.,cost["parent_fresh_total"]-sum(cost.get(k,0.) for k in ("encode_joint","panel_all_bank","bank_fusion","parent_scoring","temporal_rival","runs_folds"))),"halo_art_operational":max(0.,cost["halo_art_including_reference"]-cost["halo_reference_extra"])}
        records.append({"kind":kind,"rows":len(rows),"parent_alerts":len(d["idx"]),"cost_seconds":costs,"bank_cost_seconds":bank})
    checkpoint("profile-"+str(i+1).zfill(2),{"stage":"RAW_PROFILE","completed_groups":i+1,"total_groups":11,"target":name,"day":g["day"],"raw_tiles":sorted(hs),"parity":True,"acquisition_plus_injections_seconds":acq,"cold_setup_seconds":cold,"full_observed_wall_seconds":time.perf_counter()-t0,"full_observed_cpu_seconds":time.process_time()-p0,"rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"turbo":accel.report(),"records":records})
def subsets():
    vals=[(x,) for x in BLOCKS]+list(itertools.combinations(BLOCKS,2))+[BLOCKS]+[tuple(x for x in BLOCKS if x!=omit) for omit in BLOCKS]
    vals += [("temporal","rival"),("temporal","rival","runs"),("temporal","rival","runs","folds"),("temporal","rival","halo","artmoney"),("halo","artmoney","repeat"),("temporal","rival","repeat")]
    vals=list(dict.fromkeys(tuple(x for x in BLOCKS if x in v) for v in vals))
    return vals
def feat(d,blocks,stress=False):
    return np.column_stack([d["repeat_stress"] if k=="repeat" and stress else d[k] for k in blocks])
def support(d,blocks):
    use=np.ones(len(d["idx"]),bool)
    for k in blocks:
        if k in ("temporal","rival"):use &= d["supported"]
        if k in ("runs","folds"):use &= d["supported"]&d["fold_supported"]
        if k in ("halo","artmoney"):use &= d["halo_supported"]
        if k=="repeat":use &= d["repeat_supported"]
    return use
def mask(d,s,th,use):
    keep=np.ones(len(d["y"]),bool);keep[d["idx"]]=(~use)|(s>=th);return keep
def metrics(d,al,hit,stress_lost=None):
    aa=d["parent_alert"];hh=d["parent_hit"];jm=hh|(d["joint_rescue"]&d["correct"]);weak=d["weak"];y=d["y"];idx=d["idx"]
    stress=d["repeat_stress_supported"]&support(d,BLOCKS)&hh[idx]
    return {"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum()),"control_alerts":int((al&(y==0)).sum()),"lost_parent":int((hh&~hit).sum()),"lost_joint":int((jm&~hit).sum()),"lost_joint_weak":int((jm&~hit&weak).sum()),"false_paths_on_positive":int((al&(y==1)&~hit).sum()),"repeat_stress_cases":int(stress.sum()),"repeat_stress_lost":int((stress&~hit[idx]).sum()) if stress_lost is None else int(stress_lost)}
def aggregate(rows):
    return {kind:{k:sum(x[2][k] for x in rows if x[1]==kind) for k in rows[0][2]} for kind in KINDS}
def score_masks(d,s,ss,th,use):
    keep=mask(d,s,th,use);aa=d["parent_alert"];hh=d["parent_hit"];res=d["joint_rescue"];co=d["correct"]
    stress=d["repeat_stress_supported"]&support(d,BLOCKS)&hh[d["idx"]]
    return metrics(d,aa&keep|res,hh&keep|(res&co),(stress&use&(ss<th)).sum())
def thresholds(v):
    return {"preserve":float(np.nextafter(v.min(),-np.inf)),"q01":float(np.quantile(v,.01)),"q05":float(np.quantile(v,.05))}
def safe(m,learned):
    e=m["eval"];s=m["shared"]
    return e["lost_parent"]==e["lost_joint"]==s["lost_parent"]==s["lost_joint"]==e["repeat_stress_lost"]==s["repeat_stress_lost"]==0 and e["control_alerts"]*640<=35*e["controls"] and e["control_alerts"]*8<=7*learned["eval"]["control_alerts"] and e["repeat_stress_cases"]+s["repeat_stress_cases"]>0
def fit_tournament():
    t0=time.perf_counter();p=plan();profiles=[json.loads((PROG/("profile-"+str(i+1).zfill(2)+".json")).read_text()) for i in range(11)]
    assert all(x["parity"] and x["raw_tiles"]==g["raw_tile_hashes"] for x,g in zip(profiles,p["development"]))
    data=[(i,kind,archive_row(i,kind)) for i in range(11) for kind in KINDS]
    fit=[x for x in data if p["development"][x[0]]["fit"]];audit=[x for x in data if not p["development"][x[0]]["fit"]]
    y=np.concatenate([d["parent_hit"][d["idx"]][support(d,BLOCKS)] for _,_,d in fit])
    group=np.concatenate([np.repeat(i,int(support(d,BLOCKS).sum())) for i,_,d in fit]);assert len(set(group))==6 and y.any() and (~y).any()
    L10=source_lock();L9=a.lock_cache("n9");fixed=collections.defaultdict(list)
    legacy={}
    for i,kind,d in data:
        ll=a.read_npz(O/("legacy_"+str(i+1).zfill(2)+"_"+kind+".npz"));legacy[i,kind]=ll
        if p["development"][i]["fit"]:continue
        old=a.measure(d,L10,L9)
        for key,m in old.items():
            mm=metrics(d,d[key+"_alert"],d[key+"_hit"],m["repeat_stress_lost"])
            fixed["frozen_"+key].append((i,kind,mm))
        for key in ("r28_frozen","r40c_fused","r40f_priority"):
            fixed[key].append((i,kind,metrics(d,ll[key+"_alert"],ll[key+"_hit"])))
    # Independent legacy component scores use fixed fit-only background calibration at common budget.
    bankkeys=sorted(k[5:-6] for k in next(iter(legacy.values())) if k.startswith("bank_") and k.endswith("_score"))
    bank_thresholds={}
    for key in bankkeys:
        null=np.concatenate([legacy[i,kind]["bank_"+key+"_score"][d["y"]==0] for i,kind,d in fit if kind=="eval"])
        th=float(a.l.cut(null,35/640));bank_thresholds[key]=th
        for i,kind,d in audit:
            ll=legacy[i,kind];al=ll["bank_"+key+"_score"]>=th;hit=al&ll["bank_"+key+"_localized"]
            fixed["calibrated_bank_"+key].append((i,kind,metrics(d,al,hit)))
    result={key:{"family":"fixed_or_legacy","audit":aggregate(rows),"per_group":rows} for key,rows in fixed.items()}
    learned=result["frozen_learned05"]["audit"]
    models={};scores_wall={}
    # Same fit and audit labels and common availability for every ablation.
    for blocks in subsets():
        label="+".join(blocks)
        X=np.concatenate([feat(d,blocks)[support(d,BLOCKS)] for _,_,d in fit])
        for kind in ("logistic","tree"):
            k=label+"__"+kind;start=time.perf_counter()
            ss=np.zeros(len(y));oof_rows={cut:[] for cut in CUTS}
            for gi in sorted(set(group)):
                te=group==gi;sc,ser=a.v.model_fit(X[~te],y[~te],kind);ss[te]=sc(X[te]);cuts=thresholds(sc(X[~te][y[~te]]))
                for i,typ,d in fit:
                    if i!=gi:continue
                    v=a.score(feat(d,blocks),ser,kind);sv=a.score(feat(d,blocks,True),ser,kind)
                    for cut,th in cuts.items():oof_rows[cut].append((i,typ,score_masks(d,v,sv,th,support(d,BLOCKS))))
            auc=float(roc_auc_score(y,ss));sc,ser=a.v.model_fit(X,y,kind);cuts=thresholds(sc(X[y]))
            assert np.allclose(sc(X),a.score(X,ser,kind),atol=1e-12,rtol=1e-12),"serialized parity"
            fit_seconds=time.perf_counter()-start;models[k]={"blocks":blocks,"kind":kind,"model":ser,"cuts":cuts,"fit_seconds":fit_seconds,"oof_auc":auc}
            per={cut:[] for cut in CUTS};common={cut:[] for cut in CUTS};wall=0.
            for i,typ,d in audit:
                start=time.perf_counter();v=a.score(feat(d,blocks),ser,kind);wall+=time.perf_counter()-start
                sv=a.score(feat(d,blocks,True),ser,kind)
                for cut,th in cuts.items():
                    per[cut].append((i,typ,score_masks(d,v,sv,th,support(d,blocks))))
                    common[cut].append((i,typ,score_masks(d,v,sv,th,support(d,BLOCKS))))
            scores_wall[k]=wall
            for cut in CUTS:
                r={"family":"new_ablation","model_key":k,"blocks":blocks,"classifier":kind,"cut":cut,"oof_auc":auc,"oof":aggregate(oof_rows[cut]),"audit":aggregate(per[cut]),"audit_common_support":aggregate(common[cut]),"per_group":per[cut],"score_wall_seconds":wall}
                result[k+"__"+cut]=r
        print("R41_SUBSET",label,"models_completed",len(models),flush=True)
    # Costs on the same five audit days. Shared feature dependencies charged once.
    costbase=0.;base_parent=0.;base_normalize=0.;cost_blocks=collections.Counter();bank_cost=collections.Counter()
    for i,g in enumerate(p["development"]):
        if g["fit"]:continue
        pr=profiles[i];costbase+=pr["acquisition_plus_injections_seconds"]+pr["cold_setup_seconds"]
        for row in pr["records"]:
            co=row["cost_seconds"];base_parent+=co["parent_required"];base_normalize+=co["shared_witness_normalization_and_other"]
            for key in ("temporal_rival","runs_folds","halo_art_operational","repeat"):cost_blocks[key]+=co.get(key,0.)
            bank_cost.update(row["bank_cost_seconds"])
    def cost(r):
        bs=r["blocks"];keys=set()
        if any(x in bs for x in ("temporal","rival")):keys.add("temporal_rival")
        if any(x in bs for x in ("runs","folds")):keys.add("runs_folds")
        if any(x in bs for x in ("halo","artmoney")):keys.add("halo_art_operational")
        if "repeat" in bs:keys.add("repeat")
        return costbase+base_parent+base_normalize+sum(cost_blocks[k] for k in keys)+r["score_wall_seconds"]
    # Freeze only a new candidate satisfying BOTH operational and matched-support gates.
    # Existing anchors are reported separately; no replacement if new candidate fails.
    candidates=[]
    for key,r in result.items():
        r["gate"]=safe(r["audit"],learned)
        if r["family"]=="new_ablation":
            r["matched_gate"]=safe(r["audit_common_support"],learned);r["conservative_cost_seconds"]=cost(r)
            if r["gate"] and r["matched_gate"]:candidates.append(key)
    rank=lambda k:(result[k]["audit"]["eval"]["control_alerts"],result[k]["conservative_cost_seconds"],len(result[k]["blocks"]),k)
    candidates.sort(key=rank)
    winner=candidates[0] if candidates else None
    desc=[k for k,r in result.items() if r["family"]=="new_ablation" and all(r["audit"][typ][loss]==0 for typ in KINDS for loss in ("lost_joint","lost_parent","repeat_stress_lost"))]
    desc.sort(key=rank)
    # Compare frozen anchors as well: if an older arm passes with fewer alerts, no alleged improvement.
    older=[k for k,r in result.items() if r["family"]=="fixed_or_legacy" and r["gate"]]
    older.sort(key=lambda k:(result[k]["audit"]["eval"]["control_alerts"],k))
    if winner and older and result[older[0]]["audit"]["eval"]["control_alerts"]<result[winner]["audit"]["eval"]["control_alerts"]:winner=None
    frozen={"schema":"janus.r41.development_winner.v1","protocol_sha256":c.sha(P.read_bytes()),"source_n10_lock_sha256":S["source"]["lock_sha256"],"source_plan_sha256":p["n10_plan_sha256"],"winner":winner,"status":"DEVELOPMENT_ONLY_FROZEN_CANDIDATE" if winner else "HOLD_NO_NEW_WINNER","no_production_adoption":True}
    if winner:
        m=models[result[winner]["model_key"]];frozen.update(model=m,threshold=m["cuts"][result[winner]["cut"]],audit=result[winner]["audit"])
    save("winner_lock.json",frozen)
    # Retain all serialized models for audit, not deployment.
    save("models.json",models)
    report={"status":frozen["status"],"groups":11,"fit_groups":6,"audit_groups":5,"subsets":len(subsets()),"new_arms":len(models)*3,"existing_arms":len(fixed),"fit_common_parent_rows":len(y),"fit_correct":int(y.sum()),"winner":winner,"winner_lock_sha256":c.sha((O/"winner_lock.json").read_bytes()),"source_n10_lock_sha256":S["source"]["lock_sha256"],"qualified_new_arms":candidates,"qualified_existing_arms":older,"zero_loss_rank":desc[:20],"bank_thresholds":bank_thresholds,"cost_scope":{"audit_acquisition_and_cold_seconds":costbase,"parent_required_seconds":base_parent,"shared_normalization_and_other_seconds":base_normalize,"witness_blocks":dict(cost_blocks),"legacy_bank_cost_seconds":dict(bank_cost),"observed_wall_all_11":sum(x["full_observed_wall_seconds"] for x in profiles),"observed_cpu_all_11":sum(x["full_observed_cpu_seconds"] for x in profiles),"max_rss_kib":max(x["rss_kib"] for x in profiles),"note":"Measured full common replay, estimated conservative disjoint block sum for candidate; selected implementation requires R41B standalone measurement."},"arms":result,"tournament_seconds":time.perf_counter()-t0,"production_adoption":False,"external_validation":False}
    save("result.json",report)
    checkpoint("99-summary",{k:v for k,v in report.items() if k not in ("arms","bank_thresholds")})
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="profile":profile(int(ARGS[1]))
    elif ARGS[0]=="tournament":fit_tournament()
    else:raise ValueError(ARGS)
