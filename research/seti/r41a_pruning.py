#!/usr/bin/env python3
"""R41A: exact frozen pruning audit, no fitting, no fresh exam."""
import sys
ARGS=sys.argv[1:]
import json,pathlib,time,itertools,collections,io
import numpy as np
import r41_consolidation as t
a=t.a;c=t.c
P=pathlib.Path("research/seti/r41a_pruning_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r41a");O.mkdir(parents=True,exist_ok=True);PROG=O/"progress";PROG.mkdir(exist_ok=True);CACHE=O/"cache";CACHE.mkdir(exist_ok=True)
KINDS=("eval","shared");MODES=("operational","common")
def save(k,d):a.write(O/k,d)
def load(k):return json.loads((O/k).read_text())
def cp(k,d):
    a.write(PROG/(k+".json"),d);print("R41A_PROGRESS",json.dumps(d,allow_nan=False),flush=True)
def prepare():
    z,nn=c.artifact(S["source_r41"],32);zn,nz=c.artifact(S["source_n10"],128)
    assert c.sha(z.read(nn["winner_lock.json"]))==S["source_r41"]["winner_lock_sha256"]
    assert c.sha(zn.read(nz["lock.json"]))==S["source_n10"]["lock_sha256"]
    r=json.loads(z.read(nn["result.json"]));assert r["status"]=="HOLD_NO_NEW_WINNER" and r["winner"] is None and r["new_arms"]==246 and r["existing_arms"]==64
    models=json.loads(z.read(nn["models.json"]))
    freeze=json.loads(zn.read(nz["freeze_receipt.json"]));raw=zn.read(nz["plan.json"]);assert c.sha(raw)==freeze["plan_sha256"] and freeze["lock_sha256"]==S["source_n10"]["lock_sha256"]
    plan=json.loads(raw);groups=plan["development"];assert len(groups)==11 and sum(x["fit"] for x in groups)==6
    assert len({x["group"]["day"] for x in groups})==11
    hashes=[k for x in groups for k in x["raw_tile_hashes"]];assert len(hashes)==len(set(hashes))==176
    used={k:models[v] for k,v in S["parents"].items()}
    for i,g in enumerate(groups):
        for kind in KINDS:
            fn="dev_"+str(i+1).zfill(2)+"_"+kind+".npz";(CACHE/fn).write_bytes(zn.read(nz[fn]))
            d=a.read_npz(CACHE/fn);assert len(d["y"])==(256 if kind=="eval" else 32)
            assert np.array_equal(d["idx"],np.flatnonzero(d["parent_alert"]))
            assert set(d["source_tile_sha"].tolist())<=set(g["raw_tile_hashes"])
    save("source_report.json",r);save("models.json",used)
    save("plan.json",{"protocol_sha256":c.sha(P.read_bytes()),"source_n10_plan_sha256":freeze["plan_sha256"],"groups":groups,"models_sha256":c.sha((O/"models.json").read_bytes()),"source_report_sha256":c.sha((O/"source_report.json").read_bytes())})
    cp("00-prepare",{"stage":"PREPARE","groups":11,"fit_groups":6,"audit_groups":5,"rows":3168,"prior_tiles":176,"rules":len(S["parent_rules"])*2,"models_sha256":c.sha((O/"models.json").read_bytes()),"no_fitting":True})
def verified():
    p=load("plan.json");assert p["protocol_sha256"]==c.sha(P.read_bytes()) and p["models_sha256"]==c.sha((O/"models.json").read_bytes()) and p["source_report_sha256"]==c.sha((O/"source_report.json").read_bytes());return p
def masks(d,model,mode):
    blocks=model["blocks"];use=t.support(d,blocks if mode=="operational" else t.BLOCKS)
    th=model["cuts"]["preserve"]
    v=a.score(t.feat(d,blocks),model["model"],model["kind"])
    sv=a.score(t.feat(d,blocks,True),model["model"],model["kind"])
    return t.mask(d,v,th,use),t.mask(d,sv,th,use),use
def combine(rule,parts,shape):
    if not parts:return np.ones(shape,bool)
    if "_or_" in rule:return np.logical_or.reduce(parts)
    return np.logical_and.reduce(parts)
def group(i):
    p=verified();g=p["groups"][i];models=load("models.json");r=load("source_report.json")
    rec={"stage":"GROUP","group":i+1,"target":g["group"]["target"],"day":g["group"]["day"],"fit_viewed":g["fit"],"modes":{},"parity":True};t0=time.perf_counter()
    for mode in MODES:
        bykind={}
        for kind in KINDS:
            d=a.read_npz(CACHE/("dev_"+str(i+1).zfill(2)+"_"+kind+".npz"));aa=d["parent_alert"];hh=d["parent_hit"];res=d["joint_rescue"];co=d["correct"];y=d["y"];null=y==0;idx=d["idx"]
            assert not (aa&res).any()
            learned=res&d["learned05_alert"];stress=d["repeat_stress_supported"]&t.support(d,t.BLOCKS)&hh[idx]
            vals={};ss={};supports={};cost={}
            for name,m in models.items():
                start=time.perf_counter();vals[name],ss[name],supports[name]=masks(d,m,mode);cost[name]=time.perf_counter()-start
            out={};bits={};lostbits={}
            for rule,parts in S["parent_rules"].items():
                k=combine(rule,[vals[x] for x in parts],len(y));sk=combine(rule,[ss[x] for x in parts],len(y))
                for rt,r0 in (("joint",res),("learned05",learned)):
                    name=rule+"__"+rt
                    al=(aa&k)|r0;hit=(hh&k)|(r0&co)
                    m=t.metrics(d,al,hit,(stress&~sk[idx]).sum())
                    out[name]={**m,"parent_control_alerts":int((aa&k&null).sum()),"rescue_control_alerts":int((r0&null).sum()),"removed_parent_controls":int((aa&~k&null).sum()),"removed_rescue_controls":int((res&~r0&null).sum()),"lost_rescue_hits":int((res&co&~r0).sum()),"parent_support_rows":int(len(idx))}
                    bits[name]=((aa|res)&~al&null);lostbits[name]=(hh|(res&co))&~hit
            # Exact immutable R41 reproduction, per-group and metric-by-metric.
            if not g["fit"]:
                for name,oldname in [("identity__joint","frozen_joint_min"),("identity__learned05","frozen_learned05")]+[(k+"__joint",S["parents"][k]+"__preserve") for k in S["parents"]]:
                    arm=r["arms"][oldname]
                    if mode=="operational" or name.startswith("identity"):
                        ref=next(x[2] for x in arm["per_group"] if x[0]==i and x[1]==kind)
                        assert all(out[name][k]==v for k,v in ref.items()),"R41 decision metrics parity "+name
            mini=bits["minimal__joint"];full=bits["full__joint"];jm=aa|res
            bykind[kind]={"metrics":out,"overlap":{"minimal_removed":int(mini.sum()),"full_removed":int(full.sum()),"both_removed":int((mini&full).sum()),"only_minimal":int((mini&~full).sum()),"only_full":int((full&~mini).sum()),"either_removed":int((mini|full).sum()),"residual_controls":int((jm&~mini&null).sum()),"protected_rescue_controls":int((res&null).sum()),"residual_parent_controls":int((aa&~mini&null).sum())},"support":{k:int(v.sum()) for k,v in supports.items()},"cached_score_seconds":cost}
            if mode=="operational":
                # These labels are audit outputs, never input decisions or future calibration.
                rows=[]
                for j in np.flatnonzero((jm&~mini&null)):
                    jj=int(np.searchsorted(idx,j));row={"index":int(j),"tile_sha256":str(d["source_tile_sha"][j]),"lane":"parent" if aa[j] else "rescue","full_would_remove":bool(full[j]),"rival_would_remove":bool(bits["rival__joint"][j]),"runs_would_remove":bool(bits["runs__joint"][j]),"learned05_would_remove":bool(bits["identity__learned05"][j])}
                    if aa[j] and jj<len(idx) and idx[jj]==j:row["native_parent_path"]=d["parent_path"][j].tolist();row["support"]={k:bool(v[jj]) for k,v in supports.items()}
                    rows.append(row)
                bykind[kind]["residual_ledger"]=rows
                bykind[kind]["loss_ledger"]={k:np.flatnonzero(v).tolist() for k,v in lostbits.items() if v.any()}
        rec["modes"][mode]=bykind
    rec["wall_seconds"]=time.perf_counter()-t0;save("group-"+str(i+1).zfill(2)+".json",rec)
    cp("group-"+str(i+1).zfill(2),{"stage":"GROUP","group":i+1,"total_groups":11,"target":rec["target"],"fit_viewed":g["fit"],"parity":True,"ordinary":{k:{m:v[m] for m in ("localized","control_alerts","lost_joint","repeat_stress_lost")} for k,v in rec["modes"]["operational"]["eval"]["metrics"].items()},"overlap":rec["modes"]["operational"]["eval"]["overlap"],"seconds":rec["wall_seconds"]})
def aggregate(records,mode):
    out={}
    for name in records[0]["modes"][mode]["eval"]["metrics"]:
        out[name]={}
        for kind in KINDS:
            rows=[d["modes"][mode][kind]["metrics"][name] for d in records]
            out[name][kind]={k:sum(v[k] for v in rows) for k in rows[0]}
    return out
def summary():
    p=verified();records=[load("group-"+str(i+1).zfill(2)+".json") for i in range(11)];audit=[r for r in records if not r["fit_viewed"]];fit=[r for r in records if r["fit_viewed"]]
    assert len(audit)==5 and len(fit)==6 and all(r["parity"] for r in records)
    agg={m:aggregate(audit,m) for m in MODES};tr={m:aggregate(fit,m) for m in MODES}
    for mode in MODES:
        assert agg[mode]["identity__joint"]["eval"]["controls"]==640
    old=load("source_report.json");orig=agg["operational"]["minimal__joint"];full=agg["operational"]["full__joint"]
    assert orig["eval"]["control_alerts"]==53 and full["eval"]["control_alerts"]==54
    assert orig["eval"]["localized"]==205 and orig["shared"]["localized"]==113
    overlap={k:sum(g["modes"]["operational"]["eval"]["overlap"][k] for g in audit) for k in audit[0]["modes"]["operational"]["eval"]["overlap"]}
    safe=[];qualified=[]
    def blocks(rule):
        return set(x for k in S["parent_rules"][rule] for x in load("models.json")[k]["blocks"])
    def zero(m):return all(m[typ][k]==0 for typ in KINDS for k in ("lost_joint","lost_parent","repeat_stress_lost"))
    for name,m in agg["operational"].items():
        if zero(m) and zero(agg["common"][name]):safe.append(name)
        if name.endswith("__joint") and name!="identity__joint" and all(t.safe(agg[mode][name],agg[mode]["identity__learned05"]) for mode in MODES):qualified.append(name)
    rank=lambda name:(agg["operational"][name]["eval"]["control_alerts"],len(blocks(name.split("__")[0])),len(S["parent_rules"][name.split("__")[0]]),name)
    safe.sort(key=rank);qualified.sort(key=rank);winner=qualified[0] if qualified else None
    chosen=winner or "minimal__joint";rule=chosen.split("__")[0];ms=load("models.json")
    lock={"status":"DEVELOPMENT_ONLY_CANDIDATE" if winner else "HOLD_NO_SCIENTIFIC_WINNER__ENGINEERING_COMPARATOR_ONLY","winner":winner,"engineering_comparator":chosen,"protocol_sha256":p["protocol_sha256"],"source_models_sha256":p["models_sha256"],"source_n10_lock_sha256":S["source_n10"]["lock_sha256"],"rule":rule,"rescue":"joint","models":{k:ms[k] for k in S["parent_rules"][rule]},"no_training":True,"production_adoption":False}
    save("pruning_lock.json",lock)
    report={"status":lock["status"],"groups":11,"fit_viewed_groups":6,"audit_viewed_groups":5,"rules":len(agg["operational"]),"parity_groups":11,"winner":winner,"engineering_comparator":chosen,"lock_sha256":c.sha((O/"pruning_lock.json").read_bytes()),"audit":agg,"fit_viewed_in_sample":tr,"overlap":overlap,"zero_loss_rank":safe,"qualified":qualified,"residual_controls_to_target":max(0,orig["eval"]["control_alerts"]-35),"best_zero_loss_to_target":max(0,agg["operational"][safe[0]]["eval"]["control_alerts"]-35) if safe else None,"r41_cost_scope":old["cost_scope"],"r41_minimal_conservative_seconds":old["arms"]["rival+runs__tree__preserve"]["conservative_cost_seconds"],"r41_full_conservative_seconds":old["arms"]["temporal+rival+runs+folds+halo+artmoney+repeat__tree__preserve"]["conservative_cost_seconds"],"cached_seconds":sum(g["wall_seconds"] for g in records),"no_training":True,"external_validation":False,"production_adoption":False}
    save("result.json",report)
    cp("99-summary",{**{k:v for k,v in report.items() if k not in ("audit","fit_viewed_in_sample","r41_cost_scope")},"arms":{k:{typ:{z:v[typ][z] for z in ("localized","weak","control_alerts","lost_joint","repeat_stress_lost","parent_control_alerts","rescue_control_alerts")} for typ in KINDS} for k,v in agg["operational"].items()}})
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="group":group(int(ARGS[1]))
    elif ARGS[0]=="summary":summary()
    else:raise ValueError(ARGS)
