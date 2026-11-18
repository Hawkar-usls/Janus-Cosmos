#!/usr/bin/env python3
"""R30: KUSTO reference regimes + TOPA-style empirical joint null. Research only."""
import hashlib,json,pathlib,sys,time
import numpy as np
import r28b_router as rb
import r28_hybrid as h
import r27a_benchmark as a
import r26c_reference_alignment as r
import r24_odontocrystal as od
import r17c_frozen_holdout as hold
import r9_expanded_search as eng

P=pathlib.Path("research/seti/r30_false_alarm_protocol.json")
S=json.loads(P.read_text()); C=S["shield"]; MODES=C["modes"]
OUT=pathlib.Path("runtime/r30");OUT.mkdir(parents=True,exist_ok=True)
h.S=S;rb.S=S;r.S=S;od.S=S;r.OUT=OUT;h.OUT=OUT;rb.OUT=OUT;eng.OUT=OUT
hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))}
hold.P=P;hold.OUT=OUT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def cut(v,rate):
    z=np.sort(np.asarray(v,dtype=float))[::-1]
    assert len(z)>0 and np.isfinite(z).all() and 0<=rate<1
    return float(np.nextafter(z[int(len(z)*rate)],np.inf))

def freeze():
    # Metadata only. No selected spectrum is fetched or scored at this stage.
    raw,_=hold.base.fetch(hold.r11.S["selection"]["index"],limit=32*hold.base.MIB)
    (OUT/"target-index.json").write_bytes(raw)
    pool=set(hold.r11.names(json.loads(raw)))-set(hold.S["exclude"])
    order=sorted(pool,key=lambda t:hashlib.sha256((S["seed"]+"|"+t).encode()).hexdigest())[:C["max_queries"]]
    m={"protocol_sha256":sha(P),"index_sha256":hashlib.sha256(raw).hexdigest(),"candidate_order":order,"selected":[],"attempts":[]}
    days=set(C["exclude_days"]); urls=set(); names=set()
    nfit=C["fit_targets"];nnull=C["null_targets"];total=nfit+nnull+C["eval_targets"]
    for target in order:
        try:
            fs,q=eng.query(target);fs.sort(key=lambda f:(f["day"],f["seconds"],f["url"]))
            candidates=[]
            for i in range(len(fs)-3):
                group=fs[i:i+4];day=int(group[0]["day"])
                if [f["position"] for f in group]!=["S","R","S","R"] or len({f["day"] for f in group})!=1 or day in days:continue
                if len({f["url"] for f in group})!=4 or len({f["name"] for f in group})!=4 or any(f["url"] in urls or f["name"] in names for f in group):continue
                if all(0<v["seconds"]-u["seconds"]<=600 for u,v in zip(group,group[1:])):candidates.append(group)
            candidates.sort(key=lambda x:(x[0]["day"]<58575,x[0]["day"],x[0]["seconds"],x[0]["url"]))
            for group in candidates[:5]:
                try:
                    hs=[eng.hdr(f) for f in group]
                    assert all(eng.valid(hd,f,target) and hd["shape"][0]==16 for hd,f in zip(hs,group))
                    assert len({(hd["shape"][0],round(hd["ts"],7),round(hd["df"],12)) for hd in hs})==1
                    coords=[(hd["src_raj"],hd["src_dej"]) for hd in hs]
                    assert coords[0]==coords[2] and coords[1]==coords[3] and coords[0]!=coords[1] and all(float(c[0])!=0 or float(c[1])!=0 for c in coords)
                    i=len(m["selected"]);day=int(group[0]["day"])
                    split="component_fit" if i<nfit else "null_calibration" if i<nfit+nnull else "evaluation"
                    m["selected"].append({"target":target,"day":day,"split":split,"query":q,"files":[{**f,"header":hd,"role":role} for f,hd,role in zip(group,hs,("ON1","OFF1","ON2","OFF2"))]})
                    days.add(day);urls.update(f["url"] for f in group);names.update(f["name"] for f in group)
                    save("manifest.json",m);break
                except Exception as e:m["attempts"].append({"target":target,"day":int(group[0]["day"]),"error":type(e).__name__+": "+str(e)[:300]})
        except Exception as e:m["attempts"].append({"target":target,"error":type(e).__name__+": "+str(e)[:300]})
        if len(m["selected"])==total:break
    save("manifest.json",m)
    emit("R30_MANIFEST",{"selected":[{"target":q["target"],"day":q["day"],"split":q["split"],"files":[f["name"] for f in q["files"]]} for q in m["selected"]],"attempt_count":len(m["attempts"])})
    assert len(m["selected"])==total,"Insufficient fresh day-disjoint archive targets: leave HOLD, never relax split."

def support(Y,p):
    assert Y.ndim==3 and Y.shape[1]==16 and p.shape==Y.shape[:2]
    assert np.all((p>=0)&(p<Y.shape[2]))
    v=Y[np.arange(len(Y))[:,None],np.arange(16)[None,:],p.astype(int)]
    return np.sort(v.reshape(-1,4,4).mean(2),axis=1)[:,1]

def pack(rows,R,pre,L):
    (q,raw,reg,g),Y=pre
    t=time.perf_counter()
    z,p,regsec,_=h.fuse(q,L,"registered_hybrid")
    u,up,usec,_=h.fuse(q,L,"protected_hybrid")
    # Reference-only covariate: cannot see the ON injection or truth labels.
    feature=np.quantile(np.maximum(np.median(r.norm(R),axis=1),0),.95,axis=1)
    v=support(Y,p)
    # Distinct raw route is retained even when registered path wins.
    branch=("raw_linear","raw_odonto")
    raw_marg=np.stack([(q[k][0]-L["thresholds"][k])/L["scales"][k] for k in branch],1)
    im=raw_marg.argmax(1);ix=np.arange(len(Y));rp=np.stack([q[k][1] for k in branch],1)[ix,im]
    return {"z":z,"path":p,"protected":u,"protected_path":up,"raw_score":raw_marg[ix,im],"raw_path":rp,
            "feature":feature,"support":v,"regime":g,"shared_seconds":float(raw+reg+regsec),
            "additional_scoring_seconds":float(time.perf_counter()-t),"panel_kernel_seconds":float(raw+reg+sum(x[2] for x in q.values())),
            "protected_kernel_seconds":float(raw+reg+usec)}
def group_ids(rows,target):return np.asarray([target+"|"+x[6] for x in rows])
def labels(rows):return h.labels(rows)
def concat(parts,key):return np.concatenate([x[key] for x in parts])
def get_ctrl(part,key):return part["data"][key][labels(part["rows"])==0]
def transform(x,L):
    g=(x["feature"]>L["feature_cut"]).astype(int)
    t=np.minimum(x["z"],(x["support"]-L["support_cut"])/L["support_scale"])
    assert np.isfinite(t).all()
    return g,t
def conditional(values,groups,targets,pooled):
    stats={};th=[]
    for j in (0,1):
        z=values[groups==j];nt=len(set(targets[groups==j].tolist()))
        stats[str(j)]={"controls":int(len(z)),"targets":int(nt)}
        th.append(cut(z,C["shield_rate"]) if len(z)>=C["min_stratum_controls"] and nt>=C["min_stratum_targets"] else None)
    # If either stratum is sparse, use the same pooled null for both. No tiny-stratum extrapolation.
    fallback=any(x is None for x in th)
    if fallback:th=[pooled,pooled]
    actual=int(np.sum(values>=np.asarray(th)[groups]))
    assert actual<=int(len(values)*C["shield_rate"])
    return {"thresholds":th,"strata":stats,"fallback_to_pooled":fallback,"calibration_alerts":actual}
def tile_max(v,groups):
    return np.asarray([float(np.max(v[groups==g])) for g in sorted(set(groups.tolist()))])
def fit(fit_parts,null_parts,L):
    fc=np.concatenate([get_ctrl(p,"feature") for p in fit_parts])
    sc=np.concatenate([get_ctrl(p,"support") for p in fit_parts])
    scut=cut(sc,.05)
    F={"component_lock":L,"feature_cut":float(np.median(fc)),"support_cut":scut,
       "support_scale":max(.1,scut-float(np.quantile(sc,.8)))}
    z=[];t=[];g=[];u=[];groups=[];targets=[]
    for p in null_parts:
        x=p["data"];c=labels(p["rows"])==0;gg,tt=transform(x,F)
        z.extend(x["z"][c]);t.extend(tt[c]);g.extend(gg[c]);u.extend(x["protected"][c])
        groups.extend(group_ids(p["rows"],p["target"])[c]);targets.extend([p["target"]]*int(c.sum()))
    z,t,g,u,groups,targets=map(np.asarray,(z,t,g,u,groups,targets))
    assert len(z)==128*C["null_targets"] and len(set(groups.tolist()))==8*C["null_targets"]
    assert all(np.sum(groups==q)==16 for q in set(groups.tolist()))
    strict=cut(z,C["shield_rate"]);topa=cut(t,C["shield_rate"])
    family=tile_max(t,groups)
    F.update({"thresholds":{"baseline_full":cut(z,C["baseline_rate"]),"pooled_strict":strict,"topa":topa,
              "topa_tile_family":cut(family,C["baseline_rate"]),"protected_raw_union":cut(u,C["baseline_rate"])},
              "kusto":conditional(z,g,targets,strict),"kusto_topa":conditional(t,g,targets,topa),
              "null_controls":int(len(z)),"null_tile_blocks":int(len(family)),"null_targets":sorted(set(targets.tolist())),
              "fit_targets":[p["target"] for p in fit_parts],"null_distribution":{"hybrid":sorted(z.tolist()),"topa":sorted(t.tolist()),"topa_tile_maxima":sorted(family.tolist())}})
    assert np.sum(z>=F["thresholds"]["baseline_full"])<=int(len(z)*C["baseline_rate"])
    assert np.sum(t>=topa)<=int(len(t)*C["shield_rate"])
    assert np.sum(family>=F["thresholds"]["topa_tile_family"])<=int(len(family)*C["baseline_rate"])
    return F
def decisions(x,L):
    # No truth labels, injection parameters, target identity, or exam fitting here.
    g,t=transform(x,L);th=L["thresholds"];z=x["z"]
    scores={"baseline_full":z,"pooled_strict":z,"kusto":z,"topa":t,"kusto_topa":t,
            "topa_tile_family":t,"protected_raw_union":x["protected"]}
    thresholds={k:np.broadcast_to(float(th[k]),z.shape) for k in th}
    for k in ("kusto","kusto_topa"):thresholds[k]=np.asarray(L[k]["thresholds"])[g]
    return {k:(scores[k],thresholds[k],scores[k]>=thresholds[k],x["protected_path"] if k=="protected_raw_union" else x["path"]) for k in MODES},g
def localize(path,rows):
    truth=np.stack([x[7] for x in rows]);active=np.stack([x[8] for x in rows])
    return ((abs(path-truth)<=2)&active).sum(1)>=.75*active.sum(1)
def measure(part,L):
    rows=part["rows"];x=part["data"];target=part["target"]
    dec,g=decisions(x,L);y=labels(rows);weak=np.array([r0[4]<=1.25 for r0 in rows]);base=dec["baseline_full"][2]
    base_loc=localize(dec["baseline_full"][3],rows);out={};rec={}
    for k,(s,th,alarm,p) in dec.items():
        loc=localize(p,rows);pos=y==1;neg=y==0;hit=alarm&loc&pos;bh=base&base_loc&pos
        out[k]={"localized":int(hit.sum()),"detected":int((alarm&pos).sum()),"positives":int(pos.sum()),
          "weak_localized":int((hit&weak).sum()),"control_alerts":int((alarm&neg).sum()),"controls":int(neg.sum()),
          "held_baseline_positive":int((base&~alarm&pos).sum()),"held_baseline_control":int((base&~alarm&neg).sum()),
          "added_localized":int((hit&~bh).sum()),"lost_localized":int((bh&~hit).sum()),
          "added_weak":int((hit&~bh&weak).sum()),"lost_weak":int((bh&~hit&weak).sum()),
          "added_controls":int((alarm&~base&neg).sum()),"removed_controls":int((base&~alarm&neg).sum()),
          "alerted_control_tiles":len(set(group_ids(rows,target)[alarm&neg].tolist())),
          "regimes":{str(j):{"controls":int((neg&(g==j)).sum()),"control_alerts":int((alarm&neg&(g==j)).sum()),
                        "positives":int((pos&(g==j)).sum()),"localized":int((hit&(g==j)).sum())} for j in (0,1)}}
        rec[k]={"score":s.tolist(),"threshold":th.tolist(),"alert":alarm.tolist(),"localized":loc.tolist(),"path":p.tolist(),
                "hold_baseline":(base&~alarm).tolist()}
    record={"target":target,"day":part["day"],"labels":y.tolist(),"amplitudes":[r0[4] for r0 in rows],
        "models":[r0[2] for r0 in rows],"hard":[bool(r0[5]) for r0 in rows],"source_tile":group_ids(rows,target).tolist(),
        "reference_feature":x["feature"].tolist(),"shield_regime":g.tolist(),"reference_registration_regime":x["regime"].tolist(),
        "support":x["support"].tolist(),"raw_score":x["raw_score"].tolist(),"raw_path":x["raw_path"].tolist(),"arms":rec}
    out["timing"]={k:v for k,v in x.items() if k.endswith("seconds")}
    return out,record
def selftest():
    assert cut([5,5,0,0],.25)>5
    z=np.zeros((1,16,128));p=np.full((1,16),17);z[0,:,17]=[9]*4+[4]*4+[3]*4+[-1]*4
    assert support(z,p)[0]==3
    assert tile_max(np.array([3,7,2,9]),np.array(["a","a","b","b"])).tolist()==[7,9]
    c=conditional(np.arange(256),np.zeros(256,int),np.repeat(["a","b","c","d"],64),cut(np.arange(256),C["shield_rate"]))
    assert c["fallback_to_pooled"] and c["thresholds"][0]==c["thresholds"][1] and c["calibration_alerts"]<=int(256*C["shield_rate"])
    return {"tie_safe_cut":"PASS","block_support":"PASS","tile_max":"PASS","sparse_fallback":"PASS"}
def main():
    jit=a.init_tree();checks={**h.selftest(),**selftest()};emit("R30_SELFTEST",checks)
    m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==sha(P)
    parts=m["selected"];nf=C["fit_targets"];nc=C["null_targets"];ne=C["eval_targets"]
    assert len(parts)==nf+nc+ne and len({q["target"] for q in parts})==len(parts)
    assert len({q["day"] for q in parts})==len(parts) and not set(q["day"] for q in parts)&set(C["exclude_days"])
    assert not set(q["target"] for q in parts)&set(hold.S["exclude"])
    assert [q["split"] for q in parts]==["component_fit"]*nf+["null_calibration"]*nc+["evaluation"]*ne
    assert len({f["url"] for q in parts for f in q["files"]})==4*len(parts)
    seen=set();source_ledger={};load_seconds=0.;panel_seconds=0.
    def fetch(q):
        nonlocal load_seconds,panel_seconds
        t0=time.perf_counter();rows,R,ch,hs=r.gather(q);load_seconds+=time.perf_counter()-t0
        assert not seen&hs and len(hs)==16;seen.update(hs)
        assert len(rows)==256
        source_ledger[q["target"]]={"day":q["day"],"split":q["split"],"source_tile_hashes":sorted(hs),"files":[f["name"] for f in q["files"]]}
        save("sources.json",source_ledger)
        t0=time.perf_counter();pre=rb.panel(rows,R);panel_seconds+=time.perf_counter()-t0
        emit("R30_LOADED",{"target":q["target"],"day":q["day"],"split":q["split"],"sources":sorted(hs)})
        return {**q,"rows":rows,"R":R,"pre":pre,"challenge":ch}
    cal=[fetch(q) for q in parts[:nf]]
    L=h.fit([(p["rows"],p["pre"][0]) for p in cal])
    for p in cal:p["data"]=pack(p["rows"],p["R"],p["pre"],L)
    null=[fetch(q) for q in parts[nf:nf+nc]]
    for p in null:p["data"]=pack(p["rows"],p["R"],p["pre"],L)
    F=fit(cal,null,L);F["manifest_sha256"]=sha(OUT/"manifest.json")
    F["calibration_sources"]=sorted(seen)
    save("lock.json",F);lock=sha(OUT/"lock.json")
    emit("R30_LOCK",{"sha256":lock,"fit_targets":F["fit_targets"],"null_targets":F["null_targets"],"thresholds":F["thresholds"],"kusto":F["kusto"],"kusto_topa":F["kusto_topa"],"null_tile_blocks":F["null_tile_blocks"]})
    out={};records={};leak={};leakrec={}
    for q in parts[nf+nc:]:
        p=fetch(q);p["data"]=pack(p["rows"],p["R"],p["pre"],L)
        target=q["target"];out[target],records[target]=measure(p,F)
        lr=[v[0] for v in p["challenge"]];RR=np.stack([v[1] for v in p["challenge"]])
        t0=time.perf_counter();lp=rb.panel(lr,RR);panel_seconds+=time.perf_counter()-t0
        xx=pack(lr,RR,lp,L)
        leak[target],leakrec[target]=measure({"target":target,"day":q["day"],"rows":lr,"data":xx},F)
        save("partial.json",{"per_target":out,"leak":leak})
        emit("R30_TARGET",{"target":target,"day":q["day"],"modes":out[target],"leak":{k:leak[target][k]["localized"] for k in MODES}})
    assert sha(OUT/"lock.json")==lock
    fields=("localized","detected","positives","weak_localized","control_alerts","controls","held_baseline_positive","held_baseline_control","added_localized","lost_localized","added_weak","lost_weak","added_controls","removed_controls","alerted_control_tiles")
    total={k:{f:sum(v[k][f] for v in out.values()) for f in fields} for k in MODES}
    lt={k:{f:sum(v[k][f] for v in leak.values()) for f in ("localized","positives","weak_localized","lost_localized","added_localized")} for k in MODES}
    assert all(v["positives"]==128*ne and v["controls"]==128*ne for v in total.values())
    assert all(v["positives"]==32*ne for v in lt.values())
    b=total["baseline_full"];v=total["kusto_topa"]
    gate={"localized_not_lower":v["localized"]>=b["localized"],"weak_not_lower":v["weak_localized"]>=b["weak_localized"],
          "control_half":b["control_alerts"]>0 and v["control_alerts"]<=b["control_alerts"]//2,
          "shared_reference_not_lower":lt["kusto_topa"]["localized"]>=lt["baseline_full"]["localized"]}
    for k in MODES:
        total[k]["control_alert_fraction"]=total[k]["control_alerts"]/total[k]["controls"]
        total[k]["paired_target_deltas"]={q:{"localized":v0[k]["localized"]-v0["baseline_full"]["localized"],"controls":v0[k]["control_alerts"]-v0["baseline_full"]["control_alerts"],"weak":v0[k]["weak_localized"]-v0["baseline_full"]["weak_localized"]} for q,v0 in out.items()}
    result={"schema":"janus.cosmos.r30.false_alarm.result.v1","status":"ARCHIVE_SYNTHETIC_DEVELOPMENT_EXAM",
       "gate":"PASS_DEVELOPMENT" if all(gate.values()) else "HOLD","gate_checks":gate,"samples":256*ne,
       "selftests":checks,"lock_sha256":lock,"manifest_sha256":F["manifest_sha256"],"protocol_sha256":sha(P),"sources_sha256":sha(OUT/"sources.json"),
       "targets":[{"target":q["target"],"day":q["day"],"split":q["split"]} for q in parts],
       "modes":total,"per_target":out,"leak_challenge":lt,"fit_and_null_and_exam_panel_wall_seconds":panel_seconds,
       "archive_load_seconds":load_seconds,"external_fetch_and_jit_seconds":jit,
       "timing_scope":"Panel computes five shared kernels for ablation; actual total panel seconds and separate kernel/feature subtotals. Not an end-to-end speedup claim.",
       "limitations":["Control alerts are empirical proxies; unmodified archive backgrounds can contain physical signals.","Control patches within source tiles correlate; only four evaluation days; report paired target deltas, no independence-based p-values.","TOPA block support computed on selected full path (not an independent confirmation); empirical null calibrated after this entire selection.","Tile maximum covers sampled windows of each tile, not entire telescope band. Different target/day is not independent instrument.","BLC1/Voyager external source-labelled verification not performed here. No RFI-origin verdict."],
       "scope":S["scope"],"production_adoption":False}
    save("scores.json",records);save("leak_scores.json",leakrec);save("result.json",result)
    emit("R30_SUMMARY",{k:v for k,v in result.items() if k!="per_target"})
if __name__=="__main__":
    phase=sys.argv[1]
    if phase=="freeze":freeze()
    elif phase=="selftest":emit("R30_LOCAL_SELFTEST",selftest())
    elif phase=="scan":main()
    else:raise ValueError(phase)
