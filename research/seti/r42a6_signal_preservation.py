"""JANUS R42A6: source-verified signal-loss diagnostics and fit-only spectral guard.
This is an exploratory historical check. Audit groups have previously been observed.
"""
import hashlib,io,json,pathlib,sys
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import r42a5_combined_score as prior

P=pathlib.Path("research/seti/r42a6_signal_preservation_protocol.json")
Q=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r42a6");OUT.mkdir(parents=True,exist_ok=True)
DIMS=128;EPOCHS=16
FEATURE_NAMES=["M","N","N_supported","parent_route","ridge_median_snr","ridge_q25_snr",
               "ridge_occupancy","ridge_shoulder","ridge_jitter","path_edge_frac","path_drift"]
def sha(b):return hashlib.sha256(b).hexdigest()
def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def identity(d,i):
    return {"group":int(d["group"][i]),"kind":"shared" if d["shared"][i] else "eval",
            "row":int(d["row"][i])}
def summarize(mask,d,filter_mask):
    # Positive has to be correctly localized, not merely injected.
    return prior.metrics(filter_mask,d,mask)
def load_history():
    src=Q["provenance"]
    r42=prior.archive(src["r42a5"])
    result=json.loads(r42("result.json"));lock_raw=r42("fit_lock.json")
    lock=json.loads(lock_raw)
    assert sha(lock_raw)==src["r42a5"]["fit_lock_sha256"]
    assert result["fit_lock_sha256"]==sha(lock_raw)
    assert result["rows"]==3168 and result["noise_supported"]==3143
    assert result["arms"]["primary_min"]["audit"]["reference_hits_lost"]==6
    assert result["arms"]["secondary_logistic"]["audit"]["reference_hits_lost"]==8
    assert result["arms"]["primary_min"]["audit"]["weak"]==65
    assert result["arms"]["secondary_logistic"]["audit"]["weak"]==65
    d=prior.npz(r42("all_scores.npz"))
    assert len(d["m"])==3168 and len(d["primary_min"])==3168 and len(d["secondary_logistic"])==3168
    assert d["valid"].sum()==3143 and d["eligible"].sum()==891
    assert np.array_equal(d["group"]<=6,np.isin(d["group"],Q["split"]["fit_groups"]))
    assert np.isfinite(d["m"]).all() and np.isfinite(d["n"]).all()
    audit=d["group"]>=7
    refs=prior.metrics(d["reference"],d,audit)
    assert refs["controls"]==51 and refs["localized"]==318 and refs["weak"]==69
    for arm in ("primary_min","secondary_logistic"):
        actual=prior.metrics(d[arm]>=lock["thresholds"][arm],d,audit)
        assert actual==result["arms"][arm]["audit"],arm
    return d,result,lock
def lost_diagnostic(d,result,lock):
    audit=d["group"]>=7
    reference_hits=d["reference"]&d["y"]&d["loc"]
    reference_nulls=d["reference"]&~d["y"]
    diagnosis={"status":"PAIRED_ROW_DIAGNOSIS_ONLY","source_rows":len(d),
               "audit_reference_hits":int((reference_hits&audit).sum()),
               "audit_reference_weak":int((reference_hits&d["weak"]&audit).sum()),
               "audit_reference_false":int((reference_nulls&audit).sum()),"arms":{}}
    for arm in ("primary_min","secondary_logistic"):
        t=float(lock["thresholds"][arm])
        decision=d[arm]>=t
        lost=np.flatnonzero(audit&reference_hits&~decision)
        removed=np.flatnonzero(audit&reference_nulls&~decision)
        added=np.flatnonzero(audit&~d["reference"]&~d["y"]&decision)
        def describe(i):
            fields={**identity(d,i),"weak":bool(d["weak"][i]),"parent_route":bool(d["route"][i]),
                    "noise_supported":bool(d["valid"][i]),"frozen_noise_eligible":bool(d["eligible"][i]),
                    "minimal_margin":float(d["m"][i]),"noise_margin":float(d["n"][i]),
                    "primary_margin":float(d["primary_min"][i]),
                    "secondary_score":float(d["secondary_logistic"][i]),
                    "score_below_locked_threshold":float(t-d[arm][i])}
            return fields
        result_list=[describe(i) for i in lost]
        assert len(result_list)==result["arms"][arm]["audit"]["reference_hits_lost"]
        assert sum(x["weak"] for x in result_list)==result["arms"][arm]["audit"]["reference_weak_lost"]
        assert len(removed)==result["arms"][arm]["audit"]["reference_false_removed"]
        assert len(added)==result["arms"][arm]["audit"]["reference_false_added"]
        diagnosis["arms"][arm]={"locked_threshold":t,"lost_signals":result_list,
            "removed_nulls":[describe(i) for i in removed],
            "added_nulls":[describe(i) for i in added],
            "lost_by_lane":{"parent":sum(x["parent_route"] for x in result_list),
                            "rescue":sum(not x["parent_route"] for x in result_list)},
            "lost_by_group":{str(g):sum(x["group"]==g for x in result_list) for g in range(7,12)},
            "lost_weak":sum(x["weak"] for x in result_list)}
    save("paired_loss_diagnostic.json",diagnosis)
    print("R42A6_LOSS_DIAGNOSIS",json.dumps({"reference":{k:v for k,v in diagnosis.items() if k!="arms"},
      "arms":{arm:{"lost_by_lane":x["lost_by_lane"],"lost_by_group":x["lost_by_group"],
                   "lost_weak":x["lost_weak"],"removed_nulls":len(x["removed_nulls"])} for arm,x in diagnosis["arms"].items()}}),flush=True)
    return diagnosis
def shape_rows(V,path):
    assert V.ndim==3 and V.shape[1:]==(EPOCHS,DIMS) and np.isfinite(V).all()
    assert path.ndim==2 and path.shape==(len(V),EPOCHS)
    assert np.all(np.isfinite(path)) and np.all(path==path.astype(int))
    p=path.astype(int);assert p.min()>=0 and p.max()<DIMS
    batch=np.arange(len(V))[:,None]
    row=np.arange(EPOCHS)[None,:]
    center=V[batch,row,p].astype(np.float64)
    offs=np.array([-12,-8,-5,-3,3,5,8,12])
    inds=p[:,:,None]+offs[None,None,:]
    within=(inds>=0)&(inds<DIMS)
    side=V[batch[:,:,None],row[:,:,None],np.clip(inds,0,DIMS-1)].astype(np.float64)
    side=np.where(within,side,np.nan)
    assert (within.sum(2)>=4).all()
    bg=np.nanmedian(side,axis=2)
    md=np.nanmedian(abs(side-bg[:,:,None]),axis=2)
    sigma=np.maximum(1.4826*md,0.02)
    contrast=np.clip((center-bg)/sigma,-30.,30.)
    shoulders=np.array([-2,-1,1,2])
    is2=p[:,:,None]+shoulders[None,None,:]
    vs=V[batch[:,:,None],row[:,:,None],np.clip(is2,0,DIMS-1)].astype(np.float64)
    vs=np.where((is2>=0)&(is2<DIMS),vs,np.nan)
    shape={
       "ridge_median_snr":np.median(contrast,axis=1),
       "ridge_q25_snr":np.quantile(contrast,.25,axis=1),
       "ridge_occupancy":(contrast>1).mean(1),
       "ridge_shoulder":np.median(np.clip((center-np.nanmean(vs,axis=2))/sigma,-30,30),axis=1),
       "ridge_jitter":np.median(abs(np.diff(contrast,axis=1)),axis=1)/(1+abs(np.median(contrast,axis=1))),
       "path_edge_frac":((p<8)|(p>=DIMS-8)).mean(1),
       "path_drift":np.mean(abs(np.diff(p,axis=1)),axis=1)/DIMS}
    arr=np.stack([shape[x] for x in FEATURE_NAMES[4:]],axis=1)
    assert arr.shape==(len(V),7) and np.isfinite(arr).all()
    return arr
def make_shape_features(d):
    source=Q["provenance"]
    mr=prior.archive(source["minimal"]);cr=prior.archive(source["common"])
    receipt=json.loads(cr("receipt.json"));common={x["name"]:x for x in receipt["records"]}
    assert len(common)==22 and sum(x["rows"] for x in common.values())==3168
    all_shape=[];contributions=[];start=0
    for g in range(1,12):
        num=f"{g:02d}"
        mg=json.loads(mr("group-"+num+".json"))
        assert mg["group"]==g and mg["fit"]==(g<=6)
        for kind in ("eval","shared"):
            filename=f"{num}-{kind}"
            rec=next(x for x in mg["records"] if x["kind"]==kind)
            source_record=common[filename+".npz"]
            q=prior.npz(prior.checked(cr,filename+".npz",source_record["sha256"]))
            m=prior.npz(prior.checked(mr,filename+"-scores.npz",rec["score_file_sha256"]))
            n=len(m["score"]);stop=start+n
            assert n==rec["rows"]==source_record["rows"] and n==len(q["V"])
            assert q["V"].shape==(n,16,128)
            assert sha(np.ascontiguousarray(q["V"]).tobytes())==source_record["array_sha256"]
            assert np.array_equal(m["source_tile_sha"],q["tile_sha"])
            assert np.array_equal(m["score"],d["m"][start:stop])
            assert np.array_equal(q["y"],d["y"][start:stop])
            assert np.array_equal(q["weak"],d["weak"][start:stop])
            assert np.array_equal(m["parent_route"],d["route"][start:stop])
            assert np.array_equal(d["group"][start:stop],g)
            assert np.array_equal(d["shared"][start:stop],kind=="shared")
            assert np.array_equal(d["row"][start:stop],np.arange(n))
            shape=shape_rows(q["V"],m["path"])
            all_shape.append(shape)
            contributions.append({"group":g,"kind":kind,"rows":n,
                   "frozen_common_sha256":source_record["sha256"],
                   "V_array_sha256":source_record["array_sha256"],
                   "frozen_score_sha256":rec["score_file_sha256"]})
            start=stop
    assert start==3168
    shape=np.concatenate(all_shape,axis=0)
    assert shape.shape==(3168,7) and np.isfinite(shape).all()
    np.savez_compressed(OUT/"morphology_features.npz",features=shape)
    save("feature_sources.json",{"records":contributions,"feature_names":FEATURE_NAMES[4:],
        "formula":"Fixed path-bin median contrast to +/-3,5,8,12 sidebands, robust MAD floor 0.02, clip +/-30, and fixed path diagnostics",
        "independent_instrument":False,"source_revisions_unchanged":True})
    return shape
def model_fit(X,mask,positive):
    assert mask.sum()>0 and positive[mask].sum()>0 and (~positive[mask]).sum()>0
    scaler=StandardScaler().fit(X[mask])
    mod=LogisticRegression(C=.1,class_weight="balanced",random_state=424206,max_iter=2000)
    mod.fit(scaler.transform(X[mask]),positive[mask])
    return scaler,mod
def model_score(m,X):
    scaler,model=m
    result=model.decision_function(scaler.transform(X))
    assert len(result)==len(X) and np.isfinite(result).all()
    return result
def full_audit(d,shape):
    fit=d["group"]<=6;audit=~fit
    target=d["reference"]&d["y"]&d["loc"]
    reference=d["reference"]
    null=reference&~d["y"]
    # Exclude injected positive rows that were not successfully localized from fit.
    train=fit&(target|null)
    X=np.column_stack([d["m"],d["n"],d["valid"].astype(float),d["route"].astype(float),shape]).astype(float)
    assert X.shape==(3168,len(FEATURE_NAMES)) and np.isfinite(X).all()
    cv=[]
    for g in range(1,7):
        tr=train&(d["group"]!=g)
        model=model_fit(X,tr,target)
        raw=model_score(model,X)
        threshold=float(np.nextafter(np.min(raw[tr&target]),-np.inf))
        assert np.isfinite(threshold) and np.all(raw[tr&target]>=threshold)
        eval_group=(d["group"]==g)
        decision=reference&(raw>=threshold)
        mt=prior.metrics(decision,d,eval_group)
        cv.append({"holdout_fit_group":g,"training_reference_hits":int((tr&target).sum()),
                   "training_nulls":int((tr&null).sum()),
                   "threshold":threshold,"metrics":mt})
    fit_model=model_fit(X,train,target)
    raw=model_score(fit_model,X)
    t=float(np.nextafter(np.min(raw[train&target]),-np.inf))
    assert np.isfinite(t) and np.all(raw[train&target]>=t)
    decision=reference&(raw>=t)
    assert np.array_equal(decision[train&target],reference[train&target])
    fit_m=prior.metrics(decision,d,fit)
    audit_m=prior.metrics(decision,d,audit)
    assert fit_m["reference_hits_lost"]==0
    assert ((decision & ~reference)).sum()==0
    lost=audit&target&~decision
    removed=audit&null&~decision
    miss=np.flatnonzero(lost)
    m={"status":"PRESERVATION_FIRST_HISTORICAL_TEST__NOT_INDEPENDENT","features":FEATURE_NAMES,
        "threshold_source":"Minimum fit-group reference hit score, nextafter toward -infinity",
        "fit_lock_threshold":t,"fit":fit_m,"audit":audit_m,
        "audit_goal_met":audit_m["controls"]<=35 and audit_m["localized"]>=318 and audit_m["weak"]>=69 and audit_m["reference_hits_lost"]==0,
        "audit_reference_lost_ids":[{**identity(d,i),"weak":bool(d["weak"][i]),"score":float(raw[i])} for i in miss],
        "audit_reference_removed_control_ids":[identity(d,i) for i in np.flatnonzero(removed)],
        "fit_group_leave_one_out":cv,
        "model":{"C":0.1,"class_weight":"balanced","seed":424206,
                 "coef":fit_model[1].coef_[0].tolist(),"intercept":float(fit_model[1].intercept_[0]),
                 "scaler_mean":fit_model[0].mean_.tolist(),"scaler_scale":fit_model[0].scale_.tolist()},
        "no_audit_parameter_selection":True,"audit_previously_viewed":True,
        "production_adoption":False,"external_validation":False}
    np.savez_compressed(OUT/"guard_scores.npz",score=raw,decision=decision,
        feature_matrix=X,fit=fit,audit=audit,reference=reference)
    save("guard_fit_lock.json",{"source_protocol_sha256":sha(P.read_bytes()),"threshold":t,
        "training_rows":int(train.sum()),"positive_training":int((train&target).sum()),
        "null_training":int((train&null).sum()),"model":m["model"],"no_audit_used":True})
    print("R42A6_GUARD_FIT_LOCK",json.dumps({"threshold":t,"fit":fit_m,"cv":[{ "group":a["holdout_fit_group"],
      "lost":a["metrics"]["reference_hits_lost"],"controls":a["metrics"]["controls"]} for a in cv]}),flush=True)
    print("R42A6_GUARD_AUDIT",json.dumps({"audit":audit_m,"goal_met":m["audit_goal_met"],
      "lost_ids":m["audit_reference_lost_ids"]}),flush=True)
    return m
def main():
    d,original,lock=load_history()
    diagnosis=lost_diagnostic(d,original,lock)
    shape=make_shape_features(d)
    guard=full_audit(d,shape)
    final={"stage":"R42A6","status":"PASS_VERIFIED_FEATURE_EXTRACTION__HOLD_EXTERNAL_GENERALIZATION",
      "r42a5_primary":original["arms"]["primary_min"]["audit"],"r42a5_secondary":original["arms"]["secondary_logistic"]["audit"],
      "reference":{"controls":51,"localized":318,"weak":69},
      "forensic_counts":{arm:{"lost":len(v["lost_signals"]),"weak":v["lost_weak"],
                            "removed_false":len(v["removed_nulls"])} for arm,v in diagnosis["arms"].items()},
      "new_guard":{k:v for k,v in guard.items() if k not in ("model","audit_reference_lost_ids","audit_reference_removed_control_ids","fit_group_leave_one_out")},
      "new_guard_lost_ids":guard["audit_reference_lost_ids"],
      "production_modified":False,"external_validated":False,
      "scientific_caveat":"These are simulated injected trajectories, same registered spectra and previously inspected audit groups. No prospective radio observation or external custodian test."}
    save("result.json",final)
    print("R42A6_SUMMARY",json.dumps(final),flush=True)
if __name__=="__main__":main()
