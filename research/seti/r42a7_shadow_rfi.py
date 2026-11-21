"""R42A7A: Non-destructive instrumental RFI morphology review, fit-only CV.
No veto, no production edits, no scoring on previously inspected audit groups 7–11.
"""
import hashlib,io,json,pathlib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import r42a6_signal_preservation as a6
import r42a5_combined_score as prior

P=pathlib.Path("research/seti/r42a7_shadow_rfi_protocol.json")
Q=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r42a7");OUT.mkdir(parents=True,exist_ok=True)
FEATURES=Q["new_feature_family"]["features"]
assert len(FEATURES)==8
def sha(x):return hashlib.sha256(x).hexdigest()
def save(name,v):(OUT/name).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
def features(V,path):
    """Independent transformation of raw V; NEVER receives injected labels or true path.
    Sideband/local path ridge scores from R42A6 are NOT inputs to the new model.
    """
    assert V.ndim==3 and V.shape[1:]==(16,128) and np.isfinite(V).all()
    assert path.shape==(len(V),16) and np.isfinite(path).all()
    p=path.astype(np.int32)
    assert np.array_equal(p,path) and p.min()>=0 and p.max()<128
    n=len(p);t=np.arange(16)[None,:];sample=np.arange(n)[:,None]
    a=V.astype(np.float64)
    med=np.median(a,axis=2)
    sig=np.maximum(1.4826*np.median(abs(a-med[:,:,None]),axis=2),0.02)
    normalized=np.clip((a-med[:,:,None])/sig[:,:,None],-30,30)
    grid=np.arange(128)[None,None,:]
    off=(abs(grid-p[:,:,None])>4)
    safe=np.where(off,normalized,-1e9)
    broad3=((normalized>3)&off).sum(axis=(1,2))/np.maximum(off.sum(axis=(1,2)),1)
    broad5=((normalized>5)&off).sum(axis=(1,2))/np.maximum(off.sum(axis=(1,2)),1)
    commonmode=np.std(med,axis=1)/np.maximum(np.median(sig,axis=1),0.02)
    center=normalized[sample,t,p]
    strongest_off=safe.max(axis=2)
    ratio=(1+np.maximum(0,np.median(strongest_off,axis=1)))/(1+np.maximum(0,np.median(center,axis=1)))
    # 8 fixed off-path offsets, each with two fixed drift hypotheses; 16 comparisons total.
    offsets=(-32,-24,-16,-8,8,16,24,32)
    alternatives=[]
    for offset in offsets:
        for slope in (0,1 if offset>0 else -1):
            ix=p+offset+slope*(t-8)
            within=(ix>=0)&(ix<128)
            v=normalized[sample,t,np.clip(ix,0,127)]
            # No extrapolated spectral samples outside original 128 bins.
            v=np.where(within,v,np.nan)
            count=within.sum(axis=1)
            alt=np.nanmedian(v,axis=1)
            alt=np.where(count>=12,alt,-30.)
            alternatives.append((alt,v,within))
    strongest=np.maximum.reduce([np.nan_to_num(y[0],nan=-30) for y in alternatives])
    alternative_ratio=(1+np.maximum(0,strongest))/(1+np.maximum(0,np.median(center,axis=1)))
    max_each=np.maximum.reduce([np.where(z[2],np.nan_to_num(z[1],nan=-30),-30) for z in alternatives])
    win_fraction=np.mean(max_each>center,axis=1)
    # Far frequency windows as a fixed common-mode control, excluding 4 bins around path.
    low=np.where(off[:,:,:48],a[:,:,:48],np.nan)
    high=np.where(off[:,:,80:],a[:,:,80:],np.nan)
    lo=np.nanmedian(low,axis=2);hi=np.nanmedian(high,axis=2)
    lo=np.nan_to_num(lo,nan=0);hi=np.nan_to_num(hi,nan=0)
    ld=lo-lo.mean(axis=1,keepdims=True);hd=hi-hi.mean(axis=1,keepdims=True)
    co=(ld*hd).sum(axis=1)/np.maximum(np.sqrt((ld*ld).sum(axis=1)*(hd*hd).sum(axis=1)),1e-8)
    persistent=np.mean(np.sum((normalized>3)&off,axis=1)>=8,axis=1)
    x=np.column_stack([broad3,broad5,commonmode,np.clip(ratio,0,40),np.clip(alternative_ratio,0,40),win_fraction,
        np.clip(co,-1,1),persistent])
    assert x.shape==(n,8) and np.isfinite(x).all()
    return x

def verified_fit_features(d):
    # Only read groups 1–6. Raw audit groups 7–11 are never loaded by this stage.
    src=a6.Q["provenance"]
    mr=prior.archive(src["minimal"]);cr=prior.archive(src["common"])
    receipt=json.loads(cr("receipt.json"));common={x["name"]:x for x in receipt["records"]}
    fit=d["group"]<=6
    df={k:arr[fit] for k,arr in d.items()}
    assert len(df["m"])==1728
    chunks=[];refs=[];idx=0
    for g in range(1,7):
        num=f"{g:02d}";mg=json.loads(mr("group-"+num+".json"));assert mg["fit"] and mg["group"]==g
        for kind in ("eval","shared"):
            fn=num+"-"+kind
            a=next(x for x in mg["records"] if x["kind"]==kind)
            r=common[fn+".npz"]
            z=prior.npz(prior.checked(cr,fn+".npz",r["sha256"]))
            scores=prior.npz(prior.checked(mr,fn+"-scores.npz",a["score_file_sha256"]))
            n=len(scores["score"]);end=idx+n
            assert n==a["rows"]==r["rows"] and z["V"].shape==(n,16,128)
            assert sha(np.ascontiguousarray(z["V"]).tobytes())==r["array_sha256"]
            assert np.array_equal(scores["source_tile_sha"],z["tile_sha"])
            assert np.array_equal(scores["score"],df["m"][idx:end])
            assert np.array_equal(z["y"],df["y"][idx:end])
            assert np.array_equal(z["weak"],df["weak"][idx:end])
            assert np.array_equal(df["row"][idx:end],np.arange(n))
            assert np.all(df["group"][idx:end]==g)
            # Remove truth-derived fields BEFORE the pure feature function.
            x=features(z["V"],scores["path"])
            chunks.append(x);refs.append({"group":g,"kind":kind,"rows":n,"common_file_sha256":r["sha256"],
                "common_V_sha256":r["array_sha256"],"minimal_path_sha256":a["score_file_sha256"]})
            idx=end
    assert idx==1728
    X=np.concatenate(chunks,axis=0);assert X.shape==(1728,8) and np.isfinite(X).all()
    save("frozen_feature_receipt.json",{"input_sha256s":refs,"features":FEATURES,
        "source_rows":1728,"loaded_groups":[1,2,3,4,5,6],"audit_groups_read":False,
        "uses_truth_in_extractor":False,"independent_instrument":False})
    np.savez_compressed(OUT/"feature_archive_fit_only.npz",features=X,groups=df["group"],row=df["row"],
        original_reference_alert=df["reference"])
    return df,X
def model(X,mask,y):
    assert mask.sum()>0 and 0<y[mask].sum()<mask.sum()
    scaler=StandardScaler().fit(X[mask])
    reg=LogisticRegression(C=.1,class_weight="balanced",random_state=424207,max_iter=2000)
    reg.fit(scaler.transform(X[mask]),y[mask].astype(int))
    return scaler,reg
def score(model,X):
    s,reg=model
    z=reg.decision_function(s.transform(X))
    assert np.isfinite(z).all()
    return z
def main():
    # A6 performs hash/score parity and checks old historical losses, but this stage
    # NEVER selects, scores, fits, reports or optimizes on old audit 7–11.
    original,previous,lock=a6.load_history()
    d,X=verified_fit_features(original)
    ref=d["reference"];positive=ref&d["y"]&d["loc"];null=ref&~d["y"]
    good=positive|null
    assert np.any(positive) and np.any(null)
    assert all(int(g) in range(1,7) for g in np.unique(d["group"]))
    cv=[];queue=[];total_ctrl=0;total_hits=0;total_unknown=0
    for group in range(1,7):
        tr=(d["group"]!=group)&good
        mod=model(X,tr,positive)
        pred=score(mod,X)
        holdout=(d["group"]==group)&ref
        indices=np.flatnonzero(holdout)
        selection=indices[np.argsort(pred[indices],kind="stable")[:5]]
        audit_known=(d["group"]==group)&good
        if (audit_known&positive).sum() and (audit_known&null).sum():
            auc=float(roc_auc_score(null[audit_known].astype(int),-pred[audit_known]))
        else:auc=None
        top_false=int(null[selection].sum());top_hit=int(positive[selection].sum())
        other=int(len(selection)-top_false-top_hit)
        total_ctrl+=top_false;total_hits+=top_hit;total_unknown+=other
        cv.append({"fit_holdout_group":group,"training_groups":[i for i in range(1,7) if i!=group],
            "training_reference_hits":int((tr&positive).sum()),"training_false_alerts":int((tr&null).sum()),
            "review_queue":int(len(selection)),"review_queue_actual_controls":top_false,
            "review_queue_actual_localized_hits":top_hit,"review_queue_other_injected":other,
            "holdout_true_control_rows":int((audit_known&null).sum()),
            "holdout_true_localized_hits":int((audit_known&positive).sum()),
            "auc_false_vs_localized_research_only":auc,"suppressed_alerts":0})
        for i in selection:
            queue.append({"fit_holdout_group":group,"kind":"shared" if d["shared"][i] else "eval",
                "row":int(d["row"][i]),"shadow_false_alarm_score":float(-pred[i]),
                "label_for_post_fit_audit_only":"FALSE_CONTROL" if null[i] else "LOCALIZED_INJECTION" if positive[i] else "UNLOCALIZED_INJECTION",
                "actually_suppressed":False})
    full=model(X,good,positive)
    scaler,reg=full
    frozen={"source_protocol_sha256":sha(P.read_bytes()),"fit_groups":[1,2,3,4,5,6],
        "old_audit_groups_train_or_test":[],"features":FEATURES,"scaler_mean":scaler.mean_.tolist(),
        "scaler_scale":scaler.scale_.tolist(),"coefficients":reg.coef_[0].tolist(),
        "intercept":float(reg.intercept_[0]),"hyperparams":Q["model"],"no_threshold_selected":True,
        "hard_veto_disabled":True,"production_unchanged":True}
    save("fit_only_cv.json",{"folds":cv,"summary":{"reviewed":len(queue),"correctly_flagged_false_controls":total_ctrl,
        "also_flagged_correctly_localized_injections":total_hits,"other_injected":total_unknown,
        "valid_folds":sum(x["auc_false_vs_localized_research_only"] is not None for x in cv),
        "mean_fold_auc":float(np.mean([x["auc_false_vs_localized_research_only"] for x in cv if x["auc_false_vs_localized_research_only"] is not None])),
        "actual_alerts_suppressed":0}})
    save("shadow_review_queue.json",{"queue":queue,"review_only":True,"not_selected_for_production":True})
    save("fit_model_lock.json",frozen)
    fit_original=prior.metrics(ref,d,np.ones(len(ref),bool))
    result={"stage":"R42A7A","status":"FIT_ONLY_SHADOW_FEATURE_EXPERIMENT__HOLD_NEW_PROSPECTIVE_DATA",
        "fit_reference":{"controls":fit_original["controls"],"localized":fit_original["localized"],
            "weak":fit_original["weak"]},
        "fit_only_cv":{"ranked_30_for_review":len(queue),"false_controls_in_queue":total_ctrl,
            "localized_injections_in_queue":total_hits,"other_injected_in_queue":total_unknown},
        "reference_decisions_unchanged":True,"reference_decision_diff_count":0,
        "historical_audit_groups_scored":False,
        "historical_audit_results_unchanged":{"reference":{"controls":51,"localized":318,"weak":69},
            "prior_R42A6_research":{"controls":45,"localized":312,"weak":68}},
        "does_not_claim_new_35_FA":True,"external_validation":False,"new_holdout_required":True,
        "fit_lock_sha256":sha((OUT/"fit_model_lock.json").read_bytes()),"production_adoption":False}
    save("result.json",result)
    print("R42A7_FIT_SHADOW",json.dumps(result,ensure_ascii=False),flush=True)
    print("R42A7_CV",json.dumps(cv,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
