#!/usr/bin/env python3
"""R42A8A: frozen R42A7A research score tested without retraining on new
synthetic stress spectrograms. NOT a full JANUS detector or real telescope test.
"""
import hashlib,json,pathlib
import numpy as np
from sklearn.metrics import roc_auc_score
import r42a5_combined_score as source
import r42a7_shadow_rfi as prev

P=pathlib.Path("research/seti/r42a8_ood_stress_protocol.json")
Q=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r42a8");OUT.mkdir(parents=True,exist_ok=True)
T=16;F=128;TIX=np.arange(T)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def frozen_model():
    artifact=Q["frozen_model"]
    read=source.archive({"artifact":artifact["artifact"],"run":artifact["run"],"digest":artifact["sha256"]})
    raw=read("fit_model_lock.json")
    if sha(raw)!=artifact["fit_lock_sha256"]:raise ValueError("R42A7_FROZEN_LOCK_HASH_MISMATCH")
    lock=json.loads(raw)
    assert lock["features"]==prev.FEATURES and lock["hard_veto_disabled"] and lock["no_threshold_selected"]
    assert lock["old_audit_groups_train_or_test"]==[]
    assert lock["fit_groups"]==[1,2,3,4,5,6]
    mean=np.array(lock["scaler_mean"],dtype=float);scale=np.array(lock["scaler_scale"],dtype=float)
    coef=np.array(lock["coefficients"],dtype=float)
    assert len(mean)==len(scale)==len(coef)==8 and np.isfinite(mean).all() and (scale>0).all()
    assert np.isfinite(coef).all()
    save("model_lock_verification.json",{"source_run":artifact["run"],"source_artifact":artifact["artifact"],
      "artifact_zip_sha256":artifact["sha256"],"fit_lock_sha256":sha(raw),
      "features":prev.FEATURES,"model_parameters_unchanged":True,"no_training":True,
      "historical_7_to_11_data_access":False})
    return mean,scale,coef,float(lock["intercept"])
def add_line(v,path,amp,active,shoulder=.24):
    for i in range(T):
        if not active[i]:continue
        col=int(path[i])
        if 0<=col<F:v[i,col]+=float(amp[i])
        for j in (-1,1):
            if 0<=col+j<F:v[i,col+j]+=float(shoulder*amp[i])
def example(rng,regime,positive,family):
    # Separate from R41/R42 generator. This is a simulation sanity check only.
    v=rng.standard_normal((T,F)).astype(np.float32)
    if regime==0:
        v=(.68*v+.23*np.roll(v,1,axis=1)+.09*np.roll(v,2,axis=1)).astype(np.float32)
        cm_amp=.22
    else:
        v=(.50*v+.27*np.roll(v,2,axis=1)+.23*np.roll(v,6,axis=1)).astype(np.float32)
        cm_amp=.75
    v+=rng.normal(0,.12,(T,1)).astype(np.float32)
    v+=float(rng.uniform(.05,cm_amp))*np.sin(np.linspace(0,2*np.pi,T)+rng.uniform(-np.pi,np.pi))[:,None]
    v+=rng.normal(0,.15,F).astype(np.float32)[None,:]
    truth=np.full(T,-1,np.int32);active=np.zeros(T,bool);weak=False
    if positive:
        start=int(rng.integers(34,95));slope=float(rng.uniform(-1.10,1.10))
        curve=float(rng.uniform(-.05,.05)) if family=="curved_drift" else 0.0
        path=np.rint(start+slope*(TIX-7.5)+curve*(TIX-7.5)**2).astype(np.int32)
        active[:]=True
        if family=="intermittent_line":active=rng.random(T)>.32
        if family=="weak_drift":
            amp=rng.uniform(1.4,2.2,T);weak=True
        else:amp=rng.uniform(2.1,3.7,T)
        if family=="modulated_line":amp=amp*(1+.42*np.sin(TIX*.85+1))
        truth[:]=path
        add_line(v,path,amp,active)
    else:
        if family=="broadband_bursts":
            for q in rng.choice(T,3,replace=False):
                v[q,:]+=rng.uniform(1.2,2.9)*(rng.random(F)<.65)
        elif family=="frequency_comb":
            off=int(rng.integers(1,11))
            for q in range(off,F,11):v[:,q]+=float(rng.uniform(1.3,3.6))
        elif family=="competing_rfi":
            p0=int(rng.integers(32,100))
            for slope in (-.7,.7):
                p=np.rint(p0+slope*(TIX-7.5)).astype(int)
                add_line(v,p,rng.uniform(1.3,2.7,T),np.ones(T,bool),shoulder=.55)
        else:
            # colored_noise: deliberately no added narrowband source.
            pass
    return v,truth,active,weak
def batch(rng):
    arr=[];truth=[];active=[];labels=[];group=[];family=[];weak=[]
    posfamilies=Q["panel"]["positive_families"]
    negfamilies=Q["panel"]["null_families"]
    for regime in range(2):
        for y in (0,1):
            for k in range(256):
                fam=(posfamilies if y else negfamilies)[k%4]
                v,p,act,w=example(rng,regime,bool(y),fam)
                arr.append(v);truth.append(p);active.append(act);labels.append(y)
                group.append(regime);family.append(fam);weak.append(w)
    data={"V":np.stack(arr),"truth":np.stack(truth),"active":np.stack(active),
          "label":np.asarray(labels,dtype=bool),"group":np.asarray(group,dtype=int),
          "family":np.asarray(family,dtype="U32"),"weak":np.asarray(weak,dtype=bool)}
    assert len(arr)==1024 and sum(data["label"])==512 and np.isfinite(data["V"]).all()
    return data
def blind_line_candidate(V):
    """Propose same label-blind candidate for positives and controls, not injected truth."""
    norm=V.astype(np.float64)
    med=np.median(norm,axis=2)
    scale=np.maximum(1.4826*np.median(abs(norm-med[:,:,None]),axis=2),.02)
    norm=(norm-med[:,:,None])/scale[:,:,None]
    n=len(V);best=np.full(n,-np.inf);paths=np.zeros((n,T),np.int32)
    nidx=np.arange(n)
    starts=np.arange(16,112,dtype=int)
    for slope in (-1.5,-1.,-.5,0.,.5,1.,1.5):
        pp=np.rint(starts[:,None]+slope*(TIX[None,:]-7.5)).astype(int)
        assert pp.min()>=0 and pp.max()<128
        scores=norm[:,np.arange(T)[None,None,:],pp[None,:,:]].mean(axis=2)
        idx=np.argmax(scores,axis=1)
        val=scores[nidx,idx]
        ok=val>best
        best[ok]=val[ok]
        paths[ok]=pp[idx[ok]]
    assert np.isfinite(best).all() and paths.min()>=0 and paths.max()<128
    return paths,best
def calc_auc(label,s):
    assert len(np.unique(label))==2
    return float(roc_auc_score(label,s))
def main():
    mean,scale,coef,intercept=frozen_model()
    rng=np.random.default_rng(int(Q["panel"]["seed"]))
    d=batch(rng)
    path,ridge=blind_line_candidate(d["V"])
    # Explicit no-label dependency: only spectra+candidate paths enter feature transform.
    X=prev.features(d["V"],path)
    frozen=((X-mean)/scale)@coef+intercept
    assert len(frozen)==len(ridge)==1024 and np.isfinite(frozen).all()
    # Positive-injection localization is an evaluation-only metric.
    good=np.abs(path-d["truth"])<=2
    hit=((good & d["active"]).sum(axis=1)>=.75*np.maximum(1,d["active"].sum(axis=1)))&d["label"]
    regimes=[]
    for k in (0,1):
        m=d["group"]==k
        target=m&d["label"];negative=m&~d["label"]
        regimes.append({"name":Q["panel"]["regimes"][k],"rows":int(m.sum()),
          "injected_positives":int(target.sum()),"simulated_controls":int(negative.sum()),
          "R42A7_frozen_shadow_auc":calc_auc(d["label"][m],frozen[m]),
          "untrained_ridge_energy_auc":calc_auc(d["label"][m],ridge[m]),
          "blind_path_localized":int(hit[target].sum()),
          "blind_path_localization_fraction":float(hit[target].mean()),
          "weak_injections":int((target&d["weak"]).sum()),
          "weak_blind_path_localized":int((target&d["weak"]&hit).sum()),
          "families":[{"family":str(f),"label":("POSITIVE_SIMULATED" if d["label"][m&(d["family"]==f)][0] else "CONTROL_SIMULATED"),
              "rows":int((m&(d["family"]==f)).sum()),
              "frozen_shadow_median_score":float(np.median(frozen[m&(d["family"]==f)]))}
             for f in sorted(np.unique(d["family"][m]))]})
    # Frequency-translation diagnostic in same original observations. Never calibrate on it.
    pplus=path+8
    safe=(pplus.max(axis=1)<128)
    if safe.sum()<400:raise ValueError("FREQ_SHIFT_SAFE_SAMPLE_TOO_SMALL")
    shifted=np.roll(d["V"][safe],8,axis=2)
    Xshift=prev.features(shifted,pplus[safe])
    shifted_score=((Xshift-mean)/scale)@coef+intercept
    delta=np.abs(shifted_score-frozen[safe])
    shift={"tested":int(safe.sum()),"shift_bins":8,
           "median_frozen_score_abs_delta":float(np.median(delta)),
           "p95_frozen_score_abs_delta":float(np.quantile(delta,.95)),
           "changed_model_rank_polarity":int(((shifted_score>0)!=(frozen[safe]>0)).sum()),
           "property_is_qualification_only":True}
    save("invariance_probe.json",shift)
    potential=all(x["R42A7_frozen_shadow_auc"]>=.7 and x["R42A7_frozen_shadow_auc"]>=x["untrained_ridge_energy_auc"] for x in regimes)
    result={"stage":"R42A8A",
       "status":"PASS_AUTHENTICATED_NEW_SYNTHETIC_STRESS__GENERALIZATION_POTENTIAL_ONLY" if potential else
           "PASS_AUTHENTICATED_NEW_SYNTHETIC_STRESS__FAIL_GENERALIZATION_QUALIFICATION",
       "source":"New internal simulation of time-frequency spectra, NOT real telescope ON/OFF or independent custodian.",
       "model":"Frozen R42A7A shadow scorer (no retraining), compared to simple blind ridge-energy heuristic",
       "candidate_path":"Same blind linear exhaustive search for positive and control simulated spectra",
       "new_simulated_rows":1024,"regimes":regimes,"frequency_shift":shift,
       "frozen_shadow_qualified_on_both_regimes":bool(potential),
       "historical_fit_data_loaded":False,"historical_audit_data_loaded":False,
       "actual_signal_alerts_suppressed":0,"actual_false_alarm_reduction":0,
       "R42_reference_unchanged":{"false_alerts":51,"localized":318,"weak":69},
       "full_JANUS_detector_not_rerun":True,"no_native_radio_verification":True,
       "external_blind_validation":False,"production_adoption":False,
       "scientific_next":"Implement native and same-information ON/OFF feature lane on fresh verified GBT spectra once original antenna pointing/custodian available. New genuine groups need sealed labels; no old audit retuning."}
    np.savez_compressed(OUT/"r42a8_scores.npz",score=frozen,baseline_ridge=ridge,
          simulated_label=d["label"],simulated_family=d["family"],ood_group=d["group"],
          candidate_path=path,truth_for_evaluation_only=d["truth"],localized=hit)
    save("result.json",result)
    print("R42A8_RESULT",json.dumps(result),flush=True)
if __name__=="__main__":main()
