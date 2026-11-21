"""JANUS R41S1H3 RETINA: BLC1-inspired physical RFI-family features.
All RFI origins are simulated hypotheses; historical 38 parent alerts never enter.
No veto, no H2 threshold changes, no independent telescope classification.
"""
import hashlib,json,pathlib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

P=pathlib.Path("research/seti/r41s1h3_retina_rfi_protocol.json")
Q=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r41s1h3");OUT.mkdir(parents=True,exist_ok=True)
T,F=16,256
TI=np.arange(T);FEATS=Q["feature_spec"]["features"]
assert len(FEATS)==6

def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def add_track(v,p,amps,active,spread=.22):
    for t in range(T):
        if not active[t]:continue
        j=int(p[t])
        if 0<=j<F:
            v[t,j]+=amps[t]
            if j>0:v[t,j-1]+=spread*amps[t]
            if j+1<F:v[t,j+1]+=spread*amps[t]
def simulate(seed):
    r=np.random.default_rng(seed)
    records=[]
    pos=Q["synthetic_simulators"]["positive_families"]
    neg=Q["synthetic_simulators"]["negative_families"]
    for label in (True,False):
        for family in (pos if label else neg):
            for k in range(32):
                regime=k%2
                v=r.normal(0,1,(T,F)).astype(np.float32)
                if regime:
                    # Heavy-tailed/non-stationary electronics background.
                    v=(.58*v+.23*np.roll(v,3,axis=1)+.19*np.roll(v,7,axis=1)).astype(np.float32)
                    v+=r.standard_t(4,(T,F)).astype(np.float32)*.17
                    gain=r.uniform(.5,1.7,T)
                else:
                    v=(.78*v+.16*np.roll(v,1,axis=1)+.06*np.roll(v,2,axis=1)).astype(np.float32)
                    gain=np.ones(T)
                v*=gain[:,None]
                v+=r.normal(0,.22,T).astype(np.float32)[:,None]
                v+=r.normal(0,.08,F).astype(np.float32)[None,:]
                start=int(r.integers(72,181))
                slope=float(r.uniform(-1.15,1.15))
                phase=float(r.uniform(0,2*np.pi))
                path=np.rint(start+slope*(TI-7.5)).astype(np.int32)
                active=np.ones(T,bool)
                weak=False
                if label:
                    if family=="faint_single":amp=r.uniform(1.4,2.0,T);weak=True
                    elif family=="intermittent_single":
                        amp=r.uniform(2.8,4.1,T)
                        active=(r.random(T)>.25)
                    else:amp=r.uniform(2.4,4.2,T)
                    if family=="curved_single":
                        path=np.rint(start+slope*(TI-7.5)+.07*(TI-7.5)**2).astype(np.int32)
                    elif family=="drifting_single":
                        amp*=1+.22*np.sin(TI*.7+phase)
                    add_track(v,path,amp,active)
                elif family=="clock_sideband_comb":
                    d=int(r.choice([8,12,16,24]))
                    amp=r.uniform(2.5,3.4,T)*(1+.30*np.sin(TI*.6+phase))
                    add_track(v,path,amp,active)
                    for off in (-2*d,-d,d,2*d):
                        add_track(v,path+off,amp*r.uniform(.65,1.2),active)
                elif family=="intermod_triplet":
                    d=int(r.choice([8,12,16,24]))
                    amp=r.uniform(2.6,3.5,T)*(1+.27*np.sin(TI*.5+phase))
                    add_track(v,path,amp,active)
                    # f0 = 2*f1 - f2 when f1=f0+d and f2=f0+2d.
                    sign=int(r.choice([-1,1]))
                    add_track(v,path+sign*d,amp*.85,active)
                    add_track(v,path+sign*2*d,amp*.78,active)
                elif family=="broadband_impulses":
                    for t in r.choice(T,4,replace=False):
                        v[int(t)]+=r.uniform(1.8,3.3)*(r.random(F)<.7)
                else:
                    # Deliberately unidentifiable from a lone astrophysical-looking line.
                    amp=r.uniform(1.7,3.8,T)
                    add_track(v,path,amp,active)
                records.append((v,label,family,seed,regime,path.copy() if label else np.full(T,-1,np.int32),active.copy() if label else np.zeros(T,bool),weak))
    data={"V":np.stack([x[0] for x in records]),
          "positive":np.array([x[1] for x in records],bool),
          "family":np.array([x[2] for x in records],"U32"),
          "seed":np.array([x[3] for x in records],int),
          "regime":np.array([x[4] for x in records],int),
          "truth_eval_only":np.stack([x[5] for x in records]),
          "active_eval_only":np.stack([x[6] for x in records]),
          "weak_eval_only":np.array([x[7] for x in records],bool)}
    assert data["V"].shape==(256,T,F) and np.isfinite(data["V"]).all()
    return data

def blind_candidate(V):
    # Same search for every control and positive; no ground truth path input.
    n=len(V)
    z=V.astype(np.float64)
    med=np.median(z,axis=2)
    sig=np.maximum(1.4826*np.median(abs(z-med[:,:,None]),axis=2),.1)
    z=(z-med[:,:,None])/sig[:,:,None]
    starts=np.arange(48,208,2)
    path=np.zeros((n,T),np.int32)
    best=np.full(n,-np.inf)
    ar=np.arange(n)
    for slope in (-1.25,-.75,-.25,0.,.25,.75,1.25):
        lines=np.rint(starts[:,None]+slope*(TI-7.5)).astype(int)
        assert 0<=lines.min() and lines.max()<F
        samples=z[ar[:,None,None],TI[None,None,:],lines[None,:,:]]
        score=samples.mean(axis=2)
        ix=np.argmax(score,axis=1)
        value=score[ar,ix]
        improve=value>best
        best[improve]=value[improve]
        path[improve]=lines[ix[improve]]
    assert np.isfinite(best).all()
    return path,best

def physical_features(V,path):
    """Six label-blind candidate-relative diagnostics. No physical origin labels."""
    assert V.ndim==3 and V.shape[1:]==(T,F)
    assert path.shape==(len(V),T) and path.min()>=0 and path.max()<F
    x=V.astype(np.float64)
    med=np.median(x,axis=2)
    sig=np.maximum(1.4826*np.median(abs(x-med[:,:,None]),axis=2),.10)
    z=np.clip((x-med[:,:,None])/sig[:,:,None],-20,20)
    n=len(V);ar=np.arange(n)[:,None];ts=TI[None,:]
    def sample(off):
        p=path+off
        valid=(p>=0)&(p<F)
        return np.where(valid,z[ar,ts,np.clip(p,0,F-1)],-20.)
    cand=sample(0)
    sidevals={off:sample(off) for d in (8,12,16,24) for off in (-2*d,-d,d,2*d)}
    medvals={off:np.median(v,axis=1) for off,v in sidevals.items()}
    clock=[]
    closure=[]
    for d in (8,12,16,24):
        clock.append(np.min(np.stack([medvals[-2*d],medvals[-d],medvals[d],medvals[2*d]]),axis=0))
        closure.append(np.min(np.stack([medvals[d],medvals[2*d]]),axis=0))
        closure.append(np.min(np.stack([medvals[-d],medvals[-2*d]]),axis=0))
    sym=np.max(np.stack(clock),axis=0)
    intermod=np.max(np.stack(closure),axis=0)
    tracks=sum((m>1.2).astype(float) for m in medvals.values())
    # Raw common-envelope coherence among off-path alternatives. Correlation alone
    # is not a source-origin identification; family templates only.
    cand0=cand-cand.mean(axis=1,keepdims=True)
    correlations=[]
    for v in sidevals.values():
        v0=v-v.mean(axis=1,keepdims=True)
        num=(v0*cand0).sum(axis=1)
        den=np.maximum(np.sqrt((v0*v0).sum(axis=1)*(cand0*cand0).sum(axis=1)),1e-5)
        correlations.append(num/den)
    coherence=np.max(np.stack(correlations),axis=0)
    off=np.abs(np.arange(F)[None,None,:]-path[:,:,None])>4
    occupancy=((z>3)&off).sum(axis=(1,2))/np.maximum(off.sum(axis=(1,2)),1)
    common=np.std(med,axis=1)/np.maximum(np.median(sig,axis=1),.1)
    feat=np.column_stack([sym,intermod,tracks,coherence,occupancy,common])
    assert feat.shape==(n,6) and np.isfinite(feat).all()
    return feat

def assemble(seeds):
    chunks=[simulate(int(s)) for s in seeds]
    return {k:np.concatenate([d[k] for d in chunks],axis=0) for k in chunks[0]}

def model_fit(X,positive):
    scaler=StandardScaler().fit(X)
    model=LogisticRegression(C=.1,class_weight="balanced",random_state=411103,max_iter=2000)
    model.fit(scaler.transform(X),(~positive).astype(int))
    return scaler,model
def predict(model,X):
    scaler,logit=model
    z=logit.decision_function(scaler.transform(X))
    assert np.isfinite(z).all()
    return z
def key_metrics(d,path):
    true=d["positive"]
    active=d["active_eval_only"]
    match=np.abs(path-d["truth_eval_only"])<=2
    localized=true&((match&active).sum(axis=1)>=.75*np.maximum(1,active.sum(axis=1)))
    return localized

def main():
    dev=assemble(Q["synthetic_simulators"]["independent_seeds_dev"])
    held=assemble(Q["synthetic_simulators"]["independent_seeds_heldout"])
    first_seed=Q["synthetic_simulators"]["independent_seeds_dev"][0]
    repeated=simulate(int(first_seed))
    check=(sha(np.ascontiguousarray(dev["V"][:256]).tobytes())==sha(np.ascontiguousarray(repeated["V"]).tobytes()))
    assert check and len(dev["V"])==768 and len(held["V"])==512
    source_digests={str(first_seed):sha(repeated["V"].tobytes()),"heldout_V_sha256":sha(held["V"].tobytes())}
    dpath,denergy=blind_candidate(dev["V"])
    hpath,henergy=blind_candidate(held["V"])
    Xd=physical_features(dev["V"],dpath)
    Xh=physical_features(held["V"],hpath)
    mod=model_fit(Xd,dev["positive"])
    raw_dev=predict(mod,Xd);raw_held=predict(mod,Xh)
    lock=float(np.nextafter(np.max(raw_dev[dev["positive"]]),np.inf))
    assert (~(raw_dev[dev["positive"]]>lock)).all()
    save("model_lock.json",{"source_protocol_sha256":sha(P.read_bytes()),"scaler_mean":mod[0].mean_.tolist(),
      "scaler_scale":mod[0].scale_.tolist(),"coef":mod[1].coef_[0].tolist(),
      "intercept":float(mod[1].intercept_[0]),"rfi_suspicion_threshold_fit_only":lock,
      "threshold_source":"above highest dev positive RFI score to preserve all dev positives",
      "no_historical_labels":True,"shadow_only":True,"no_production_veto":True})
    folds=[]
    for seed in Q["synthetic_simulators"]["independent_seeds_heldout"]:
        m=held["seed"]==seed
        gt=(~held["positive"][m]).astype(int)
        baseline=Xh[m,0]
        hypothetical=(raw_held[m]>lock)
        folds.append({"seed":seed,"rows":int(m.sum()),"positive_simulated":int(held["positive"][m].sum()),
          "null_simulated":int((~held["positive"][m]).sum()),
          "physical_rfi_auc":float(roc_auc_score(gt,raw_held[m])),
          "simple_clock_strength_auc":float(roc_auc_score(gt,baseline)),
          "hypothetical_false_control_removals":int((hypothetical&(~held["positive"][m])).sum()),
          "hypothetical_real_simulated_positive_losses":int((hypothetical&held["positive"][m]).sum()),
          "actual_suppressed":0})
    hit=key_metrics(held,hpath)
    family=[]
    for f in list(Q["synthetic_simulators"]["positive_families"])+list(Q["synthetic_simulators"]["negative_families"]):
        m=held["family"]==f
        family.append({"family":f,"simulated_positive":bool(held["positive"][m][0]),"rows":int(m.sum()),
          "median_rfi_score":float(np.median(raw_held[m])),"hypothetical_rejections":int((raw_held[m]>lock).sum()),
          "blind_linear_path_localized":int(hit[m].sum()) if held["positive"][m][0] else None})
    shifts=[]
    for delta in (-16,-8,8,16):
        # Noncircular translation: new edge bins = epochwise median, never wrap RFI.
        moved=np.roll(held["V"],delta,axis=2)
        if delta>0:
            moved[:,:,:delta]=np.median(held["V"],axis=2)[:,:,None]
        else:
            moved[:,:,delta:]=np.median(held["V"],axis=2)[:,:,None]
        new_path=hpath+delta
        assert new_path.min()>=0 and new_path.max()<F
        transformed=physical_features(moved,new_path)
        fresh=predict(mod,transformed)
        error=abs(fresh-raw_held)
        cross=((fresh>lock)!=(raw_held>lock))
        shifts.append({"delta_bins":delta,"rows":len(fresh),"median_abs_change":float(np.median(error)),
          "p95_abs_change":float(np.quantile(error,.95)),
          "hypothetical_threshold_crossings":int(cross.sum()),
          "positive_injection_crossings":int((cross&held["positive"]).sum())})
    weak=(held["positive"]&held["weak_eval_only"])
    weak_hit=int((hit&weak).sum())
    success_auc=all(f["physical_rfi_auc"]>=.7 for f in folds)
    success_zero=all(f["hypothetical_real_simulated_positive_losses"]==0 for f in folds)
    shift_cross=sum(f["hypothetical_threshold_crossings"] for f in shifts)
    shift_rate=shift_cross/sum(f["rows"] for f in shifts)
    shift_p95=max(f["p95_abs_change"] for f in shifts)
    success_shift=shift_rate<=.01 and shift_p95<=.2
    m={"stage":"R41S1H3","status":"PASS_SHADOW_CODE_AND_NEW_SYNTHETIC_EXAM__HOLD_REAL_R42B" if (success_auc and success_zero and success_shift) else "PASS_SHADOW_CODE__FAIL_ONE_OR_MORE_PREDECLARED_GATES__HOLD_REAL_R42B",
      "simulation_source":"Synthetic clock/intermod perturbations, NOT real RFI origin evidence.",
      "development_rows":len(dev["V"]),"new_heldout_rows":len(held["V"]),
      "generated_seed_byte_reproducible":check,"seed_archive_hashes":source_digests,
      "heldout_by_seed":folds,"family_metrics":family,
      "blind_candidate_path_localized_simulated_positive":int(hit.sum()),
      "weak_simulated_count":int(weak.sum()),"weak_blind_path_localized":weak_hit,
      "frequency_shifts":shifts,"shift_threshold_flip_fraction":shift_rate,
      "max_p95_abs_score_shift":shift_p95,
      "preregistered_gates":{"auc_at_least_0.7_both":success_auc,"zero_hypothetical_positive_loss":success_zero,"shift_flip_at_most_1_percent_and_p95_below_point_2":success_shift},
      "historical_38_parent_alerts_scored":False,"historical_audit_used":False,
      "external_telescope_test":False,"tcs_custodian_R42B_status":"HOLD",
      "literal_f0_double_harmonic_available":False,
      "actual_suppressed_alerts":0,"actual_false_alarm_reduction":0,
      "previous_reference":{"false_alerts":51,"localized":318,"weak":69},
      "production_modified":False,"hawking_transient":"HOLD",
      "scientific_next":"Physical Hz normalization and known oscillator frequency controls require qualified native ON/OFF and actual telescope pointing; no hard veto and do not retune from failed heldout."}
    save("feature_contract.json",{"features":FEATS,"explicit_2f0_harmonics_unsupported":True,
      "clock_offset_bins":[8,12,16,24],"native_Hz_clock_frequency_known":False,
      "no_label_dependency":True,"no_ONOFF_emulation":True,
      "single_electronic_carrier_cannot_be_attributed_as_RFI_by_one_view":True})
    save("family_matrix.json",{"families":family,"heldout":folds})
    np.savez_compressed(OUT/"per_example_scores.npz",heldout_scores=raw_held,
        heldout_features=Xh,heldout_labels_evaluation_only=held["positive"],
        heldout_family=held["family"],heldout_seed=held["seed"],candidate_path=hpath,
        ridge_score=henergy,simulated_truth_evaluation_only=held["truth_eval_only"])
    save("result.json",m)
    print("R41S1H3_RESULT",json.dumps(m),flush=True)
if __name__=="__main__":main()
