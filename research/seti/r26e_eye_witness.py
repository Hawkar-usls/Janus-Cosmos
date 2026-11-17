#!/usr/bin/env python3
import hashlib,json,pathlib,itertools
import numpy as np
import r26d_family_discriminator as fam
import donor_eye_r4_4 as eye
S=json.loads(pathlib.Path("research/seti/r26e_eye_protocol.json").read_text())
OUT=pathlib.Path("runtime/r26e");OUT.mkdir(parents=True,exist_ok=True)
fam.OUT=OUT;fam.P["seed"]=S["seed"]
CL=fam.P["classes"];MOD=fam.P["modes"]
OUTCOMES=["".join(x) for x in itertools.product("01",repeat=3)]
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def load_sample(data,target):
    X,R,l,t=fam.sample(data[target],target);q,p,_=fam.scores(X,R)
    return l,q,np.mean(abs(p-t)<=2,axis=1)>=.75
def outcome(q,thresholds):
    return np.array(["".join(str(int(q[k][i]>=thresholds[k])) for k in MOD) for i in range(len(q[MOD[0]]))])
def likelihood(labels,outs):
    res={}
    for c in CL:
        z=outs[labels==c];a={o:float(np.sum(z==o))+S["smoothing"] for o in OUTCOMES};n=sum(a.values());res[c]={o:v/n for o,v in a.items()}
    return res
def main():
    assert fam.test()=="PASS"
    data=fam.restore();cal={t:load_sample(data,t) for t in fam.P["calibration"]}
    th={k:fam.cut(np.concatenate([q[k][np.isin(l,fam.P["calibration_negative"])] for l,q,_ in cal.values()])) for k in MOD}
    likes={t:likelihood(l,outcome(q,th)) for t,(l,q,_) in cal.items()}
    causes=[{"id":c,"target_class":S["target_classes"][c],"prior":1/6} for c in CL]
    model={"id":"JANUS-COSMOS-R26E-SYNTHETIC-PATTERN_ONLY","cause_classes":causes,"requirements":{"confidence_threshold":S["confidence_threshold"],"require_failure_domain_disjoint_path":True,"robust_objective":"minimax_expected_cost"},"tests":[{"id":"joint_family_panel","cost":1,"calibrated":True,"failure_domains":["shared_receiver_and_spectrum"],"likelihood_by_cause":next(iter(likes.values()))},{"id":"second_receiver_repeat","available":False,"decision_usable":False,"calibrated":False,"cost":10},{"id":"fullband_clock_relation","available":False,"decision_usable":False,"calibrated":False,"cost":5}],"model_set":[]}
    for t,L in likes.items():
        for k,p in S["prior_models"].items():
            model["model_set"].append({"id":t+"_"+k,"provenance":"synthetic calibration at "+t+"; "+k+" declared prior sensitivity","prior_by_cause":dict(zip(CL,p)),"likelihood_overrides":{"joint_family_panel":L}})
    save("lock.json",{"thresholds":th,"model":model})
    receipt=eye.solve(model);save("eye_policy.json",receipt)
    causes,tests,models,acq,r=eye.parse(model);ids=tuple(sorted(models));prior=tuple(eye.canon_belief(models[m]["priors"],causes,r["belief_round_digits"]) for m in ids)
    decisions={}
    for o in OUTCOMES:
        post=tuple(eye.model_posterior(prior[i],"joint_family_panel",o,m,causes,models,r) for i,m in enumerate(ids))
        decisions[o]=eye.robust_terminal(post,ids,causes,r,1)
    save("outcome_decisions.json",decisions)
    counts={c:{"examples":0,"anchor_localized":0,"abstain":0,"correct":0,"wrong":0,"outcome_counts":{o:0 for o in OUTCOMES}} for c in CL}
    for t in fam.P["evaluation"]:
        labels,q,h=load_sample(data,t);outs=outcome(q,th)
        for c,o,ok in zip(labels,outs,h):
            z=counts[c];z["examples"]+=1;z["anchor_localized"]+=int(ok);z["outcome_counts"][o]+=1
            decision=decisions[o]["identified_target_class"]
            if not ok or decision is None:z["abstain"]+=1
            elif decision==S["target_classes"][c]:z["correct"]+=1
            else:z["wrong"]+=1
    # Real origin model remains uncalibrated: shared receiver descriptors are advisory only.
    real={"id":"JANUS-COSMOS-R26E-REAL-ORIGIN-HOLD","cause_classes":[{"id":"terrestrial_rfi"},{"id":"human_spacecraft"},{"id":"other_sky_or_instrument"}],"tests":[{"id":"joint_family_panel","available":True,"calibrated":False,"decision_usable":False,"cost":1},{"id":"independent_receiver","available":False,"calibrated":False,"decision_usable":False,"cost":10}],"model_set":[{"id":"balanced_sensitivity","provenance":"declared prior only; no measured origin probabilities"},{"id":"rfi_heavy_sensitivity","provenance":"declared prior sensitivity only","prior_by_cause":{"terrestrial_rfi":.8,"human_spacecraft":.1,"other_sky_or_instrument":.1}}],"requirements":{"confidence_threshold":.95}}
    real_result=eye.solve(real);assert real_result["status"]=="ROBUST_NON_IDENTIFIABLE_UNDER_DECLARED_MODEL_SET"
    save("real_origin_policy.json",real_result)
    n=sum(z["examples"] for z in counts.values());correct=sum(z["correct"] for z in counts.values());wrong=sum(z["wrong"] for z in counts.values());abstain=sum(z["abstain"] for z in counts.values());assert n==768 and correct+wrong+abstain==n
    d={"status":"DEVELOPMENT_ONLY","calibration_examples":sum(len(l) for l,q,h in cal.values()),"evaluation_examples":n,"model_count":len(models),"joint_outcome_count":len(OUTCOMES),"policy_status":receipt["status"],"real_origin_status":real_result["status"],"correct":correct,"wrong":wrong,"abstain":abstain,"counts":counts,"outcomes":{o:{"target":z["identified_target_class"],"minimum_model_confidence":z["minimum_model_confidence"]} for o,z in decisions.items()},"scope":S["scope"],"production_adoption":False}
    save("result.json",d);print("R26E_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
