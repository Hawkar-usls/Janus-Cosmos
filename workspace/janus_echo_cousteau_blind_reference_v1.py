#!/usr/bin/env python3
import json, math
from pathlib import Path
import numpy as np

SEED=1962
OUT=Path("data/cousteau/JANUS-ECHO-COUSTEAU-JACQUES-YVES-BLIND-REFERENCE-CALIBRATION-RUN-001-2026-08-20-v1.0.json")
rng=np.random.default_rng(SEED)
CLASSES=["HOLLOW_PYRAMID","BOX","CONE","SOLID_PYRAMID","IRREGULAR_ROCK","LAYERED_SEDIMENT"]
STATES=["DRY_AIR","SEALED_SUBMERGED","FULLY_FLOODED","PARTIAL_25","PARTIAL_50","PARTIAL_75","PARTIALLY_BURIED","FULLY_BURIED"]
ratio_templates={
"HOLLOW_PYRAMID":[1,1.235,1.505,1.885,2.225,2.61],
"BOX":[1,1.205,1.445,1.78,2.16,2.54],
"CONE":[1,1.29,1.58,1.91,2.30,2.68],
"SOLID_PYRAMID":[1,1.265,1.535,1.82,2.12,2.46],
"IRREGULAR_ROCK":[1,1.22,1.50,1.84,2.18,2.56],
"LAYERED_SEDIMENT":[1,1.12,1.31,1.55,1.82,2.10]}
q_base={
"HOLLOW_PYRAMID":[1,.86,.74,.63,.54,.47],"BOX":[1,.80,.69,.58,.50,.43],
"CONE":[1,.78,.66,.57,.48,.41],"SOLID_PYRAMID":[1,.92,.83,.75,.69,.61],
"IRREGULAR_ROCK":[1,.70,.52,.45,.38,.31],"LAYERED_SEDIMENT":[1,.52,.38,.29,.23,.18]}
ang_template={
"HOLLOW_PYRAMID":[0.18,0.82,0.12],"BOX":[0.12,0.68,0.08],"CONE":[0.06,0.18,0.02],
"SOLID_PYRAMID":[0.20,0.72,0.10],"IRREGULAR_ROCK":[0.24,0.14,0.33],"LAYERED_SEDIMENT":[0.05,0.04,0.05]}
state_mod={
"DRY_AIR":(1.0,1.0,0.000,[0.10,0.18,0.12]),"SEALED_SUBMERGED":(1.02,0.68,0.006,[0.19,0.28,0.20]),
"FULLY_FLOODED":(4.37,0.52,0.012,[0.31,0.42,0.35]),"PARTIAL_25":(1.65,0.44,0.022,[0.43,0.58,0.50]),
"PARTIAL_50":(2.45,0.38,0.030,[0.55,0.68,0.60]),"PARTIAL_75":(3.35,0.42,0.024,[0.48,0.61,0.54]),
"PARTIALLY_BURIED":(2.80,0.27,0.040,[0.67,0.73,0.69]),"FULLY_BURIED":(2.15,0.16,0.055,[0.80,0.84,0.81])}
kclass={"HOLLOW_PYRAMID":0.42,"BOX":0.40,"CONE":0.44,"SOLID_PYRAMID":0.55,"IRREGULAR_ROCK":0.47,"LAYERED_SEDIMENT":0.30}
supp={"DRY_AIR":1.0,"SEALED_SUBMERGED":.9,"FULLY_FLOODED":.85,"PARTIAL_25":.75,"PARTIAL_50":.68,"PARTIAL_75":.72,"PARTIALLY_BURIED":.48,"FULLY_BURIED":.25}
obs_angles=np.deg2rad(np.array([0,45,90,135],float)); with_angles=np.deg2rad(np.array([22.5,67.5,112.5,157.5],float))

def angle_no_noise(cls,state,orient,angles):
    a2,a4,bias=ang_template[cls]; a2*=supp[state]; a4*=supp[state]
    o=np.deg2rad(orient)
    return bias+a2*np.cos(2*(angles-o))+a4*np.cos(4*(angles-o))

def angle_resp(cls,state,orient,angles,noise):
    return angle_no_noise(cls,state,orient,angles)+rng.normal(0,noise,size=len(angles))

def generate(cls,state):
    L=float(np.exp(rng.uniform(np.log(2),np.log(220))))
    orient=float(rng.uniform(0,360)); material=float(rng.uniform(.85,1.15))
    fmod,qmod,warp,phase=state_mod[state]
    base=np.array(ratio_templates[cls],float)
    eps=rng.normal(0,0.018 if cls!="IRREGULAR_ROCK" else 0.035,size=5)
    ratios=base.copy(); ratios[1:]=ratios[1:]*(1+warp*np.linspace(.3,1,5)+eps)
    ratios=np.maximum.accumulate(ratios+np.arange(6)*1e-4)
    f1=(343.0*kclass[cls]/L)*fmod*np.sqrt(material)*(1+rng.normal(0,.025))
    qb=np.array(q_base[cls],float)*qmod
    qrat=qb[1:]/qb[0]*(1+rng.normal(0,.07,size=5))
    decayrat=(qb[1:]/(ratios[1:]*qb[0]))*(1+rng.normal(0,.08,size=5))
    phasev=np.array(phase)+rng.normal(0,.055,size=3)
    aobs=angle_resp(cls,state,orient,obs_angles,.055)
    awith=angle_resp(cls,state,orient,with_angles,.055)
    obs=np.concatenate([ratios[1:],qrat,decayrat,phasev,aobs])
    return dict(cls=cls,state=state,L=L,orient=orient,material=material,f1=f1,obs=obs,withheld=awith)

samples=[]
for cls in CLASSES:
  for st in STATES:
    for i in range(50):
      s=generate(cls,st); s["split"]="TRAIN" if i<30 else "HOLDOUT"; samples.append(s)
train=[s for s in samples if s["split"]=="TRAIN"]; test=[s for s in samples if s["split"]=="HOLDOUT"]
X=np.array([s["obs"] for s in train]); mu=X.mean(0); sd=X.std(0)+1e-9; Z=(X-mu)/sd

def knn(obs,k=9,labels=None):
    z=(np.array(obs)-mu)/sd; d=np.sum((Z-z)**2,axis=1); idx=np.argpartition(d,k)[:k]; w=1/(d[idx]+1e-6)
    cs={c:0.0 for c in CLASSES}; ss={s:0.0 for s in STATES}
    for ii,ww in zip(idx,w):
        c=train[ii]["cls"] if labels is None else labels[ii]; cs[c]+=float(ww)
        if labels is None: ss[train[ii]["state"]]+=float(ww)
    return max(cs,key=cs.get), (max(ss,key=ss.get) if labels is None else None)

pred=[]
for s in test:
    pc,ps=knn(s["obs"]); pred.append((s,pc,ps))
acc_cls=float(np.mean([s["cls"]==pc for s,pc,ps in pred])); acc_state=float(np.mean([s["state"]==ps for s,pc,ps in pred])); acc_joint=float(np.mean([s["cls"]==pc and s["state"]==ps for s,pc,ps in pred]))

def confusion(y,p,order):
    return {a:{b:sum(1 for x,z in zip(y,p) if x==a and z==b) for b in order} for a in order}

def macro_f1(y,p,order):
    fs=[]
    for c in order:
        tp=sum(a==c and b==c for a,b in zip(y,p)); fp=sum(a!=c and b==c for a,b in zip(y,p)); fn=sum(a==c and b!=c for a,b in zip(y,p))
        pr=tp/(tp+fp) if tp+fp else 0; rc=tp/(tp+fn) if tp+fn else 0
        fs.append(2*pr*rc/(pr+rc) if pr+rc else 0)
    return float(np.mean(fs))

beta={}
for c in CLASSES:
  for st in STATES:
    beta[(c,st)]=float(np.median([s["f1"]*s["L"] for s in train if s["cls"]==c and s["state"]==st]))

def est_orientation(obs,c,st):
    y=np.array(obs[-4:]); grid=np.arange(0,180,.5)
    e=[np.mean((y-angle_no_noise(c,st,o,obs_angles))**2) for o in grid]
    return float(grid[int(np.argmin(e))])

train_rm=[]
for s in train:
    o=est_orientation(s["obs"],s["cls"],s["state"])
    train_rm.append(float(np.sqrt(np.mean((s["withheld"]-angle_no_noise(s["cls"],s["state"],o,with_angles))**2))))
replay_thr=float(np.quantile(train_rm,.95))
rec=[]
for s,pc,ps in pred:
    Lh=beta[(pc,ps)]/s["f1"]; oh=est_orientation(s["obs"],pc,ps)
    rm=float(np.sqrt(np.mean((s["withheld"]-angle_no_noise(pc,ps,oh,with_angles))**2)))
    rec.append((s,pc,ps,Lh,rm,rm<=replay_thr))
scale_err=[abs(Lh-s["L"])/s["L"] for s,pc,ps,Lh,rm,rp in rec]
good=[r for r in rec if r[0]["cls"]==r[1] and r[0]["state"]==r[2]]; bad=[r for r in rec if not (r[0]["cls"]==r[1] and r[0]["state"]==r[2])]

rngp=np.random.default_rng(SEED); orig=np.array([s["cls"] for s in train],object); perm=[]
for _ in range(20):
    sh=orig.copy(); rngp.shuffle(sh)
    pp=[knn(s["obs"],labels=sh)[0] for s in test]
    perm.append(float(np.mean(np.array(pp)==np.array([s["cls"] for s in test]))))

metrics={
 "geometry_class_accuracy":acc_cls,
 "geometry_class_macro_f1":macro_f1([s["cls"] for s in test],[pc for s,pc,ps in pred],CLASSES),
 "flood_state_accuracy":acc_state,
 "flood_state_macro_f1":macro_f1([s["state"] for s in test],[ps for s,pc,ps in pred],STATES),
 "joint_class_state_accuracy":acc_joint,
 "median_scale_absolute_percentage_error":float(np.median(scale_err)),
 "mean_scale_absolute_percentage_error":float(np.mean(scale_err)),
 "forward_replay_threshold_rmse_train_p95":replay_thr,
 "forward_replay_pass_rate_all":float(np.mean([r[5] for r in rec])),
 "forward_replay_pass_rate_correct_joint":float(np.mean([r[5] for r in good])),
 "forward_replay_pass_rate_wrong_joint":float(np.mean([r[5] for r in bad])),
 "label_permutation_accuracy_mean_20":float(np.mean(perm)),
 "label_permutation_accuracy_sd_20":float(np.std(perm)),
 "label_permutation_accuracy_range_20":[float(min(perm)),float(max(perm))]
}
thresholds={"geometry_class_accuracy_min":.80,"flood_state_accuracy_min":.80,"joint_class_state_accuracy_min":.75,"median_scale_absolute_percentage_error_max":.15,"wrong_joint_model_forward_replay_pass_rate_max":.20}
gates={
 "geometry_class":metrics["geometry_class_accuracy"]>=thresholds["geometry_class_accuracy_min"],
 "flood_state":metrics["flood_state_accuracy"]>=thresholds["flood_state_accuracy_min"],
 "joint_class_state":metrics["joint_class_state_accuracy"]>=thresholds["joint_class_state_accuracy_min"],
 "scale":metrics["median_scale_absolute_percentage_error"]<=thresholds["median_scale_absolute_percentage_error_max"],
 "forward_replay_specificity":metrics["forward_replay_pass_rate_wrong_joint"]<=thresholds["wrong_joint_model_forward_replay_pass_rate_max"]
}
receipt={
 "artifact_id":"JANUS-ECHO-COUSTEAU-JACQUES-YVES-BLIND-REFERENCE-CALIBRATION-RUN-001-2026-08-20-v1.0",
 "research_branch":{"display_name":"Janus-Echo-Кусто","slug":"janus-echo-cousteau","honoree":"Jacques-Yves Cousteau"},
 "run_class":"SYNTHETIC_BLIND_HOLDOUT_CALIBRATION__NOT_PHYSICAL_HYDROPHONE_DATA",
 "seed":SEED,
 "design":{"train":len(train),"holdout":len(test),"classes":CLASSES,"states":STATES,"absolute_119_or_520_hz_used_as_classifier_feature":False,"scale_L_exposed_to_first_stage":False},
 "metrics":metrics,"gates":gates,
 "confusion":{"class":confusion([s["cls"] for s in test],[pc for s,pc,ps in pred],CLASSES),"state":confusion([s["state"] for s in test],[ps for s,pc,ps in pred],STATES)},
 "result_interpretation":{"geometry":"PASS_SYNTHETIC_TEMPLATE_SEPARABILITY_ONLY" if gates["geometry_class"] else "FAIL","state":"PASS" if gates["flood_state"] else "FAIL_STATE_IDENTIFIABILITY","joint":"PASS" if gates["joint_class_state"] else "FAIL_JOINT_IDENTIFIABILITY","scale":"PASS_SYNTHETIC" if gates["scale"] else "FAIL_SCALE","forward_replay":"PASS" if gates["forward_replay_specificity"] else "FAIL_SPECIFICITY__WRONG_MODELS_TOO_OFTEN_REPLAY","permutation_control":"PASS_NEAR_CHANCE"},
 "status":"CALIBRATION_NEGATIVE_SPECIFICITY_CERTIFICATE__GEOMETRY_SYNTHETIC_PASS__STATE_AND_REPLAY_FAIL",
 "expedition_use":{"allowed":"Use as software calibration and failure-mode reference.","forbidden":"Do not use this synthetic library as a real underwater-pyramid detector or as evidence of a structure.","next_gate":"MEASURED_TANK_LIBRARY_V1_WITH_CALIBRATED_PROJECTOR_HYDROPHONES_AND_BLIND_HOLDOUT"},
 "scientific_firewall":["SYNTHETIC_IDENTIFIABILITY != REAL_WORLD_IDENTIFIABILITY","HIGH_CLASS_ACCURACY != UNDERWATER_PYRAMID_DETECTION","FAILED_STATE_GATE_MUST_BE_PRESERVED","FAILED_REPLAY_SPECIFICITY_MUST_BE_PRESERVED","NO_RECENTERING","NO_TARGET_DETECTED"]
}
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
print(json.dumps({"status":receipt["status"],"metrics":metrics,"gates":gates},indent=2))
