#!/usr/bin/env python3
"""Fresh archive paired prospective comparison of retained spectral detectors."""
import hashlib,json,pathlib,time,statistics,math,sys
import numpy as np
import torch
import juxtapose_r18a as j
import juxtapose_undina_r18bc as old
import juxtapose_undina_r18b2_vector as vec
import r19a_physarum_curves as curved
import r17c_frozen_holdout as hold
import r17b_hybrid as hybrid
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r19e_integrated_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r19e");OUT.mkdir(parents=True,exist_ok=True)
hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))}
hold.OUT=OUT;hold.P=P;old.OUT=OUT;old.S={**old.S,"seed":S["seed"]};eng.OUT=OUT
R=tuple(S["rank"])
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def emit(k,d):print(k,json.dumps(d),flush=True)
def gather(q):
    data=old.gather(q);rows=[]
    rng=np.random.default_rng(int.from_bytes(hashlib.sha256((S["seed"]+"|"+q["target"]).encode()).digest()[:8],"big"))
    for k,(raw,label,hard,sha) in enumerate(v for v in data if v[1]==0):
        n,w=raw.shape;t=np.arange(n);mode=k%4;curve=bool(mode&1);inter=bool(mode&2)
        amp=S["amplitudes"][(k//4)%4];center=int(rng.integers(42,w-42))
        if curve:
            slope=int(rng.integers(-1,2));a=float(rng.choice(S["curve_accel"]))
            path=center+np.rint(slope*t+a*(t-(n-1)/2)**2).astype(int)
        else:
            shift=int(rng.integers(-32,33));path=center+np.rint(t*shift/(n-1)).astype(int)
        assert path.min()>=0 and path.max()<w and np.max(np.abs(np.diff(path)))<=3
        active=np.ones(n,bool)
        if inter:active[(np.arange(n)+int(rng.integers(0,5)))%5==0]=False
        scale=1.4826*np.median(np.abs(raw-np.median(raw,axis=1)[:,None]),axis=1)
        for y in (0,1):
            x=raw.copy()
            if y:x[t[active],path[active]]+=amp*scale[active]
            rows.append((x,y,mode,amp,hard,sha))
    assert len(rows)==256
    return rows
def dp_score(v):
    p=curved.dp(v)
    return float(np.sum(v[np.arange(len(v)),p],dtype=np.float64)/np.sqrt(len(v)))
def hybrid_score(v,model,mu,sd):
    fe=np.asarray(hybrid.track_features(v),np.float32)
    with torch.no_grad():return float(model(torch.from_numpy(np.clip((fe-mu)/sd,-20,20))[None]).softmax(1)[0,1])
def cutoff(x):
    a=sorted(map(float,x),reverse=True)
    return float(np.nextafter(a[math.floor(len(a)*S["calibration_control_fraction"])],np.inf))
def calibrate(rows):
    start=time.perf_counter()
    V=vec.normalize([r[0] for r in rows]);y=np.array([r[1] for r in rows])
    a=np.array([hold.classical(v) for v in V]);b=np.array([dp_score(v) for v in V])
    ct=cutoff(a[y==0]);dt=cutoff(b[y==0])
    # Scale each score by its control threshold; calibrate their union on the same joint control budget.
    assert ct>0 and dt>0
    joint=cutoff(np.maximum(a[y==0]/ct,b[y==0]/dt))
    lock={"straight":ct,"dp":dt,"joint_straight":ct*joint,"joint_dp":dt*joint,"hybrid":S["hybrid_threshold"],"band":S["selective_band"],"calibration_seconds":time.perf_counter()-start,"calibration_controls":int(np.sum(y==0)),"calibration_alerts":{"straight":int(np.sum(a[y==0]>=ct)),"dp":int(np.sum(b[y==0]>=dt)),"joint":int(np.sum((a[y==0]>=ct*joint)|(b[y==0]>=dt*joint)))}}
    return lock
def score_run(rows,lock,model,mu,sd):
    X=[r[0] for r in rows];y=np.array([r[1] for r in rows]);mode=np.array([r[2] for r in rows]);hard=np.array([r[4] for r in rows])
    def vector(V,t):
        # Underlying vector kernel sums float32 rows; near-threshold exact fallback.
        vec.TH=t
        return vec.vector_route(V,R,S["wave_width"])[0]
    def run(k):
        V=vec.normalize(X);calls=0;dp_calls=0
        if k=="scalar_straight":pred=np.asarray([hold.classical(v)>=lock["straight"] for v in V])
        elif k=="vector_straight":pred=vector(V,lock["straight"])
        elif k=="curved_dp":
            pred=np.asarray([dp_score(v)>=lock["dp"] for v in V]);dp_calls=len(V)
        elif k=="joint_matched_control":
            pred=vector(V,lock["joint_straight"])
            for i in np.flatnonzero(~pred):pred[i]=dp_score(V[i])>=lock["joint_dp"];dp_calls+=1
        elif k=="full_hybrid":
            pred=np.asarray([hybrid_score(v,model,mu,sd)>=lock["hybrid"] for v in V]);calls=len(V)
        elif k=="selective_hybrid":
            # Exact full score needed to identify band, no hypothetical speed credit.
            a=np.asarray([hold.classical(v) for v in V]);pred=a>=lock["straight"]
            for i in np.flatnonzero((a>=lock["straight"]-lock["band"])&~pred):
                pred[i]=hybrid_score(V[i],model,mu,sd)>=lock["hybrid"];calls+=1
        else:raise ValueError(k)
        return pred,calls,dp_calls
    modes=S["modes"];out={k:[] for k in modes};first={}
    for rep in range(S["repeats"]):
        keys=list(modes);np.random.default_rng(19050+rep).shuffle(keys)
        for k in keys:
            t=time.perf_counter();cpu=time.process_time();pred,calls,dp_calls=run(k);dt=time.perf_counter()-t;dc=time.process_time()-cpu
            if k not in first:first[k]=pred.copy()
            else:assert np.array_equal(first[k],pred),k
            out[k].append({"wall":dt,"cpu":dc,"ai_calls":calls,"dp_calls":dp_calls})
    assert np.array_equal(first["scalar_straight"],first["vector_straight"])
    assert np.all(first["selective_hybrid"]|~first["scalar_straight"])
    def subset(p,mask):
        return {"detected":int(np.sum(p&mask&(y==1))),"positives":int(np.sum(mask&(y==1))),"control_alerts":int(np.sum(p&mask&(y==0))),"controls":int(np.sum(mask&(y==0)))}
    return {k:{**subset(first[k],np.ones(len(y),bool)),"subtypes":{str(i):subset(first[k],mode==i) for i in range(4)},"hard":subset(first[k],hard),"wall":statistics.median(r["wall"] for r in out[k]),"cpu":statistics.median(r["cpu"] for r in out[k]),"ai_calls":out[k][0]["ai_calls"],"dp_calls":out[k][0]["dp_calls"],"repeats":out[k]} for k in modes}
def main():
    torch.set_num_threads(1)
    m=json.loads((OUT/"manifest.json").read_text())
    assert len(m["selected"])==6 and m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest()
    assert len({q["target"] for q in m["selected"]})==6 and all(q["target"] not in hold.S["exclude"] for q in m["selected"])
    t=time.perf_counter();model,mu,sd,_=hold.load_model();model_load=time.perf_counter()-t
    t=time.perf_counter();cal=[]
    for q in m["selected"][:2]:cal+=gather(q)
    cal_load=time.perf_counter()-t
    lock=calibrate(cal);lock["calibration_targets"]=[q["target"] for q in m["selected"][:2]]
    lock["calibration_load_seconds"]=cal_load;save("lock.json",lock);emit("R19E_LOCK",lock)
    rows={};load=0;calsha={r[5] for r in cal}
    for q in m["selected"][2:]:
        t=time.perf_counter();data=gather(q);load+=time.perf_counter()-t
        assert not calsha&{r[5] for r in data}
        rows[q["target"]]=score_run(data,lock,model,mu,sd);save("partial.json",rows)
        emit("R19E_TARGET",{"target":q["target"],"metrics":rows[q["target"]]})
    keys=("detected","positives","control_alerts","controls","wall","cpu","ai_calls","dp_calls")
    total={k:{a:sum(x[k][a] for x in rows.values()) for a in keys} for k in S["modes"]}
    for k in total:
        total[k]["subtypes"]={str(i):{a:sum(x[k]["subtypes"][str(i)][a] for x in rows.values()) for a in keys[:4]} for i in range(4)}
        total[k]["hard"]={a:sum(x[k]["hard"][a] for x in rows.values()) for a in keys[:4]}
    summary={"targets":list(rows),"samples":1024,"vector_mismatches":0,"lock":lock,"modes":total,"per_target":rows,"model_load_seconds":model_load,"evaluation_load_seconds":load,"scope":"Fresh archive synthetic injection test; balanced curve and gaps, fixed weights; joint detector calibrated to same control budget, selective hybrid no veto; times include normalization and actual scoring, no I/O, separate load times. No origin inference.","production_adoption":False}
    save("result.json",summary);emit("R19E_SUMMARY",{k:v for k,v in summary.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="scan":main()
    else:raise ValueError("phase")
