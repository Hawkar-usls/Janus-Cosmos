#!/usr/bin/env python3
"""Curved-track pilot with original Physarum swarm; development archive replay."""
import hashlib,json,pathlib,random,time,math
import numpy as np
import juxtapose_r18a as j
from donor_physarum_engine import PhysarumSwarm
P=pathlib.Path("research/seti/r19a_physarum_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r19a");OUT.mkdir(parents=True,exist_ok=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def dp(y):
    n,w=y.shape;radius=S["curve_max_step"];cur=y[0].copy();trace=[]
    for r in range(1,n):
        pad=np.pad(cur,(radius,radius),constant_values=-np.inf)
        cand=np.stack([pad[k:k+w] for k in range(2*radius+1)])
        k=np.argmax(cand,axis=0);trace.append(k-radius)
        cur=y[r]+np.max(cand,axis=0)
    path=np.zeros(n,dtype=int);path[-1]=int(np.argmax(cur))
    for r in range(n-1,0,-1):path[r-1]=path[r]+trace[r-1][path[r]]
    assert np.max(np.abs(np.diff(path)))<=radius and np.all((path>=0)&(path<w))
    return path
def straight(y):
    a,b,s=max((*j.score(y,s),s) for s in j.ASC)
    return b+np.rint(np.arange(len(y))*s/(len(y)-1)).astype(int)
def smooth(y):
    p=np.pad(y,((0,0),(1,1)),mode="edge")
    return (p[:,:-2]+2*p[:,1:-1]+p[:,2:])/4
def swarm(y,index):
    random.seed(S["seed"]+index);n,w=y.shape
    z=PhysarumSwarm(w,n,n_particles=S["particles"],sensor_distance=2.0)
    env=np.clip(y,0,10)
    for _ in range(S["steps"]):z.step(env)
    trail=z.get_trail_map()
    return dp(y+trail/max(float(trail.max()),1)*2)
def dataset(a,seed):
    X,Y,meta,targets,_=j.build(a);rng=np.random.default_rng(seed)
    ix=rng.choice(np.flatnonzero(Y==0),size=S["samples_per_split"]//2,replace=False);rows=[]
    for k,i in enumerate(ix):
        raw=X[i].copy();n,w=raw.shape
        # Generate bounded curved path, with an intermittent active mask; no truth supplied to proposer.
        t=np.arange(n);c=int(rng.integers(42,w-42));d=int(rng.integers(-1,2));a2=float(rng.choice([-.09,.09,.14,-.14]))
        path=c+np.rint(d*t+a2*(t-(n-1)/2)**2).astype(int)
        assert np.max(np.abs(np.diff(path)))<=3 and path.min()>=0 and path.max()<w
        active=(t%5!=0) if k%2 else np.ones(n,bool);amp=S["amplitudes"][k%4]
        for label in (0,1):
            x=raw.copy()
            if label:
                sc=1.4826*np.median(np.abs(raw-np.median(raw,axis=1)[:,None]),axis=1)
                x[t[active],path[active]]+=amp*sc[active]
            rows.append((x,label,path,active,amp))
    return rows,sorted(targets)
def evaluate(rows):
    m={k:[] for k in S["modes"]};timing={k:0.0 for k in m}
    for i,(x,label,truth,active,amp) in enumerate(rows):
        y=j.normalized(x);n,w=y.shape
        methods={"straight":lambda:straight(y),"exact_dp":lambda:dp(y),"smooth_dp":lambda:dp(smooth(y)),"physarum_trail_dp":lambda:swarm(y,i)}
        scores={}
        for k,fn in methods.items():
            start=time.perf_counter();path=fn();score=float(np.sum(y[np.arange(n),path],dtype=np.float64)/np.sqrt(n));timing[k]+=time.perf_counter()-start
            scores[k]=score
            m[k].append({"label":label,"score":score,"amp":amp,"intermittent":not bool(np.all(active)),"path_accuracy":float(np.mean(np.abs(path[active]-truth[active])<=2)) if label else None})
        assert scores["exact_dp"]+1e-5>=max(scores["smooth_dp"],scores["physarum_trail_dp"])
    return m,timing
def threshold(rows):
    a=sorted((r["score"] for r in rows if not r["label"]),reverse=True)
    return float(np.nextafter(a[math.floor(len(a)*S["threshold_target_control_max"])],np.inf))
def metrics(rows,t,seconds):
    pos=[r for r in rows if r["label"]];neg=[r for r in rows if not r["label"]]
    return {"detected":sum(r["score"]>=t for r in pos),"positives":len(pos),"control_alerts":sum(r["score"]>=t for r in neg),"controls":len(neg),"mean_path_accuracy":float(np.mean([r["path_accuracy"] for r in pos])),"intermittent_detected":sum(r["score"]>=t for r in pos if r["intermittent"]),"intermittent_total":sum(r["intermittent"] for r in pos),"seconds":seconds}
def main():
    raw=pathlib.Path("research/seti/donor_physarum_engine.py").read_bytes()
    assert hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()==S["source"]["blob"]
    t=time.perf_counter();cal,ct=dataset(j.S["source_artifacts"][0],S["seed"]);m,sec=evaluate(cal)
    th={k:threshold(v) for k,v in m.items()}
    lock={"thresholds":th,"targets":ct,"metrics":{k:metrics(v,th[k],sec[k]) for k,v in m.items()},"calibration_seconds":time.perf_counter()-t}
    save("lock.json",lock);print("R19A_LOCK",json.dumps(lock),flush=True)
    t=time.perf_counter();test,tt=dataset(j.S["source_artifacts"][1],S["seed"]+1);assert not set(ct)&set(tt)
    m,sec=evaluate(test)
    summary={"targets":tt,"samples":len(test),"metrics":{k:metrics(v,th[k],sec[k]) for k,v in m.items()},"lock":lock,"evaluation_total_seconds":time.perf_counter()-t,"scope":"Development only; known archive, synthetic bounded curves, swarm proposal scored on unchanged raw normalized power, DP is reference within specified path domain; no origin evidence.","production_adoption":False}
    save("result.json",summary);save("scores.json",m);print("R19A_SUMMARY",json.dumps(summary),flush=True)
if __name__=="__main__":main()
