#!/usr/bin/env python3
"""R21E: small neural ranker trained only on previously inspected weak-signal spectra."""
import hashlib,json,pathlib,statistics,time
import numpy as np
import torch
from torch import nn
import r20c_verified_tile_cache as cache
import r21a_weak_template as prior
import r20d_weak_joint as reg
import r20a_vector_dp as dp
import r17c_frozen_holdout as hold
import juxtapose_undina_r18b2_vector as vec
import juxtapose_undina_r18bc as old
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r21e_neural_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r21e");OUT.mkdir(parents=True,exist_ok=True)
cache.S=S;cache.OUT=OUT;cache.CACHE=OUT/"cache";cache.CACHE.mkdir(exist_ok=True)
cache.COLD=OUT/"cold"
def emit(k,d):print(k,json.dumps(d),flush=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def features(raw):
    V=vec.normalize(raw).astype(np.float64)
    straight=np.asarray([hold.classical(v) for v in V]);ds,_=dp.sliced(V,256,True)
    rs,p=reg.reg_sliced(V,.7,True);ts,_=prior.template_sliced(V,None)
    n=V.shape[1];w=V.shape[2];r=np.arange(len(V))[:,None];t=np.arange(n)[None,:]
    v=V[r,t,p];med=np.median(v,axis=1);iqr=np.quantile(v,.75,axis=1)-np.quantile(v,.25,axis=1)
    bg=(V[r,t,np.clip(p-4,0,w-1)]+V[r,t,np.clip(p+4,0,w-1)])/2
    share=np.max(np.abs(v),axis=1)/(1e-8+np.sum(np.abs(v),axis=1))
    a=np.column_stack((straight,ds,rs,ts,v.mean(1),med,iqr,v.reshape(len(V),4,4).mean(2).min(1),(v>1).mean(1),np.abs(np.diff(p,n=2,axis=1)).sum(1),(v-bg).mean(1),share))
    assert a.shape==(len(V),12) and np.isfinite(a).all()
    return a.astype(np.float32)
def net(kind):
    if kind=="logistic":return nn.Linear(12,1)
    if kind=="mlp":return nn.Sequential(nn.Linear(12,24),nn.ReLU(),nn.Linear(24,8),nn.ReLU(),nn.Linear(8,1))
    raise ValueError(kind)
def fit(x,y,kind):
    torch.manual_seed(S["seed"]+(kind=="mlp"));m=net(kind)
    opt=torch.optim.AdamW(m.parameters(),lr=S["model"]["learning_rate"],weight_decay=S["model"]["weight_decay"])
    a=torch.from_numpy(x);b=torch.tensor(y,dtype=torch.float32)
    rng=np.random.default_rng(S["seed"]);start=time.perf_counter()
    for epoch in range(S["model"]["epochs"]):
        m.train()
        for i in np.array_split(rng.permutation(len(x)),12):
            loss=nn.functional.binary_cross_entropy_with_logits(m(a[i]).reshape(-1),b[i]);opt.zero_grad();loss.backward();opt.step()
    return m.eval(),time.perf_counter()-start
def predict(m,x):
    with torch.no_grad():return m(torch.from_numpy(x)).reshape(-1).numpy().astype(float)
def cutoff(x):
    a=sorted(map(float,x),reverse=True)
    return float(np.nextafter(a[int(len(a)*S["calibration_control_fraction"])],np.inf))
def met(p,y,m):return {"detected":int(np.sum(p&(y==1)&m)),"positives":int(np.sum((y==1)&m)),"control_alerts":int(np.sum(p&(y==0)&m)),"controls":int(np.sum((y==0)&m))}
def joint(a,b,y):
    # Convert scores to additive normalized excess; scale from controls only, supports negative logits.
    a0=cutoff(a[y==0]);b0=cutoff(b[y==0]);sa=max(.001,float(np.std(a[y==0])));sb=max(.001,float(np.std(b[y==0])))
    c=np.maximum((a-a0)/sa,(b-b0)/sb);th=cutoff(c[y==0])
    return {"a0":a0,"b0":b0,"sa":sa,"sb":sb,"threshold":th}
def apply(a,b,l):return np.maximum((a-l["a0"])/l["sa"],(b-l["b0"])/l["sb"])>=l["threshold"]
def main():
    torch.set_num_threads(1);torch.use_deterministic_algorithms(True)
    sources,by,m,restore=cache.restore();eng.get_tiles=lambda f:cache.cached_get(f,sources,by);old.OUT=cache.CACHE
    assert m["protocol_sha256"]==hashlib.sha256(prior.P.read_bytes()).hexdigest()
    rows={};groups={};seen=set();times={}
    for q in m["selected"]:
        t=time.perf_counter();r=prior.gather(q);h={a[6] for a in r};assert not seen&h;seen|=h
        rows[q["target"]]=r;f=features([a[0] for a in r]);times[q["target"]]=time.perf_counter()-t;groups[q["target"]]=f
    assert set(groups)==set(sum(S["roles"].values(),[]))
    def split(role):
        names=S["roles"][role];return np.concatenate([groups[k] for k in names]),sum([rows[k] for k in names],[])
    x,train=split("train");y=np.asarray([r[1] for r in train]);mu=x.mean(0);sd=np.maximum(x.std(0),.001)
    models={};training={}
    for k in ("logistic","mlp"):models[k],training[k]=fit(np.clip((x-mu)/sd,-20,20),y,k)
    weights={"schema":"janus.cosmos.weak_neural.v1","features":S["features"],"roles":S["roles"],"mu":mu.tolist(),"scale":sd.tolist(),"models":{k:{n:v.detach().cpu().tolist() for n,v in m.state_dict().items()} for k,m in models.items()}}
    save("weights.json",weights);digest=hashlib.sha256((OUT/"weights.json").read_bytes()).hexdigest()
    X,cal=split("calibration");Y=np.asarray([r[1] for r in cal]);z=np.clip((X-mu)/sd,-20,20)
    scores={k:predict(m,z) for k,m in models.items()}
    # Recreate frozen R21A math joint ratio calibration.
    A=X[:,0].astype(float);B=X[:,2].astype(float);a0=cutoff(A[Y==0]);b0=cutoff(B[Y==0]);math_score=np.maximum(A/a0,B/b0);math_t=cutoff(math_score[Y==0])
    lock={"math_a0":a0,"math_b0":b0,"math_t":math_t,"logistic":cutoff(scores["logistic"][Y==0]),"mlp":cutoff(scores["mlp"][Y==0]),"math_mlp_union":joint(math_score,scores["mlp"],Y),"weights_sha256":digest}
    lock["calibration"]={k:met(p,Y,np.ones(len(Y),bool)) for k,p in {"math_joint":math_score>=math_t,"logistic":scores["logistic"]>=lock["logistic"],"mlp":scores["mlp"]>=lock["mlp"],"math_mlp_union":apply(math_score,scores["mlp"],lock["math_mlp_union"])}.items()}
    save("lock.json",lock);emit("R21E_LOCK",lock)
    X,rr=split("development_check");Y=np.asarray([r[1] for r in rr]);amp=np.asarray([r[4] for r in rr]);mod=np.asarray([r[2] for r in rr])
    Z=np.clip((X-mu)/sd,-20,20);a=np.maximum(X[:,0].astype(float)/a0,X[:,2].astype(float)/b0);l=predict(models["logistic"],Z);m=predict(models["mlp"],Z)
    pred={"math_joint":a>=math_t,"logistic":l>=lock["logistic"],"mlp":m>=lock["mlp"],"math_mlp_union":apply(a,m,lock["math_mlp_union"])}
    base=pred["math_joint"];out={}
    for k,p in pred.items():
        out[k]={**met(p,Y,np.ones(len(Y),bool)),"amplitudes":{str(v):met(p,Y,amp==v) for v in prior.S["amplitudes"]},"models":{v:met(p,Y,mod==v) for v in prior.S["injection_models"]},"paired_vs_math":{"added":int(np.sum(p&~base&(Y==1))),"lost":int(np.sum(~p&base&(Y==1))),"added_controls":int(np.sum(p&~base&(Y==0))),"removed_controls":int(np.sum(~p&base&(Y==0)))}}
    # Time actual full feature extraction + neural inference versus math-only computation.
    timing={k:[] for k in S["modes"]};raw=[r[0] for r in rr]
    for rep in range(3):
        for k in S["modes"]:
            t=time.perf_counter()
            if k=="math_joint":
                v=vec.normalize(raw);s=np.asarray([hold.classical(q) for q in v]);r,_=reg.reg_sliced(v,.7,True);p=np.maximum(s/a0,r/b0)>=math_t
            else:
                f=features(raw);z=np.clip((f-mu)/sd,-20,20);a=np.maximum(f[:,0].astype(float)/a0,f[:,2].astype(float)/b0)
                q=predict(models["logistic" if k=="logistic" else "mlp"],z)
                p=apply(a,q,lock[k]) if k=="math_mlp_union" else q>=lock[k]
            timing[k].append(time.perf_counter()-t)
            assert np.array_equal(p,pred[k]),k
    for k in out:out[k]["wall_seconds"]=statistics.median(timing[k])
    d={"status":"DEVELOPMENT_ONLY","train_samples":len(train),"calibration_samples":len(cal),"development_samples":len(rr),"weights_sha256":digest,"training_seconds":training,"feature_and_patch_seconds":times,"cache_restore":restore,"lock":lock,"modes":out,"scope":"All targets previously viewed; training, calibration and development target disjoint. No fresh-test, independent validation or production claim. New archive exam must use exact frozen feature and weight hashes.","production_adoption":False}
    save("result.json",d);emit("R21E_SUMMARY",d)
if __name__=="__main__":main()
