#!/usr/bin/env python3
"""R28B: select a bounded linear-first policy on calibration data only."""
import pathlib,json,hashlib,sys,time,itertools
import numpy as np
import r28_hybrid as h
import r27a_benchmark as a
import r26c_reference_alignment as r
import r24_odontocrystal as od
import r17c_frozen_holdout as hold
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r28b_router_protocol.json");S=json.loads(P.read_text());OUT=pathlib.Path("runtime/r28b");OUT.mkdir(parents=True,exist_ok=True)
h.S={**S,"modes":S["modes"][:7]};r.S=S;r.OUT=OUT;od.S=S
hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))};hold.OUT=OUT;hold.P=P;eng.OUT=OUT
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def emit(k,d):print(k,json.dumps(d),flush=True)
def inputs(rows,R):
    X=r.norm(np.stack([x[0] for x in rows]));Y,g,_,_=r.register(X,r.norm(R))
    return X.astype(np.float32).astype(np.float64),Y.astype(np.float32).astype(np.float64),g
def panel(rows,R):
    t=time.perf_counter();X=r.norm(np.stack([x[0] for x in rows]));raw=time.perf_counter()-t
    t=time.perf_counter();Y,g,_,_=r.register(X,r.norm(R));X=X.astype(np.float32).astype(np.float64);Y=Y.astype(np.float32).astype(np.float64);reg=time.perf_counter()-t
    q={}
    for k in h.BASE:
        V=X if k.startswith("raw_") else Y
        fn=a.straight if k.endswith("linear") else a.tree if k.endswith("taylor") else lambda v:od.sliced(v,0,None)[:2]
        t=time.perf_counter();s,p=fn(V);q[k]=(s,p,time.perf_counter()-t)
    return (q,raw,reg,g),Y
def marg(s,L,k):return (s-L["thresholds"][k])/L["scales"][k]
def cut(v,rate):
    z=np.sort(np.asarray(v))[::-1];assert len(z)>=64
    return float(np.nextafter(z[int(len(z)*rate)],np.inf))
def route(q,L,policy):
    sl,pl,_=q["registered_linear"];sd,pd,_=q["registered_odonto"]
    zl=marg(sl,L,"registered_linear");zd=marg(sd,L,"registered_odonto")+policy["bias"]
    run=(zl>=policy["low"])&(zl<=policy["high"])
    use=run&(zd>zl)
    return np.where(use,zd,zl),np.where(use[:,None],pd,pl),run,use
def loc(p,rows):
    t=np.stack([x[7] for x in rows]);m=np.stack([x[8] for x in rows])
    return ((abs(p-t)<=2)&m).sum(1)>=.75*m.sum(1)
def fit_router(cal,L):
    y=np.concatenate([h.labels(rows) for rows,pack in cal]);grid=S["routing"]["candidate_grid"];trials=[]
    for low,high,bias in itertools.product(grid["low"],grid["high"],grid["odonto_bias"]):
        p={"low":low,"high":high,"bias":bias}
        zz=[route(pack[0],L,p) for rows,pack in cal]
        threshold=cut(np.concatenate([v[0] for v in zz])[y==0],.04)
        counts=[];allloc=0;calls=0;alerts=0
        for (rows,pack),(z,paths,run,use) in zip(cal,zz):
            yy=h.labels(rows);weak=np.array([x[4]<=1.25 for x in rows]);hit=(z>=threshold)&loc(paths,rows)&(yy==1)
            counts.append(int(np.sum(hit&weak)));allloc+=int(hit.sum());calls+=int(run.sum());alerts+=int(np.sum((z>=threshold)&(yy==0)))
        trials.append({**p,"threshold":threshold,"weak_by_target":counts,"weak":sum(counts),"localized":allloc,"dp_calls":calls,"controls":alerts})
    # Predeclared deterministic order for ties.
    pick=max(trials,key=lambda t:(min(t["weak_by_target"]),t["weak"],t["localized"],-t["dp_calls"],-t["controls"]))
    assert pick["controls"]<=int(np.sum(y==0)*.04)
    cons=cut(np.concatenate([h.fuse(pack[0],L,"registered_hybrid")[0] for rows,pack in cal])[y==0],.04)
    return {"policy":pick,"trials":trials,"conservative_threshold":cons}
def actual(rows,R,L,p):
    t=time.perf_counter();X,Y,g=inputs(rows,R)
    s,paths=a.straight(Y);z=marg(s,L,"registered_linear")
    run=(z>=p["low"])&(z<=p["high"]);use=np.zeros(len(z),bool)
    if run.any():
        ds,dp,_=od.sliced(Y[run],0,None);v=marg(ds,L,"registered_odonto")+p["bias"];idx=np.flatnonzero(run);take=v>z[run];use[idx[take]]=True
        z[idx[take]]=v[take];paths[idx[take]]=dp[take]
    return z,paths,run,use,time.perf_counter()-t
def measure(rows,pack,R,L,RL):
    q,raw,reg,g=pack;out,record=h.met(rows,pack,L)
    y=h.labels(rows)
    z,p,run,use,dt=actual(rows,R,L,RL["policy"]);z0,p0,r0,u0=route(q,L,RL["policy"])
    assert np.allclose(z,z0,atol=1e-10) and np.array_equal(p,p0) and np.array_equal(run,r0)
    v,w,t,choice=h.fuse(q,L,"registered_hybrid")
    for k,s,paths,seconds,threshold in (("adaptive_router",z,p,dt,RL["policy"]["threshold"]),("conservative_hybrid",v,w,raw+reg+t,RL["conservative_threshold"])):
        pred=s>=threshold;local=loc(paths,rows)
        out[k]={**od.measure(pred,paths,np.zeros_like(paths,bool),rows),"seconds":seconds}
        record["arms"][k]={"scores":s.tolist(),"localized":local.tolist(),"alert":pred.tolist(),"paths":paths.tolist()}
        out[k]["paired"]={}
        for base in ("registered_hybrid","registered_linear","registered_odonto","registered_taylor"):
            b=record["arms"][base];bp=np.asarray(b["alert"]);bl=np.asarray(b["localized"]);bb=bp&bl&(y==1);aa=pred&local&(y==1)
            out[k]["paired"][base]={"added_localized":int(np.sum(aa&~bb)),"lost_localized":int(np.sum(bb&~aa)),"added_controls":int(np.sum(pred&~bp&(y==0))),"removed_controls":int(np.sum(bp&~pred&(y==0)))}
    out["adaptive_router"]["dp_calls"]=int(run.sum());out["adaptive_router"]["odonto_selected"]=int(use.sum())
    return out,record
def export(rows,Y,group,target,record=None):
    d={"V":Y.astype(np.float32),"y":h.labels(rows),"truth":np.stack([x[7] for x in rows]),"active":np.stack([x[8] for x in rows]),"amp":np.array([x[4] for x in rows]),"model":np.array([x[2] for x in rows])}
    if record:
        for k,v in record["arms"].items():d[k+"_alert"]=np.asarray(v["alert"]);d[k+"_localized"]=np.asarray(v["localized"])
    p=OUT/"bliss_inputs";p.mkdir(exist_ok=True);np.savez_compressed(p/(group+"_"+target+".npz"),**d)
def main():
    jit=a.init_tree();checks=h.selftest();emit("R28B_SELFTEST",checks)
    M=json.loads((OUT/"manifest.json").read_text());assert M["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest() and len(M["selected"])==6
    names=[q["target"] for q in M["selected"]];assert len(set(names))==6 and not set(names)&set(hold.S["exclude"])
    seen=set();cal=[];load=0
    for q in M["selected"][:2]:
        t=time.perf_counter();rows,R,_,hs=r.gather(q);load+=time.perf_counter()-t;assert not seen&hs;seen|=hs
        pack,Y=panel(rows,R);cal.append((rows,pack));export(rows,Y,"cal",q["target"])
    L=h.fit(cal);RL=fit_router(cal,L);lock={"baseline":L,"router":RL,"targets":names[:2]};save("lock.json",lock);locksha=hashlib.sha256((OUT/"lock.json").read_bytes()).hexdigest();emit("R28B_LOCK",{"policy":RL["policy"],"conservative_threshold":RL["conservative_threshold"],"targets":names[:2],"sha256":locksha})
    out={};rec={};leak={}
    for q in M["selected"][2:]:
        t=time.perf_counter();rows,R,ch,hs=r.gather(q);load+=time.perf_counter()-t;assert not seen&hs;seen|=hs
        pack,Y=panel(rows,R);target=q["target"];out[target],rec[target]=measure(rows,pack,R,L,RL);export(rows,Y,"eval",target,rec[target])
        lr=[x[0] for x in ch];rr=np.stack([x[1] for x in ch]);lp,_=panel(lr,rr);leak[target]=measure(lr,lp,rr,L,RL)[0]
        save("partial.json",out);emit("R28B_TARGET",{"target":target,"metrics":out[target]})
    assert hashlib.sha256((OUT/"lock.json").read_bytes()).hexdigest()==locksha
    keys=("detected","localized","positives","control_alerts","controls","seconds")
    total={k:{a:sum(v[k][a] for v in out.values()) for a in keys} for k in S["modes"]}
    for k in total:
        for f,vs in (("amplitudes",S["amplitudes"]),("gaps",S["gap_levels"]),("models",S["injection_models"])):
            total[k][f]={str(z):{a:sum(v[k][f][str(z)][a] for v in out.values()) for a in keys[:5]} for z in vs}
        total[k]["weak_localized"]=sum(total[k]["amplitudes"][str(z)]["localized"] for z in (.75,1.25))
        if k in ("adaptive_router","conservative_hybrid"):
            total[k]["paired"]={b:{a:sum(x[k]["paired"][b][a] for x in out.values()) for a in ("added_localized","lost_localized","added_controls","removed_controls")} for b in ("registered_hybrid","registered_linear","registered_odonto","registered_taylor")}
    total["adaptive_router"]["dp_calls"]=sum(v["adaptive_router"]["dp_calls"] for v in out.values())
    x=total["adaptive_router"];b=total["registered_hybrid"];gate=x["localized"]>=b["localized"] and x["control_alerts"]<=min(25,b["control_alerts"]) and x["weak_localized"]>b["weak_localized"]
    d={"status":"FRESH_ARCHIVE_SYNTHETIC_EXAM","gate":"PASS_DEVELOPMENT" if gate else "HOLD","samples":1024,"targets":names[2:],"selftests":checks,"lock_sha256":locksha,"policy":RL["policy"],"modes":total,"per_target":out,"leak_challenge":{k:{a:sum(v[k][a] for v in leak.values()) for a in ("localized","detected","positives")} for k in total},"archive_load_seconds":load,"external_fetch_and_jit_seconds":jit,"scope":S["scope"],"production_adoption":False}
    save("scores.json",rec);save("result.json",d);emit("R28B_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    else:main()
