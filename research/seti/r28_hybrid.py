#!/usr/bin/env python3
"""R28: frozen common control budget, linear + nonlinear, preserve raw evidence."""
import hashlib,json,pathlib,sys,time
import numpy as np
import r27a_benchmark as a
import r26c_reference_alignment as r
import r24_odontocrystal as od
import r17c_frozen_holdout as hold
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r28_hybrid_protocol.json");S=json.loads(P.read_text());OUT=pathlib.Path("runtime/r28");OUT.mkdir(parents=True,exist_ok=True)
r.S=S;r.OUT=OUT;od.S=S;hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))};hold.OUT=OUT;hold.P=P;eng.OUT=OUT
BASE=S["modes"][:5]
GROUPS={"registered_hybrid":["registered_linear","registered_odonto"],"protected_hybrid":["raw_linear","raw_odonto","registered_linear","registered_odonto"]}
def emit(k,d):print(k,json.dumps(d),flush=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def prep(rows,R):
    t=time.perf_counter();X=r.norm(np.stack([x[0] for x in rows]));raw=time.perf_counter()-t
    t=time.perf_counter();Y,g,sh,conf=r.register(X,r.norm(R));reg=time.perf_counter()-t
    res={}
    for k in BASE:
        V=X if k.startswith("raw_") else Y
        fn=a.straight if k.endswith("linear") else a.tree if k.endswith("taylor") else lambda v:od.sliced(v,0,None)[:2]
        t=time.perf_counter();s,p=fn(V);elapsed=time.perf_counter()-t
        assert np.isfinite(s).all() and p.shape==(len(rows),16) and np.all((p>=0)&(p<V.shape[2]))
        res[k]=(s,p,elapsed)
    return res,raw,reg,g
def labels(rows):return np.asarray([x[1] for x in rows])
def cut(x):
    z=np.sort(np.asarray(x))[::-1];assert len(z)>=64
    return float(np.nextafter(z[int(len(z)*S["calibration_control_fraction"])],np.inf))
def fuse(q,L,k):
    members=GROUPS[k];t=time.perf_counter()
    z=np.stack([(q[m][0]-L["thresholds"][m])/L["scales"][m] for m in members],axis=1)
    choice=z.argmax(1);i=np.arange(len(z));p=np.stack([q[m][1] for m in members],axis=1)[i,choice]
    return z[i,choice],p,sum(q[m][2] for m in members)+time.perf_counter()-t,choice
def fit(cal):
    y=np.concatenate([labels(rows) for rows,q in cal])
    L={"thresholds":{},"scales":{},"calibration_controls":int((y==0).sum())}
    for k in BASE:
        x=np.concatenate([q[0][k][0] for rows,q in cal])[y==0];th=cut(x)
        L["thresholds"][k]=th;L["scales"][k]=max(.1,th-float(np.quantile(x,.8)))
    for k in GROUPS:
        s=np.concatenate([fuse(q[0],L,k)[0] for rows,q in cal]);L["thresholds"][k]=cut(s[y==0])
    L["actual_control_alerts"]={}
    for k in S["modes"]:
        v=np.concatenate([q[0][k][0] if k in BASE else fuse(q[0],L,k)[0] for rows,q in cal])
        L["actual_control_alerts"][k]=int(np.sum((v>=L["thresholds"][k])&(y==0)))
        assert L["actual_control_alerts"][k]<=int((y==0).sum()*.05)
    return L
def met(rows,pack,L):
    q,raw,reg,g=pack;out={};detail={};record={}
    y=labels(rows);truth=np.stack([x[7] for x in rows]);active=np.stack([x[8] for x in rows])
    for k in S["modes"]:
        if k in BASE:s,p,t=q[k];choice=None
        else:s,p,t,choice=fuse(q,L,k)
        pred=s>=L["thresholds"][k];loc=((abs(p-truth)<=2)&active).sum(1)>=.75*active.sum(1)
        out[k]={**od.measure(pred,p,np.zeros_like(p,bool),rows),"seconds":raw+t+(0 if k.startswith("raw_") else reg)}
        if choice is not None:
            out[k]["selected_branch_alerts"]={m:int(np.sum(pred&(choice==i))) for i,m in enumerate(GROUPS[k])}
        detail[k]=(pred,loc)
        record[k]={"scores":s.tolist(),"localized":loc.tolist(),"alert":pred.tolist(),"paths":p.tolist()}
    for k in GROUPS:
        pred,loc=detail[k];v=pred&loc&(y==1)
        out[k]["paired"]={}
        for b in ("registered_linear","registered_odonto","registered_taylor"):
            bp,bl=detail[b];u=bp&bl&(y==1);out[k]["paired"][b]={"added_localized":int(np.sum(v&~u)),"lost_localized":int(np.sum(u&~v)),"added_controls":int(np.sum(pred&~bp&(y==0))),"removed_controls":int(np.sum(bp&~pred&(y==0)))}
    # Preserve raw candidate paths and scores (including beneath threshold); no destructive veto.
    return out,{"labels":y.tolist(),"amplitude":[x[4] for x in rows],"model":[x[2] for x in rows],"regime":g.tolist(),"arms":record}
def selftest():
    rng=np.random.default_rng(2801);n=128;rows=[(None,int(i%2),"stationary",0,.75,False,"synthetic",np.full(16,10),np.ones(16,bool)) for i in range(n)]
    q={k:(rng.normal(size=n),np.full((n,16),8+i),0.) for i,k in enumerate(BASE)}
    L=fit([(rows,(q,0.,0.,np.zeros(n,int)))])
    for k in GROUPS:
        z,p,t,c=fuse(q,L,k)
        for i,m in enumerate(GROUPS[k]):
            mask=c==i;assert np.array_equal(p[mask],q[m][1][mask])
        assert np.all(np.isfinite(z))
    return {"joint_calibration":"PASS","selected_witness":"PASS",**a.selftest(),**r.selftest()}
def main():
    jit=a.init_tree();checks=selftest();emit("R28_SELFTEST",checks)
    m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest() and len(m["selected"])==6
    names=[q["target"] for q in m["selected"]];assert len(set(names))==6 and not set(names)&set(hold.S["exclude"])
    seen=set();cal=[];load=0
    for q in m["selected"][:2]:
        t=time.perf_counter();rows,R,_,h=r.gather(q);load+=time.perf_counter()-t;assert not seen&h;seen|=h;cal.append((rows,prep(rows,R)))
    L=fit(cal);L["targets"]=names[:2];save("lock.json",L);emit("R28_LOCK",L)
    hashes=hashlib.sha256((OUT/"lock.json").read_bytes()).hexdigest()
    out={};records={};leak={};leakrec={}
    for q in m["selected"][2:]:
        t=time.perf_counter();rows,R,ch,h=r.gather(q);load+=time.perf_counter()-t;assert not seen&h;seen|=h
        out[q["target"]],records[q["target"]]=met(rows,prep(rows,R),L)
        lr=[x[0] for x in ch];rr=np.stack([x[1] for x in ch]);leak[q["target"]],leakrec[q["target"]]=met(lr,prep(lr,rr),L)
        save("partial.json",out);emit("R28_TARGET",{"target":q["target"],"metrics":out[q["target"]]})
    assert hashes==hashlib.sha256((OUT/"lock.json").read_bytes()).hexdigest()
    keys=("detected","localized","positives","control_alerts","controls","seconds");total={k:{a:sum(x[k][a] for x in out.values()) for a in keys} for k in S["modes"]}
    for k in total:
        for f,vals in (("amplitudes",S["amplitudes"]),("gaps",S["gap_levels"]),("models",S["injection_models"])):
            total[k][f]={str(v):{a:sum(x[k][f][str(v)][a] for x in out.values()) for a in keys[:5]} for v in vals}
        total[k]["weak_localized"]=sum(total[k]["amplitudes"][str(v)]["localized"] for v in (.75,1.25))
        if k in GROUPS:
            total[k]["paired"]={b:{a:sum(x[k]["paired"][b][a] for x in out.values()) for a in ("added_localized","lost_localized","added_controls","removed_controls")} for b in ("registered_linear","registered_odonto","registered_taylor")}
            total[k]["selected_branch_alerts"]={b:sum(x[k]["selected_branch_alerts"][b] for x in out.values()) for b in GROUPS[k]}
    x=total["registered_hybrid"];bs=[total[k] for k in ("registered_linear","registered_odonto")];gate=x["localized"]>max(b["localized"] for b in bs) and x["control_alerts"]<=min(b["control_alerts"] for b in bs)
    d={"status":"FRESH_ARCHIVE_SYNTHETIC_EXAM","gate":"PASS_DEVELOPMENT" if gate else "HOLD_NO_JOINT_GAIN_AT_CONTROL_BUDGET","targets":names[2:],"selftests":checks,"lock":L,"lock_sha256":hashes,"samples":1024,"modes":total,"per_target":out,"leak_challenge":{k:{a:sum(x[k][a] for x in leak.values()) for a in ("localized","detected","positives")} for k in total},"external_fetch_and_jit_seconds":jit,"load_seconds_all_targets":load,"scope":S["scope"],"production_adoption":False}
    save("scores.json",records);save("leak_scores.json",leakrec);save("result.json",d);emit("R28_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    elif sys.argv[1]=="selftest":a.init_tree();emit("R28_SELFTEST",selftest())
    else:main()
