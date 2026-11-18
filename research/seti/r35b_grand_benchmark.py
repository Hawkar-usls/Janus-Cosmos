#!/usr/bin/env python3
"""Fresh archive paired benchmark; template hypotheses and exact JUXTAPOSE decision ablation."""
import hashlib,json,pathlib,resource,sys,time,itertools
import numpy as np
import r30_false_alarm_shield as sh
import r28_hybrid as h
import r28b_router as rb
import r26c_reference_alignment as ref
import r27a_benchmark as taylor
import r24_odontocrystal as od
import juxtapose_r18a as j
P=pathlib.Path("research/seti/r35b_grand_benchmark_protocol.json");S=json.loads(P.read_text());B=S["r35"];OUT=pathlib.Path("runtime/r35b");OUT.mkdir(parents=True,exist_ok=True)
sh.P=P;sh.S=S;sh.C=S["shield"];sh.OUT=OUT;h.S=S;rb.S=S;ref.S=S;od.S=S
ref.OUT=OUT;sh.eng.OUT=OUT
sh.hold.S={**sh.hold.S,"seed":S["seed"],"exclude":sorted(set(sh.hold.S["exclude"]+S["exclude_extra"]))};sh.hold.P=P;sh.hold.OUT=OUT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(k,d):print(k,json.dumps(d,allow_nan=False),flush=True)
def bank():
    t=np.arange(16);u=t-7.5;b=B["bank"]
    arr=[np.rint(t*s/15).astype(int) for s in range(-32,33)]
    arr += [np.rint(s*u+a*u*u).astype(int) for s,a in itertools.product(b["slopes"],b["curvature"])]
    arr += [np.rint(s*u+v*np.sin(2*np.pi*t/16+k*2*np.pi/b["phases"])).astype(int) for s,v,k in itertools.product(b["sinusoid_slopes"],b["sinusoid_amplitude"],range(b["phases"]))]
    paths=[x-x[0] for x in arr if np.max(abs(np.diff(x)))<=b["max_step"]]
    return np.asarray(sorted(set(map(tuple,paths))),int)
D=bank()
def template(V,kind):
    x=V
    if kind=="clip":x=np.clip(x,-3,3)
    if kind=="neighborhood":
        # Fixed weighted adjacent-bin combination, no circular wrap.
        x=.5*x.copy();x[:,:,1:]+=.25*V[:,:,:-1];x[:,:,:-1]+=.25*V[:,:,1:]
    N,n,w=x.shape;best=np.full(N,-np.inf);path=np.zeros((N,n),int)
    for d in D:
        lo=max(0,-int(d.min()));hi=min(w,w-int(d.max()))
        if kind=="neighborhood":lo+=1;hi-=1
        if hi<=lo:continue
        z=sum(x[:,i,lo+k:hi+k] for i,k in enumerate(d))/np.sqrt(n)
        pos=z.argmax(1);v=z[np.arange(N),pos];m=v>best
        best[m]=v[m];path[m]=lo+pos[m,None]+d
    assert np.isfinite(best).all()
    replay=x[np.arange(N)[:,None],np.arange(n)[None,:],path].sum(1)/np.sqrt(n)
    assert np.allclose(replay,best,atol=1e-8,rtol=1e-10)
    return best,path
def panel(rows,R):
    pack,Y=rb.panel(rows,R);q,raw,reg,g=pack
    X=ref.norm(np.stack([v[0] for v in rows])).astype(np.float32).astype(np.float64)
    ext={}
    for name,V,kind in (("raw_smooth_bank",X,"raw"),("registered_smooth_bank",Y,"raw"),("registered_smooth_clipped",Y,"clip"),("registered_smooth_neighborhood",Y,"neighborhood")):
        start=time.perf_counter();s,p=template(V,kind);ext[name]=(s,p,time.perf_counter()-start)
    return {"rows":rows,"pack":pack,"Y":Y,"ext":ext,"rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
def base_arms(d,L):
    q,raw,reg,_=d["pack"];z,p,t,_=h.fuse(q,L,"registered_hybrid")
    k={"hybrid_r28":(z,p,raw+reg+t)}
    for label,key in (("native_taylor_core","registered_taylor"),("registered_linear","registered_linear"),("registered_odonto","registered_odonto")):
        s,p,t=q[key];k[label]=(s,p,raw+reg+t)
    for name,(s,p,t) in d["ext"].items():k[name]=(s,p,raw+t+(0 if name.startswith("raw_") else reg))
    return k
def labels(d):return np.asarray([r[1] for r in d["rows"]])
def controls(parts,arm,L):
    return np.concatenate([base_arms(v,L)[arm][0][labels(v)==0] for v in parts])
def cut(v):return sh.cut(v,.05)
def fused(d,L,marg):
    a=base_arms(d,L);keys=("hybrid_r28","registered_smooth_bank","registered_smooth_clipped","registered_smooth_neighborhood","raw_smooth_bank")
    z=np.stack([(a[k][0]-marg[k][0])/marg[k][1] for k in keys],axis=1);i=z.argmax(1);ix=np.arange(len(i));paths=np.stack([a[k][1] for k in keys],axis=1)
    raw,reg=d["pack"][1:3]
    cost=a["hybrid_r28"][2]+sum(d["ext"][k][2] for k in keys[1:])
    return z[ix,i],paths[ix,i],cost
def allarms(d,L,marg):
    a=base_arms(d,L);a["hybrid_plus_smooth"]=fused(d,L,marg);return a
def jrank(fit,threshold):
    # Original exact JUXTAPOSE shift groups, only fit examples. No evaluation tuning.
    hits=[]
    for group in j.GROUPS:
        n=0
        for part in fit:
            for v in part["Y"]:
                n+=int(any(j.score(v,s)[0]>=threshold for s in group))
        hits.append(n)
    return sorted(range(8),key=lambda i:(-hits[i]/len(j.GROUPS[i]),i)),hits
def jtest(V,threshold,rank):
    def run(order,early):
        start=time.perf_counter();pred=[];cost=0
        for v in V:
            ok=False
            for s in order:
                a,b=j.score(v,s);cost+=1
                if a>=threshold:
                    assert abs(j.witness(v,s,b)-a)<=max(1e-7,abs(a)*1e-9)
                    ok=True
                    if early:break
            pred.append(ok)
        return np.asarray(pred),cost,time.perf_counter()-start
    e,c,t=run(list(range(-32,33)),False);a,cc,tt=run([s for i in rank for s in j.GROUPS[i]],True)
    assert np.array_equal(e,a)
    return {"samples":len(V),"alerts":int(a.sum()),"exhaustive_checks":c,"early_checks":cc,"exhaustive_seconds":t,"early_seconds":tt,"decision_equivalent":True}
def metrics(d,a,threshold):
    rows=d["rows"];s,p,t=a;y=labels(d);truth=np.stack([r[7] for r in rows]);active=np.stack([r[8] for r in rows]);amp=np.asarray([r[4] for r in rows]);pred=s>=threshold;loc=((abs(p-truth)<=2)&active).sum(1)>=.75*active.sum(1);z=pred&loc&(y==1)
    out={"localized":int(z.sum()),"weak":int((z&(amp<=1.25)).sum()),"very_weak":int((z&(amp==.75)).sum()),"control_alerts":int((pred&(y==0)).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum()),"seconds":float(t),"models":{k:int((z&np.array([r[2]==k for r in rows])).sum()) for k in S["injection_models"]}}
    return out,{"alert":pred,"localized_positive":z,"paths":p,"scores":s,"weak":amp<=1.25}
def selftest():
    rng=np.random.default_rng(3511);v=rng.normal(size=(4,16,128));p=48+D[len(D)//2];v[0,np.arange(16),p]+=20;s,q=template(v,"raw");assert np.array_equal(q[0],p)
    assert len(D)>65 and len(set(tuple(x) for x in D))==len(D)
    r=jtest(v,12.,list(range(8)));assert r["decision_equivalent"]
    return {"template_exact_witness":"PASS","bank_unique":"PASS","juxtapose_exact_decision":"PASS"}
def main():
    init=taylor.init_tree();checks={**h.selftest(),**selftest()}
    mf=json.loads((OUT/"manifest.json").read_text());assert mf["protocol_sha256"]==sha(P)
    selected=mf["selected"];assert len(selected)==8 and len({q["day"] for q in selected})==8
    assert [q["split"] for q in selected]==["component_fit"]*2+["null_calibration"]*2+["evaluation"]*4
    seen=set();total_load=0
    def acquire(q):
        nonlocal total_load
        st=time.perf_counter();rows,R,ch,hs=ref.gather(q);total_load+=time.perf_counter()-st
        assert not seen&hs;seen.update(hs)
        return panel(rows,R),ch
    fit=[acquire(q)[0] for q in selected[:2]];L=h.fit([(d["rows"],d["pack"]) for d in fit])
    marg={}
    for k in B["arms"][:-1]:
        v=controls(fit,k,L);th=cut(v);marg[k]=[th,max(.1,th-float(np.quantile(v,.8)))]
    null=[acquire(q)[0] for q in selected[2:4]]
    thresholds={}
    for k in B["arms"]:
        zz=np.concatenate([allarms(d,L,marg)[k][0][labels(d)==0] for d in null]);thresholds[k]=cut(zz)
    rank,hits=jrank(fit,thresholds["registered_linear"])
    lock={"component_fit":L,"margins":marg,"thresholds":thresholds,"juxtapose_rank":rank,"fit_group_hits":hits,"bank":D.tolist(),"manifest_sha256":sha(OUT/"manifest.json")}
    save("lock.json",lock);locksha=sha(OUT/"lock.json");emit("R35B_LOCK",{"sha256":locksha,"thresholds":thresholds,"bank_size":len(D),"juxtapose_rank":rank})
    per={};leak={};record={};js=[]
    for q in selected[4:]:
        d,ch=acquire(q);name=q["target"];dd=panel([r[0] for r in ch],np.stack([r[1] for r in ch]));per[name]={};leak[name]={};rr={}
        for k,a in allarms(d,L,marg).items():per[name][k],rr[k]=metrics(d,a,thresholds[k])
        lr={}
        for k,a in allarms(dd,L,marg).items():leak[name][k],lr[k]=metrics(dd,a,thresholds[k])
        for k in B["arms"]:
            x=rr[k]["localized_positive"];base=rr["hybrid_r28"]["localized_positive"]
            per[name][k].update(added=int((x&~base).sum()),lost=int((base&~x).sum()),lost_weak=int((base&~x&rr[k]["weak"]).sum()),lost_shared=int((lr["hybrid_r28"]["localized_positive"]&~lr[k]["localized_positive"]).sum()))
        js.append(jtest(d["Y"],thresholds["registered_linear"],rank))
        record[name]={k:{z:v.tolist() for z,v in rr[k].items()} for k in rr}
        # Export common registered data for native pipeline followup, no hidden labels used in detector.
        path=OUT/"common_inputs";path.mkdir(exist_ok=True)
        rb.OUT=OUT
        np.savez_compressed(path/("eval_"+name+".npz"),V=d["Y"].astype(np.float32),y=labels(d),truth=np.stack([r[7] for r in d["rows"]]),active=np.stack([r[8] for r in d["rows"]]))
        save("partial.json",per);emit("R35B_TARGET",{"target":name,"metrics":per[name]})
    assert sha(OUT/"lock.json")==locksha
    total={k:{f:sum(v[k][f] for v in per.values()) for f in ("localized","weak","very_weak","control_alerts","positives","controls","seconds","added","lost","lost_weak","lost_shared")} for k in B["arms"]}
    for k in total:
        total[k]["models"]={m:sum(v[k]["models"][m] for v in per.values()) for m in S["injection_models"]};total[k]["shared_localized"]=sum(v[k]["localized"] for v in leak.values())
    assert all(v["positives"]==v["controls"]==512 for v in total.values())
    b=total["hybrid_r28"];v=total["hybrid_plus_smooth"]
    gate={"overall_gain":v["localized"]>b["localized"],"weak_gain":v["weak"]>b["weak"],"controls_not_higher":v["control_alerts"]<=b["control_alerts"],"shared_not_lower":v["shared_localized"]>=b["shared_localized"],"zero_paired_weak_loss":v["lost_weak"]==0,"zero_shared_loss":v["lost_shared"]==0}
    jtot={k:sum(x[k] for x in js) for k in ("samples","alerts","exhaustive_checks","early_checks","exhaustive_seconds","early_seconds")}
    result={"status":"FRESH_ARCHIVE_SYNTHETIC_DEVELOPMENT","gate":"PASS_DEVELOPMENT" if all(gate.values()) else "HOLD","gate_checks":gate,"modes":total,"per_target":per,"juxtapose":{**jtot,"decision_equivalent":True},"bank_size":len(D),"selftests":checks,"peak_rss_kib_process":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"archive_load_seconds":total_load,"upstream_compile_seconds":init,"lock_sha256":locksha,"scope":S["scope"],"production_adoption":False}
    save("scores.json",record);save("result.json",result);emit("R35B_SUMMARY",{k:v for k,v in result.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":sh.freeze()
    elif sys.argv[1]=="selftest":emit("R35B_SELFTEST",selftest())
    else:main()
