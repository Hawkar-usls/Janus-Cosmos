#!/usr/bin/env python3
"""R29: Star Forge-inspired five-lane scheduling, no exam label feedback."""
import json,pathlib,sys,math,hashlib,time
import numpy as np
import r28b_router as rb
import r21a_weak_template as tmpl
import r28_hybrid as h
import r27a_benchmark as a
import r26c_reference_alignment as r
import r24_odontocrystal as od
import r17c_frozen_holdout as hold
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r29_starforge_protocol.json");S=json.loads(P.read_text());OUT=pathlib.Path("runtime/r29");OUT.mkdir(parents=True,exist_ok=True)
rb.S=S;h.S={**S,"modes":S["modes"][:7]};r.S=S;r.OUT=OUT;od.S=S
hold.S={**hold.S,"seed":S["seed"],"exclude":sorted(set(hold.S["exclude"]+S["exclude_extra"]))};hold.OUT=OUT;hold.P=P;eng.OUT=OUT
LANES=["DAWN","ORBIT","KYBER","MERCY","GROGU"]
MODES=S["star_forge"]["compare"]
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def emit(k,d):print(k,json.dumps(d),flush=True)
def perm(n,tag):
    seed=int.from_bytes(hashlib.sha256((S["seed"]+"|"+tag).encode()).digest()[:8],"big")
    stride=1+seed%max(1,n-1)
    while math.gcd(stride,n)!=1:stride+=1
    p=(seed%n+np.arange(n)*stride)%n
    assert len(np.unique(p))==n
    return p
def kernel(k,X,Y):
    if k=="DAWN":return a.straight(X)
    if k=="ORBIT":return od.sliced(Y,0,None)[:2]
    if k=="KYBER":return tmpl.template_sliced(Y,None)
    if k=="MERCY":return od.sliced(Y,4,None)[:2]
    if k=="GROGU":return od.sliced(X,0,None)[:2]
    raise ValueError(k)
def margin(s,L,k):return (s-L["thresholds"][k])/L["scales"][k]
def pack(rows,R,full):
    t=time.perf_counter();X,Y,g=rb.inputs(rows,R);prep=time.perf_counter()-t
    t=time.perf_counter();s,p=a.straight(Y);dt=time.perf_counter()-t
    q={"registered_linear":(s,p,dt)}
    if full:
        for k in LANES:
            t=time.perf_counter();s,p=kernel(k,X,Y);q[k]=(s,p,time.perf_counter()-t)
        q["raw_linear"]=q["DAWN"];q["registered_odonto"]=q["ORBIT"];q["raw_odonto"]=q["GROGU"]
        t=time.perf_counter();s,p=a.tree(Y);q["registered_taylor"]=(s,p,time.perf_counter()-t)
    return q,X,Y,g,prep
def context(q,g,L):
    z=margin(q["registered_linear"][0],L,"registered_linear")
    return (4*g+np.digitize(z,[-2,0,1])).astype(int)
def merged(z,p,v,w):
    use=v>z;return np.maximum(z,v),np.where(use[:,None],w,p)
def reward(rows,z,p):
    y=h.labels(rows);weak=np.array([x[4]<=1.25 for x in rows])
    return (z>=0)*np.where(y==0,-2,rb.loc(p,rows)*np.where(weak,2,1))
def lane_stats(values,cost):
    v=np.asarray(values,float);n=len(v);ema=0.
    for x in v:ema=(31*ema+x)/32
    mean=float(v.mean());positive=float(np.mean(v>0));peak=float(v.max());hits=int((v>0).sum())
    # Original score anatomy adapted to signed, verified utility, normalized per cost.
    score=(3*ema+2*mean+.15*peak+.5*positive+.05*min(hits,8)+.2/(n+1))/cost
    return {"count":n,"ema":ema,"mean":mean,"peak":peak,"positive_rate":positive,"useful_hits":hits,"score":float(score)}
def fit(cal):
    hcal=[(rows,(q,0.,0.,g)) for rows,(q,X,Y,g,t) in cal]
    L=h.fit(hcal);RL=rb.fit_router(hcal,L)
    y=np.concatenate([h.labels(rows) for rows,_ in cal])
    costs0={k:float(np.median([pack[0][k][2]/len(rows) for rows,pack in cal])) for k in LANES}
    costs={k:max(.05,costs0[k]/costs0["ORBIT"]) for k in LANES}
    for k in LANES:
        x=np.concatenate([q[0][k][0] for _,q in cal])[y==0]
        L["thresholds"][k]=rb.cut(x,.05);L["scales"][k]=max(.1,L["thresholds"][k]-float(np.quantile(x,.8)))
    contexts=[];rewards={k:[] for k in LANES}
    for idx,(rows,(q,X,Y,g,t)) in enumerate(cal):
        order=perm(len(rows),"cal"+str(idx));base=q["registered_linear"];z=margin(base[0],L,"registered_linear")
        base_reward=reward(rows,z,base[1]);contexts.append(context(q,g,L)[order])
        for k in LANES:
            v,w=merged(z,base[1],margin(q[k][0],L,k),q[k][1]);rr=reward(rows,v,w)-base_reward
            rewards[k].extend(rr[order].tolist())
    cx=np.concatenate(contexts);stats={"global":{},"context":{}}
    for k in LANES:stats["global"][k]=lane_stats(rewards[k],costs[k])
    for j in range(12):
        mask=cx==j
        if mask.sum()>=24:stats["context"][str(j)]={k:lane_stats(np.array(rewards[k])[mask],costs[k]) for k in LANES}
    F={"baseline":L,"router":RL,"costs":costs,"cal_seconds_per_patch":costs0,"budget_per_patch":RL["policy"]["dp_calls"]/len(y),"stats":stats,"thresholds":{}}
    for k in MODES[2:]:
        zz=[]
        for idx,(rows,pk) in enumerate(cal):
            out=select(pk,F,k,"cal"+str(idx),True);zz.append(out[0])
        z=np.concatenate(zz);F["thresholds"][k]=rb.cut(z[y==0],.04)
        assert np.sum(z[y==0]>=F["thresholds"][k])<=int(.04*np.sum(y==0))
    return F
def allocate(q,g,F,mode,tag):
    n=len(g);L=F["baseline"];z=margin(q["registered_linear"][0],L,"registered_linear");cx=context(q,g,L);order=perm(n,tag)
    choices=np.full(n,-1,int);budget=F["budget_per_patch"]*n;used=0.
    proposals=[]
    for rank,i in enumerate(order):
        for j,k in enumerate(LANES):
            if mode=="orbit_budget" and k!="ORBIT":continue
            table=F["stats"]["context"].get(str(cx[i]),F["stats"]["global"]) if mode=="star_forge" else F["stats"]["global"]
            # For orbit-only baseline prioritize same R28B uncertainty band, then margin proximity.
            score=float(table[k]["score"]) if mode!="orbit_budget" else (100 if F["router"]["policy"]["low"]<=z[i]<=F["router"]["policy"]["high"] else 0)-abs(float(z[i]))
            proposals.append((-score,rank,j,int(i)))
    for _,rank,j,i in sorted(proposals):
        if choices[i]>=0:continue
        cost=F["costs"][LANES[j]]
        if used+cost<=budget+1e-9:
            choices[i]=j;used+=cost
    assert used<=budget+1e-8
    return choices,used,budget
def select(pk,F,mode,tag,cached=False):
    q,X,Y,g,prep=pk;L=F["baseline"];t=time.perf_counter()
    z=margin(q["registered_linear"][0],L,"registered_linear").copy();p=q["registered_linear"][1].copy()
    choices,used,budget=allocate(q,g,F,mode,tag);calls={}
    for j,k in enumerate(LANES):
        idx=np.flatnonzero(choices==j);calls[k]=int(len(idx))
        if not len(idx):continue
        if cached:s,w,_=q[k];s=s[idx];w=w[idx]
        else:s,w=kernel(k,X[idx],Y[idx])
        z[idx],p[idx]=merged(z[idx],p[idx],margin(s,L,k),w)
    assert sum(calls.values())==int(np.sum(choices>=0))
    return z,p,choices,{"lane_calls":calls,"unique_optional":sum(calls.values()),"possible_optional":5*len(g),"proxy_cost":used,"proxy_budget":budget,"seconds":prep+q["registered_linear"][2]+time.perf_counter()-t}
def exam(rows,R,F,tag):
    pk=pack(rows,R,False);q,X,Y,g,prep=pk;L=F["baseline"];RL=F["router"];n=len(rows);y=h.labels(rows)
    base=q["registered_linear"];z=margin(base[0],L,"registered_linear")
    t=time.perf_counter();ds,dp=kernel("ORBIT",X,Y);v=margin(ds,L,"registered_odonto");zz,pp=merged(z,base[1],v,dp)
    fulltime=prep+base[2]+time.perf_counter()-t
    z1,p1,run,use,dt=rb.actual(rows,R,L,RL["policy"])
    expected,ep,er,eu=rb.route({"registered_linear":base,"registered_odonto":(ds,dp,0)},L,RL["policy"])
    assert np.allclose(z1,expected) and np.array_equal(p1,ep) and np.array_equal(run,er)
    def cv(calls,seconds):return {"lane_calls":{"DAWN":0,"ORBIT":calls,"KYBER":0,"MERCY":0,"GROGU":0},"unique_optional":calls,"possible_optional":5*n,"proxy_cost":calls,"proxy_budget":F["budget_per_patch"]*n,"seconds":seconds}
    vals={"r28_full":(zz,pp,L["thresholds"]["registered_hybrid"],cv(n,fulltime)),"r28b_fixed":(z1,p1,RL["policy"]["threshold"],cv(int(run.sum()),dt))}
    assignments={}
    for mode in MODES[2:]:
        s,p,ch,c=select(pk,F,mode,tag);vals[mode]=(s,p,F["thresholds"][mode],c);assignments[mode]=ch.tolist()
    out={};records={}
    for k,(s,p,th,c) in vals.items():
        alert=s>=th;located=rb.loc(p,rows);out[k]={**od.measure(alert,p,np.zeros_like(p,bool),rows),**c}
        out[k]["weak_localized"]=int(np.sum(alert&located&(y==1)&np.array([x[4]<=1.25 for x in rows])))
        records[k]={"scores":s.tolist(),"paths":p.tolist(),"alert":alert.tolist(),"localized":located.tolist(),"assigned":assignments.get(k)}
    for k in MODES[1:]:
        out[k]["paired"]={}
        for b in ("r28_full","r28b_fixed","fixed_five"):
            if b==k:continue
            x=records[k];v=records[b];aa=np.array(x["alert"]);bb=np.array(v["alert"]);al=aa&np.array(x["localized"])&(y==1);bl=bb&np.array(v["localized"])&(y==1)
            out[k]["paired"][b]={"added_localized":int(np.sum(al&~bl)),"lost_localized":int(np.sum(bl&~al)),"added_controls":int(np.sum(aa&~bb&(y==0))),"removed_controls":int(np.sum(bb&~aa&(y==0)))}
    return out,{"labels":y.tolist(),"amp":[x[4] for x in rows],"arms":records}
def main():
    jit=a.init_tree();checks=h.selftest()
    for n in (1,2,7,16,127,256,257):assert sorted(perm(n,"test").tolist())==list(range(n))
    checks["coprime_coverage"]="PASS"
    m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest() and len(m["selected"])==6
    names=[q["target"] for q in m["selected"]];assert len(set(names))==6 and not set(names)&set(hold.S["exclude"])
    cal=[];seen=set()
    for q in m["selected"][:2]:
        rows,R,ch,hs=r.gather(q);assert not seen&hs;seen|=hs;cal.append((rows,pack(rows,R,True)))
    F=fit(cal);F["targets"]=names[:2];save("lock.json",F);lock=hashlib.sha256((OUT/"lock.json").read_bytes()).hexdigest()
    emit("R29_LOCK",{"targets":names[:2],"lock_sha256":lock,"costs":F["costs"],"budget":F["budget_per_patch"],"policy":F["router"]["policy"],"context_count":len(F["stats"]["context"]),"thresholds":F["thresholds"]})
    out={};records={};leak={}
    for q in m["selected"][2:]:
        rows,R,ch,hs=r.gather(q);assert not seen&hs;seen|=hs;target=q["target"]
        out[target],records[target]=exam(rows,R,F,target)
        lr=[x[0] for x in ch];rr=np.stack([x[1] for x in ch]);leak[target]=exam(lr,rr,F,target+"-leak")[0]
        save("partial.json",out);emit("R29_TARGET",{"target":target,"metrics":out[target]})
    assert hashlib.sha256((OUT/"lock.json").read_bytes()).hexdigest()==lock
    keys=("localized","detected","positives","control_alerts","controls","weak_localized","seconds","unique_optional","possible_optional","proxy_cost","proxy_budget")
    total={k:{a:sum(x[k][a] for x in out.values()) for a in keys} for k in MODES}
    for k in total:
        total[k]["lane_calls"]={j:sum(x[k]["lane_calls"][j] for x in out.values()) for j in LANES}
        total[k]["coverage"]=total[k]["unique_optional"]/total[k]["possible_optional"]
        if k!="r28_full":
            total[k]["paired"]={b:{a:sum(x[k]["paired"][b][a] for x in out.values()) for a in ("added_localized","lost_localized","added_controls","removed_controls")} for b in out[names[2]][k]["paired"]}
        total[k]["amplitudes"]={str(v):{a:sum(x[k]["amplitudes"][str(v)][a] for x in out.values()) for a in keys[:5]} for v in S["amplitudes"]}
    x=total["star_forge"];gate=all(x["localized"]>=total[k]["localized"] and x["control_alerts"]<=total[k]["control_alerts"] for k in ("r28_full","r28b_fixed","fixed_five"))
    d={"status":"FRESH_ARCHIVE_SYNTHETIC_EXAM","gate":"PASS_DEVELOPMENT" if gate else "HOLD","samples":1024,"targets":names[2:],"selftests":checks,"lock_sha256":lock,"costs":F["costs"],"budget_per_patch":F["budget_per_patch"],"router_policy":F["router"]["policy"],"modes":total,"per_target":out,"leak_challenge":{k:{a:sum(v[k][a] for v in leak.values()) for a in ("localized","positives")} for k in total},"external_fetch_and_jit_seconds":jit,"scope":S["scope"],"production_adoption":False}
    save("scores.json",records);save("result.json",d);emit("R29_SUMMARY",{k:v for k,v in d.items() if k!="per_target"})
if __name__=="__main__":
    if sys.argv[1]=="freeze":hold.freeze()
    else:main()
