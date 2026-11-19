#!/usr/bin/env python3
"""R41P: Penrose-inspired phase-space scheduling of exact native rival search."""
import sys
ARGS=sys.argv[1:]
import json,time,pathlib,collections,math,statistics
import numpy as np
from numba import njit
import r41b_minimal_hybrid as base
a=base.a;n=base.n;bankmod=base.b;c=base.c
P=pathlib.Path("research/seti/r41p_penrose_billiard_protocol.json");S=json.loads(P.read_text());G=S["geometry"]
O=pathlib.Path("runtime/r41p");O.mkdir(parents=True,exist_ok=True);PROG=O/"progress";PROG.mkdir(exist_ok=True)
METHODS=tuple(S["methods"]);KINDS=("eval","shared")
def save(k,d):a.write(O/k,d)
def load(k):return json.loads((O/k).read_text())
def cp(k,d):
    a.write(PROG/(k+".json"),d);print("R41P_PROGRESS",json.dumps(d,allow_nan=False),flush=True)
def native_bank():
    paths=[];offset=[];counts=[];templates=[]
    for delta in bankmod.D:
        delta=np.asarray(delta,dtype=np.int64);lo=5-int(delta.min());hi=123-int(delta.max())
        if hi<=lo:continue
        offset.append(len(paths));counts.append(hi-lo);templates.append(delta)
        paths.extend(np.arange(lo,hi)[:,None]+delta)
    paths=np.ascontiguousarray(paths,dtype=np.int64);counts=np.array(counts,dtype=int);offset=np.array(offset,dtype=int)
    # Preserve original enumeration for tie-breaking; geometry uses only a separate rank.
    rank=sorted(range(len(templates)),key=lambda k:(int(templates[k][-1]-templates[k][0]),int(np.abs(np.diff(templates[k],n=2)).sum()),tuple(templates[k])))
    sector=np.empty(len(paths),dtype=np.int64)
    for j,k in enumerate(rank):
        v=min(7,8*j//len(rank))
        for i in range(counts[k]):sector[offset[k]+i]=8*v+min(7,8*i//counts[k])
    assert paths.ndim==2 and paths.shape[1]==16 and np.all((paths>=5)&(paths<123)) and len(paths)>=512
    return paths,offset,counts,np.array(rank),sector
def billiard(seed,offset,counts,rank):
    rng=np.random.default_rng(seed);A=G["a"];B=G["b"];scale=np.array([A,B]);angle=rng.uniform(0,2*np.pi)
    q=np.array([A*np.cos(angle),B*np.sin(angle)])
    normal=q/(scale*scale);normal/=np.linalg.norm(normal);tangent=np.array([-normal[1],normal[0]])
    theta=rng.uniform(-1.32,1.32);v=-np.cos(theta)*normal+np.sin(theta)*tangent
    found=[];seen=set();err=0.;velerr=0.
    for _ in range(G["bounces_per_orbit"]):
        dt=-2*np.dot(q/scale,v/scale)/np.dot(v/scale,v/scale)
        assert dt>0 and np.isfinite(dt)
        end=q+dt*v;err=max(err,abs(float(np.dot(end/scale,end/scale))-1))
        for frac in G["chord_samples"]:
            u=q+frac*(end-q);yy=np.clip(u[1]/B,-1+1e-12,1-1e-12)
            xx=np.clip(u[0]/(A*np.sqrt(1-yy*yy)),-1,1)
            k=rank[min(len(rank)-1,int((yy+1)*.5*len(rank)))]
            idx=int(offset[k]+min(counts[k]-1,int((xx+1)*.5*counts[k])))
            if idx not in seen:found.append(idx);seen.add(idx)
        # Analytic specular reflection; remove numerical drift onto ellipse only.
        q=end/np.linalg.norm(end/scale);normal=q/(scale*scale);normal/=np.linalg.norm(normal)
        old=v.copy();v=v-2*np.dot(v,normal)*normal
        velerr=max(velerr,abs(float(np.linalg.norm(v))-1))
        assert np.isclose(np.dot(old,normal),-np.dot(v,normal),rtol=0,atol=1e-12)
        v/=np.linalg.norm(v)
    assert err<1e-10 and velerr<1e-10
    return found,{"unique":len(found),"max_boundary_error":err,"max_speed_error":velerr}
def schedules(paths,offset,counts,rank,sector):
    N=len(paths);rng=np.random.default_rng(G["seed"]);rnd=rng.permutation(N)
    queues=[collections.deque(int(i) for i in rnd if sector[i]==s) for s in range(64)]
    strat=[];sq=[q.copy() for q in queues]
    while len(strat)<N:
        for q in sq:
            if q:strat.append(q.popleft())
    orbits=[];reports=[]
    for j in range(G["orbits"]):
        x,d=billiard(G["seed"]+j*9973,offset,counts,rank);orbits.append(x);reports.append(d)
    multi=[];seen=set()
    for k in range(max(map(len,orbits))):
        for x in orbits:
            if k<len(x) and x[k] not in seen:multi.append(x[k]);seen.add(x[k])
    def complete(x):
        ss=set(x);return np.array(list(x)+[i for i in range(N) if i not in ss],dtype=np.int64)
    dark=[];used=set();vis=np.zeros(64,int);q0=collections.deque(multi);forced=0;native=0
    while len(dark)<N:
        chosen=None;force=(len(dark)+1)%G["darkfill_every"]==0 or not q0
        if not force:
            while q0 and q0[0] in used:q0.popleft()
            if q0:chosen=q0.popleft();native+=1
        if chosen is None:
            for s in sorted(range(64),key=lambda s:(vis[s],s)):
                q=queues[s]
                while q and q[0] in used:q.popleft()
                if q:chosen=q.popleft();forced+=1;break
        assert chosen is not None and chosen not in used
        used.add(chosen);dark.append(chosen);vis[sector[chosen]]+=1
    orders={"canonical":np.arange(N,dtype=np.int64),"random":np.asarray(rnd,dtype=np.int64),"stratified":np.array(strat,dtype=np.int64),"ellipse_single":complete(orbits[0]),"ellipse_multi":complete(multi),"ellipse_multi_darkfill":np.array(dark,dtype=np.int64)}
    budgets=sorted(set(min(N,x) for x in G["budgets"]+[N]))
    coverage={}
    for key,order in orders.items():
        assert np.array_equal(np.sort(order),np.arange(N))
        coverage[key]={str(B):{"covered_sectors":int(len(set(sector[order[:B]]))),"sector_counts":np.bincount(sector[order[:B]],minlength=64).tolist()} for B in budgets}
    return orders,budgets,{"bank_size":N,"templates":len(rank),"sectors_nonempty":int(len(set(sector))),"single_orbit_unique":len(orbits[0]),"multi_orbit_unique":len(multi),"darkfill_billiard_entries":native,"darkfill_forced_entries":forced,"orbit_checks":reports,"coverage":coverage}
@njit(cache=False,fastmath=False)
def scan(U,V,p,paths,order,budgets):
    N=len(p);K=len(budgets)
    best=np.full((K,N,2),-np.inf);ids=np.full((K,N,2),-1,np.int64)
    # counts: all geometric attempts, legal candidates, score terms, pruned sides, full sides.
    counts=np.zeros((K,5),np.int64)
    for i in range(N):
        suffix=np.zeros((2,17),np.float64)
        for j in range(15,-1,-1):
            u=-np.inf;v=-np.inf
            for f in range(5,123):
                u=max(u,U[i,j,f]);v=max(v,V[i,j,f])
            suffix[0,j]=suffix[0,j+1]+u;suffix[1,j]=suffix[1,j+1]+v
        bs=np.full(2,-np.inf);bi=np.full(2,-1,np.int64);cts=np.zeros(5,np.int64);k=0
        for q0 in range(len(order)):
            j=order[q0];cts[0]+=1;distinct=0
            for row in range(16):
                if abs(paths[j,row]-p[i,row])>2:distinct+=1
            if distinct>=5:
                cts[1]+=1
                for side in range(2):
                    total=0.;pruned=False
                    for row in range(16):
                        total+=U[i,row,paths[j,row]] if side==0 else V[i,row,paths[j,row]]
                        cts[2]+=1
                        # Strict upper bound, huge guard versus roundoff of <=16 terms in [-4,4].
                        # Never prune a tied incumbent; canonical ID selects ties regardless of visit order.
                        if row<15 and total+suffix[side,row+1]<bs[side]-1e-10:
                            cts[3]+=1;pruned=True;break
                    if not pruned:
                        cts[4]+=1
                        if total>bs[side] or (total==bs[side] and (bi[side]<0 or j<bi[side])):
                            bs[side]=total;bi[side]=j
            if q0+1==budgets[k]:
                best[k,i,:]=bs/16.;ids[k,i,:]=bi;counts[k]+=cts;k+=1
                if k==K:break
        assert k==K
    return best,ids,counts
def witness_from(U,V,p,paths,best,ids):
    use=(ids>=0).all(1);safe=np.maximum(ids,0);rp=paths[safe]
    b0=np.where(use[:,None],best,0.);ii=np.arange(len(p))[:,None];rr=np.arange(16)[None,:]
    valid=(p>=5)&(p<123);pp=np.clip(p,5,122)
    x=np.where(valid,U[ii,rr,pp],0.);y=np.where(valid,V[ii,rr,pp],0.)
    den=np.maximum(1,valid.sum(1))
    rv=np.column_stack([b0,x.sum(1)/den-b0[:,0],y.sum(1)/den-b0[:,1],b0[:,0]-b0[:,1],np.mean(abs(rp[:,0]-rp[:,1]),axis=1)])
    return rv,rp,use
def synthetic(paths,orders,budgets):
    rng=np.random.default_rng(41209);X=rng.normal(size=(3,16,128));D=rng.normal(size=(3,16,128));X[1]=0;D[1]=0
    p=np.tile(np.arange(16)+40,(3,1));p[2,:8]=0
    U=n.normalized_contrast(X);V=n.normalized_contrast(D)
    ref,supported,rp,bs=base.rival_only(X,D,p);start=time.perf_counter()
    for key in METHODS:
        scores,ix,cts=scan(U,V,p,paths,orders[key],np.array(budgets,np.int64));rv,path,use=witness_from(U,V,p,paths,scores[-1],ix[-1]);assert use.all()
        assert np.array_equal(scores[-1],bs) and np.array_equal(path,rp) and np.array_equal(rv,ref),"synthetic full parity"
        # Independent brute oracle at every prefix, sequential sums and canonical tie IDs.
        for k,B in enumerate(budgets):
            for j in range(len(p)):
                cand=orders[key][:B];cand=cand[(abs(paths[cand]-p[j])>2).sum(1)>=5];cand=np.sort(cand)
                if not len(cand):assert (ix[k,j]<0).all();continue
                for side,Z in enumerate((U,V)):
                    pp=paths[cand];val=sum(Z[j,row,pp[:,row]] for row in range(16))/16
                    at=int(val.argmax());assert scores[k,j,side]==val[at] and ix[k,j,side]==cand[at],"prefix oracle"
    return {"synthetic_rows":3,"methods":len(METHODS),"prefix_oracle":True,"full_tie_parity":True,"seconds":time.perf_counter()-start}
def prepare():
    start=time.perf_counter();base.prepare();paths,offset,counts,rank,sector=native_bank()
    t0=time.perf_counter();orders,budgets,report=schedules(paths,offset,counts,rank,sector);elapsed=time.perf_counter()-t0
    report["schedule_seconds"]=elapsed;report["synthetic"]=synthetic(paths,orders,budgets)
    np.savez_compressed(O/"schedule.npz",paths=paths,budgets=np.array(budgets),sector=sector,**orders)
    save("schedule_receipt.json",{"protocol_sha256":c.sha(P.read_bytes()),"schedule_sha256":c.sha((O/"schedule.npz").read_bytes()),"geometry":report})
    cp("00-prepare",{"stage":"PREPARE","groups":11,"rows":3168,"bank_size":len(paths),"budgets":budgets,"coverage":{k:{B:v["covered_sectors"] for B,v in d.items()} for k,d in report["coverage"].items()},"single_orbit_unique":report["single_orbit_unique"],"multi_orbit_unique":report["multi_orbit_unique"],"synthetic":report["synthetic"],"seconds":time.perf_counter()-start})
def group(i):
    start=time.perf_counter();p0=base.plan();g=p0["groups"][i];num=str(i+1).zfill(2)
    rec=load("schedule_receipt.json");assert rec["protocol_sha256"]==c.sha(P.read_bytes()) and rec["schedule_sha256"]==c.sha((O/"schedule.npz").read_bytes())
    data=a.read_npz(O/"schedule.npz");paths=data["paths"];budgets=data["budgets"];m,L,old=base.locks("minimal",p0)
    seed={"n6":a.v.C["seed"],"n7":base.w.C["seed"],"n9":base.old9.C["seed"]}[g["source"]];n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);bankmod.taylor.init_tree()
    t0=time.perf_counter();rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r41p_"+num);acq=time.perf_counter()-t0
    assert sorted(hs)==g["raw_tile_hashes"]
    t0=time.perf_counter();scan(np.zeros((1,16,128)),np.zeros((1,16,128)),np.full((1,16),64,np.int64),paths,data["canonical"],budgets);jit=time.perf_counter()-t0
    out={"stage":"GROUP","group":i+1,"target":g["group"]["target"],"fit_viewed":g["fit"],"parity":True,"acquisition_seconds":acq,"jit_plus_warmup_seconds":jit,"kinds":{}}
    for kind,rows,refs,mm in (("eval",rr,R,meta),("shared",[q[0] for q in ch],np.stack([q[1] for q in ch]),cm)):
        t0=time.perf_counter();d=base.infer("minimal",rows,refs,mm,pre,L,old,m,None,kind);reference_seconds=time.perf_counter()-t0
        ar=a.read_npz(base.C/("dev_"+num+"_"+kind+".npz"));a.verify_rows(ar,rows);reference=base.check(d,ar,m)
        idx=d["idx"];pp=d["parent_path"][idx].astype(np.int64)
        X=n.ref.norm(np.stack([rows[q][0] for q in idx])) if len(idx) else np.empty((0,16,128));D=n.ref.norm(refs[idx]) if len(idx) else X.copy()
        t0=time.perf_counter();U=n.normalized_contrast(X);V=n.normalized_contrast(D);norm_seconds=time.perf_counter()-t0
        t0=time.perf_counter()
        if len(idx):frv,fsup,frp,fbs=base.rival_only(X,D,pp);assert np.array_equal(frv,d["rival"]) and np.array_equal(frp,d["rival_paths"])
        reference_rival_seconds=time.perf_counter()-t0
        aa=d["parent_alert"];hh=d["parent_hit"];res=d["joint_rescue"];co=d["correct"];null=d["y"]==0;weak=d["weak"];jm=hh|(res&co)
        stress=ar["repeat_stress_supported"]&ar["supported"]&ar["fold_supported"]&ar["halo_supported"]&ar["repeat_supported"]&hh[idx]
        stats={}
        # Counterbalance order across groups for timing, no outcome-dependent schedule.
        order=METHODS[i%len(METHODS):]+METHODS[:i%len(METHODS)]
        for key in order:
            t0=time.perf_counter();bs,ids,cts=scan(U,V,pp,paths,data[key],budgets);sec=time.perf_counter()-t0
            t0=time.perf_counter();bb,ii,cc=scan(U,V,pp,paths,data[key],budgets);sec2=time.perf_counter()-t0
            assert np.array_equal(bs,bb) and np.array_equal(ids,ii) and np.array_equal(cts,cc)
            rv,rp,use=witness_from(U,V,pp,paths,bs[-1],ids[-1]);assert use.all()
            assert np.array_equal(bs[-1],d["rival_scores"]) and np.array_equal(rp,d["rival_paths"]) and np.array_equal(rv,d["rival"]),"actual exact rival parity "+key
            bb={}
            for k,B in enumerate(budgets):
                rv,rp,available=witness_from(U,V,pp,paths,bs[k],ids[k])
                score=a.score(np.column_stack([rv,d["runs"]]),m["model"],m["kind"])
                available=available&d["supported"]&d["fold_supported"];keep=np.ones(len(aa),bool);keep[idx]=(~available)|(score>=m["cuts"]["preserve"])
                alert=aa&keep|res;hit=hh&keep|(res&co)
                gap=np.where(np.isfinite(bs[k]),bs[-1]-np.where(np.isfinite(bs[k]),bs[k],0),0)
                # Score/path tests include only real legal candidates; no reconstructed radiodata.
                bb[str(int(B))]={"parents":len(idx),"available_parents":int(available.sum()),"both_winners":int((ids[k]==ids[-1]).all(1).sum()),"winner_sides":int((ids[k]==ids[-1]).sum()),"score_gap_sum":float(gap.sum()),"score_gap_count":int(np.isfinite(bs[k]).sum()),"score_gap_max":float(gap.max()) if gap.size else 0.,"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"controls":int((alert&null).sum()),"lost_joint":int((jm&~hit).sum()),"lost_vs_minimal":int((d["hit"]&~hit).sum()),"changed_decisions":int((alert!=d["alert"]).sum()),"repeat_stress_losses":int((stress&~keep[idx]).sum()),"attempted":int(cts[k,0]),"legal":int(cts[k,1]),"terms":int(cts[k,2]),"pruned_sides":int(cts[k,3]),"fully_scored_sides":int(cts[k,4])}
                if k==len(budgets)-1:
                    assert np.array_equal(score,d["score"]) and np.array_equal(keep,d["keep"]) and np.array_equal(alert,d["alert"]) and np.array_equal(hit,d["hit"]),"exact detector parity "+key
            stats[key]={"kernel_wall_seconds":min(sec,sec2),"paired_repeat_wall_seconds":[sec,sec2],"budgets":bb}
        out["kinds"][kind]={"reference":reference,"reference_full_infer_seconds":reference_seconds,"reference_rival_seconds":reference_rival_seconds,"contrast_seconds":norm_seconds,"methods":stats}
    out["seconds"]=time.perf_counter()-start;save("group-"+num+".json",out)
    cp("group-"+num,{"stage":"GROUP","group":i+1,"total_groups":11,"target":out["target"],"fit_viewed":g["fit"],"parity":True,"acq_seconds":acq,"seconds":out["seconds"],"methods":{key:{"exact_terms":sum(out["kinds"][k]["methods"][key]["budgets"][str(int(budgets[-1]))]["terms"] for k in KINDS),"exact_kernel_seconds":sum(out["kinds"][k]["methods"][key]["kernel_wall_seconds"] for k in KINDS)} for key in METHODS}})
def summary():
    rows=[load("group-"+str(i+1).zfill(2)+".json") for i in range(11)];assert all(d["parity"] for d in rows);audit=[d for d in rows if not d["fit_viewed"]];assert len(audit)==5
    sr=load("schedule_receipt.json");budgets=[str(v) for v in sorted(set(min(sr["geometry"]["bank_size"],x) for x in G["budgets"]+[sr["geometry"]["bank_size"]]))];last=budgets[-1]
    sums={}
    for key in METHODS:
        sums[key]={}
        for kind in KINDS:
            r0=[d["kinds"][kind]["methods"][key] for d in audit]
            agg={}
            for B in budgets:
                rr=[v["budgets"][B] for v in r0];agg[B]={k:max(v[k] for v in rr) if k=="score_gap_max" else sum(v[k] for v in rr) for k in rr[0]}
            sums[key][kind]={"budgets":agg,"audit_kernel_seconds":sum(v["kernel_wall_seconds"] for v in r0),"all_11_kernel_seconds":sum(d["kinds"][kind]["methods"][key]["kernel_wall_seconds"] for d in rows)}
    baseline=sums["canonical"]["eval"]["budgets"][last];sh=sums["canonical"]["shared"]["budgets"][last]
    assert (baseline["localized"],sh["localized"],baseline["controls"],baseline["lost_joint"],sh["lost_joint"])==(205,113,53,0,0)
    refsec=sum(d["kinds"][k]["reference_rival_seconds"] for d in rows for k in KINDS)
    normsec=sum(d["kinds"][k]["contrast_seconds"] for d in rows for k in KINDS)
    can=sum(sums["canonical"][k]["all_11_kernel_seconds"] for k in KINDS);exact={}
    for key in METHODS:
        e=sums[key]["eval"]["budgets"][last];sh=sums[key]["shared"]["budgets"][last];assert e["changed_decisions"]+sh["changed_decisions"]==0 and e["localized"]==205 and sh["localized"]==113 and e["controls"]==53
        term=e["terms"]+sh["terms"];canon=baseline["terms"]+sums["canonical"]["shared"]["budgets"][last]["terms"]
        sec=sum(sums[key][k]["all_11_kernel_seconds"] for k in KINDS)
        exact[key]={"audit_terms":term,"audit_terms_vs_canonical":term/canon,"all_11_kernel_seconds":sec,"canonical_over_kernel_speedup":can/sec,"original_rival_over_new_with_contrast_speedup":refsec/(sec+normsec)}
    result={"status":"PASS_EXACT_BILLIARD_ORDERING__DEVELOPMENT_ONLY","groups":11,"audit_groups":5,"rows":3168,"geometry":sr["geometry"],"budgets":budgets,"exact":exact,"results":sums,"original_rival_seconds":refsec,"contrast_seconds":normsec,"acquisition_seconds":sum(d["acquisition_seconds"] for d in rows),"schedule_receipt_sha256":c.sha((O/"schedule_receipt.json").read_bytes()),"no_training":True,"production_adoption":False,"external_validation":False,"scope":"Exact schedule cannot improve scientific control counts; approximate prefix results are descriptive only."}
    save("result.json",result)
    cp("99-summary",{k:v for k,v in result.items() if k not in ("geometry","results")}|{"prefix_audit":{m:{k:{B:{z:v[z] for z in ("parents","both_winners","winner_sides","score_gap_max","localized","weak","controls","lost_joint","lost_vs_minimal","changed_decisions","repeat_stress_losses","terms")} for B,v in d["budgets"].items()} for k,d in v.items()} for m,v in sums.items()},"coverage":{m:{B:d["covered_sectors"] for B,d in v.items()} for m,v in sr["geometry"]["coverage"].items()}})
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="group":group(int(ARGS[1]))
    elif ARGS[0]=="summary":summary()
    else:raise ValueError(ARGS)
