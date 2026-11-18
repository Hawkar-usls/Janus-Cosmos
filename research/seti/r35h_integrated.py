#!/usr/bin/env python3
"""R35H exact JUXTAPOSE maximum, ODONTO baseline, protected raw rescue and verified-only HANNAH advisory."""
import hashlib,itertools,json,math,pathlib,resource,sys,time
import numpy as np
from numba import njit
import r35b_grand_benchmark as b
import r31_hannah_paired as prev
import r30c_hannah_memory as mem
P=pathlib.Path("research/seti/r35h_integrated_protocol.json");S=json.loads(P.read_text());H=S["r35h"];OUT=pathlib.Path("runtime/r35h");OUT.mkdir(parents=True,exist_ok=True)
b.P=P;b.S=S;b.OUT=OUT;prev.S=S;prev.OUT=OUT;mem.OUT=OUT
b.sh.P=P;b.sh.S=S;b.sh.C=S["shield"];b.sh.OUT=OUT;b.h.S=S;b.rb.S=S;b.ref.S=S;b.od.S=S;b.ref.OUT=OUT;b.sh.eng.OUT=OUT
b.sh.hold.S={**b.sh.hold.S,"seed":S["seed"],"exclude":sorted(set(b.sh.hold.S["exclude"]+S["exclude_extra"]))};b.sh.hold.P=P;b.sh.hold.OUT=OUT
D=np.asarray([np.rint(np.arange(16)*s/15).astype(int) for s in range(-32,33)],np.int64)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,v):mem.save(n,v)
def emit(k,v):mem.emit(k,v)
@njit(fastmath=False,cache=False)
def ordered(V,D,order):
    B,n,w=V.shape;best=np.full(B,-np.inf);paths=np.zeros((B,n),np.int64);full=0;pruned=0;row_ops=0
    for i in range(B):
        suffix=np.zeros(n+1);mxabs=0.
        for t in range(n-1,-1,-1):
            mx=-np.inf
            for k in range(w):
                v=V[i,t,k]
                if v>mx:mx=v
                if abs(v)>mxabs:mxabs=abs(v)
            suffix[t]=suffix[t+1]+mx
        # Strict slack generously exceeds 16-term float64 roundoff; ties never pruned.
        eps=1e-10*(1.+n*mxabs)
        bk=65;bj=w
        for k in order:
            d=D[k];lo=max(0,-d.min());hi=min(w,w-d.max())
            for j in range(lo,hi):
                full+=1;s=0.;stop=False
                for t in range(n):
                    s+=V[i,t,j+d[t]];row_ops+=1
                    if t<n-1 and s+suffix[t+1]+eps<best[i]*4.:
                        stop=True;pruned+=1;break
                if not stop:
                    z=s/4.
                    if z>best[i] or (z==best[i] and (k<bk or (k==bk and j<bj))):
                        best[i]=z;bk=k;bj=j
        for t in range(n):paths[i,t]=bj+D[bk,t]
    return best,paths,full,pruned,row_ops
def order_for(rank):
    r=np.asarray([s+32 for k in rank for s in b.j.GROUPS[k]],np.int64)
    assert sorted(r.tolist())==list(range(65));return r
def verify(V,rank):
    t=time.perf_counter();s,p=b.taylor.straight(V);old=time.perf_counter()-t
    t=time.perf_counter();q,w,n,k,ops=ordered(np.ascontiguousarray(V),D,order_for(rank));new=time.perf_counter()-t
    assert np.array_equal(s,q) and np.array_equal(p,w),"Exact maximum or tie witness mismatch"
    return {"samples":len(V),"full_paths":int(n),"pruned_paths":int(k),"row_operations":int(ops),"exhaustive_row_operations":int(16*n),"exhaustive_seconds":old,"ordered_seconds":new,"score_path_exact":True}
def selftest():
    rng=np.random.default_rng(35081);x=rng.normal(size=(8,16,128));x[0]=0.;x[1]=1.;x[2]=-1.;x[3]*=1e8
    t=time.perf_counter();ordered(x[:2],D,order_for(list(range(8))));compile=time.perf_counter()-t
    for rank in (list(range(8)),list(reversed(range(8)))):verify(x,rank)
    return {"exact_maximum_ties_and_large_scale":"PASS","compile_seconds":compile}
def loadmemory():
    def read(k,hh):
        p=pathlib.Path(H["memory"][k]);assert sha(p)==H["memory"][hh];return json.loads(p.read_text())
    lib,_=mem.append(read("g4","g4_sha256"),[]);assert len(lib["entries"])==169
    known=[e for e in lib["entries"] if e["label"]=="VERIFIED_RFI"];assert len(known)==1
    ext=read("sidecar","sidecar_sha256");assert len(ext["entries"])==6
    for e in ext["entries"]:
        assert e["origin_label"]=="VERIFIED_RFI_EXTERNAL_SOURCE_REPLAY";p=e["raw_evidence"]["passport"];mem.validate(p);known.append(p)
    drift=read("drift_candidates","drift_sha256");unconfirmed=[]
    for e in drift["entries"]:
        for p in [e["passport"]]+[o["passport"] for o in e["off_matches"]]:
            mem.validate(p);assert p["label"]=="UNKNOWN";unconfirmed.append(p)
    assert len(unconfirmed)==4
    old={e["provenance"]["source_tile_sha256"] for e in lib["entries"]+known+unconfirmed}
    return known,unconfirmed,old
class Guard:
    def __init__(self,old):self.old=old;self.seen=set()
    def ingest(self,q,f,t,x):
        h=t["slice_sha256"];assert hashlib.sha256(x.tobytes()).hexdigest()==h and h not in self.old and h not in self.seen
        self.seen.add(h)
def off_support(raw,path):
    try:y=mem.norm(raw.astype(float))
    except ValueError:return -1e6
    t=np.arange(len(path));best=-1e6
    for start,end in itertools.product(range(-2,3),(-2,0,2)):
        p=path+start+np.rint(t*end/15).astype(int)
        if p.min()<0 or p.max()>=y.shape[1]:continue
        best=max(best,float(np.sort(y[t,p].reshape(4,4).mean(1))[1]))
    return best
def evidence(row,R,path,known,unknown):
    m=row[9];out={"shape":0.,"reference":None,"review_shape":0.,"review_reference":None,"off":-1e6,"compatible":0,"valid":False}
    try:f=mem.features(row[0],path,m["df_hz"],m["ts"])
    except ValueError:return out
    out["valid"]=True;out["off"]=off_support(R,path);freq=m["channel0_mhz"]*1e6+int(path[0])*m["df_hz"]
    for refs,key,idkey in ((known,"shape","reference"),(unknown,"review_shape","review_reference")):
        for e in refs:
            u=e["features"]
            if not (math.isclose(f["channel_width_hz"],u["channel_width_hz"],rel_tol=.02) and math.isclose(f["sample_interval_s"],u["sample_interval_s"],rel_tol=.02)):continue
            if key=="shape":out["compatible"]+=1
            if abs(freq-e["frequency_hz"])>8.4 or abs(f["drift_hz_s"]-u["drift_hz_s"])>.011:continue
            score=float(np.exp(-np.sqrt(np.mean((np.asarray(f["vector"])-np.asarray(u["vector"]))**2))))
            if score>out[key]:out[key]=score;out[idkey]=e["id"]
    return out
def prepare(q,guard):
    rows,R,ch,hs=prev.gather(q,guard);assert len(hs)==16
    def part(rows,R):
        d=b.panel(rows,R);d["R"]=R;return d
    return part(rows,R),part([v[0] for v in ch],np.stack([v[1] for v in ch]))
def extract(d,L,threshold,known,unknown):
    a=b.base_arms(d,L);z,p,sec=a["hybrid_r28"];r,rp,rc=a["raw_smooth_bank"];y=b.labels(d);truth=np.stack([v[7] for v in d["rows"]]);act=np.stack([v[8] for v in d["rows"]]);amp=np.array([v[4] for v in d["rows"]]);alarm=z>=threshold
    loc=lambda paths:((abs(paths-truth)<=2)&act).sum(1)>=.75*act.sum(1)
    t=time.perf_counter();ev=[evidence(row,R,path,known,unknown) if al else {"shape":0.,"reference":None,"review_shape":0.,"review_reference":None,"off":-1e6,"compatible":0,"valid":False} for row,R,path,al in zip(d["rows"],d["R"],p,alarm)]
    dt=time.perf_counter()-t
    return {"y":y,"weak":amp<=1.25,"very_weak":amp==.75,"alarm":alarm,"loc":loc(p),"raw_loc":loc(rp),"raw":r,"p":p,"rp":rp,"z":z,"ev":ev,"time":sec+rc+dt,"baseline_time":sec,"hannah_seconds":dt,"tile":[v[6] for v in d["rows"]]}
def apply(d,policy):
    a=d["alarm"];hold=np.zeros(len(a),bool);rescue=np.zeros(len(a),bool)
    if policy.get("shape") is not None:
        hold=a&np.array([v["shape"]>=policy["shape"] and v["off"]>=policy["off"] for v in d["ev"]])
    if policy.get("raw_threshold") is not None:rescue=(~a)&(d["raw"]>=policy["raw_threshold"])
    keep=(a&~hold)|rescue;loc=((a&~hold)&d["loc"])|(rescue&d["raw_loc"])
    return keep,loc,hold,rescue
def metric(d,policy):
    keep,loc,hold,rescue=apply(d,policy);y=d["y"];hit=keep&loc&(y==1);base=d["alarm"]&d["loc"]&(y==1)
    return {"localized":int(hit.sum()),"weak":int((hit&d["weak"]).sum()),"very_weak":int((hit&d["very_weak"]).sum()),"controls":int((y==0).sum()),"positives":int((y==1).sum()),"control_alerts":int((keep&(y==0)).sum()),"held_positives":int((hold&(y==1)).sum()),"lost":int((base&~hit).sum()),"lost_weak":int((base&~hit&d["weak"]).sum()),"added":int((hit&~base).sum()),"rescued":int((rescue&(y==1)&loc).sum()),"verified_matches":sum(v["reference"] is not None for v in d["ev"]),"compatible_queries":sum(v["compatible"]>0 for v in d["ev"]),"review_matches":sum(v["review_reference"] is not None for v in d["ev"])}
def sum_metrics(ds,p):
    ms=[metric(d,p) for d in ds];return {k:sum(m[k] for m in ms) for k in ms[0]}
def fit_policy(cal,shared,rawcuts):
    base=sum_metrics(cal,{});ps=[{}]+[{"raw_threshold":v} for v in rawcuts]+[{"shape":s,"off":o} for s,o in itertools.product(H["policy"]["shape"],H["policy"]["off"])]+[{"raw_threshold":v,"shape":s,"off":o} for v,s,o in itertools.product(rawcuts,H["policy"]["shape"],H["policy"]["off"])]
    trials=[]
    for p in ps:
        a=sum_metrics(cal,p);l=sum_metrics(shared,p)
        admissible=a["held_positives"]==l["held_positives"]==a["lost"]==l["lost"]==0 and a["control_alerts"]<=base["control_alerts"]
        trials.append({"policy":p,"cal":a,"shared":l,"admissible":admissible})
    good=[x for x in trials if x["admissible"]];assert good
    winner=max(good,key=lambda v:(v["cal"]["weak"],v["cal"]["localized"],-v["cal"]["control_alerts"],v["shared"]["localized"],-len(v["policy"])))
    return winner["policy"],trials
def main():
    checks=selftest();b.taylor.init_tree();checks.update(b.h.selftest());known,unknown,old=loadmemory();guard=Guard(old)
    mf=json.loads((OUT/"manifest.json").read_text());assert mf["protocol_sha256"]==sha(P);ss=mf["selected"];assert len(ss)==8 and len({v["day"] for v in ss})==8
    assert not set(v["day"] for v in ss)&set(S["shield"]["exclude_days"]) and not set(v["target"] for v in ss)&set(S["exclude_extra"])
    assert [v["split"] for v in ss]==["component_fit"]*2+["null_calibration"]*2+["evaluation"]*4
    fit=[prepare(q,guard)[0] for q in ss[:2]];L=b.h.fit([(v["rows"],v["pack"]) for v in fit])
    pairs=[prepare(q,guard) for q in ss[2:4]]
    def control(k):
        return np.concatenate([b.base_arms(d,L)[k][0][b.labels(d)==0] for d,l in pairs])
    th=b.sh.cut(control("hybrid_r28"),.05);rawcuts=[b.sh.cut(control("raw_smooth_bank"),r) for r in H["policy"]["raw_rates"]]
    linear_th=b.sh.cut(control("registered_linear"),.05);rank,_=b.jrank(fit,linear_th)
    def ex(d):return extract(d,L,th,known,unknown)
    cal=[ex(d) for d,l in pairs];cl=[ex(l) for d,l in pairs];policy,trials=fit_policy(cal,cl,rawcuts)
    arms={"baseline":{},"protected_raw_1pct":{"raw_threshold":rawcuts[0]},"protected_raw_5pct":{"raw_threshold":rawcuts[-1]},"verified_hannah":{"shape":.95,"off":6},"joint_calibrated":policy}
    save("lock.json",{"component_fit":L,"baseline_threshold":th,"raw_thresholds":rawcuts,"rank":rank,"policy":policy,"trials":trials,"arms":arms,"manifest_sha256":sha(OUT/"manifest.json")});lock=sha(OUT/"lock.json");emit("R35H_LOCK",{"sha256":lock,"policy":policy,"cal_trials":len(trials),"admissible":sum(t["admissible"] for t in trials),"rank":rank})
    ev=[];shared=[];cost=[];records={}
    for q in ss[4:]:
        d,l=prepare(q,guard);cost.append(verify(d["Y"],rank));cost.append(verify(l["Y"],rank));a=ex(d);z=ex(l);ev.append(a);shared.append(z)
        records[q["target"]]={"eval":{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in a.items()},"shared":{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in z.items()}}
        save("partial.json",{"eval":{k:sum_metrics(ev,p) for k,p in arms.items()},"shared":{k:sum_metrics(shared,p) for k,p in arms.items()},"exact_search":cost});emit("R35H_TARGET",{"target":q["target"],"arms":{k:metric(a,p) for k,p in arms.items()},"shared":{k:metric(z,p) for k,p in arms.items()},"exact_search":cost[-2:]})
    assert sha(OUT/"lock.json")==lock and len(guard.seen)==128
    e={k:sum_metrics(ev,p) for k,p in arms.items()};l={k:sum_metrics(shared,p) for k,p in arms.items()};m=e["joint_calibrated"];a=e["baseline"]
    assert all(v["positives"]==v["controls"]==512 for v in e.values()) and all(v["positives"]==128 for v in l.values())
    speed={k:sum(v[k] for v in cost) for k in ("samples","full_paths","pruned_paths","row_operations","exhaustive_row_operations","exhaustive_seconds","ordered_seconds")};speed["score_path_exact"]=True
    gate={"overall_gain":m["localized"]>a["localized"],"weak_gain":m["weak"]>a["weak"],"fewer_controls":m["control_alerts"]<a["control_alerts"],"no_paired_loss":m["lost"]==m["held_positives"]==l["joint_calibrated"]["lost"]==l["joint_calibrated"]["held_positives"]==0,"ordered_max_faster_with_compile":speed["ordered_seconds"]+checks["compile_seconds"]<speed["exhaustive_seconds"]}
    res={"status":"INTEGRATED_PAIRED_DEVELOPMENT","gate":"PASS_DEVELOPMENT" if all(gate.values()) else "HOLD","gate_checks":gate,"eval":e,"shared":l,"policy":policy,"search":speed,"checks":checks,"known_references":len(known),"review_only_references":len(unknown),"distinct_measured_tiles":len(guard.seen),"lock_sha256":lock,"peak_rss_kib_process":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"scope":S["scope"],"production_adoption":False}
    save("records.json",records);save("result.json",res);emit("R35H_SUMMARY",res)
if __name__=="__main__":
    if sys.argv[1]=="freeze":b.sh.freeze()
    elif sys.argv[1]=="selftest":emit("R35H_SELFTEST",selftest())
    else:main()
