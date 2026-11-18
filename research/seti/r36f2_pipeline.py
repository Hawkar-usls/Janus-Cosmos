#!/usr/bin/env python3
"""R36F2 new source-group common pipeline. No destructive filtering, no production changes."""
import concurrent.futures,hashlib,io,json,os,pathlib,platform,resource,sys,time,zipfile
import numpy as np
import requests
import r35h_integrated as h
import r36_galaxy_swarm as g
import r36f1_hybrid as f
P=pathlib.Path("research/seti/r36f2_pipeline_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r36f2");OUT.mkdir(parents=True,exist_ok=True)
IN=OUT/"inputs";IN.mkdir(exist_ok=True)
b=h.b;ref=b.ref
def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,v):(OUT/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def emit(k,v):print(k,json.dumps(v,allow_nan=False),flush=True)
def setup(ex=None):
    h.P=P;h.S=S;h.H=S["r35h"];h.OUT=OUT;h.prev.S=S;h.prev.OUT=OUT;h.mem.OUT=OUT
    b.P=P;b.S=S;b.OUT=OUT;b.sh.P=P;b.sh.S=S;b.sh.C=S["shield"];b.sh.OUT=OUT
    b.h.S=S;b.rb.S=S;ref.S=S;ref.OUT=OUT;b.od.S=S;b.sh.eng.OUT=OUT
    targets=set(S["exclude_extra"]);days=set(S["shield"]["exclude_days"])
    if ex:
        targets.update(ex["targets"]);days.update(ex["days"])
    b.sh.hold.S={**b.sh.hold.S,"seed":S["seed"],"exclude":sorted(set(b.sh.hold.S["exclude"])|targets)}
    b.sh.hold.P=P;b.sh.hold.OUT=OUT
    # Keep original protocol bytes unchanged; dynamic exclusions are separately hash-bound.
    b.sh.C={**S["shield"],"exclude_days":sorted(days)}
    return targets,days
def prior():
    targets=set();days=set();names=set();urls=set();receipts=[]
    for a in S["r36f2"]["prior_manifests"]:
        u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip"
        r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90);r.raise_for_status()
        raw=r.content;assert len(raw)<16*1024**2 and sha(raw)==a["digest"]
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            n=[v for v in z.infolist() if pathlib.PurePosixPath(v.filename).name=="manifest.json"];assert len(n)==1 and n[0].file_size<8*1024**2
            x=z.read(n[0])
        m=json.loads(x);ss=m["selected"];assert len(ss)==(6 if a["stage"]=="R28B" else 8)
        for q in ss:
            targets.add(q["target"]);days.add(int(q["day"]))
            for v in q["files"]:names.add(v["name"]);urls.add(v["url"])
        receipts.append({**a,"manifest_sha256":sha(x),"groups":[{"target":q["target"],"day":q["day"]} for q in ss]})
    return {"targets":sorted(targets),"days":sorted(days),"names":sorted(names),"urls":sorted(urls),"sources":receipts}
def audit(m,ex):
    targets,days=setup(ex)
    ss=m["selected"];assert len(ss)==8 and len({q["day"] for q in ss})==8
    assert len({q["target"] for q in ss})==8 and not {q["target"] for q in ss}&targets and not {q["day"] for q in ss}&days
    assert [q["split"] for q in ss]==["component_fit"]*2+["null_calibration"]*2+["evaluation"]*4
    ns=[v["name"] for q in ss for v in q["files"]];us=[v["url"] for q in ss for v in q["files"]]
    assert len(set(ns))==len(set(us))==32 and not set(ns)&set(ex["names"]) and not set(us)&set(ex["urls"])
    return {"status":"PASS_KNOWN_GROUP_DISJOINTNESS","targets_excluded":len(targets),"days_excluded":len(days),"prior_scan_names":len(ex["names"]),"selected":[{k:q[k] for k in ("target","day","split")} for q in ss],"scope":"Known manifests and memory only, not external custody"}
def freeze():
    ex=prior();save("prior_exclusion.json",ex);setup(ex)
    start=time.perf_counter();b.sh.freeze()
    m=json.loads((OUT/"manifest.json").read_text());assert m["protocol_sha256"]==sha(P.read_bytes())
    a=audit(m,ex);a["prior_exclusion_sha256"]=sha((OUT/"prior_exclusion.json").read_bytes());a["manifest_sha256"]=sha((OUT/"manifest.json").read_bytes());a["freeze_wall_seconds"]=time.perf_counter()-start
    save("source_audit.json",a);emit("R36F2_FREEZE",a)
def clock(fn):
    t=time.perf_counter();c=time.process_time();v=fn()
    return v,{"wall":time.perf_counter()-t,"cpu":time.process_time()-c}
def add(*xs):return {k:sum(x[k] for x in xs) for k in ("wall","cpu")}
def prep(rows,R):
    def fn():
        X=ref.norm(np.stack([x[0] for x in rows]));Y=ref.register(X,ref.norm(R))[0]
        return np.ascontiguousarray(X.astype(np.float32),dtype=np.float64),np.ascontiguousarray(Y.astype(np.float32),dtype=np.float64)
    (X,Y),t=clock(fn)
    return {"rows":rows,"R":R,"X":X,"Y":Y,"preprocess":t}
def gather(q,guard):
    (rows,R,ch,hs),t=clock(lambda:h.prev.gather(q,guard));assert len(hs)==16
    return prep(rows,R),prep([v[0] for v in ch],np.stack([v[1] for v in ch])),t
def labels(d):return np.array([r[1] for r in d["rows"]])
def fit(parts):
    y=np.concatenate([labels(d) for d in parts]);a=[(f.lin(d["Y"]),f.od(d["Y"])) for d in parts];L={}
    for k,i in (("linear",0),("odonto",1)):
        z=np.concatenate([v[i][0] for v in a])[y==0];th=f.cut(z)
        L[k+"_threshold"]=th;L[k+"_scale"]=max(.1,th-float(np.quantile(z,.8)))
    return L
def compute(d,L,pool):
    Y=d["Y"];X=d["X"]
    def base():return f.fuse(f.lin(Y),f.od(Y),L)
    def twin():
        aa=pool.submit(g.scan,Y,g.D,g.R,g.L,False);bb=pool.submit(f.od,Y)
        return f.fuse(f.get(aa.result(),0),bb.result(),L)
    # Raw route is independent and never replaces an already detected baseline path.
    def full():
        z,p=twin();r,rp=f.get(g.scan(X,g.D,g.R,g.T,False),1)
        return z,p,r,rp
    (z,p),t0=clock(base);(u,v),t1=clock(twin);(s,sp,r,rp),t2=clock(full)
    f.same((z,p),(u,v));f.same((z,p),(s,sp))
    return {"z":z,"p":p,"raw":r,"rp":rp,"timing":{"reference":t0,"twin":t1,"galaxy_search":t2}}
def evidence(d,c,th,known,unknown):
    y=labels(d);truth=np.stack([r[7] for r in d["rows"]]);act=np.stack([r[8] for r in d["rows"]]);amp=np.array([r[4] for r in d["rows"]])
    alarm=c["z"]>=th;loc=lambda p:((abs(p-truth)<=2)&act).sum(1)>=.75*act.sum(1)
    # Raw rescue is not silently vetoed by the reference recognizer. Preserve all source candidates.
    def ev():
        return [h.evidence(row,R,path,known,unknown) if a else {"shape":0.,"reference":None,"review_shape":0.,"review_reference":None,"off":-1e6,"compatible":0,"valid":False} for row,R,path,a in zip(d["rows"],d["R"],c["p"],alarm)]
    vv,t=clock(ev)
    return {**c,"y":y,"weak":amp<=1.25,"very_weak":amp==.75,"alarm":alarm,"loc":loc(c["p"]),"raw_loc":loc(c["rp"]),"ev":vv,"hannah_time":t}
def encode(d,e,policy,label,target):
    rows=d["rows"];keep,loc,hold,rescue=h.apply(e,policy)
    z=e["z"];yp=e["y"];p=np.where(rescue[:,None],e["rp"],e["p"]);base=e["alarm"]&e["loc"]&(yp==1);hit=keep&loc&(yp==1)
    arr={"V":d["Y"].astype(np.float32),"y":yp,"truth":np.stack([r[7] for r in rows]),"active":np.stack([r[8] for r in rows]),"amp":np.array([r[4] for r in rows]),"model":np.array([r[2] for r in rows],dtype="U32"),"tile_sha":np.array([r[6] for r in rows],dtype="U128"),"hybrid_r28_localized":base,"galaxy_localized":hit,"hybrid_r28_alert":e["alarm"],"galaxy_alert":keep,"hybrid_r28_paths":e["p"],"galaxy_paths":p,"hybrid_r28_score":z,"raw_score":e["raw"],"raw_path":e["rp"],"hold_advisory":hold,"raw_rescue":rescue}
    path=IN/(label+"_"+target+".npz");np.savez_compressed(path,**arr)
    return {"name":path.name,"sha256":sha(path.read_bytes()),"array_sha256":sha(arr["V"].tobytes()),"examples":len(rows),"metrics":{"reference":h.metric(e,{}),"galaxy":h.metric(e,policy)},"preprocess":d["preprocess"],"timing":{**e["timing"],"hannah":e["hannah_time"]}}
def main():
    ex=json.loads((OUT/"prior_exclusion.json").read_text());m=json.loads((OUT/"manifest.json").read_text());a=json.loads((OUT/"source_audit.json").read_text())
    assert a["prior_exclusion_sha256"]==sha((OUT/"prior_exclusion.json").read_bytes()) and a["manifest_sha256"]==sha((OUT/"manifest.json").read_bytes())
    audit(m,ex);ss=m["selected"];started=time.perf_counter();known,unknown,old=h.loadmemory();guard=h.Guard(old)
    def init():
        x=np.zeros((2,16,128));f.same(f.lin(x),f.get(g.scan(x,g.D,g.R,g.L,False),0));g.scan(x,g.D,g.R,g.T,False)
    _,compile=clock(init)
    times=[];records=[];all_ev=[];all_shared=[];detail={};calstart=time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        fits=[]
        for q in ss[:2]:
            d,l,t=gather(q,guard);fits.append(d);times.append({"split":"fit","target":q["target"],"acquire":t,"preprocess":d["preprocess"]})
        L,fit_time=clock(lambda:fit(fits))
        pairs=[]
        for q in ss[2:4]:
            d,l,t=gather(q,guard);times.append({"split":"cal","target":q["target"],"acquire":t})
            pairs.append((d,l,compute(d,L,pool),compute(l,L,pool)))
        controls=lambda k:np.concatenate([c[k][labels(d)==0] for d,l,c,s in pairs])
        th=f.cut(controls("z"));rawcuts=[b.sh.cut(controls("raw"),r) for r in h.H["policy"]["raw_rates"]]
        ee=[evidence(d,c,th,known,unknown) for d,l,c,s in pairs];ll=[evidence(l,s,th,known,unknown) for d,l,c,s in pairs]
        (policy,trials),policy_time=clock(lambda:h.fit_policy(ee,ll,rawcuts))
        save("lock.json",{"component_fit":L,"baseline_threshold":th,"raw_thresholds":rawcuts,"policy":policy,"trials":trials,"source_audit_sha256":sha((OUT/"source_audit.json").read_bytes())});lock=sha((OUT/"lock.json").read_bytes())
        for q,(d,l,c,s),e in zip(ss[2:4],pairs,ee):records.append(encode(d,e,policy,"cal",q["target"]))
        assert sum(r["metrics"]["galaxy"]["control_alerts"] for r in records)<=12
        emit("R36F2_LOCK",{"sha256":lock,"policy":policy,"calibration_trials":len(trials),"calibration_seconds":time.perf_counter()-calstart})
        for q in ss[4:]:
            assert sha((OUT/"lock.json").read_bytes())==lock
            d,l,t=gather(q,guard);times.append({"split":"eval","target":q["target"],"acquire":t})
            e=evidence(d,compute(d,L,pool),th,known,unknown);s=evidence(l,compute(l,L,pool),th,known,unknown)
            all_ev.append(e);all_shared.append(s)
            records.extend([encode(d,e,policy,"eval",q["target"]),encode(l,s,policy,"shared",q["target"])])
            detail[q["target"]]={"eval":e["ev"],"shared":s["ev"]}
            save("partial.json",records);emit("R36F2_TARGET",{"target":q["target"],"eval":records[-2],"shared":records[-1]})
    assert len(guard.seen)==128 and len(records)==10 and sha((OUT/"lock.json").read_bytes())==lock
    def metrics(parts):return {"reference":h.sum_metrics(parts,{}),"galaxy":h.sum_metrics(parts,policy)}
    e=metrics(all_ev);s=metrics(all_shared);assert e["reference"]["positives"]==e["reference"]["controls"]==512 and s["reference"]["positives"]==128
    def cost(prefix):
        rr=[r for r in records if r["name"].startswith(prefix)]
        pre=add(*(r["preprocess"] for r in rr))
        out={k:add(*(r["timing"][k] for r in rr)) for k in ("reference","twin","galaxy_search","hannah")}
        out["shared_preprocess"]=pre;out["reference_pipeline"]=add(pre,out["reference"]);out["galaxy_pipeline"]=add(pre,out["galaxy_search"],out["hannah"]);return out
    cost_eval=cost("eval_");cost_cal=cost("cal_");cost_shared=cost("shared_")
    ga=e["galaxy"];base=e["reference"];gs=s["galaxy"]
    gate={"overall":ga["localized"]>base["localized"],"weak":ga["weak"]>base["weak"],"fewer_controls":ga["control_alerts"]<base["control_alerts"],"zero_losses":ga["lost"]==ga["held_positives"]==gs["lost"]==gs["held_positives"]==0,"warm_compute_wall":cost_eval["galaxy_pipeline"]["wall"]<cost_eval["reference_pipeline"]["wall"],"warm_compute_cpu":cost_eval["galaxy_pipeline"]["cpu"]<cost_eval["reference_pipeline"]["cpu"]}
    result={"status":"PASS_COMMON_REPLAY","gate":"PASS_DEVELOPMENT" if all(gate.values()) else "HOLD","gate_checks":gate,"records":records,"lock_sha256":lock,"freshness_status":"KNOWN_GROUP_DISJOINT_NOT_EXTERNAL_BLIND","r32_overlap":[],"source_audit":a,"metrics":e,"shared":s,"policy":policy,"timing":{"evaluation":cost_eval,"calibration":cost_cal,"shared":cost_shared,"startup_compile":compile,"component_fit":fit_time,"policy_fit":policy_time,"acquisition":times,"total_experiment_wall":time.perf_counter()-started},"known_references":len(known),"review_only_references":len(unknown),"distinct_raw_tiles":len(guard.seen),"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"cpu":platform.processor(),"production_adoption":False,"origin_validation":"HOLD_REAL_PROTECTED_AND_EXTERNAL_CUSTODY"}
    save("advisory.json",detail);save("receipt.json",result)
    emit("R36F2_SUMMARY",{k:v for k,v in result.items() if k not in ("records","source_audit")})
if __name__=="__main__":
    {"freeze":freeze,"scan":main}[sys.argv[1]]()
