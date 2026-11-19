#!/usr/bin/env python3
"""R41B: independently executable fixed minimal pipeline, paired real acquisition."""
import sys,time
BOOT=time.perf_counter();BOOT_CPU=time.process_time();ARGS=sys.argv[1:]
import json,pathlib,resource,subprocess,collections,hashlib,statistics
import numpy as np
import r40n10_artmoney_turbo as a
n=a.n;c=a.c;b=a.b;h=a.h;w=a.w;old9=a.old9;halo=a.halo;turbo=a.turbo;l=a.l;eng=a.eng
P=pathlib.Path("research/seti/r41b_minimal_hybrid_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r41b");O.mkdir(parents=True,exist_ok=True);C=O/"cache";C.mkdir(exist_ok=True);PROG=O/"progress";PROG.mkdir(exist_ok=True)
a.O=O;a.CACHE=C;a.PROG=PROG;n.O=O;n.k.O=O;n.l.O=O;eng.OUT=O
def write(path,d):a.write(path,d)
def save(k,d):write(O/k,d)
def load(k):return json.loads((O/k).read_text())
def cp(k,d):
    write(PROG/(k+".json"),d);print("R41B_PROGRESS",json.dumps(d,allow_nan=False),flush=True)
def sha(raw):return c.sha(raw)
def r41a_archive(spec):
    # This archive intentionally contains group-NN.json both at its root and in progress/.
    # Authenticate the whole ZIP first, then resolve required files relative to the unique
    # pruning lock directory; never collapse paths to basenames or extract untrusted files.
    import io,os,zipfile,requests
    u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{spec['id']}/zip"
    r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=120)
    r.raise_for_status();raw=r.content
    assert len(raw)<32*1024**2 and sha(raw)==spec["sha256"]
    z=zipfile.ZipFile(io.BytesIO(raw));entries=[i for i in z.infolist() if not i.is_dir()]
    assert sum(i.file_size for i in entries)<512*1024**2
    names=[i.filename for i in entries];assert len(names)==len(set(names))
    assert all(not pathlib.PurePosixPath(s).is_absolute() and ".." not in pathlib.PurePosixPath(s).parts for s in names)
    locks=[s for s in names if pathlib.PurePosixPath(s).name=="pruning_lock.json"]
    assert len(locks)==1 and sha(z.read(locks[0]))==spec["lock_sha256"]
    root=pathlib.PurePosixPath(locks[0]).parent
    required=["pruning_lock.json"]+["group-"+str(i+1).zfill(2)+".json" for i in range(11)]
    n={s:str(root/s) for s in required};assert all(v in names for v in n.values())
    # Verify group payload schema rather than silently using a short checkpoint.
    for s in required[1:]:
        d=json.loads(z.read(n[s]));assert "modes" in d and "operational" in d["modes"] and d["stage"]=="GROUP"
    print("R41B_ARCHIVE_NAMESPACE",json.dumps({"root":str(root),"entries":len(names),"resolved":len(n)}),flush=True)
    return z,n
def prepare():
    za,na=r41a_archive(S["source_r41a"]);zr,nr=c.artifact(S["source_r41"],32);z,nn=c.artifact(S["source_n10"],128)
    raw=za.read(na["pruning_lock.json"]);assert sha(raw)==S["source_r41a"]["lock_sha256"]
    lock=json.loads(raw);assert lock["engineering_comparator"]=="minimal__joint" and lock["winner"] is None and lock["rule"]=="minimal" and lock["rescue"]=="joint"
    (C/"pruning_lock.json").write_bytes(raw)
    assert sha(zr.read(nr["winner_lock.json"]))==S["source_r41"]["winner_lock_sha256"]
    m=json.loads(zr.read(nr["models.json"]));mini=m["rival+runs__tree"];full=m["temporal+rival+runs+folds+halo+artmoney+repeat__tree"]
    assert mini==lock["models"]["minimal"] and mini["blocks"]==["rival","runs"] and mini["kind"]==full["kind"]=="tree"
    write(C/"full_model.json",full);assert sha(z.read(nn["lock.json"]))==S["source_n10"]["lock_sha256"]
    rec=json.loads(z.read(nn["freeze_receipt.json"]));raw=z.read(nn["plan.json"]);assert rec["plan_sha256"]==sha(raw) and rec["lock_sha256"]==S["source_n10"]["lock_sha256"]
    gs=json.loads(raw)["development"];assert len(gs)==11 and sum(g["fit"] for g in gs)==6 and len({g["group"]["day"] for g in gs})==11
    hs=[v for g in gs for v in g["raw_tile_hashes"]];assert len(hs)==len(set(hs))==176
    for i,g in enumerate(gs):
        num=str(i+1).zfill(2)
        old=json.loads(za.read(na["group-"+num+".json"]));assert old["target"]==g["group"]["target"] and old["fit_viewed"]==g["fit"];write(C/("r41a-"+num+".json"),old)
        for kind in ("eval","shared"):
            fn="dev_"+num+"_"+kind+".npz";(C/fn).write_bytes(z.read(nn[fn]))
    h.O=C;L=h.readlocks();write(C/"parent_locks.json",L)
    zz,nn0=n.archive("r40l");raw=zz.read(nn0["lock.json"]);assert sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];(C/"joint_lock.json").write_bytes(raw)
    save("plan.json",{"protocol_sha256":sha(P.read_bytes()),"groups":gs,"n10_plan_sha256":rec["plan_sha256"],"pruning_lock_sha256":S["source_r41a"]["lock_sha256"],"full_model_sha256":sha((C/"full_model.json").read_bytes()),"parent_locks_sha256":sha((C/"parent_locks.json").read_bytes())})
    cp("00-prepare",{"stage":"PREPARE","groups":11,"paired_processes":22,"old_rows":3168,"raw_tiles":176,"no_fitting":True})
def plan():
    p=load("plan.json");assert p["protocol_sha256"]==sha(P.read_bytes());return p
def locks(mode,p):
    raw=(C/"pruning_lock.json").read_bytes();assert sha(raw)==p["pruning_lock_sha256"]
    if mode=="minimal":m=json.loads(raw)["models"]["minimal"]
    else:
        raw=(C/"full_model.json").read_bytes();assert sha(raw)==p["full_model_sha256"];m=json.loads(raw)
    raw=(C/"parent_locks.json").read_bytes();assert sha(raw)==p["parent_locks_sha256"];L=json.loads(raw)
    raw=(C/"joint_lock.json").read_bytes();assert sha(raw)==n.S["sources"]["r40l"]["lock_sha256"];return m,L,json.loads(raw)
def rival_only(X,D,p):
    # Only measured, legal bins; no wrapping or artificial edge values as evidence.
    assert X.shape==D.shape and X.shape[1:]==(16,128) and p.shape==(len(X),16)
    rows=np.arange(16)[None,:];ii=np.arange(len(X))[:,None]
    valid=(p>=5)&(p<123);count=valid.sum(1);supported=count>=12
    pp=np.clip(p,5,122);U=n.normalized_contrast(X);V=n.normalized_contrast(D)
    a=U[ii,rows,pp];d=V[ii,rows,pp]
    a=np.where(valid,a,0.);d=np.where(valid,d,0.)
    def mean(z,mask):
        return (z*mask).sum(1)/np.maximum(1,mask.sum(1))
    # Recalculate raw ON and OFF separately on candidate and alternatives.
    best=np.full((len(X),2),-np.inf);rpaths=np.zeros((len(X),2,16),int);ix=np.arange(len(X))
    for delta in b.D:
        low=5-int(delta.min());high=123-int(delta.max())
        if high<=low:continue
        starts=np.arange(low,high);paths=starts[:,None]+delta
        distinct=(abs(paths[None,:,:]-p[:,None,:])>2).sum(2)>=5
        for k,Z in enumerate((U,V)):
            sc=sum(Z[:,j,low+dd:high+dd] for j,dd in enumerate(delta))/16
            sc=np.where(distinct,sc,-np.inf);arg=sc.argmax(1);vv=sc[ix,arg];m=vv>best[:,k]
            best[m,k]=vv[m];rpaths[m,k]=paths[arg[m]]
    assert np.isfinite(best).all()
    for k,Z in enumerate((U,V)):
        path=rpaths[:,k];assert ((abs(path-p)>2).sum(1)>=5).all()
        assert np.allclose(Z[ii,rows,path].mean(1),best[:,k],atol=1e-12,rtol=1e-12)
    rv=np.column_stack([best,mean(a,valid)-best[:,0],mean(d,valid)-best[:,1],best[:,0]-best[:,1],np.mean(abs(rpaths[:,0]-rpaths[:,1]),axis=1)])
    assert np.isfinite(rv).all()
    return rv,supported,rpaths,best

def runs_only(X,D,p):
    N=len(X);ii=np.arange(N)[:,None];rr=np.arange(16)[None,:]
    valid=(p>=5)&(p<123);pp=np.clip(p,5,122);ok=(valid.sum(1)>=12)
    def mean(z,mask):return (z*mask).sum(1)/np.maximum(1,mask.sum(1))
    U=n.normalized_contrast(X);V=n.normalized_contrast(D);run=[]
    for Z in (U,V):
        sel=Z[ii,rr,pp];run+=old9.runstats(sel,valid)
        # Native ±2 residuals; strongest local bin relative to selected path.
        # Restrict all candidate neighbours to legal contrast support.
        nei=pp[:,:,None]+np.arange(-2,3);nv=(nei>=5)&(nei<123)
        z=Z[ii[:,:,None],rr[:,:,None],np.clip(nei,5,122)];z=np.where(nv,z,-np.inf)
        best=z.max(2);offset=z.argmax(2)-2;run += [mean(abs(offset),valid),mean(offset==0,valid),mean(best-sel,valid),mean(abs(np.diff(offset,axis=1)),valid[:,:-1]&valid[:,1:])]
    for tr,te in old9.FOLDS:
        vt=valid[:,tr];ve=valid[:,te];ok &= (vt.sum(1)>=4)&(ve.sum(1)>=4)
    R=np.column_stack(run)
    assert R.shape==(N,44) and np.isfinite(R).all()
    return R,ok

def core(rows,refs,L,old):
    # Identical common upstream path in both arms; use only required joint score/path.
    # labels/amp are carried solely for localization metrics and prior API shape.
    X=n.ref.norm(np.stack([v[0] for v in rows]));D=n.ref.norm(refs);Y,_,_,_=n.ref.register(X,D)
    q={"V":Y.astype(np.float32),"y":np.array([v[1] for v in rows]),"truth":np.stack([v[7] for v in rows]),"active":np.stack([v[8] for v in rows]),"amp":np.array([v[4] for v in rows])}
    s,p,sup,reg,F=l.evidence(rows,refs,q)
    pan=b.panel(rows,refs);arms=b.base_arms(pan,L["r40c"]["original"]["component_fit"])
    inp={**q,"V":pan["Y"].astype(np.float32),"tile_sha":np.array([v[6] for v in rows])}
    for path in a.t.S["paths"]:inp[path+"_score"],inp[path+"_paths"]=arms[path][:2]
    aa,pp,hh=l.parent(inp,L);aa=aa.astype(bool);hh=hh.astype(bool);res=(~aa)&(s>=old["parent_joint_min_threshold"])
    return {"parent_alert":aa,"parent_hit":hh,"parent_path":pp,"joint_rescue":res,"correct":l.correct(q,p),"p":p,"joint_score":s,"joint_support":sup,"joint_regime":reg,"idx":np.flatnonzero(aa),"y":q["y"],"truth":q["truth"],"active":q["active"],"amp":q["amp"],"weak":np.isin(q["amp"],a.t.S["weak_amplitudes"]),"source_tile_sha":inp["tile_sha"]}
def infer(mode,rows,refs,meta,pre,L,old,m,paths,kind):
    d=core(rows,refs,L,old);idx=d["idx"];pp=d["parent_path"][idx].astype(int)
    if len(idx):X=n.ref.norm(np.stack([rows[i][0] for i in idx]));D=n.ref.norm(refs[idx])
    else:X=np.empty((0,16,128));D=X.copy()
    if mode=="minimal":
        if len(idx):rv,supported,rp,bs=rival_only(X,D,pp);runs,ok=runs_only(X,D,pp)
        else:rv=np.empty((0,6));supported=np.empty(0,bool);rp=np.empty((0,2,16),int);bs=np.empty((0,2));runs=np.empty((0,44));ok=np.empty(0,bool)
        d.update(rival=rv,runs=runs,supported=supported,fold_supported=ok,rival_paths=rp,rival_scores=bs)
        use=supported&ok
    else:
        if len(idx):st,rv,supported,_,_,rp,bs=w.witness(X,D,pp);runs,folds,ok,fp,fs=old9.extra(X,D,pp)
        else:st=np.empty((0,36));rv=np.empty((0,6));supported=np.empty(0,bool);rp=np.empty((0,2,16),int);bs=np.empty((0,2));runs=np.empty((0,44));folds=np.empty((0,56));ok=np.empty(0,bool)
        d.update(temporal=st,rival=rv,runs=runs,folds=folds,supported=supported,fold_supported=ok,rival_paths=rp,rival_scores=bs)
        d.update(halo.extract(rows,refs,meta,pre,idx,pp,n.ref,c.sha,paths))
        d.update(turbo.repeat_witness(rows,meta,pre,idx,pp,n.ref,kind))
        use=supported&ok&d["halo_supported"]&d["repeat_supported"]
    x=np.column_stack([d[k] for k in m["blocks"]]);sc=a.score(x,m["model"],m["kind"]);keep=np.ones(len(d["y"]),bool);keep[idx]=(~use)|(sc>=m["cuts"]["preserve"])
    aa=d["parent_alert"];res=d["joint_rescue"];d.update(score=sc,keep=keep,use=use,alert=(aa&keep)|res,hit=(d["parent_hit"]&keep)|(res&d["correct"]))
    return d
def check(d,old,m):
    keys=["parent_alert","parent_hit","parent_path","joint_rescue","correct","p","idx","y","truth","active","amp","weak","source_tile_sha","rival","runs","supported","fold_supported","rival_paths","rival_scores"]+([k for k in ("temporal","folds","halo","artmoney","halo_supported","repeat","repeat_stress","repeat_supported","repeat_stress_supported") if k in d])
    for k in keys:assert np.array_equal(d[k],old[k]),"R41B parity "+k
    x=np.column_stack([old[k] for k in m["blocks"]]);s=a.score(x,m["model"],m["kind"]);assert np.array_equal(s,d["score"]),"score parity"
    use=old["supported"]&old["fold_supported"]
    if "halo" in m["blocks"]:use=use&old["halo_supported"]&old["repeat_supported"]
    keep=np.ones(len(old["y"]),bool);keep[old["idx"]]=(~use)|(s>=m["cuts"]["preserve"])
    assert np.array_equal(use,d["use"]) and np.array_equal(keep,d["keep"])
    assert np.array_equal(d["alert"],(old["parent_alert"]&keep)|old["joint_rescue"]) and np.array_equal(d["hit"],old["parent_hit"]&keep|(old["joint_rescue"]&old["correct"]))
    stress=old["repeat_stress_supported"]&old["supported"]&old["fold_supported"]&old["halo_supported"]&old["repeat_supported"]&old["parent_hit"][old["idx"]]
    xs=np.column_stack([old["repeat_stress"] if k=="repeat" else old[k] for k in m["blocks"]]);ss=a.score(xs,m["model"],m["kind"])
    r=d["joint_rescue"];aa=d["parent_alert"];hh=d["parent_hit"];null=d["y"]==0;hit=d["hit"];weak=d["weak"];jm=hh|(r&d["correct"])
    return {"rows":len(null),"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"controls":int((d["alert"]&null).sum()),"parent_controls":int((aa&d["keep"]&null).sum()),"rescue_controls":int((r&null).sum()),"lost_joint":int((jm&~hit).sum()),"stress_cases":int(stress.sum()),"stress_loss":int((stress&use&(ss<m["cuts"]["preserve"])).sum()),"matched_keys":len(keys)}
def residual(d,rows,refs,meta,pre,m):
    ids=np.flatnonzero(d["alert"]&(d["y"]==0));out=[]
    if len(ids)==0:return out
    pa=d["parent_alert"][ids];pp=np.where(pa[:,None],d["parent_path"][ids],d["p"][ids]).astype(int)
    X=n.ref.norm(np.stack([rows[i][0] for i in ids]));D=n.ref.norm(refs[ids]);rv,valid,rp,bs=rival_only(X,D,pp)
    r=turbo.repeat_witness(rows,meta,pre,ids,pp,n.ref,"eval")
    for j,i in enumerate(ids):
        _,wi,f0,df=meta[i];v=(pp[j]>=5)&(pp[j]<123);p0=np.clip(pp[j],5,122)
        on=n.normalized_contrast(X[j:j+1])[0,np.arange(16),p0];off=n.normalized_contrast(D[j:j+1])[0,np.arange(16),p0]
        k=int(np.searchsorted(d["idx"],i))
        margin=float(d["score"][k]-m["cuts"]["preserve"]) if pa[j] else float(d["joint_score"][i])
        out.append({"index":int(i),"lane":"parent" if pa[j] else "rescue","tile_sha256":str(d["source_tile_sha"][i]),"window_index":int(wi),"frequency_start_mhz":float(f0+pp[j,0]*df/1e6),"df_hz":float(df),"path":pp[j].tolist(),"legal_rows":int(v.sum()),"rival_supported":bool(valid[j]),"rival":rv[j].tolist(),"rival_paths":rp[j].tolist(),"on_contrast":on[v].tolist(),"off_contrast":off[v].tolist(),"other_scan_supported":bool(r["repeat_supported"][j]),"other_scan_features":r["repeat"][j].tolist() if r["repeat_supported"][j] else None,"score_or_margin":margin,"score_kind":"tree_margin" if pa[j] else "joint_score","fail_open":bool(pa[j] and not d["use"][k]),"hard_control":bool(rows[i][5])})
    return out
def worker(i,mode):
    assert mode in ("minimal","full");p=plan();g=p["groups"][i];name=g["group"]["target"];num=str(i+1).zfill(2)
    m,L,old=locks(mode,p);seed={"n6":a.v.C["seed"],"n7":w.C["seed"],"n9":old9.C["seed"]}[g["source"]]
    n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);b.taylor.init_tree()
    paths=None
    if mode=="full":
        paths=halo.candidates(b.D);halo.context_scores=turbo.Turbo(paths,0)
    setup_wall=time.perf_counter()-BOOT
    net={"payload_bytes":0,"requests":0,"files":0};orig=eng.base.RangeFile
    class Counted(orig):
        def close(self):
            if not self.closed and hasattr(self,"total") and hasattr(self,"size"):
                net["payload_bytes"]+=int(self.total);net["requests"]+=int(self.requests_count);net["files"]+=1
            super().close()
    eng.base.RangeFile=Counted
    start=time.perf_counter();cpu=time.process_time()
    try:rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r41b_"+mode+"_"+num)
    finally:eng.base.RangeFile=orig
    acq=time.perf_counter()-start;acqcpu=time.process_time()-cpu
    assert sorted(hs)==g["raw_tile_hashes"] and net["files"]==4 and net["payload_bytes"]>0
    pairs=[("eval",rr,R,meta),("shared",[x[0] for x in ch],np.stack([x[1] for x in ch]),cm)]
    # Two actual calls, first cold including lazy JIT, then warm on identical arrays.
    cold_start=time.perf_counter();cold_cpu=time.process_time()
    first={kind:infer(mode,rows,refs,mm,pre,L,old,m,paths,kind) for kind,rows,refs,mm in pairs}
    perf={"setup_wall":setup_wall,"acquisition_wall":acq,"acquisition_cpu":acqcpu,**net,"cold_compute_wall":time.perf_counter()-cold_start,"cold_compute_cpu":time.process_time()-cold_cpu,"cold_total_wall":time.perf_counter()-BOOT,"cold_total_cpu":time.process_time()-BOOT_CPU,"cold_peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    warm_start=time.perf_counter();warm_cpu=time.process_time()
    # Keep only small decision/path receipts from cold run so warm peak is not inflated by duplicate full feature arrays.
    first={k:{z:v[z] for z in ("alert","hit","parent_path","p","score")} for k,v in first.items()}
    second={kind:infer(mode,rows,refs,mm,pre,L,old,m,paths,kind) for kind,rows,refs,mm in pairs}
    perf.update(warm_compute_wall=time.perf_counter()-warm_start,warm_compute_cpu=time.process_time()-warm_cpu,warm_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    metrics={};ledger=[];refold=load("cache/r41a-"+num+".json")
    for kind,rows,refs,mm in pairs:
        d=second[kind];oldarr=a.read_npz(C/("dev_"+num+"_"+kind+".npz"));a.verify_rows(oldarr,rows)
        for z,v in first[kind].items():assert np.array_equal(d[z],v),"cold warm parity "+z
        q=check(d,oldarr,m);ref=refold["modes"]["operational"][kind]["metrics"][mode+"__joint"]
        for k,v in {"localized":"localized","weak":"weak","controls":"control_alerts","parent_controls":"parent_control_alerts","rescue_controls":"rescue_control_alerts","lost_joint":"lost_joint","stress_cases":"repeat_stress_cases","stress_loss":"repeat_stress_lost"}.items():assert q[k]==ref[v],"R41A metric parity "+k
        metrics[kind]=q
        if mode=="minimal" and kind=="eval" and not g["fit"]:
            ledger=residual(d,rows,refs,mm,pre,m)
            expected=refold["modes"]["operational"]["eval"]["residual_ledger"]
            assert [(v["index"],v["lane"],v["tile_sha256"]) for v in ledger]==[(v["index"],v["lane"],v["tile_sha256"]) for v in expected],"R41A residual identity"
    result={"group":i+1,"mode":mode,"target":name,"day":g["group"]["day"],"fit_viewed":g["fit"],"parity":True,"raw_tiles":sorted(hs),"metrics":metrics,"performance":perf,"residual":ledger}
    save(mode+"-"+num+".json",result)
    print("R41B_WORKER",json.dumps({k:v for k,v in result.items() if k!="residual"},allow_nan=False),flush=True)
def pair(i):
    modes=["minimal","full"] if i%2==0 else ["full","minimal"];num=str(i+1).zfill(2)
    for mode in modes:subprocess.run([sys.executable,__file__,"worker",str(i),mode],check=True)
    v={k:load(k+"-"+num+".json") for k in modes};assert v["minimal"]["raw_tiles"]==v["full"]["raw_tiles"]
    cp("group-"+num,{"stage":"PAIRED_GROUP","group":i+1,"total_groups":11,"target":v["minimal"]["target"],"fit_viewed":v["minimal"]["fit_viewed"],"order":modes,"parity":True,"minimal":v["minimal"]["performance"],"full":v["full"]["performance"],"metrics":{k:x["metrics"] for k,x in v.items()}})
def summary():
    allr={k:[load(k+"-"+str(i+1).zfill(2)+".json") for i in range(11)] for k in ("minimal","full")}
    assert all(x["parity"] for v in allr.values() for x in v)
    metrics={};performance={}
    for k,rr in allr.items():
        audit=[x for x in rr if not x["fit_viewed"]]
        metrics[k]={typ:{z:sum(x["metrics"][typ][z] for x in audit) for z in audit[0]["metrics"][typ]} for typ in ("eval","shared")}
        fields=rr[0]["performance"];performance[k]={z:(max(x["performance"][z] for x in rr) if "rss" in z else sum(x["performance"][z] for x in rr)) for z in fields}
    e=metrics["minimal"]["eval"];sh=metrics["minimal"]["shared"];assert (e["localized"],sh["localized"],e["weak"],sh["weak"],e["controls"],e["parent_controls"],e["rescue_controls"])==(205,113,31,38,53,38,15)
    assert e["lost_joint"]+sh["lost_joint"]+e["stress_loss"]+sh["stress_loss"]==0 and e["stress_cases"]+sh["stress_cases"]==119 and metrics["full"]["eval"]["controls"]==54
    ledger=[{"target":g["target"],"day":g["day"],**x} for g in allr["minimal"] if not g["fit_viewed"] for x in g["residual"]];assert len(ledger)==53
    lanes={}
    for lane in ("parent","rescue"):
        r=[x for x in ledger if x["lane"]==lane];lanes[lane]={"rows":len(r),"hard_controls":sum(x["hard_control"] for x in r),"fail_open":sum(x["fail_open"] for x in r),"legal_rival_support":sum(x["rival_supported"] for x in r),"other_scan_support":sum(x["other_scan_supported"] for x in r),"by_target":dict(collections.Counter(x["target"] for x in r)),"median_on_contrast":float(np.median([np.mean(x["on_contrast"]) for x in r if x["legal_rows"]])),"median_off_contrast":float(np.median([np.mean(x["off_contrast"]) for x in r if x["legal_rows"]]))}
    assert lanes["parent"]["rows"]==38 and lanes["rescue"]["rows"]==15
    ratio={}
    for key in ("cold_total_wall","cold_total_cpu","cold_compute_wall","warm_compute_wall","warm_compute_cpu","cold_peak_rss_kib","payload_bytes"):
        vals=[f["performance"][key]/m["performance"][key] for m,f in zip(allr["minimal"],allr["full"]) if m["performance"][key]>0]
        ratio[key]={"aggregate_full_over_minimal":performance["full"][key]/performance["minimal"][key],"median_paired":statistics.median(vals),"min":min(vals),"max":max(vals)}
    save("residual_ledger.json",ledger)
    result={"status":"PASS_STANDALONE_EXACT_PARITY__HOLD_SCIENTIFIC_CONTROL_BUDGET","groups":11,"paired_processes":22,"rows":3168,"metrics":metrics,"performance":performance,"speedup_full_over_minimal":ratio,"residual":lanes,"residual_ledger_sha256":sha((O/"residual_ledger.json").read_bytes()),"source_pruning_lock_sha256":S["source_r41a"]["lock_sha256"],"no_training":True,"external_validation":False,"production_adoption":False}
    save("result.json",result);cp("99-summary",result)
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="worker":worker(int(ARGS[1]),ARGS[2])
    elif ARGS[0]=="pair":pair(int(ARGS[1]))
    elif ARGS[0]=="summary":summary()
    else:raise ValueError(ARGS)
