#!/usr/bin/env python3
"""R41T: frozen Minimal raw archival replay. No injection, tuning, or origin claim."""
import sys,time
ARGS=sys.argv[1:];BOOT=time.perf_counter();CPU=time.process_time()
import json,pathlib,os,io,zipfile,hashlib,resource,collections
import requests
import numpy as np
import r41b_minimal_hybrid as base
a=base.a;n=base.n
O=pathlib.Path("runtime/r41t");O.mkdir(parents=True,exist_ok=True);PROG=O/"progress";PROG.mkdir(exist_ok=True)
P=pathlib.Path("research/seti/r41t_raw_technosignature_protocol.json");S=json.loads(P.read_text())
ROLES=S["scan"]["roles"];W=S["scan"]["patch"];STRIDE=S["scan"]["stride"];CHUNK=S["scan"]["chunk"]
def sha(x):return hashlib.sha256(x).hexdigest()
def save(name,d):a.write(O/name,d)
def read(name):return json.loads((O/name).read_text())
def cp(name,d):
    a.write(PROG/(name+".json"),d);print("R41T_PROGRESS",json.dumps(d,allow_nan=False),flush=True)
def prepare():
    q=S["source"]["r41b"];headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"};root="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions"
    def get(url):
        rr=requests.get(url,headers=headers,timeout=120);rr.raise_for_status();return rr
    run=get(f"{root}/runs/{q['run']}").json();assert run["status"]=="completed" and run["conclusion"]=="success" and run["head_sha"]==q["head_sha"]
    md=get(f"{root}/artifacts/{q['artifact']}").json();assert md["workflow_run"]["id"]==q["run"] and md["digest"]=="sha256:"+q["archive_sha256"]
    raw=get(f"{root}/artifacts/{q['artifact']}/zip").content;assert sha(raw)==q["archive_sha256"] and len(raw)<64*1024**2
    z=zipfile.ZipFile(io.BytesIO(raw));names=[x.filename for x in z.infolist() if not x.is_dir()]
    assert len(names)==len(set(names)) and sum(x.file_size for x in z.infolist())<512*1024**2 and all(not pathlib.PurePosixPath(x).is_absolute() and ".." not in pathlib.PurePosixPath(x).parts for x in names)
    f=[x for x in names if pathlib.PurePosixPath(x).name=="result.json"];assert len(f)==1
    raw=z.read(f[0]);assert sha(raw)==q["result_sha256"];d=json.loads(raw)
    assert d["status"]=="PASS_STANDALONE_EXACT_PARITY__HOLD_SCIENTIFIC_CONTROL_BUDGET" and d["groups"]==11 and d["paired_processes"]==22
    e=d["metrics"]["minimal"]["eval"];sh=d["metrics"]["minimal"]["shared"]
    assert (e["localized"],sh["localized"],e["weak"],sh["weak"],e["controls"])==(205,113,31,38,53)
    assert d["source_pruning_lock_sha256"]==S["source"]["pruning_lock_sha256"]
    base.prepare();p=base.plan();assert p["pruning_lock_sha256"]==S["source"]["pruning_lock_sha256"]
    selected=[{"index":i,"target":g["group"]["target"],"day":g["group"]["day"],"hashes":g["raw_tile_hashes"]} for i,g in enumerate(p["groups"]) if not g["fit"]]
    assert len(selected)==5 and len({h for g in selected for h in g["hashes"]})==80
    save("plan.json",{"groups":selected,"protocol_sha256":sha(P.read_bytes()),"code_sha256":sha(pathlib.Path(__file__).read_bytes()),"parent":q,"locked_reference":{"localized":318,"weak":69,"controls":53},"no_injection":True})
    cp("00-prepare",{"stage":"PREPARE","groups":5,"native_tiles":80,"status":"PASS_AUTHENTICATED_FROZEN_MINIMAL","targets":[g["target"] for g in selected]})
def row(x,src,alternate=False):
    # Fields 1,4,7,8 are legacy evaluation API placeholders, not real labels.
    return (x,int(alternate),"unknown",0,1. if alternate else 0.,False,src,np.full(16,40 if alternate else 64,int),np.ones(16,bool))
def paths(d):
    parent=d["parent_alert"];joint=parent|d["joint_rescue"];pp=np.where(parent[:,None],d["parent_path"],d["p"]).astype(int)
    return {"r40g":parent,"joint_min":joint,"minimal":d["alert"],"lane":np.where(parent,1,np.where(d["joint_rescue"],2,0)),"path":pp,"joint_score":d["joint_score"]}
def scan(A,X,B,Y,role,window,m,LL,old,test):
    xx,yy,lo=n.ref.common(A,B,X,Y);assert xx.shape==yy.shape and xx.shape[0]==16
    starts=np.asarray(sorted(set(list(range(0,xx.shape[1]-W+1,STRIDE))+[xx.shape[1]-W])),int);df=float(A["df_hz"]);f0=float(A["channel0_mhz"])*1e6+lo*df
    acc={k:[] for k in ("r40g","joint_min","minimal","lane","path","joint_score")};acc["tree_margin"]=[]
    for pos in range(0,len(starts),CHUNK):
        ss=starts[pos:pos+CHUNK];rs=[row(xx[:,k:k+W].astype(float),A["slice_sha256"]) for k in ss];ref=np.stack([yy[:,k:k+W].astype(float) for k in ss])
        d=base.infer("minimal",rs,ref,None,None,LL,old,m,None,"eval");v=paths(d)
        if test and pos==0:
            alt=[row(x[0],x[6],True) for x in rs];chk=base.infer("minimal",alt,ref,None,None,LL,old,m,None,"eval");w=paths(chk)
            assert all(np.array_equal(v[k],w[k]) for k in v) and np.array_equal(d["score"],chk["score"]) and np.array_equal(d["keep"],chk["keep"]),"evaluation placeholder changed decisions"
        for k,x in v.items():acc[k].append(x)
        margin=np.zeros(len(rs));margin[d["idx"]]=d["score"]-m["cuts"]["preserve"];acc["tree_margin"].append(margin)
    z={k:np.concatenate(v) for k,v in acc.items()};assert np.all(~z["minimal"]|z["joint_min"]) and np.all(~z["r40g"]|z["joint_min"])
    assert z["path"].shape==(len(starts),16) and np.all((z["path"]>=0)&(z["path"]<128))
    fn=f"{role}-{window:02d}.npz";return {"role":role,"window_index":window,"center_mhz":float(A["center_mhz"]),"f0_hz":f0,"df_hz":df,"ts":float(A["ts"]),"source_sha256":A["slice_sha256"],"reference_sha256":B["slice_sha256"],"patches":len(starts),"channels":int(xx.shape[1]),"bandwidth_hz":float(xx.shape[1]*abs(df)),"counts":{k:int(z[k].sum()) for k in ("r40g","joint_min","minimal")},"file":fn},starts,z
def freq(part,mode):
    rec,starts,z=part;ids=np.flatnonzero(z[mode]);return ids,rec["f0_hz"]+(starts[ids,None]+z["path"][ids])*rec["df_hz"]
def matches(f,part,mode):
    r=part[0];ix,ff=freq(part,mode);tol=S["association"]["tolerance_bins"]*abs(r["df_hz"]);needs=S["association"]["min_rows"]
    if not len(ix):return {"count":0,"best":None}
    delta=np.abs(ff-f[None,:]);support=(delta<=tol+1e-6).sum(1);j=np.flatnonzero(support>=needs)
    if not len(j):return {"count":0,"best":None}
    # First highest support, then minimum median error, canonical patch index.
    k=min(j,key=lambda i:(-int(support[i]),float(np.median(delta[i])),int(ix[i])))
    return {"count":len(j),"best":{"patch":int(ix[k]),"rows_within_tolerance":int(support[k]),"median_abs_error_hz":float(np.median(delta[k])),"path_hz":ff[k].tolist()}}
def direct(f,anchor_start_hz,A,X):
    df=float(A["df_hz"]);x=(anchor_start_hz-float(A["channel0_mhz"])*1e6)/df;k=int(round(x));p=(f-anchor_start_hz)/df;pp=np.rint(p).astype(int)
    if abs(x-k)>.001 or np.max(abs(p-pp))>.001 or k<0 or k+128>X.shape[1] or not np.all((pp>=5)&(pp<123)):return None
    v=n.ref.norm(X[:,k:k+128].astype(float)[None,:,:]);u=n.normalized_contrast(v)[0,np.arange(16),pp]
    return {"source_sha256":A["slice_sha256"],"start_native_bin":k,"mean_contrast":float(u.mean()),"halves":[float(u[:8].mean()),float(u[8:].mean())],"positive_rows":int((u>0).sum())}
def group(i):
    plan=read("plan.json");assert plan["protocol_sha256"]==sha(P.read_bytes()) and plan["code_sha256"]==sha(pathlib.Path(__file__).read_bytes())
    g0=plan["groups"][i];p=base.plan();g=p["groups"][g0["index"]];num=f"{i+1:02d}";out=O/("group-"+num);out.mkdir(exist_ok=True)
    m,LL,old=base.locks("minimal",p);seed={"n6":a.v.C["seed"],"n7":base.w.C["seed"],"n9":base.old9.C["seed"]}[g["source"]]
    n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);base.b.taylor.init_tree()
    fs=g["group"]["files"];assert [f["role"] for f in fs]==ROLES and len({f["url"] for f in fs})==4
    rawdir=O/("raw-"+num);rawdir.mkdir(exist_ok=True);n.ref.OUT=rawdir;base.b.sh.eng.OUT=rawdir
    net={"payload_bytes":0,"requests":0,"files":0};orig=base.eng.base.RangeFile
    class Counted(orig):
        def close(self):
            if not self.closed and hasattr(self,"total") and hasattr(self,"size"):
                net["payload_bytes"]+=int(self.total);net["requests"]+=int(self.requests_count);net["files"]+=1
            super().close()
    base.eng.base.RangeFile=Counted;t0=time.perf_counter()
    try:pre=[n.ref.tiles(f) for f in fs]
    finally:base.eng.base.RangeFile=orig
    acq=time.perf_counter()-t0;hs=sorted({A["slice_sha256"] for scan0 in pre for A,X in scan0});assert hs==g0["hashes"] and len(hs)==16 and net["files"]==4 and all(len(x)==4 for x in pre)
    records=[];review=[];t0=time.perf_counter()
    for wi in range(4):
        parts={}
        for j,role in enumerate(ROLES):
            A,X=pre[j][wi];B,Y=pre[j^1][wi];assert A["center_mhz"]==B["center_mhz"]
            rec,starts,z=scan(A,X,B,Y,role,wi,m,LL,old,j==0 and wi==0)
            dest=out/rec["file"];np.savez_compressed(dest,starts=starts,**z);rec["sha256"]=sha(dest.read_bytes());parts[role]=(rec,starts,z);records.append(rec)
        ref=parts["ON1"];r,starts,z=ref;ix,ff=freq(ref,"minimal");used=[]
        for j,f in zip(ix,ff):
            if any(np.count_nonzero(abs(f-v)<=4*abs(r["df_hz"])+1e-6)>=12 for v in used):continue
            used.append(f.copy());other={role:{k:matches(f,parts[role],k) for k in ("minimal","joint_min")} for role in ("ON2","OFF1","OFF2")};anchor_start=r["f0_hz"]+starts[j]*r["df_hz"]
            witness={role:direct(f,anchor_start,*pre[k][wi]) for k,role in enumerate(ROLES)}
            repeat=other["ON2"]["minimal"]["count"]>0;off=any(other[v]["joint_min"]["count"]>0 for v in ("OFF1","OFF2"))
            review.append({"window_index":wi,"center_mhz":r["center_mhz"],"patch":int(j),"lane":"R40G" if z["lane"][j]==1 else "Joint_Min_rescue","path_hz":f.tolist(),"start_frequency_hz":float(f[0]),"end_frequency_hz":float(f[-1]),"drift_hz_s":float((f[-1]-f[0])/(15*r["ts"])),"source_sha256":r["source_sha256"],"matches":other,"direct_witness":witness,"repeat_without_matched_off":bool(repeat and not off),"status":"UNVERIFIED_REVIEW_ONLY"})
    compute=time.perf_counter()-t0;counts={k:sum(x["counts"][k] for x in records) for k in ("r40g","joint_min","minimal")}
    result={"stage":"GROUP","group":i+1,"target":g0["target"],"day":g0["day"],"raw_tiles":hs,"records":records,"patches":sum(x["patches"] for x in records),"counts":counts,"on1_anchors":len(review),"on2_matched":sum(x["matches"]["ON2"]["minimal"]["count"]>0 for x in review),"off_matched":sum(any(x["matches"][v]["joint_min"]["count"]>0 for v in ("OFF1","OFF2")) for x in review),"repeat_without_matched_off":sum(x["repeat_without_matched_off"] for x in review),"review":review,"injection":False,"placeholder_invariance":True,"acquisition_wall":acq,"compute_and_receipts_wall":compute,"full_wall":time.perf_counter()-BOOT,"cpu":time.process_time()-CPU,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,**net}
    save("group-"+num+".json",result);cp("group-"+num,{k:v for k,v in result.items() if k not in ("records","review","raw_tiles")})
def summary():
    groups=[read("group-"+f"{i+1:02d}"+".json") for i in range(5)];assert all(not g["injection"] and g["placeholder_invariance"] for g in groups) and len({h for g in groups for h in g["raw_tiles"]})==80
    for g in groups:
        for r in g["records"]:assert sha((O/("group-"+f"{g['group']:02d}")/r["file"]).read_bytes())==r["sha256"]
    d={"status":"PASS_BOUNDED_RAW_REPLAY__REVIEW_ONLY","groups":5,"raw_tiles":80,"totals":{k:sum(g[k] for g in groups) for k in ("patches","on1_anchors","on2_matched","off_matched","repeat_without_matched_off","payload_bytes","requests","files","acquisition_wall","compute_and_receipts_wall","full_wall","cpu")},"counts":{k:sum(g["counts"][k] for g in groups) for k in ("r40g","joint_min","minimal")},"targets":[{k:g[k] for k in ("target","day","patches","counts","on1_anchors","on2_matched","off_matched","repeat_without_matched_off")} for g in groups],"review_repeat_without_matched_off":[{"target":g["target"],"day":g["day"],**x} for g in groups for x in g["review"] if x["repeat_without_matched_off"]],"no_injection":True,"global_false_alarm_calibrated":False,"confirmed_technosignatures":0,"external_validation":False,"production_adoption":False}
    save("result.json",d);cp("99-summary",{k:v for k,v in d.items() if k!="review_repeat_without_matched_off"})
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="group":group(int(ARGS[1]))
    elif ARGS[0]=="summary":summary()
    else:raise ValueError(ARGS)
