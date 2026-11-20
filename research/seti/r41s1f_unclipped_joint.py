#!/usr/bin/env python3
"""Paired native shoulder contrast with lane-separated frozen retrospective ablation."""
import sys,time
ARGS=sys.argv[1:];BOOT=time.perf_counter();CPU=time.process_time()
import json,pathlib,os,io,zipfile,resource
import numpy as np
import requests
from scipy.signal import savgol_filter
from scipy.ndimage import median_filter
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import r41s1d_retina_rubik as retina
parent=retina.parent;audit=parent.audit;s=parent.s;base=s.base;a=s.a;n=s.n;HALO=parent.HALO
O=pathlib.Path("runtime/r41s1f");O.mkdir(parents=True,exist_ok=True)
PC=retina.PC;C=parent.C;PROG=O/"progress";PROG.mkdir(exist_ok=True)
P=pathlib.Path("research/seti/r41s1f_unclipped_joint_protocol.json");L=json.loads(P.read_text())
def save(k,d):a.write(O/k,d)
def read(k):return json.loads((O/k).read_text())
def cp(k,d):
    a.write(PROG/(k+".json"),d);print("R41S1F_PROGRESS",json.dumps(d,allow_nan=False),flush=True)
def mean(x):return float(np.mean(x))
def rms(x):return float(np.sqrt(np.mean(x*x)))
def mad(x):return float(1.4826*np.median(abs(x-np.median(x))))
def corr(x,y):
    x=x-x.mean();y=y-y.mean();d=np.linalg.norm(x)*np.linalg.norm(y);return float(np.clip(x@y/d,-1,1)) if d>1e-12 else None
def div(a,b):return float(a/b) if b>1e-12 else None
def noise(x):
    med=x-median_filter(x,size=3,mode="nearest");sg=x-savgol_filter(x,5,2,mode="interp");y=x-x.mean();e=abs(np.fft.rfft(y))**2;e[1:8]*=2
    return {"mean":mean(x),"std":float(x.std()),"median_rms":rms(med),"savgol_rms":rms(sg),"mad":mad(med),"hf":div(e[4:].sum(),e[1:].sum()),"haar":div(((y[::2]-y[1::2])**2).sum()/2,(y*y).sum()),"lag1":corr(sg[:-1],sg[1:]),"diff_rms":rms(np.diff(x))}
SIDE=retina.SIDE
def rawcore(x):return x[...,24]-np.median(x[...,SIDE],axis=-1)
def lay(x):
    v=x-np.median(x[...,SIDE],axis=-1,keepdims=True);occ=(v>0).sum(axis=-3)>=3;av=v.mean(axis=-3);mi=v.min(axis=-3);return v,occ,av,mi,occ[...,24].mean(axis=-1),rawcore(mi).mean(axis=-1)
def raw_diagnose(cube,p):
    # Original power arrays, not normalized or clipped. Per-scan pooled shoulder scale.
    assert cube.shape==(4,16,208) and np.all((p>=24)&(p<=183)) and np.isfinite(cube).all()
    strip=cube[:,np.arange(16)[:,None],p[:,None]+np.arange(-24,25)[None,:]].astype(float)
    shoulders=strip[:,:,SIDE];bg=np.median(shoulders,axis=-1,keepdims=True)
    scale=np.array([mad(z) for z in (shoulders-bg)])
    if not np.all(np.isfinite(scale)&(scale>0)):return None
    z=(strip-bg)/scale[:,None,None];feat={};details={};cores=[]
    for i,k in enumerate(retina.FACES):
        x=z[i];core=rawcore(x);cores.append(core);layer=x.reshape(4,4,49);v,occ,av,mi,ov,mc=lay(layer)
        idx=(np.arange(49)[None,None,None,:]-retina.SHIFTS[:,i,:,None,None])%49
        shifted=np.take_along_axis(np.broadcast_to(layer,(32,4,4,49)),np.broadcast_to(idx,(32,4,4,49)),axis=-1)
        _,_,_,_,no,nm=lay(shifted);pos=np.maximum(av,0)
        lf={"mean":mean(rawcore(av)),"min":float(mc),"median":mean(rawcore(np.median(v,axis=0))),"overlap":float(ov),"dispersion":mean(np.std(v[:,:,24],axis=0)),"leave_min":min(mean(rawcore(np.delete(v,j,axis=0).mean(0))) for j in range(4)),"held_min":float(core.reshape(4,4).mean(1).min()),"concentration":div(pos[:,23:26].sum(),pos[:,12:37].sum()),"null_overlap_excess":float(ov-np.quantile(no,.9)),"null_min_excess":float(mc-np.quantile(nm,.9))}
        feat.update({"retina_"+k+"_"+j:w for j,w in lf.items()});feat.update({"noise_"+k+"_"+j:w for j,w in noise(core).items()})
        norm=s.norm(cube[i]);r=np.arange(16);pc=norm[r,p]-np.median(norm[r[:,None],p[:,None]+np.array([-5,-4,-3,3,4,5])],axis=1);cl=np.clip(pc,-4,4)
        assert np.array_equal(cl,s.contrast(norm,p))
        feat.update({"norm_"+k+"_preclip_mean":mean(pc),"norm_"+k+"_preclip_std":float(pc.std()),"norm_"+k+"_clip_mean":mean(cl),"norm_"+k+"_exceeds4":mean(abs(pc)>4),"norm_"+k+"_equals4_after":mean(abs(cl)>=4)})
        details[k]={"raw_power_contrast":(strip[i,:,24]-bg[i,:,0]).tolist(),"raw_shoulder_scale":float(scale[i]),"unclipped_shoulder_z":core.tolist(),"normalized_preclip":pc.tolist(),"normalized_clipped":cl.tolist(),"input_sha256":s.sha(np.ascontiguousarray(cube[i]).tobytes())}
    u,v,w,z=cores
    feat.update(noise_on_repeat_disagreement=rms(u-w),noise_off_repeat_disagreement=rms(v-z),noise_on_repeat_corr=corr(u,w),noise_off_repeat_corr=corr(v,z),retina_primary_on_off=mean(u)-mean(v),retina_repeat_on_off=mean(w)-mean(z),retina_shared_min=min(mean(u),mean(v)))
    assert all(x is None or np.isfinite(x) for x in feat.values())
    return {"features":feat,"details":details,"scope":"fixed native receiver path, shoulder normalized without clipping; not calibrated physical noise likelihood"}
def download(spec):
    root="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions";h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    def get(url):
        r=requests.get(url,headers=h,timeout=120);r.raise_for_status();return r
    md=get(f"{root}/artifacts/{spec['artifact']}").json();assert md["workflow_run"]["id"]==spec["run"] and md["digest"]=="sha256:"+spec["sha256"]
    raw=get(f"{root}/artifacts/{spec['artifact']}/zip").content;assert s.sha(raw)==spec["sha256"] and len(raw)<256*1024**2
    z=zipfile.ZipFile(io.BytesIO(raw));names=z.namelist();assert len(names)==len(set(names)) and sum(x.file_size for x in z.infolist())<1024*1024**2
    root=pathlib.PurePosixPath(next(x for x in names if pathlib.PurePosixPath(x).name=="result.json")).parent
    def rd(k):return z.read(str(root/k))
    d=rd("result.json");assert s.sha(d)==spec["result_sha256"];return json.loads(d),rd
def prepare():
    retina.prepare()
    d,rd=download(L["source"]["retina"]);e,re=download(L["source"]["noise"]);assert d["qualified"]==e["qualified"]==891 and d["decision_changes"]==e["decision_changes"]==0
    for k,h in d["trace_sha256"].items():assert s.sha(rd(k))==h
    assert len(d["trace_sha256"])==22
    raw=re("records.json");assert s.sha(raw)==e["records_sha256"];nr=json.loads(raw)
    for i in range(1,12):
        k=f"{i:02d}";(O/f"retina-{k}.json").write_bytes(rd(f"group-{k}.json"));save(f"noise-{k}.json",[x for x in nr if x["group"]==i])
    # Explicit contrast scale test: no saturation at >4 and normalized clipping parity.
    rng=np.random.default_rng(413401);p=np.full(16,80);x=rng.normal(size=(4,16,208));x[0,np.arange(16),p]+=20
    q=raw_diagnose(x,p);assert q is not None and max(q["details"]["primary_on"]["unclipped_shoulder_z"])>4 and max(q["details"]["primary_on"]["normalized_clipped"])<=4
    save("prereq_receipt.json",{"operator_sha256":s.sha(pathlib.Path(__file__).read_bytes()),"protocol_sha256":s.sha(P.read_bytes()),"bank_sha256":s.sha((parent.O/"bank.npy").read_bytes()),"archive_traces":22,"noise_records":891,"unclipped_selftest":True});cp("00-prepare",read("prereq_receipt.json"))
def group(i):
    p0=read("prereq_receipt.json");assert p0["operator_sha256"]==s.sha(pathlib.Path(__file__).read_bytes()) and p0["protocol_sha256"]==s.sha(P.read_bytes()) and p0["bank_sha256"]==s.sha((parent.O/"bank.npy").read_bytes())
    bank=np.load(parent.O/"bank.npy",allow_pickle=False);oldbank=np.load(audit.O/"bank.npy",allow_pickle=False);g=base.plan()["groups"][i];num=str(i+1).zfill(2)
    prior=json.loads((C/("group-"+num+".json")).read_text());old=json.loads((audit.C/("C-"+num+".json")).read_text());ref=json.loads((audit.C/("parent-"+num+".json")).read_text())
    m,LL,old0=base.locks("minimal",base.plan());seed={"n6":a.v.C["seed"],"n7":base.w.C["seed"],"n9":base.old9.C["seed"]}[g["source"]];n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);base.b.taylor.init_tree()
    fs=g["group"]["files"];net={"payload_bytes":0,"requests":0,"files":0};orig=base.eng.base.RangeFile
    class Counted(orig):
        def close(self):
            if not self.closed and hasattr(self,"total") and hasattr(self,"size"):
                net["payload_bytes"]+=int(self.total);net["requests"]+=int(self.requests_count);net["files"]+=1
            super().close()
    base.eng.base.RangeFile=Counted;t0=time.perf_counter()
    try:rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r41s1f_"+num)
    finally:base.eng.base.RangeFile=orig
    acq=time.perf_counter()-t0;assert sorted(hs)==g["raw_tile_hashes"]==ref["raw_tiles"] and net["files"]==4
    out={"stage":"GROUP","group":i+1,"target":g["group"]["target"],"fit_viewed":g["fit"],"parity":True,"records":[],"acquisition_wall":acq,"baseline_wall":0.,"parent_replay_wall":0.,"diagnostic_wall":0.,**net}
    for kind,rows,refs,mm in (("eval",rr,R,meta),("shared",[x[0] for x in ch],np.stack([x[1] for x in ch]),cm)):
        t=time.perf_counter();d=base.infer("minimal",rows,refs,mm,pre,LL,old0,m,None,kind);out["baseline_wall"]+=time.perf_counter()-t
        oldarr=a.read_npz(base.C/("dev_"+num+"_"+kind+".npz"));a.verify_rows(oldarr,rows);assert base.check(d,oldarr,m)==ref["metrics"][kind]
        mask=s.sha(np.ascontiguousarray(d["alert"]).tobytes()+np.ascontiguousarray(d["hit"]).tobytes());assert mask==old["baseline_masks"][kind];before=d["alert"].copy();hits=d["hit"].copy()
        original=json.loads((audit.C/("trace-"+num+"-"+kind+".json")).read_text())["normal"]
        archived=json.loads((C/("trace-"+num+"-"+kind+".json")).read_text());assert archived["baseline_mask"]==mask
        records=archived["records"];assert [x["row"] for x in records]==[x["row"] for x in original]==np.flatnonzero(before).tolist()
        sealed=json.loads((PC/("trace-"+num+"-"+kind+".json")).read_text());assert sealed["baseline_mask"]==mask and [v["row"] for v in sealed["records"]]==[v["row"] for v in records]
        oldd=read("retina-"+num+".json"); oldd=[v for v in oldd["records"] if v["kind"]==kind]; olde=read("noise-"+num+".json"); olde={(v["kind"],v["row"]):v for v in olde}
        assert [v["row"] for v in oldd]==[v["row"] for v in records]
        advtrace=[]
        for x,y0,previous,dr in zip(records,original,sealed["records"],oldd):
            j=x["row"];lane="parent" if d["parent_alert"][j] else "rescue";path=d["parent_path"][j] if lane=="parent" else d["p"][j];src=str(d["source_tile_sha"][j]);assert lane==x["lane"]==y0["lane"] and src==x["source_tile_sha256"]==y0["source_tile_sha256"]
            t=time.perf_counter();tr=s.spiral(rows[j][0],refs[j],path,lane,mm[j],src,fs,pre,oldbank);assert tr==y0["trace"] and not tr["veto"]
            ad=audit.advisory(tr,lane,mm[j],src,fs,pre);assert ad==x["advisory"] and not ad["veto"];out["parent_replay_wall"]+=time.perf_counter()-t
            h=list(ad["states"]);features=None;clipped=None;raw=None
            if ad["advanced"]:
                t=time.perf_counter();pv=s.provenance(fs,pre,mm[j],src);assert pv is not None
                blocks,receipts=pv;slot=mm[j][0];oth=2-slot
                au=blocks[slot].copy();af=blocks[slot+1].copy();au[:,HALO:HALO+128]=rows[j][0];af[:,HALO:HALO+128]=refs[j]
                z=parent.diagnose(au,af,blocks[oth],blocks[oth+1],np.asarray(path,int)+HALO,float(mm[j][3]),float(receipts[slot]["ts"]),bank)
                e={**h[-1]["snapshot"]["evidence"],"differentiation":z,"parent_measurement_hashes":[s.sha(v.tobytes()) for v in (au,af,blocks[oth],blocks[oth+1])]}
                assert s.state(e,h,["physical_origin","independent_confirmation"]) and h[:-1]==ad["states"];audit.chain(h)
                assert h==previous["states"] and previous["advanced"] and previous["source_tile_sha256"]==src and not previous["decision_changed"],"R41S1C sealed replay"
                out["parent_replay_wall"]+=time.perf_counter()-t
                t=time.perf_counter();cube=np.stack([au,af,blocks[oth],blocks[oth+1]]);p=np.asarray(path,int)+HALO
                z2=retina.diagnose(cube,p);assert z2["features"]==dr["features"] and dr["supported"] and dr["lane"]==lane
                ne=olde[(kind,j)];assert ne["parent_state_sha256"]==previous["states"][-1]["state_sha256"] and ne["source_tile_sha256"]==src and ne["lane"]==lane and ne["label_evaluation_only"]==int(d["y"][j]);clipped={**{"retina_"+k:v for k,v in z2["features"].items()},**{"noise_"+k:v for k,v in ne["features"].items()}}
                raw=raw_diagnose(cube,p);features=raw["features"] if raw is not None else None
                e2={**h[-1]["snapshot"]["evidence"],"retina_rubik":z2}
                assert s.state(e2,h,["physical_origin","independent_confirmation"]) and h[:-1]==previous["states"];audit.chain(h)
                out["diagnostic_wall"]+=time.perf_counter()-t
            else:
                assert h==previous["states"] and not previous["advanced"] and not previous["decision_changed"] and not dr["supported"]
            out["records"].append({"row":j,"kind":kind,"lane":lane,"supported":ad["advanced"],"features":features,"clipped":clipped,"weak":bool(d["weak"][j]),"label_evaluation_only":int(d["y"][j]),"localized_evaluation_only":bool(d["hit"][j])})
            advtrace.append({"row":j,"lane":lane,"source_tile_sha256":src,"parent_state_sha256":h[-1]["state_sha256"],"raw":raw,"advanced":ad["advanced"],"decision_changed":False})
        assert np.array_equal(before,d["alert"]) and np.array_equal(hits,d["hit"])
        save("trace-"+num+"-"+kind+".json",{"records":advtrace,"baseline_mask":mask,"decision_changes":0})
    count=sum(x["supported"] for x in out["records"]);assert count==prior["summary"]["qualified"] and len(out["records"])==prior["summary"]["triggered"]
    out.update(triggered=len(out["records"]),qualified=count,diagnostic_advanced=count,decision_changes=0,full_wall=time.perf_counter()-BOOT,cpu=time.process_time()-CPU,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    save("group-"+num+".json",out);cp("group-"+num,{k:v for k,v in out.items() if k!="records"})

def arm(row,name):
    f=row["features"]
    if name=="clipped_combined":return row["clipped"]
    if f is None:return None
    return {k:v for k,v in f.items() if k.startswith("retina_") or k.startswith("noise_")} if name=="raw_combined" else {k:v for k,v in f.items() if k.startswith(name.split("_")[1]+"_")}
def fit_arm(train,test,name):
    tr=[r for r in train if arm(r,name) is not None];te=[r for r in test if arm(r,name) is not None]
    if len(set(r["label_evaluation_only"] for r in tr))!=2:return {"status":"HOLD_TRAIN_CLASSES","predictions":[]}
    keys=sorted(arm(tr[0],name));assert all(sorted(arm(r,name))==keys for r in tr+te)
    def mat(rows):return np.array([[np.nan if arm(r,name)[k] is None else arm(r,name)[k] for k in keys] for r in rows],float)
    X=mat(tr);missing=np.all(~np.isfinite(X),axis=0);im=np.zeros(len(keys));im[~missing]=np.nanmedian(X[:,~missing],axis=0)
    def imp(x):return np.where(np.isfinite(x),x,im)[:,~missing]
    X=imp(X);Y=np.array([r["label_evaluation_only"] for r in tr]);sc=StandardScaler().fit(X);model=LogisticRegression(C=.1,class_weight="balanced",max_iter=2000,random_state=413401).fit(sc.transform(X),Y)
    train_scores=model.decision_function(sc.transform(X));cut=float(train_scores[Y==1].min());score=model.decision_function(sc.transform(imp(mat(te)))) if te else []
    pred=[{"group":r["group"],"kind":r["kind"],"row":r["row"],"score":float(v),"would_keep":bool(v>=cut)} for r,v in zip(te,score)]
    yy=[r["label_evaluation_only"] for r in te];auc=float(roc_auc_score(yy,score)) if len(set(yy))==2 else None
    return {"status":"HISTORICAL_ONLY","features":keys,"train_positive_threshold":cut,"train_n":len(tr),"test_n":len(te),"audit_auc_positive":auc,"coefficients":model.coef_[0].tolist(),"scaler_mean":sc.mean_.tolist(),"scaler_scale":sc.scale_.tolist(),"intercept":float(model.intercept_[0]),"imputation":im.tolist(),"dropped_all_missing":missing.tolist(),"predictions":pred}
def summary():
    gr=[read(f"group-{i:02d}.json") for i in range(1,12)];rows=[{**v,"group":g["group"],"fit_viewed":g["fit_viewed"]} for g in gr for v in g["records"]];assert sum(r["supported"] for r in rows)==891 and all(x["decision_changes"]==0 for x in gr)
    train=[v for v in rows if v["fit_viewed"]];test=[v for v in rows if not v["fit_viewed"]]
    def counts(rr):return {"controls":sum(v["label_evaluation_only"]==0 and v["kind"]=="eval" for v in rr),"localized":sum(v["localized_evaluation_only"] for v in rr),"weak":sum(v["localized_evaluation_only"] and v["weak"] for v in rr),"ordinary":sum(v["localized_evaluation_only"] and v["kind"]=="eval" for v in rr),"shared":sum(v["localized_evaluation_only"] and v["kind"]=="shared" for v in rr)}
    before=counts(test);assert before["controls"]==53 and before["localized"]==318 and before["weak"]==69
    models={};tt=time.perf_counter()
    for name in ("raw_retina","raw_noise","raw_combined","clipped_combined"):
        models[name]={}
        for lane in ("parent","rescue"):models[name][lane]=fit_arm([v for v in train if v["lane"]==lane],[v for v in test if v["lane"]==lane],name)
        keep={(p["group"],p["kind"],p["row"]):p["would_keep"] for lane in ("parent","rescue") for p in models[name][lane]["predictions"]}
        retained=[v for v in test if keep.get((v["group"],v["kind"],v["row"]),True)];models[name]["hypothetical"]=counts(retained);models[name]["lost"]={k:before[k]-models[name]["hypothetical"][k] for k in before};models[name]["lanes"]={lane:{"before":counts([v for v in test if v["lane"]==lane]),"retained":counts([v for v in retained if v["lane"]==lane])} for lane in ("parent","rescue")}
    result={"status":"PASS_PAIRED_UNCLIPPED_ABLATION__HOLD_FRESH_VALIDATION","qualified":891,"raw_supported":sum(v["features"] is not None for v in rows),"baseline":before,"decision_changes":0,"models":models,"fit_wall_s":time.perf_counter()-tt,"source":L["source"],"historical_only":True,"production_adoption":False,"external_validation":False}
    save("result.json",result);cp("99-summary",{k:v for k,v in result.items() if k!="models"}|{"arms":{k:{"hypothetical":v["hypothetical"],"lost":v["lost"],"lanes":v["lanes"],"auc":{l:v[l].get("audit_auc_positive") for l in ("parent","rescue")}} for k,v in models.items()}})
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="group":group(int(ARGS[1]))
    elif ARGS[0]=="summary":summary()
