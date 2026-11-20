#!/usr/bin/env python3
"""R41S1C: morphology / moving native paths / instrument proxies; advisory only."""
import sys,time
ARGS=sys.argv[1:];BOOT=time.perf_counter();CPU=time.process_time()
import json,pathlib,os,io,zipfile,resource,collections
import numpy as np
import requests
from sklearn.metrics import roc_auc_score
import r41s1b_transition_audit as audit
s=audit.s;base=s.base;a=s.a;n=s.n
O=pathlib.Path("runtime/r41s1c");O.mkdir(parents=True,exist_ok=True)
C=O/"parent";C.mkdir(exist_ok=True);PROG=O/"progress";PROG.mkdir(exist_ok=True)
P=pathlib.Path("research/seti/r41s1c_evidence_differentiation_protocol.json");L=json.loads(P.read_text());T=L["limits"];HALO=s.T["halo_bins"];K=("eval","shared")
FEATURES=("on_shape_cosine","off_shape_cosine","on_time_pearson","off_time_pearson","on_width_delta","off_width_delta","on_repeat_width","off_repeat_width","on_repeat_core_shoulder","off_repeat_core_shoulder","on_heldout_gain","off_heldout_gain","on_heldout_score","off_heldout_score","on_consensus_fraction","off_consensus_fraction","on_abs_offset_hz","off_abs_offset_hz","on_slope_delta_abs_hz_s","off_slope_delta_abs_hz_s","on_off_path_agreement","on_off_persistent_jaccard","repeat_on_occupancy","repeat_off_occupancy","repeat_on_persistent_fraction","repeat_off_persistent_fraction","repeat_on_broad_row_mean","repeat_off_broad_row_mean","repeat_on_fixed_contrast","repeat_off_fixed_contrast","repeat_off_minus_on_fixed","repeat_on_active","repeat_off_active")
FEATURES=FEATURES+tuple(k+"_template_"+v for k in ("on","off") for v in ("heldout_gain","heldout_score","abs_offset_hz","shape_cosine","width","consensus_fraction"))
def save(name,d):a.write(O/name,d)
def read(name):return json.loads((O/name).read_text())
def cp(name,d):
    a.write(PROG/(name+".json"),d);print("R41S1C_PROGRESS",json.dumps(d,allow_nan=False),flush=True)
def mean(x):return float(np.mean(x))
def makebank():
    paths=[]
    for d in base.b.D:
        d=np.asarray(d,np.int64);lo=5-int(d.min());hi=203-int(d.max())
        if hi>lo:paths.extend(np.arange(lo,hi)[:,None]+d)
    bank=np.asarray(paths,np.int64)
    assert bank.ndim==2 and bank.shape[1]==16 and np.all((bank>=5)&(bank<203))
    return bank
def native_contrast(z):
    assert z.shape==(16,208)
    idx=np.arange(5,203);off=np.array([-5,-4,-3,3,4,5]);v=np.zeros_like(z)
    v[:,idx]=np.clip(z[:,idx]-np.median(z[:,idx[:,None]+off],axis=2),-4.,4.)
    return v
def search(z,p,bank,df,ts):
    u=native_contrast(z);r=np.arange(16);vals=u[r[None,:],bank]
    i=int(vals[:,:8].sum(1).argmax());j=int(vals[:,8:].sum(1).argmax())
    f=bank[i];b=bank[j];held=np.r_[vals[j,:8],vals[i,8:]];fixed=u[r,p]
    hp=np.r_[b[:8],f[8:]] # two independent heldout fragments, not a continuous physical track
    stats=lambda x:{"start_offset_hz":float((x[0]-p[0])*df),"end_offset_hz":float((x[-1]-p[-1])*df),"slope_hz_s":float((x[-1]-x[0])*df/(15*ts)),"curvature_hz_per_sample2":mean(np.abs(np.diff(x,2)))*abs(df),"edge":bool(x.min()<12 or x.max()>195)}
    result={"forward_fit_path":f.tolist(),"backward_fit_path":b.tolist(),"heldout_fragments":hp.tolist(),"forward_fit_id":i,"backward_fit_id":j,"searched_paths":len(bank),"heldout_blocks":s.halves(held),"fixed_blocks":s.halves(fixed),"heldout_gain_blocks":s.halves(held-fixed),"consensus_fraction":mean(abs(f-b)<=2),"forward_geometry":stats(f),"backward_geometry":stats(b),"mean_abs_offset_hz":mean(abs(hp-p))*abs(df),"slope_delta_abs_hz_s":mean([abs((x[-1]-x[0])-(p[-1]-p[0]))*abs(df)/(15*ts) for x in (f,b)])}
    return result,hp,held,fixed
def profile(z,p,values=None):
    p=np.asarray(p,int);v=s.contrast(z,p) if values is None else values
    out={"trace":v.tolist(),"halves":s.halves(v),"active":mean(v>0),"width":None,"core_shoulder":None,"profile":None}
    if not np.all((p>=12)&(p<=195)):return out
    band=z[np.arange(16)[:,None],p[:,None]+np.arange(-12,13)[None,:]]
    floor=np.median(np.c_[band[:,:5],band[:,-5:]],axis=1)
    excess=np.maximum(band-floor[:,None],0.);prof=excess.mean(0);den=prof.sum()
    out.update(profile=prof.tolist(),width=float(np.sum(prof*np.abs(np.arange(-12,13)))/den) if den>0 else None,core_shoulder=float(prof[11:14].mean()-np.r_[prof[4:9],prof[16:21]].mean()))
    return out
def cosine(x,y):
    if x is None or y is None:return None
    x=np.asarray(x,float);y=np.asarray(y,float);d=np.linalg.norm(x)*np.linalg.norm(y)
    return float(np.clip(x@y/d,-1,1)) if d>1e-12 else None
def pearson(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    return cosine(x-x.mean(),y-y.mean())
def occupancy(z):
    high=z>T["normalized_occupancy_threshold"];persist=high.mean(0)>T["channel_persistence_fraction"]
    return {"occupancy":mean(high),"persistent_fraction":mean(persist),"broad_row_mean":mean(high.mean(1)**2)},persist
def diagnose(primary_on,primary_off,repeat_on,repeat_off,p,df,ts,bank):
    assert df!=0 and ts>0 and all(x.shape==(16,208) and np.isfinite(x).all() for x in (primary_on,primary_off,repeat_on,repeat_off))
    z=[s.norm(x) for x in (primary_on,primary_off,repeat_on,repeat_off)]
    pa,pb=profile(z[0],p),profile(z[1],p);result={};features={};hps={};occ={}
    shift_bank=p[None,:]+np.arange(5-int(p.min()),203-int(p.max()))[:,None]
    assert len(shift_bank)>0 and np.all((shift_bank>=5)&(shift_bank<203))
    for k,pri,rep in (("on",pa,z[2]),("off",pb,z[3])):
        q,hp,held,fixed=search(rep,p,bank,df,ts);pr=profile(rep,hp,held);oc,persist=occupancy(rep);hps[k]=hp;occ[k]=persist
        result[k]={"primary_morphology":pri,"repeat_crossfit_morphology":pr,"trajectory":q,"instrument_proxies":oc}
        features.update({k+"_shape_cosine":cosine(pri["profile"],pr["profile"]),k+"_time_pearson":pearson(pri["trace"],pr["trace"]),k+"_width_delta":None if pri["width"] is None or pr["width"] is None else pr["width"]-pri["width"],k+"_repeat_width":pr["width"],k+"_repeat_core_shoulder":pr["core_shoulder"],k+"_heldout_gain":mean(held-fixed),k+"_heldout_score":mean(held),k+"_consensus_fraction":q["consensus_fraction"],k+"_abs_offset_hz":q["mean_abs_offset_hz"],k+"_slope_delta_abs_hz_s":q["slope_delta_abs_hz_s"],"repeat_"+k+"_occupancy":oc["occupancy"],"repeat_"+k+"_persistent_fraction":oc["persistent_fraction"],"repeat_"+k+"_broad_row_mean":oc["broad_row_mean"],"repeat_"+k+"_fixed_contrast":mean(fixed),"repeat_"+k+"_active":pr["active"]})
        tq,thp,theld,tfixed=search(rep,p,shift_bank,df,ts);tp=profile(rep,thp,theld)
        result[k]["template_shift"]={"trajectory":tq,"crossfit_morphology":tp,"family":"translation_only"}
        features.update({k+"_template_heldout_gain":mean(theld-tfixed),k+"_template_heldout_score":mean(theld),k+"_template_abs_offset_hz":tq["mean_abs_offset_hz"],k+"_template_shape_cosine":cosine(pri["profile"],tp["profile"]),k+"_template_width":tp["width"],k+"_template_consensus_fraction":tq["consensus_fraction"]})
    union=(occ["on"]|occ["off"]).sum()
    features.update(on_off_path_agreement=mean(abs(hps["on"]-hps["off"])<=2),on_off_persistent_jaccard=float((occ["on"]&occ["off"]).sum()/union) if union else None,repeat_off_minus_on_fixed=features["repeat_off_fixed_contrast"]-features["repeat_on_fixed_contrast"])
    assert set(features)==set(FEATURES) and all(v is None or np.isfinite(v) for v in features.values())
    return {"features":features,"details":result,"scope":"fixed receiver-native crossfit fragments; not extrapolated physical trajectory or instrument identity"}
def synthetic(bank):
    rng=np.random.default_rng(T["synthetic_seed"]);rows=[];r=np.arange(16);classes=("noise","narrow_drifting","broad_burst","shared_stationary_proxy")
    def add(x,p,mask,offs):
        x=x.copy();scale=1.4826*np.median(abs(x-np.median(x,axis=1)[:,None]),axis=1)
        for d in offs:x[r[mask],p[mask]+d]+=T["synthetic_amplitude"]*scale[mask]
        return x
    for i in range(T["synthetic_n_per_class"]):
        raw=[rng.normal(size=(16,208)) for _ in range(4)]
        for cl in classes:
            o,f,u,v=[x.copy() for x in raw];p=80+r if cl=="narrow_drifting" else np.full(16,80);full=np.ones(16,bool);burst=(r>=4)&(r<12)
            if cl=="narrow_drifting":
                o=add(o,p,full,[0]);u=add(u,p+T["shift_bins"],full,[0])
            elif cl=="broad_burst":
                off=range(-8,9);o=add(o,p,burst,off);u=add(u,p+T["shift_bins"],burst,off)
            elif cl=="shared_stationary_proxy":
                o,f,u,v=[add(x,p,full,[0]) for x in (o,f,u,v)]
            d=diagnose(o,f,u,v,p,2.8,1.,bank);rows.append({"pair":i,"class":cl,"features":d["features"]})
    def med(cl,k):
        v=[x["features"][k] for x in rows if x["class"]==cl and x["features"][k] is not None]
        return float(np.median(v)) if v else None
    medians={cl:{k:med(cl,k) for k in FEATURES} for cl in classes}
    full_bank_check=med("narrow_drifting","on_heldout_gain")>0
    d={"shifted_narrow_template_positive_heldout_gain":med("narrow_drifting","on_template_heldout_gain")>0,"broad_more_occupied_than_narrow":med("broad_burst","repeat_on_occupancy")>med("narrow_drifting","repeat_on_occupancy"),"shared_off_fixed_stronger_than_narrow":med("shared_stationary_proxy","repeat_off_fixed_contrast")>med("narrow_drifting","repeat_off_fixed_contrast")}
    save("synthetic_controls.json",{"seed":T["synthetic_seed"],"synthetic_only":True,"rows":rows,"medians":medians,"direction_checks":d,"full_bank_direction_check_diagnostic_only":full_bank_check});assert all(d.values()),d
    return d
def prepare():
    spec=L["source"];head={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"};root="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions"
    def get(url):
        rr=requests.get(url,headers=head,timeout=120);rr.raise_for_status();return rr
    run=get(f"{root}/runs/{spec['run']}").json();assert run["status"]=="completed" and run["conclusion"]=="success" and run["head_sha"]==spec["head_sha"]
    md=get(f"{root}/artifacts/{spec['artifact']}").json();assert md["workflow_run"]["id"]==spec["run"] and md["digest"]=="sha256:"+spec["archive_sha256"]
    raw=get(f"{root}/artifacts/{spec['artifact']}/zip").content;assert s.sha(raw)==spec["archive_sha256"] and len(raw)<128*1024**2
    z=zipfile.ZipFile(io.BytesIO(raw));names=[x.filename for x in z.infolist() if not x.is_dir()];assert len(names)==len(set(names)) and sum(x.file_size for x in z.infolist())<512*1024**2 and all(not pathlib.PurePosixPath(x).is_absolute() and ".." not in pathlib.PurePosixPath(x).parts for x in names)
    assert s.sha(z.read("result.json"))==spec["result_sha256"];d=json.loads(z.read("result.json"));assert d["groups"]==11 and d["decision_changes"]==0 and d["status"]=="PASS_TRANSITION_DECOUPLED_ADVISORY__HOLD_SCIENTIFIC_SCREEN" and len(d["trace_sha256"])==22
    for k,h in d["trace_sha256"].items():
        assert s.sha(z.read(k))==h;(C/k).write_bytes(z.read(k))
    for i in range(1,12):
        k=f"group-{i:02d}.json";(C/k).write_bytes(z.read(k))
    # Authenticate sealed original C too; restore frozen baseline locks and original bank.
    audit.prepare();bank=makebank();np.save(O/"bank.npy",bank,allow_pickle=False)
    checks=synthetic(bank)
    receipt={"stage":"PREPARE","source":spec,"source_traces_verified":22,"bank_paths":len(bank),"operator_sha256":s.sha(pathlib.Path(__file__).read_bytes()),"protocol_sha256":s.sha(P.read_bytes()),"bank_sha256":s.sha((O/"bank.npy").read_bytes()),"synthetic_checks":checks}
    save("prereq_receipt.json",receipt);cp("00-prepare",receipt)
def group(i):
    p0=read("prereq_receipt.json");assert p0["operator_sha256"]==s.sha(pathlib.Path(__file__).read_bytes()) and p0["protocol_sha256"]==s.sha(P.read_bytes()) and p0["bank_sha256"]==s.sha((O/"bank.npy").read_bytes())
    bank=np.load(O/"bank.npy",allow_pickle=False);oldbank=np.load(audit.O/"bank.npy",allow_pickle=False);g=base.plan()["groups"][i];num=str(i+1).zfill(2)
    prior=json.loads((C/("group-"+num+".json")).read_text());old=json.loads((audit.C/("C-"+num+".json")).read_text());ref=json.loads((audit.C/("parent-"+num+".json")).read_text())
    m,LL,old0=base.locks("minimal",base.plan());seed={"n6":a.v.C["seed"],"n7":base.w.C["seed"],"n9":base.old9.C["seed"]}[g["source"]];n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);base.b.taylor.init_tree()
    fs=g["group"]["files"];net={"payload_bytes":0,"requests":0,"files":0};orig=base.eng.base.RangeFile
    class Counted(orig):
        def close(self):
            if not self.closed and hasattr(self,"total") and hasattr(self,"size"):
                net["payload_bytes"]+=int(self.total);net["requests"]+=int(self.requests_count);net["files"]+=1
            super().close()
    base.eng.base.RangeFile=Counted;t0=time.perf_counter()
    try:rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r41s1c_"+num)
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
        advtrace=[]
        for x,y0 in zip(records,original):
            j=x["row"];lane="parent" if d["parent_alert"][j] else "rescue";path=d["parent_path"][j] if lane=="parent" else d["p"][j];src=str(d["source_tile_sha"][j]);assert lane==x["lane"]==y0["lane"] and src==x["source_tile_sha256"]==y0["source_tile_sha256"]
            t=time.perf_counter();tr=s.spiral(rows[j][0],refs[j],path,lane,mm[j],src,fs,pre,oldbank);assert tr==y0["trace"] and not tr["veto"]
            ad=audit.advisory(tr,lane,mm[j],src,fs,pre);assert ad==x["advisory"] and not ad["veto"];out["parent_replay_wall"]+=time.perf_counter()-t
            h=list(ad["states"]);features=None
            if ad["advanced"]:
                t=time.perf_counter();pv=s.provenance(fs,pre,mm[j],src);assert pv is not None
                blocks,receipts=pv;slot=mm[j][0];oth=2-slot
                au=blocks[slot].copy();af=blocks[slot+1].copy();au[:,HALO:HALO+128]=rows[j][0];af[:,HALO:HALO+128]=refs[j]
                z=diagnose(au,af,blocks[oth],blocks[oth+1],np.asarray(path,int)+HALO,float(mm[j][3]),float(receipts[slot]["ts"]),bank);features=z["features"]
                e={**h[-1]["snapshot"]["evidence"],"differentiation":z,"parent_measurement_hashes":[s.sha(v.tobytes()) for v in (au,af,blocks[oth],blocks[oth+1])]}
                assert s.state(e,h,["physical_origin","independent_confirmation"]) and h[:-1]==ad["states"];audit.chain(h)
                out["diagnostic_wall"]+=time.perf_counter()-t
            out["records"].append({"row":j,"kind":kind,"lane":lane,"supported":ad["advanced"],"features":features,"label_evaluation_only":int(d["y"][j]),"localized_evaluation_only":bool(d["hit"][j])})
            advtrace.append({"row":j,"lane":lane,"source_tile_sha256":src,"states":h,"advanced":ad["advanced"],"decision_changed":False})
        assert np.array_equal(before,d["alert"]) and np.array_equal(hits,d["hit"])
        save("trace-"+num+"-"+kind+".json",{"records":advtrace,"baseline_mask":mask,"decision_changes":0})
    count=sum(x["supported"] for x in out["records"]);assert count==prior["summary"]["qualified"] and len(out["records"])==prior["summary"]["triggered"]
    out.update(triggered=len(out["records"]),qualified=count,diagnostic_advanced=count,decision_changes=0,full_wall=time.perf_counter()-BOOT,cpu=time.process_time()-CPU,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    save("group-"+num+".json",out);cp("group-"+num,{k:v for k,v in out.items() if k!="records"})
def stats(rows):
    res={}
    for k in FEATURES:
        v=[x["features"][k] for x in rows if x["features"] is not None and x["features"][k] is not None]
        res[k]={"n":len(v),"missing":len(rows)-len(v),"median":float(np.median(v)) if v else None,"q25":float(np.quantile(v,.25)) if v else None,"q75":float(np.quantile(v,.75)) if v else None}
    return res
def auc(rows,k):
    v=[(x["features"][k],1-x["label_evaluation_only"]) for x in rows if x["features"] is not None and x["features"][k] is not None]
    return float(roc_auc_score([x[1] for x in v],[x[0] for x in v])) if v and len({x[1] for x in v})==2 else None
def summary():
    gr=[read("group-"+str(i+1).zfill(2)+".json") for i in range(11)];assert all(x["parity"] and x["decision_changes"]==0 for x in gr)
    tables={}
    for split in ("fit","audit"):
        gg=[x for x in gr if x["fit_viewed"]==(split=="fit")];rr=[x for g in gg for x in g["records"]];ss=[x for x in rr if x["supported"]]
        tab={"triggered":len(rr),"qualified":len(ss),"all":stats(rr),"posthoc_control":stats([v for v in rr if v["label_evaluation_only"]==0]),"posthoc_positive":stats([v for v in rr if v["label_evaluation_only"]==1]),"by_lane":{k:stats([v for v in rr if v["lane"]==k]) for k in ("parent","rescue")}}
        tab["descriptive_auc_control_positive"]={k:{"pooled":auc(rr,k),"group_median":float(np.median(v)) if (v:=[val for g in gg if (val:=auc(g["records"],k)) is not None]) else None,"groups_with_both":len(v)} for k in FEATURES}
        tables[split]=tab
    assert tables["fit"]["qualified"]==482 and tables["audit"]["qualified"]==409 and sum(x["qualified"] for x in gr)==891
    cost={k:max(g[k] for g in gr) if "rss" in k else sum(g[k] for g in gr) for k in ("acquisition_wall","baseline_wall","parent_replay_wall","diagnostic_wall","full_wall","cpu","peak_rss_kib","payload_bytes","requests","files")}
    ss=read("synthetic_controls.json");assert all(ss["direction_checks"].values())
    d={"status":"PASS_EVIDENCE_DIFFERENTIATION_DIAGNOSTIC__HOLD_SCIENTIFIC_SCREEN","groups":11,"rows":3168,"qualified":891,"decision_changes":0,"baseline":{"localized":318,"weak":69,"controls":53},"tables":tables,"cost":cost,"synthetic_direction_checks":ss["direction_checks"],"full_bank_direction_check_diagnostic_only":ss["full_bank_direction_check_diagnostic_only"],"synthetic_sha256":s.sha((O/"synthetic_controls.json").read_bytes()),"source":L["source"],"trace_sha256":{p.name:s.sha(p.read_bytes()) for p in sorted(O.glob("trace-*.json"))},"production_adoption":False,"external_validation":False}
    save("result.json",d);cp("99-summary",{"status":d["status"],"groups":11,"qualified":891,"decision_changes":0,"baseline":d["baseline"],"audit_auc":tables["audit"]["descriptive_auc_control_positive"],"synthetic_direction_checks":ss["direction_checks"],"full_bank_direction_check_diagnostic_only":ss["full_bank_direction_check_diagnostic_only"],"cost":cost})
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="group":group(int(ARGS[1]))
    elif ARGS[0]=="summary":summary()
    else:raise ValueError(ARGS)
