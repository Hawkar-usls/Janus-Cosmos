#!/usr/bin/env python3
"""Retina-Rubik registered pixel layer overlap, advisory only."""
import sys,time
ARGS=sys.argv[1:];BOOT=time.perf_counter();CPU=time.process_time()
import json,pathlib,os,io,zipfile,resource
import numpy as np
import requests
from sklearn.metrics import roc_auc_score
import r41s1c_evidence_differentiation as parent
audit=parent.audit;s=parent.s;base=s.base;a=s.a;n=s.n;HALO=parent.HALO
O=pathlib.Path("runtime/r41s1d");O.mkdir(parents=True,exist_ok=True)
PC=O/"sealed_parent";PC.mkdir(exist_ok=True);C=parent.C;PROG=O/"progress";PROG.mkdir(exist_ok=True)
P=pathlib.Path("research/seti/r41s1d_retina_rubik_protocol.json");L=json.loads(P.read_text())
FACES=("primary_on","primary_off","repeat_on","repeat_off")
KEYS=("mean_core","median_core","min_core","positive_overlap","overlap_core_minus_side","layer_dispersion","leave_one_min","held_layer_min","concentration","native_persistence","null_overlap_excess","null_min_excess","null_overlap_q90","null_min_q90")
FEATURES=tuple(k+"_"+v for k in FACES for v in KEYS)+("primary_on_off_delta","repeat_on_off_delta","primary_shared_min","repeat_shared_min","on_repeat_min","off_repeat_min","on_overlap_min","off_overlap_min")
SHIFTS=np.random.default_rng(413201).choice(np.r_[np.arange(-18,-5),np.arange(6,19)],size=(32,4,4),replace=True)
def save(name,d):a.write(O/name,d)
def read(name):return json.loads((O/name).read_text())
def cp(name,d):
    a.write(PROG/(name+".json"),d);print("R41S1D_PROGRESS",json.dumps(d,allow_nan=False),flush=True)
def mean(v):return float(np.mean(v))
def geometry(cube,p):
    assert cube.shape==(4,16,208) and p.shape==(16,) and np.all((p>=24)&(p<=183)) and np.isfinite(cube).all()
    r=np.rot90(cube,axes=(1,2));assert np.array_equal(np.rot90(r,k=-1,axes=(1,2)),cube)
    t=cube.transpose(2,0,1);assert np.array_equal(t.transpose(1,2,0),cube)
    for i in range(4):assert np.array_equal(r[i,207-p,np.arange(16)],cube[i,np.arange(16),p])
    normalized=np.stack([s.norm(x) for x in cube])
    return normalized[:,np.arange(16)[:,None],p[:,None]+np.arange(-24,25)[None,:]]
SIDE=np.r_[np.arange(12,17),np.arange(32,37)]
def contrast(x):
    return np.clip(x[...,24]-np.median(x[...,SIDE],axis=-1),-4.,4.)
def overlay(layer):
    # layer shape [..., four temporal blocks, four phase rows, 49 native channels].
    v=layer-np.median(layer[...,SIDE],axis=-1,keepdims=True)
    occ=(v>0).sum(axis=-3)>=3
    meanmap=v.mean(axis=-3);medmap=np.median(v,axis=-3);minmap=v.min(axis=-3)
    ov=occ[...,24].mean(axis=-1)
    mn=contrast(minmap).mean(axis=-1)
    return v,occ,meanmap,medmap,minmap,ov,mn
def face(x,i):
    layer=x.reshape(4,4,49);v,occ,avg,med,mn,ov,mc=overlay(layer)
    idx=(np.arange(49)[None,None,None,:]-SHIFTS[:,i,:,None,None])%49
    shifted=np.take_along_axis(np.broadcast_to(layer,(32,4,4,49)),np.broadcast_to(idx,(32,4,4,49)),axis=-1)
    inv=(np.arange(49)[None,None,None,:]+SHIFTS[:,i,:,None,None])%49
    recovered=np.take_along_axis(shifted,np.broadcast_to(inv,shifted.shape),axis=-1)
    assert np.array_equal(recovered,np.broadcast_to(layer,shifted.shape))
    _,_,_,_,_,no,nm=overlay(shifted);oq=float(np.quantile(no,.9));mq=float(np.quantile(nm,.9))
    layer_core=contrast(layer).mean(1)
    leave=[mean(contrast(np.delete(v,j,axis=0).mean(0))) for j in range(4)]
    pos=np.maximum(avg,0);den=pos[:,12:37].sum()
    out={"mean_core":mean(contrast(avg)),"median_core":mean(contrast(med)),"min_core":float(mc),"positive_overlap":float(ov),"overlap_core_minus_side":float(ov-occ[:,SIDE].mean()),"layer_dispersion":mean(np.std(v[:,:,24],axis=0)),"leave_one_min":min(leave),"held_layer_min":float(layer_core.min()),"concentration":float(pos[:,23:26].sum()/den) if den>0 else None,"native_persistence":mean(contrast(x)>0),"null_overlap_excess":float(ov-oq),"null_min_excess":float(mc-mq),"null_overlap_q90":oq,"null_min_q90":mq}
    return out,{"mean_plane":avg.tolist(),"median_plane":med.tolist(),"minimum_plane":mn.tolist(),"positive_overlap_plane":occ.astype(int).tolist(),"core_by_layer":layer_core.tolist(),"leave_one_layer":leave,"null_overlap":no.tolist(),"null_minimum":nm.tolist(),"strip_sha256":s.sha(np.ascontiguousarray(x).tobytes())}
def diagnose(cube,p):
    x=geometry(cube,p);feat={};detail={}
    for i,k in enumerate(FACES):
        v,d=face(x[i],i);detail[k]=d;feat.update({k+"_"+j:v[j] for j in KEYS})
    m=[feat[k+"_mean_core"] for k in FACES];q=[feat[k+"_positive_overlap"] for k in FACES]
    feat.update(primary_on_off_delta=m[0]-m[1],repeat_on_off_delta=m[2]-m[3],primary_shared_min=min(m[:2]),repeat_shared_min=min(m[2:]),on_repeat_min=min(m[0],m[2]),off_repeat_min=min(m[1],m[3]),on_overlap_min=min(q[0],q[2]),off_overlap_min=min(q[1],q[3]))
    assert set(feat)==set(FEATURES) and all(v is None or np.isfinite(v) for v in feat.values())
    return {"features":feat,"faces":detail,"shape":list(cube.shape),"rotation_inverse_exact":True,"axis_inverse_exact":True,"slab_inverse_exact":True,"permutations":32,"scope":"registered overlapping data layers; conditioned slab-shift contrast not p-value; not physical identity"}
def synthetic():
    rng=np.random.default_rng(413202);rows=[];r=np.arange(16);p=80+r//2
    def add(c,p,amp,width,burst=False,faces=(0,)):
        c=c.copy()
        for k in faces:
            x=c[k];scale=1.4826*np.median(abs(x-np.median(x,axis=1)[:,None]),axis=1)
            mask=(r>=4)&(r<12) if burst else np.ones(16,bool)
            for d in range(-(width//2),width//2+1):x[r[mask],p[mask]+d]+=amp*scale[mask]
        return c
    for i in range(32):
        c=rng.normal(size=(4,16,208));strong=add(c,p,6,1)
        weak=add(c,p,1.5,1);broad=add(c,p,6,17,True)
        fixed=np.full(16,84);shared=add(c,fixed,6,1,faces=(0,1,2,3))
        scramble=strong.copy()
        for j,k in enumerate((11,-13,16,-9)):
            scramble[0,j*4:j*4+4]=np.roll(strong[0,j*4:j*4+4],k,axis=1)
        assert np.array_equal(np.sort(scramble[0],axis=1),np.sort(strong[0],axis=1))
        z=s.norm(c[0]);bank=p[None,:]+np.arange(-30,31)[:,None]
        sel=bank[np.argmax(parent.native_contrast(z)[r[None,:],bank].mean(1))]
        for name,v,path in (("noise",c,sel),("coherent",strong,p),("weak",weak,p),("scrambled",scramble,p),("broad",broad,p),("shared",shared,fixed)):
            rows.append({"pair":i,"class":name,"features":diagnose(v,path)["features"]})
    classes=sorted({x["class"] for x in rows})
    med={k:{j:float(np.median(v)) if (v:=[x["features"][j] for x in rows if x["class"]==k and x["features"][j] is not None]) else None for j in FEATURES} for k in classes}
    checks={"coherent_overlap_gt_scrambled":med["coherent"]["primary_on_positive_overlap"]>med["scrambled"]["primary_on_positive_overlap"],"coherent_null_excess_gt_scrambled":med["coherent"]["primary_on_null_min_excess"]>med["scrambled"]["primary_on_null_min_excess"],"narrow_concentration_gt_broad":med["coherent"]["primary_on_concentration"]>med["broad"]["primary_on_concentration"],"shared_off_gt_target_only":med["shared"]["primary_off_mean_core"]>med["coherent"]["primary_off_mean_core"]}
    save("synthetic_controls.json",{"scope":"simulated; no actual origin evidence","seed":413202,"rows":rows,"medians":med,"direction_checks":checks})
    return checks
def prepare():
    spec=L["parent"];head={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"};root="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions"
    def get(url):
        rr=requests.get(url,headers=head,timeout=120);rr.raise_for_status();return rr
    run=get(f"{root}/runs/{spec['run']}").json();assert run["status"]=="completed" and run["conclusion"]=="success" and run["head_sha"]==spec["head_sha"]
    md=get(f"{root}/artifacts/{spec['artifact']}").json();assert md["workflow_run"]["id"]==spec["run"] and md["digest"]=="sha256:"+spec["archive_sha256"]
    raw=get(f"{root}/artifacts/{spec['artifact']}/zip").content;assert s.sha(raw)==spec["archive_sha256"] and len(raw)<128*1024**2
    z=zipfile.ZipFile(io.BytesIO(raw));names=[x.filename for x in z.infolist() if not x.is_dir()];assert len(names)==len(set(names)) and sum(x.file_size for x in z.infolist())<512*1024**2 and all(not pathlib.PurePosixPath(x).is_absolute() and ".." not in pathlib.PurePosixPath(x).parts for x in names)
    result=[k for k in names if pathlib.PurePosixPath(k).name=="result.json"];assert len(result)==1
    root=pathlib.PurePosixPath(result[0]).parent
    assert s.sha(z.read(result[0]))==spec["result_sha256"];d=json.loads(z.read(result[0]))
    assert d["groups"]==11 and d["qualified"]==891 and d["decision_changes"]==0 and len(d["trace_sha256"])==22
    for k,h in d["trace_sha256"].items():
        assert s.sha(z.read(str(root/k)))==h;(PC/k).write_bytes(z.read(str(root/k)))
    parent.prepare();checks=synthetic()
    receipt={"stage":"PREPARE","archive_traces":22,"operator_sha256":s.sha(pathlib.Path(__file__).read_bytes()),"protocol_sha256":s.sha(P.read_bytes()),"bank_sha256":s.sha((parent.O/"bank.npy").read_bytes()),"synthetic_checks":checks,"features":len(FEATURES)}
    save("prereq_receipt.json",receipt);cp("00-prepare",receipt)
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
    try:rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r41s1d_"+num)
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
        advtrace=[]
        for x,y0,previous in zip(records,original,sealed["records"]):
            j=x["row"];lane="parent" if d["parent_alert"][j] else "rescue";path=d["parent_path"][j] if lane=="parent" else d["p"][j];src=str(d["source_tile_sha"][j]);assert lane==x["lane"]==y0["lane"] and src==x["source_tile_sha256"]==y0["source_tile_sha256"]
            t=time.perf_counter();tr=s.spiral(rows[j][0],refs[j],path,lane,mm[j],src,fs,pre,oldbank);assert tr==y0["trace"] and not tr["veto"]
            ad=audit.advisory(tr,lane,mm[j],src,fs,pre);assert ad==x["advisory"] and not ad["veto"];out["parent_replay_wall"]+=time.perf_counter()-t
            h=list(ad["states"]);features=None
            if ad["advanced"]:
                t=time.perf_counter();pv=s.provenance(fs,pre,mm[j],src);assert pv is not None
                blocks,receipts=pv;slot=mm[j][0];oth=2-slot
                au=blocks[slot].copy();af=blocks[slot+1].copy();au[:,HALO:HALO+128]=rows[j][0];af[:,HALO:HALO+128]=refs[j]
                z=parent.diagnose(au,af,blocks[oth],blocks[oth+1],np.asarray(path,int)+HALO,float(mm[j][3]),float(receipts[slot]["ts"]),bank)
                e={**h[-1]["snapshot"]["evidence"],"differentiation":z,"parent_measurement_hashes":[s.sha(v.tobytes()) for v in (au,af,blocks[oth],blocks[oth+1])]}
                assert s.state(e,h,["physical_origin","independent_confirmation"]) and h[:-1]==ad["states"];audit.chain(h)
                assert h==previous["states"] and previous["advanced"] and previous["source_tile_sha256"]==src and not previous["decision_changed"],"R41S1C sealed replay"
                out["parent_replay_wall"]+=time.perf_counter()-t
                t=time.perf_counter();z2=diagnose(np.stack([au,af,blocks[oth],blocks[oth+1]]),np.asarray(path,int)+HALO);features=z2["features"]
                e2={**h[-1]["snapshot"]["evidence"],"retina_rubik":z2}
                assert s.state(e2,h,["physical_origin","independent_confirmation"]) and h[:-1]==previous["states"];audit.chain(h)
                out["diagnostic_wall"]+=time.perf_counter()-t
            else:
                assert h==previous["states"] and not previous["advanced"] and not previous["decision_changed"]
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
    gr=[read("group-"+str(i+1).zfill(2)+".json") for i in range(11)]
    assert all(x["parity"] and x["decision_changes"]==0 for x in gr)
    tables={}
    for split in ("fit","audit"):
        gg=[x for x in gr if x["fit_viewed"]==(split=="fit")];rr=[x for g in gg for x in g["records"]];ss=[x for x in rr if x["supported"]]
        tab={"triggered":len(rr),"qualified":len(ss),"all":stats(rr),"control":stats([v for v in rr if v["label_evaluation_only"]==0]),"positive":stats([v for v in rr if v["label_evaluation_only"]==1]),"by_lane":{k:stats([v for v in rr if v["lane"]==k]) for k in ("parent","rescue")}}
        tab["descriptive_auc"]={k:{"pooled":auc(rr,k),"group_median":float(np.median(v)) if (v:=[val for g in gg if (val:=auc(g["records"],k)) is not None]) else None,"groups":len(v)} for k in FEATURES}
        tables[split]=tab
    assert tables["fit"]["qualified"]==482 and tables["audit"]["qualified"]==409
    cost={k:max(g[k] for g in gr) if "rss" in k else sum(g[k] for g in gr) for k in ("acquisition_wall","baseline_wall","parent_replay_wall","diagnostic_wall","full_wall","cpu","peak_rss_kib","payload_bytes","requests","files")}
    ss=read("synthetic_controls.json")
    d={"status":"PASS_RETINA_RUBIK_LAYER_DIAGNOSTIC__HOLD_SCIENTIFIC_SCREEN","groups":11,"qualified":891,"decision_changes":0,"features":len(FEATURES),"baseline":{"localized":318,"weak":69,"controls":53},"synthetic_direction_checks":ss["direction_checks"],"tables":tables,"cost":cost,"source":L["parent"],"trace_sha256":{p.name:s.sha(p.read_bytes()) for p in sorted(O.glob("trace-*.json"))},"production_adoption":False,"external_validation":False}
    save("result.json",d);cp("99-summary",{k:v for k,v in d.items() if k not in ("tables","trace_sha256")}|{"audit_auc":tables["audit"]["descriptive_auc"]})
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="group":group(int(ARGS[1]))
    elif ARGS[0]=="summary":summary()
    else:raise ValueError(ARGS)
