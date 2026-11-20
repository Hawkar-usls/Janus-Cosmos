#!/usr/bin/env python3
"""R41S: label-blind local native-evidence spiral with fail-open, frozen veto."""
import sys,time
BOOT=time.perf_counter();BOOTCPU=time.process_time();ARGS=sys.argv[1:]
import json,pathlib,hashlib,io,os,zipfile,subprocess,resource,statistics,collections
import numpy as np
import requests
import r41b_minimal_hybrid as base
a=base.a;n=base.n;c=base.c
O=pathlib.Path("runtime/r41s");O.mkdir(parents=True,exist_ok=True);PROG=O/"progress";PROG.mkdir(exist_ok=True)
P=pathlib.Path("research/seti/r41s_5d_spiral_protocol.json");L=pathlib.Path("research/seti/r41s_numeric_operator_lock.json");S=json.loads(P.read_text());LOCK=json.loads(L.read_text());T=LOCK["thresholds"]
K=("eval","shared")
def sha(x):return hashlib.sha256(x).hexdigest()
def digest(d):return sha(json.dumps(d,sort_keys=True,separators=(",",":"),allow_nan=False).encode())
def save(name,d):a.write(O/name,d)
def read(name):return json.loads((O/name).read_text())
def cp(name,d):
    a.write(PROG/(name+".json"),d);print("R41S_PROGRESS",json.dumps(d,allow_nan=False),flush=True)
def archive(spec):
    url=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{spec['id']}/zip";headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    run=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/runs/{spec['run']}",headers=headers,timeout=60);run.raise_for_status();run=run.json()
    assert run["status"]=="completed" and run["conclusion"]=="success" and run["head_sha"]==spec["head_sha"]
    meta=requests.get(url[:-4],headers=headers,timeout=60);meta.raise_for_status();meta=meta.json()
    assert meta["workflow_run"]["id"]==spec["run"] and meta["digest"]=="sha256:"+spec["sha256"]
    r=requests.get(url,headers=headers,timeout=120);r.raise_for_status();assert len(r.content)<64*1024**2 and sha(r.content)==spec["sha256"]
    z=zipfile.ZipFile(io.BytesIO(r.content));entries=[x for x in z.infolist() if not x.is_dir()];names=[x.filename for x in entries]
    assert len(names)==len(set(names)) and sum(x.file_size for x in entries)<512*1024**2
    assert all(not pathlib.PurePosixPath(v).is_absolute() and ".." not in pathlib.PurePosixPath(v).parts for v in names)
    result=[v for v in names if pathlib.PurePosixPath(v).name=="result.json"];assert len(result)==1
    root=pathlib.PurePosixPath(result[0]).parent
    raw=z.read(result[0]);assert sha(raw)==spec["result_sha256"];d=json.loads(raw)
    assert d["status"]==spec["status"] and d["groups"]==11 and d["rows"]==3168 and not d["production_adoption"] and not d["external_validation"]
    return z,root,d
def makebank():
    p=[]
    for d in base.b.D:
        d=np.asarray(d,dtype=np.int64);lo=5-int(d.min());hi=123-int(d.max())
        if hi>lo:p.extend(np.arange(lo,hi)[:,None]+d)
    p=np.asarray(p,np.int64);assert p.shape[1]==16 and np.all((p>=5)&(p<123))
    return p
BANK=None
def norm(x):return n.ref.norm(np.asarray(x,float)[None,:,:])[0]
def contrast(z,p,offsets=(-5,-4,-3,3,4,5)):
    p=np.asarray(p,int);r=np.arange(16);off=np.asarray(offsets,int)
    assert np.all((p+off.min()>=0)&(p+off.max()<z.shape[1]))
    return np.clip(z[r,p]-np.median(z[r[:,None],p[:,None]+off],axis=1),-4.,4.)
def halves(x):return [float(np.mean(x[:8])),float(np.mean(x[8:]))]
def fit(z,p,bank):
    # Distinct hypothesis selection on one half only. Full-bank canonical ties.
    u=n.normalized_contrast(z[None,:,:])[0];inds=np.arange(len(bank));legal=(abs(bank-p)>T["distinct_bin_distance"]).sum(1)>=T["distinct_rows"];ids=inds[legal]
    if not len(ids):return None
    paths=bank[ids];vals=u[np.arange(16)[None,:],paths];pc=contrast(z,p)
    i=int(np.argmax(vals[:,:8].sum(1)));j=int(np.argmax(vals[:,8:].sum(1)));p1=paths[i];p2=paths[j]
    pred=[float(vals[j,:8].mean()),float(vals[i,8:].mean())];pc0=halves(pc)
    return {"forward_path":p1.tolist(),"backward_path":p2.tolist(),"crossfit":pred,"candidate":pc0,"margins":[pred[k]-pc0[k] for k in range(2)],"consensus":bool(np.max(abs(p1-p2))<=T["fit_consensus_max_bin_difference"]),"stationary_agreement":bool((p1==p1[0]).all() and (p2==p2[0]).all() and p1[0]==p2[0])}
def narrow(z,p):
    return {"near":halves(contrast(z,p)),"wide":halves(contrast(z,p,tuple(range(-12,-7))+tuple(range(8,13))))}
def ge(x,v):return all(a>=v for a in x)
def le(x,v):return all(a<=v for a in x)
def minus(x,y):return [x[i]-y[i] for i in range(2)]
def limits(lane):
    return (T["rescue_on_ceiling"],T["rescue_off_floor"],T["rescue_margin"]) if lane=="rescue" else (T["on_ceiling"],T["off_floor"],T["block_margin"])
def provenance(fs,pre,meta,src):
    # No use of injection labels or truth. Physical scan records + independently measured slices.
    slot,wi,f0,df=meta
    if slot not in (0,2) or len(fs)!=4 or [x["role"] for x in fs]!=["ON1","OFF1","ON2","OFF2"]:return None
    if len({x["url"] for x in fs})!=4 or len({x["name"] for x in fs})!=4:return None
    try:
        hh=[x["header"] for x in fs];ts=[float(x["ts"]) for x in hh];mjd=[float(x["mjd"]) for x in hh]
        if not (np.isfinite(ts).all() and np.isfinite(mjd).all() and min(ts)>0 and all(np.isclose(t,ts[0],rtol=T["ts_rtol"],atol=0) for t in ts)):return None
        if any((mjd[k+1]-mjd[k])*86400<16*ts[k] for k in range(3)):return None
        if any(x["shape"][0]!=16 or not np.isclose(float(x["df"])*1e6,df,rtol=T["df_rtol"],atol=0) for x in hh):return None
        aa=[];receipt=[];halo=T["halo_bins"]
        for k in range(4):
            at,X=pre[k][wi];vdf=float(at["df_hz"]);x=(float(f0)-float(at["channel0_mhz"]))*1e6/vdf;j=int(round(x));h=str(at["slice_sha256"])
            if not np.isclose(vdf,df,rtol=T["df_rtol"],atol=0) or abs(x-j)>T["grid_error_bins"] or X.shape[0]!=16 or j<halo or j+128+halo>X.shape[1]:return None
            if not len(h)==64:return None
            raw=X[:,j-halo:j+128+halo].astype(float).copy()
            if not np.isfinite(raw).all():return None
            aa.append(raw);receipt.append({"role":fs[k]["role"],"url":fs[k]["url"],"mjd":mjd[k],"ts":ts[k],"slice_sha256":h,"mapped_sha256":sha(raw.tobytes()),"start_bin":j,"df_hz":vdf})
        if receipt[slot]["slice_sha256"]!=src or len({v["slice_sha256"] for v in receipt})!=4:return None
        return aa,receipt
    except (KeyError,IndexError,TypeError,ValueError,ZeroDivisionError,OverflowError):return None
def state(e,history,unresolved):
    # Semantic content only, not clock or generation. Monotone evidence graph.
    v={"evidence":e,"unresolved":unresolved};h=digest(v)
    if history and h==history[-1]["state_sha256"]:return False
    assert h not in {z["state_sha256"] for z in history}
    history.append({"parent_state_sha256":history[-1]["state_sha256"] if history else None,"state_sha256":h,"snapshot":v})
    return True
def spiral(x,y,p,lane,meta,src,fs,pre,bank,repeat_override=None):
    h=[];e={};calls=0;reason="HOLD_NO_INFORMATION_GAIN";veto=False
    # The only inputs are spectra, legal path, lane, physical metadata and fixed bank.
    # No labels, true paths, injected amplitude, target name, or prior residual ids.
    src=str(src);p=np.asarray(p,int);origin={"candidate":p.tolist(),"lane":lane,"source_sha256":src,"on_sha256":sha(np.ascontiguousarray(x).tobytes()),"off_sha256":sha(np.ascontiguousarray(y).tobytes())}
    state({"origin":origin},h,["native_support","competitor","repeat"])
    def finish():
        assert calls<=T["max_calls"] and len(h)-1<=T["max_turns"]
        return {"veto":veto,"reason":reason,"calls":calls,"turns":len(h)-1,"states":h,"last_sha256":h[-1]["state_sha256"]}
    if p.shape!=(16,) or not np.all((p>=5)&(p<123)) or not np.isfinite(x).all() or not np.isfinite(y).all():
        reason="HOLD_UNSUPPORTED_NATIVE_PATH";return finish()
    z=norm(x);off=norm(y);ce,fl,ma=limits(lane)
    f=fit(z,p,bank);calls+=1
    if f is None:reason="HOLD_NO_RIVAL";return finish()
    e["D1"]=f
    on=halves(contrast(z,p));of=halves(contrast(off,p));gap=minus(of,on)
    e["D2"]={"on":on,"off":of,"off_minus_on":gap};calls+=1
    # Native wide-band evidence needs actual measured support. No artificial padding.
    station=bool((p==p[0]).all());safe=bool(np.all((p>=12)&(p<=115)))
    e["D3"]={"stationary":station,"full_neighborhood":safe,"off_native":narrow(off,p) if safe else None};calls+=1
    primary=bool(station and safe and f["consensus"] and f["stationary_agreement"] and ge(f["crossfit"],T["alternative_on_floor"]) and ge(f["margins"],ma) and le(on,ce) and ge(of,fl) and ge(gap,ma))
    # Two radial scales check actual narrow feature contrast; broad context may not manufacture it.
    if primary:
        nw=e["D3"]["off_native"];primary=ge(nw["near"],fl) and ge(nw["wide"],fl) and ge(minus(nw["wide"],on),T["wide_narrow_margin"])
    e["D4"]={"primary_conjunction":primary,"hypotheses":["candidate","coherent_rival","common_frequency_contaminant","unresolved"],"contradiction":not primary};calls+=1
    prov=provenance(fs,pre,meta,src);e["D5"]={"qualified":prov is not None,"scans":prov[1] if prov is not None else []};calls+=1
    state({"origin":origin,"pass1":e},h,["repeat","physical_origin"] if primary and prov else ["unresolved_candidate","physical_origin"])
    if not primary or prov is None:
        reason="HOLD_PRIMARY_CONFLICT_NOT_CORROBORATED" if not primary else "HOLD_PROVENANCE";return finish()
    blocks,receipt=prov;slot=meta[0];oth=2-slot;halo=T["halo_bins"]
    # Repeat evidence physically separate; stress augmentation is built by evaluation caller only.
    r_on=blocks[oth] if repeat_override is None else repeat_override[0]
    r_off=blocks[oth+1] if repeat_override is None else repeat_override[1]
    if not np.isfinite(r_on).all() or not np.isfinite(r_off).all():reason="HOLD_INVALID_REPEAT";return finish()
    Z=norm(r_on);D=norm(r_off);pp=p+halo;alt=np.asarray(f["forward_path"])+halo
    repeat_on=halves(contrast(Z,pp));repeat_alt=halves(contrast(Z,alt));repeat_off=halves(contrast(D,pp))
    rr={"D1":{"repeat_on":repeat_on,"repeat_held_path":alt.tolist(),"repeat_alt":repeat_alt,"alt_margin":minus(repeat_alt,repeat_on)}};calls+=1
    rr["D2"]={"repeat_off":repeat_off,"off_minus_on":minus(repeat_off,repeat_on)};calls+=1
    rr["D3"]={"repeat_off_native":narrow(D,pp),"repeat_on_native":narrow(Z,pp),"independent_halo":True};calls+=1
    # D4: evidence gate uses independent scan conjunction; no weak evidence is a veto.
    nn=rr["D3"]["repeat_off_native"];ok=le(repeat_on,ce) and ge(repeat_off,fl) and ge(minus(repeat_off,repeat_on),ma) and ge(repeat_alt,T["alternative_on_floor"]) and ge(minus(repeat_alt,repeat_on),ma) and ge(nn["near"],fl) and ge(nn["wide"],fl) and ge(minus(nn["wide"],repeat_on),T["wide_narrow_margin"])
    rr["D4"]={"independent_conjunction":bool(ok),"retained_alternatives":["physical_origin_not_established"],"veto_hypothesis":"common_receiver_frequency_contamination"};calls+=1
    rr["D5"]={"source_sha256":[receipt[oth]["slice_sha256"],receipt[oth+1]["slice_sha256"]],"measurement_hashes":[sha(np.ascontiguousarray(r_on).tobytes()),sha(np.ascontiguousarray(r_off).tobytes())],"evaluation_perturbation":repeat_override is not None};calls+=1
    new=state({"origin":origin,"pass1":e,"pass2":rr},h,["physical_origin"] if ok else ["unresolved_candidate","physical_origin"])
    veto=bool(new and ok);reason="EXPERIMENTAL_VETO" if veto else ("HOLD_NO_INFORMATION_GAIN" if not new else "HOLD_REPEAT_DISAGREEMENT")
    return finish()
def selftest(bank):
    # Fixed synthetic data; scientific threshold not selected on outcome.
    rng=np.random.default_rng(41231);x=rng.normal(0,0.04,(16,128));p=np.full(16,64);z=x.copy();z[:,40]+=10
    small=np.array([[20]*16,[40]*16,[64]*16],np.int64)
    f=fit(norm(z),p,small);assert f["stationary_agreement"] and f["forward_path"]==[40]*16
    # Full bank test: each reported rival matches an independent direct fit objective.
    u=n.normalized_contrast(norm(z)[None,:,:])[0];legal=bank[(abs(bank-p)>2).sum(1)>=5];v=u[np.arange(16)[None,:],legal];f0=fit(norm(z),p,bank)
    assert f0["forward_path"]==legal[int(v[:,:8].sum(1).argmax())].tolist() and f0["backward_path"]==legal[int(v[:,8:].sum(1).argmax())].tolist()
    hist=[];assert state({"x":1},hist,["a"]);assert not state({"x":1},hist,["a"]);assert state({"x":2},hist,["a"]);assert hist[1]["parent_state_sha256"]==hist[0]["state_sha256"]
    q=spiral(z,z,p,"parent",(0,0,3000.,2.8),"0"*64,[],[],bank);assert not q["veto"]
    q=spiral(z,z,np.full(16,0),"rescue",(0,0,3000.,2.8),"0"*64,[],[],bank);assert not q["veto"] and q["calls"]==0
    # Adversarial leakage and explicit source gate tests (missing/duplicate/misordered provenance).
    fs=[{"url":f"u{k}","name":f"f{k}","role":role,"header":{"ts":1.,"mjd":60000.+k*40/86400,"shape":[16,1,256],"df":2.8e-6}} for k,role in enumerate(("ON1","OFF1","ON2","OFF2"))]
    # Fabricate distinct slice sources and measured halo for a synthetic witness only.
    pp=[]
    for k in range(4):
        raw=rng.normal(0,.04,(16,256));raw[:,104 if k%2 else 80]+=10
        if k%2==0:raw[:,104]-=10
        pp.append([({"df_hz":2.8,"channel0_mhz":3000.-40*2.8/1e6,"slice_sha256":str(k)*64},raw)])
    assert provenance(fs,pp,(0,0,3000.,2.8),"0"*64) is not None
    assert provenance(fs,pp,(0,0,3000.,2.8),"f"*64) is None
    bad=[dict(q) for q in fs];bad[2]["url"]=bad[0]["url"];assert provenance(bad,pp,(0,0,3000.,2.8),"0"*64) is None
    on=pp[0][0][1][:,40:168];off=pp[1][0][1][:,40:168]
    rr=spiral(on,off,p,"parent",(0,0,3000.,2.8),"0"*64,fs,pp,small)
    assert rr["veto"] and rr["turns"]==2 and rr["calls"]==10
    # Corroboration removed: preserve the candidate despite the same primary spectrum.
    repeat_bad=(pp[2][0][1][:,:208].copy(),pp[3][0][1][:,:208].copy());repeat_bad[0][:,104]+=20
    assert not spiral(on,off,p,"parent",(0,0,3000.,2.8),"0"*64,fs,pp,small,repeat_bad)["veto"]
    broken=spiral(on,off,p,"parent",(0,0,3000.,2.8),"f"*64,fs,pp,small);assert not broken["veto"]
    assert not spiral(np.full_like(on,np.nan),off,p,"parent",(0,0,3000.,2.8),"0"*64,fs,pp,bank)["veto"]
    return {"bank_size":len(bank),"disjoint_crossfit":True,"content_hash_no_gain":True,"source_invalidation":True,"unsupported_fail_open":True,"synthetic_example_veto":rr["veto"]}
def prepare():
    bzip,root,B=archive(LOCK["parents"]["r41b"]);pzip,pr,P0=archive(LOCK["parents"]["r41p"])
    assert B["paired_processes"]==22 and B["source_pruning_lock_sha256"]==S["frozen_data"]["source_r41a"]["lock_sha256"]
    assert P0["status"]=="PASS_EXACT_BILLIARD_ORDERING__DEVELOPMENT_ONLY" and len(P0["exact"])==6
    for i in range(1,12):
        k=str(i).zfill(2);br=json.loads(bzip.read(str(root/("minimal-"+k+".json"))));rr=json.loads(pzip.read(str(pr/("group-"+k+".json"))))
        assert br["parity"] and rr["parity"]
        save("parent-"+k+".json",{"group":i,"raw_tiles":br["raw_tiles"],"metrics":br["metrics"],"target":br["target"],"fit_viewed":br["fit_viewed"]})
    e=B["metrics"]["minimal"]["eval"];sh=B["metrics"]["minimal"]["shared"]
    assert (e["localized"],sh["localized"],e["weak"],sh["weak"],e["controls"],e["parent_controls"],e["rescue_controls"],e["stress_cases"]+sh["stress_cases"],e["lost_joint"]+sh["lost_joint"],e["stress_loss"]+sh["stress_loss"])==(205,113,31,38,53,38,15,119,0,0)
    base.prepare();p0=base.plan();assert p0["pruning_lock_sha256"]==B["source_pruning_lock_sha256"]
    bank=makebank();tests=selftest(bank);np.save(O/"bank.npy",bank,allow_pickle=False)
    out={"stage":"PREPARE","status":"PASS_AUTHENTICATED_PARENTS_AND_OPERATOR_SELFTEST","protocol_sha256":sha(P.read_bytes()),"operator_lock_sha256":sha(L.read_bytes()),"code_sha256":sha(pathlib.Path(__file__).read_bytes()),"bank_sha256":sha((O/"bank.npy").read_bytes()),"parent_receipts":LOCK["parents"],"selftests":tests,"groups":11}
    save("prereq_receipt.json",out);cp("00-prepare",out)
def stress_override(row,kind,meta,src,fs,pre):
    # Evaluation-only data construction. Label-blind spiral receives perturbed arrays, never truth.
    pv=provenance(fs,pre,meta,src)
    if pv is None:return None
    blocks,_=pv;oth=2-meta[0];A,B=blocks[oth].copy(),blocks[oth+1].copy();halo=T["halo_bins"]
    active=np.asarray(row[8],bool);truth=np.asarray(row[7],int)+halo
    def inject(X):
        scale=1.4826*np.median(abs(X[:,halo:halo+128]-np.median(X[:,halo:halo+128],axis=1)[:,None]),axis=1)
        assert np.all(scale>0);X[np.arange(16)[active],truth[active]]+=row[4]*scale[active];return X
    return inject(A),inject(B) if kind=="shared" else B
def metrics(d,removed,oldarr,stress_removed):
    alert=d["alert"]&~removed;hit=d["hit"]&~removed;y=d["y"];weak=d["weak"];null=y==0
    jm=d["parent_hit"]|(d["joint_rescue"]&d["correct"]);ix=oldarr["idx"];ss=oldarr["repeat_stress_supported"]&oldarr["supported"]&oldarr["fold_supported"]&oldarr["halo_supported"]&oldarr["repeat_supported"]&oldarr["parent_hit"][ix]
    return {"rows":len(y),"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"controls":int((alert&null).sum()),"parent_controls":int((alert&null&d["parent_alert"]).sum()),"rescue_controls":int((alert&null&d["joint_rescue"]).sum()),"lost_vs_A":int((d["hit"]&~hit).sum()),"lost_joint":int((jm&~hit).sum()),"lost_r40g":int((d["parent_hit"]&~hit).sum()),"stress_cases":int(ss.sum()),"stress_loss":int((ss&stress_removed[ix]).sum()),"removed":int(removed.sum()),"removed_controls":int((removed&null).sum()),"removed_positives":int((removed&~null).sum())}
def worker(i,mode):
    rec=read("prereq_receipt.json");assert rec["protocol_sha256"]==sha(P.read_bytes()) and rec["operator_lock_sha256"]==sha(L.read_bytes()) and rec["code_sha256"]==sha(pathlib.Path(__file__).read_bytes()) and rec["bank_sha256"]==sha((O/"bank.npy").read_bytes())
    bank=np.load(O/"bank.npy",allow_pickle=False);assert mode in ("A","C");p0=base.plan();g=p0["groups"][i];num=str(i+1).zfill(2);ref=read("parent-"+num+".json")
    m,LL,old=base.locks("minimal",p0);seed={"n6":a.v.C["seed"],"n7":base.w.C["seed"],"n9":base.old9.C["seed"]}[g["source"]];n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);base.b.taylor.init_tree()
    fs=g["group"]["files"];net={"payload_bytes":0,"requests":0,"files":0};orig=base.eng.base.RangeFile
    class Counted(orig):
        def close(self):
            if not self.closed and hasattr(self,"total") and hasattr(self,"size"):
                net["payload_bytes"]+=int(self.total);net["requests"]+=int(self.requests_count);net["files"]+=1
            super().close()
    base.eng.base.RangeFile=Counted;t0=time.perf_counter();c0=time.process_time()
    try:rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r41s_"+mode+"_"+num)
    finally:base.eng.base.RangeFile=orig
    acq=time.perf_counter()-t0;acqcpu=time.process_time()-c0
    assert sorted(hs)==g["raw_tile_hashes"]==ref["raw_tiles"] and net["files"]==4
    pairs=[("eval",rr,R,meta),("shared",[x[0] for x in ch],np.stack([x[1] for x in ch]),cm)]
    def infer():
        dd={};led={};rm={}
        for kind,rows,refs,mm in pairs:
            d=base.infer("minimal",rows,refs,mm,pre,LL,old,m,None,kind);dd[kind]=d;remove=np.zeros(len(rows),bool);ledger=[]
            if mode=="C":
                for j in np.flatnonzero(d["alert"]):
                    lane="parent" if d["parent_alert"][j] else "rescue";path=d["parent_path"][j] if lane=="parent" else d["p"][j]
                    v=spiral(rows[j][0],refs[j],path,lane,mm[j],d["source_tile_sha"][j],fs,pre,bank)
                    remove[j]=v["veto"];ledger.append({"row":int(j),"lane":lane,"source_tile_sha256":str(d["source_tile_sha"][j]),"trace":v})
            led[kind]=ledger;rm[kind]=remove
        return dd,led,rm
    t0=time.perf_counter();c0=time.process_time();d0,led0,rm0=infer();cold=time.perf_counter()-t0;coldcpu=time.process_time()-c0;total=time.perf_counter()-BOOT;totalcpu=time.process_time()-BOOTCPU;peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    masks={k:sha(np.ascontiguousarray(d0[k]["alert"]).tobytes()+np.ascontiguousarray(d0[k]["hit"]).tobytes()) for k in K};trace={k:digest(led0[k]) for k in K};del d0,led0
    t0=time.perf_counter();c0=time.process_time();dd,led,rm=infer();warm=time.perf_counter()-t0;warmcpu=time.process_time()-c0
    assert all(np.array_equal(rm[k],rm0[k]) and digest(led[k])==trace[k] for k in K)
    perf={"cold_total_wall":total,"cold_total_cpu":totalcpu,"cold_compute_wall":cold,"cold_compute_cpu":coldcpu,"warm_compute_wall":warm,"warm_compute_cpu":warmcpu,"cold_peak_rss_kib":peak,"warm_peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"acquisition_wall":acq,"acquisition_cpu":acqcpu,**net}
    result={"group":i+1,"mode":mode,"target":ref["target"],"fit_viewed":g["fit"],"parity":True,"raw_tiles":sorted(hs),"metrics":{},"baseline_metrics":{},"baseline_masks":masks,"performance":perf,"diagnostics":{}}
    t0=time.perf_counter()
    for kind,rows,refs,mm in pairs:
        d=dd[kind];oldarr=a.read_npz(base.C/("dev_"+num+"_"+kind+".npz"));a.verify_rows(oldarr,rows);q=base.check(d,oldarr,m);assert q==ref["metrics"][kind]
        assert masks[kind]==sha(np.ascontiguousarray(d["alert"]).tobytes()+np.ascontiguousarray(d["hit"]).tobytes())
        stress_removed=np.zeros(len(rows),bool);stress_receipts=[]
        if mode=="C":
            ix=oldarr["idx"];st=oldarr["repeat_stress_supported"]&oldarr["supported"]&oldarr["fold_supported"]&oldarr["halo_supported"]&oldarr["repeat_supported"]&oldarr["parent_hit"][ix]
            for j in ix[st]:
                if not d["alert"][j]:continue
                path=d["parent_path"][j];src=d["source_tile_sha"][j];override=stress_override(rows[j],kind,mm[j],src,fs,pre)
                if override is None:v={"veto":False,"reason":"HOLD_STRESS_PROVENANCE"}
                else:v=spiral(rows[j][0],refs[j],path,"parent",mm[j],src,fs,pre,bank,override)
                stress_removed[j]=v["veto"];stress_receipts.append({"row":int(j),"trace":v})
        result["metrics"][kind]=metrics(d,rm[kind],oldarr,stress_removed);result["baseline_metrics"][kind]=metrics(d,np.zeros(len(rows),bool),oldarr,np.zeros(len(rows),bool))
        if mode=="C":
            assert [v["row"] for v in led[kind]]==np.flatnonzero(d["alert"]).tolist()
            result["diagnostics"][kind]={"triggered":len(led[kind]),"turns":sum(x["trace"]["turns"] for x in led[kind]),"calls":sum(x["trace"]["calls"] for x in led[kind]),"two_turns":sum(x["trace"]["turns"]==2 for x in led[kind]),"qualified_sources":sum(any(s["snapshot"]["evidence"].get("pass1",{}).get("D5",{}).get("qualified",False) for s in x["trace"]["states"]) for x in led[kind]),"reasons":dict(collections.Counter(x["trace"]["reason"] for x in led[kind])),"stress_replayed":len(stress_receipts)}
            save("trace-"+num+"-"+kind+".json",{"normal":led[kind],"stress":stress_receipts,"decision_mask_sha256":sha(rm[kind].tobytes())})
    result["evaluation_stress_and_validation_seconds"]=time.perf_counter()-t0;save(mode+"-"+num+".json",result)
def pair(i):
    num=str(i+1).zfill(2);modes=("A","C") if i%2==0 else ("C","A")
    for mode in modes:subprocess.run([sys.executable,__file__,"worker",str(i),mode],check=True)
    A=read("A-"+num+".json");C=read("C-"+num+".json");assert A["raw_tiles"]==C["raw_tiles"] and A["baseline_masks"]==C["baseline_masks"] and A["metrics"]==C["baseline_metrics"]
    cp("group-"+num,{"stage":"GROUP","group":i+1,"target":A["target"],"fit_viewed":A["fit_viewed"],"parity":True,"order":modes,"A":A["metrics"],"C0":C["baseline_metrics"],"C":C["metrics"],"diagnostics":C["diagnostics"],"performance":{"A":A["performance"],"C":C["performance"]}})
def summary():
    allr={m:[read(m+"-"+str(i+1).zfill(2)+".json") for i in range(11)] for m in ("A","C")};agg={};cost={}
    for mode,rows in allr.items():
        audit=[v for v in rows if not v["fit_viewed"]];assert len(audit)==5
        agg[mode]={k:{f:sum(v["metrics"][k][f] for v in audit) for f in audit[0]["metrics"][k]} for k in K}
        cost[mode]={f:max(v["performance"][f] for v in rows) if "rss" in f else sum(v["performance"][f] for v in rows) for f in rows[0]["performance"]}
    agg["C0"]=agg["A"];A=agg["A"];C=agg["C"];e=C["eval"];sh=C["shared"]
    assert (A["eval"]["localized"],A["shared"]["localized"],A["eval"]["controls"],A["eval"]["weak"],A["shared"]["weak"])==(205,113,53,31,38)
    assert A["eval"]["stress_cases"]+A["shared"]["stress_cases"]==119
    qualifies=(e["localized"],sh["localized"],e["weak"],sh["weak"])==(205,113,31,38) and e["controls"]<=35 and all(e[k]+sh[k]==0 for k in ("lost_vs_A","lost_joint","lost_r40g","stress_loss"))
    dg={k:{"reasons":dict(sum((collections.Counter(v["diagnostics"][k]["reasons"]) for v in allr["C"] if not v["fit_viewed"]),collections.Counter())),**{f:sum(v["diagnostics"][k][f] for v in allr["C"] if not v["fit_viewed"]) for f in ("triggered","turns","calls","two_turns","qualified_sources","stress_replayed")}} for k in K}
    seal={"protocol_sha256":sha(P.read_bytes()),"operator_lock_sha256":sha(L.read_bytes()),"code_sha256":sha(pathlib.Path(__file__).read_bytes()),"group_results_sha256":{v.name:sha(v.read_bytes()) for v in sorted(O.glob("C-*.json"))},"trace_sha256":{v.name:sha(v.read_bytes()) for v in sorted(O.glob("trace-*.json"))},"gate_pass":bool(qualifies)}
    save("C_sealed_lock.json",seal)
    result={"status":"PASS_HISTORICAL_SCREEN__HOLD_FRESH_EXTERNAL" if qualifies else "PASS_OPERATOR_PARITY__HOLD_HISTORICAL_SCREEN","groups":11,"audit_groups":5,"rows":3168,"paired_processes":22,"metrics":agg,"cost":cost,"diagnostics":dg,"historical_screen":bool(qualifies),"C_sealed_lock_sha256":sha((O/"C_sealed_lock.json").read_bytes()),"prereq_receipt_sha256":sha((O/"prereq_receipt.json").read_bytes()),"C0_cost_scope":"C0/C share actual C operator execution; not separate independent run","fresh_external_validation":False,"production_adoption":False,"penrose_executed":False}
    save("result.json",result);cp("99-summary",result)
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="worker":worker(int(ARGS[1]),ARGS[2])
    elif ARGS[0]=="pair":pair(int(ARGS[1]))
    elif ARGS[0]=="summary":summary()
    else:raise ValueError(ARGS)
