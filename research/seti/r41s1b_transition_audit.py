#!/usr/bin/env python3
"""R41S1B: replay frozen first turns, separate advisory evidence acquisition from veto."""
import sys,time
ARGS=sys.argv[1:];BOOT=time.perf_counter();CPU=time.process_time()
import json,pathlib,os,io,zipfile,collections,itertools,resource,hashlib
import numpy as np
import requests
import r41s_5d_spiral as s
base=s.base;a=s.a;n=s.n
O=pathlib.Path("runtime/r41s1b");O.mkdir(parents=True,exist_ok=True);C=O/"parent";C.mkdir(exist_ok=True);PROG=O/"progress";PROG.mkdir(exist_ok=True)
P=pathlib.Path("research/seti/r41s1b_transition_audit_protocol.json");L=json.loads(P.read_text());T=s.T
K=("eval","shared")
NAMES=("stationary","native_wide","rival_consensus","rival_stationary_agreement","crossfit_floor","crossfit_margin","on_ceiling","off_floor","off_minus_on","near_floor","wide_floor","wide_margin","source")
def save(name,d):a.write(O/name,d)
def read(name):return json.loads((O/name).read_text())
def cp(name,d):
    a.write(PROG/(name+".json"),d);print("R41S1B_PROGRESS",json.dumps(d,allow_nan=False),flush=True)
def tests(trace,lane):
    h=trace["states"];p=h[-1]["snapshot"]["evidence"].get("pass1")
    if p is None:return None
    f=p["D1"];d=p["D2"];nw=p["D3"]["off_native"];ce,fl,ma=s.limits(lane)
    t={"stationary":bool(p["D3"]["stationary"]),"native_wide":bool(p["D3"]["full_neighborhood"]),"rival_consensus":bool(f["consensus"]),"rival_stationary_agreement":bool(f["stationary_agreement"]),"crossfit_floor":s.ge(f["crossfit"],T["alternative_on_floor"]),"crossfit_margin":s.ge(f["margins"],ma),"on_ceiling":s.le(d["on"],ce),"off_floor":s.ge(d["off"],fl),"off_minus_on":s.ge(d["off_minus_on"],ma),"near_floor":None if nw is None else s.ge(nw["near"],fl),"wide_floor":None if nw is None else s.ge(nw["wide"],fl),"wide_margin":None if nw is None else s.ge(s.minus(nw["wide"],d["on"]),T["wide_narrow_margin"]),"source":bool(p["D5"]["qualified"])}
    assert all(t[k] is True for k in NAMES[:-1])==p["D4"]["primary_conjunction"]
    return t
def chain(h):
    seen=set();prev=None
    for q in h:
        assert q["parent_state_sha256"]==prev and q["state_sha256"]==s.digest(q["snapshot"]) and q["state_sha256"] not in seen
        seen.add(q["state_sha256"]);prev=q["state_sha256"]
def advisory(trace,lane,meta,src,fs,pre):
    # No new label access, no selection on repeat data, original decision never changed.
    h=list(trace["states"]);chain(h);old=h[-1]["snapshot"]["evidence"];p1=old.get("pass1")
    if p1 is None or not p1["D5"]["qualified"]:return {"status":"HOLD_UNQUALIFIED","advanced":False,"veto":False,"states":h,"calls":0}
    prov=s.provenance(fs,pre,meta,src);assert prov is not None and prov[1]==p1["D5"]["scans"]
    blocks,receipts=prov;oth=2-meta[0];Z=s.norm(blocks[oth]);D=s.norm(blocks[oth+1]);halo=T["halo_bins"]
    paths={"candidate":old["origin"]["candidate"],"forward_rival":p1["D1"]["forward_path"],"backward_rival":p1["D1"]["backward_path"]}
    out={};ce,fl,ma=s.limits(lane)
    for k,path in paths.items():
        p=np.asarray(path,int)+halo
        assert p.shape==(16,) and (p>=12).all() and (p+12<Z.shape[1]).all()
        on=s.narrow(Z,p);off=s.narrow(D,p);out[k]={"path":path,"stationary":bool((p==p[0]).all()),"on":on,"off":off,"off_minus_on":s.minus(off["near"],on["near"]),"on_block_delta":on["near"][1]-on["near"][0],"off_block_delta":off["near"][1]-off["near"][0]}
    ca=out["candidate"];alt={k:s.minus(out[k]["on"]["near"],ca["on"]["near"]) for k in ("forward_rival","backward_rival")}
    ev={"D1":{"fixed_paths":{k:v["path"] for k,v in out.items()},"candidate_primary_blocks":p1["D2"]["on"],"candidate_repeat_blocks":ca["on"]["near"],"alternate_repeat_margins":alt},"D2":{"on_off":{k:{"on":v["on"]["near"],"off":v["off"]["near"],"off_minus_on":v["off_minus_on"]} for k,v in out.items()}},"D3":{"native_wide":{k:{"on":v["on"]["wide"],"off":v["off"]["wide"]} for k,v in out.items()},"all_within_actual_halo":True},"D4":{"hypotheses_retained":["candidate","coherent_forward_rival","coherent_backward_rival","common_receiver_frequency_contaminant","unresolved_physical_origin"],"candidate_repeat_off_floor":s.ge(ca["off"]["near"],fl),"candidate_repeat_on_ceiling":s.le(ca["on"]["near"],ce),"candidate_repeat_off_margin":s.ge(ca["off_minus_on"],ma),"both_rivals_repeat_margin":all(s.ge(x,ma) for x in alt.values()),"candidate_stationary":ca["stationary"],"receiver_frame_shape_only":True},"D5":{"repeat_scans":[receipts[oth],receipts[oth+1]],"independent_scan_pair":True,"measurement_hashes":[s.sha(np.ascontiguousarray(blocks[j]).tobytes()) for j in (oth,oth+1)],"scope":"same native receiver channels; drifting shape replay is not a prediction of physical frequency evolution"}}
    # Keep everything (including failed primary conjunction) in the append-only state.
    e={**old,"advisory_repeat":ev};ok=s.state(e,h,["physical_origin","hypothesis_resolution"])
    assert ok and len(h)==len(trace["states"])+1 and h[:-1]==trace["states"];chain(h)
    return {"status":"ADVISORY_SECOND_TURN","advanced":True,"veto":False,"states":h,"calls":5}
def prepare():
    q=L["parent"];hh={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"};root="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions"
    def get(url):
        r=requests.get(url,headers=hh,timeout=120);r.raise_for_status();return r
    run=get(f"{root}/runs/{q['run']}").json();assert run["status"]=="completed" and run["conclusion"]=="success" and run["head_sha"]==q["head_sha"]
    md=get(f"{root}/artifacts/{q['id']}").json();assert md["workflow_run"]["id"]==q["run"] and md["digest"]=="sha256:"+q["sha256"]
    raw=get(f"{root}/artifacts/{q['id']}/zip").content;assert s.sha(raw)==q["sha256"] and len(raw)<128*1024**2
    z=zipfile.ZipFile(io.BytesIO(raw));names=[x.filename for x in z.infolist() if not x.is_dir()];assert len(names)==len(set(names)) and sum(x.file_size for x in z.infolist())<512*1024**2 and all(not pathlib.PurePosixPath(x).is_absolute() and ".." not in pathlib.PurePosixPath(x).parts for x in names)
    name=[x for x in names if pathlib.PurePosixPath(x).name=="result.json"];assert len(name)==1;root=pathlib.PurePosixPath(name[0]).parent
    def rd(k):return z.read(str(root/k))
    assert s.sha(rd("result.json"))==q["result_sha256"] and s.sha(rd("C_sealed_lock.json"))==q["seal_sha256"]
    d=json.loads(rd("result.json"));seal=json.loads(rd("C_sealed_lock.json"));assert d["status"]=="PASS_OPERATOR_PARITY__HOLD_HISTORICAL_SCREEN" and d["groups"]==11 and d["paired_processes"]==22 and d["metrics"]["A"]==d["metrics"]["C"]==d["metrics"]["C0"] and not d["historical_screen"]
    assert seal["code_sha256"]==L["frozen_original"]["operator_sha256"]==s.sha(pathlib.Path(s.__file__).read_bytes()) and seal["operator_lock_sha256"]==L["frozen_original"]["numeric_lock_sha256"]==s.sha(s.L.read_bytes())
    for k,h in {**seal["group_results_sha256"],**seal["trace_sha256"]}.items():assert s.sha(rd(k))==h;(C/k).write_bytes(rd(k))
    for i in range(1,12):
        k=str(i).zfill(2)
        for f in ("A-"+k+".json","parent-"+k+".json"):(C/f).write_bytes(rd(f))
    base.prepare();bank=s.makebank();np.save(O/"bank.npy",bank,allow_pickle=False)
    st=s.selftest(bank);assert st["synthetic_example_veto"]
    # Existing synthetic positive and negative gates; advisory source invalidation handled by primary replay.
    save("prereq_receipt.json",{"stage":"PREPARE","parent":q,"protocol_sha256":s.sha(P.read_bytes()),"code_sha256":s.sha(pathlib.Path(__file__).read_bytes()),"bank_sha256":s.sha((O/"bank.npy").read_bytes()),"selftest":st,"sealed_original":True})
    cp("00-prepare",read("prereq_receipt.json"))
def summarize(records):
    c={"triggered":len(records),"pass1":sum(x["tests"] is not None for x in records),"advisory_advanced":sum(x["advisory_advanced"] for x in records),"qualified":sum(x["tests"] is not None and x["tests"]["source"] for x in records),"additional_calls":sum(x["additional_calls"] for x in records),"decision_changes":0,"original_open":sum(x["original_open"] for x in records)}
    c["gates"]={k:{"pass":sum(x["tests"] is not None and x["tests"][k] is True for x in records),"fail":sum(x["tests"] is not None and x["tests"][k] is False for x in records),"na":sum(x["tests"] is None or x["tests"][k] is None for x in records),"leave_one_out_pass":sum(x["tests"] is not None and all(x["tests"][j] is True for j in NAMES if j!=k) for x in records),"cumulative_pass":sum(x["tests"] is not None and all(x["tests"][j] is True for j in NAMES[:i+1]) for x in records)} for i,k in enumerate(NAMES)}
    c["cofail"]={k+"|"+j:sum(x["tests"] is not None and x["tests"][k] is False and x["tests"][j] is False for x in records) for k,j in itertools.combinations(NAMES,2)}
    c["repeat_diagnostics"]={k:sum(bool(x.get("repeat",{}).get(k,False)) for x in records) for k in ("candidate_repeat_off_floor","candidate_repeat_on_ceiling","candidate_repeat_off_margin","both_rivals_repeat_margin","candidate_stationary")}
    assert c["advisory_advanced"]==c["qualified"] and c["decision_changes"]==0
    return c
def group(i):
    p=read("prereq_receipt.json");assert p["protocol_sha256"]==s.sha(P.read_bytes()) and p["code_sha256"]==s.sha(pathlib.Path(__file__).read_bytes()) and p["bank_sha256"]==s.sha((O/"bank.npy").read_bytes())
    bank=np.load(O/"bank.npy",allow_pickle=False);p0=base.plan();g=p0["groups"][i];num=str(i+1).zfill(2)
    old=json.loads((C/("C-"+num+".json")).read_text());ref=json.loads((C/("parent-"+num+".json")).read_text())
    m,LL,old0=base.locks("minimal",p0);seed={"n6":a.v.C["seed"],"n7":base.w.C["seed"],"n9":base.old9.C["seed"]}[g["source"]];n.set_config(seed,a.C["windows_mhz"],a.C["half_hz"],a.C["max_mib_per_scan"]);base.b.taylor.init_tree()
    fs=g["group"]["files"];net={"payload_bytes":0,"requests":0,"files":0};orig=base.eng.base.RangeFile
    class Counted(orig):
        def close(self):
            if not self.closed and hasattr(self,"total") and hasattr(self,"size"):
                net["payload_bytes"]+=int(self.total);net["requests"]+=int(self.requests_count);net["files"]+=1
            super().close()
    base.eng.base.RangeFile=Counted;t0=time.perf_counter()
    try:rr,R,ch,hs,meta,cm,pre=n.gather(g["group"],"r41s1b_"+num)
    finally:base.eng.base.RangeFile=orig
    acq=time.perf_counter()-t0;assert sorted(hs)==g["raw_tile_hashes"]==ref["raw_tiles"] and net["files"]==4
    out={"stage":"GROUP","group":i+1,"target":g["group"]["target"],"fit_viewed":g["fit"],"parity":True,"records":[],"acquisition_wall":acq,"baseline_wall":0.,"original_trace_wall":0.,"advisory_wall":0.,**net}
    for kind,rows,refs,mm in (("eval",rr,R,meta),("shared",[x[0] for x in ch],np.stack([x[1] for x in ch]),cm)):
        t=time.perf_counter();d=base.infer("minimal",rows,refs,mm,pre,LL,old0,m,None,kind);out["baseline_wall"]+=time.perf_counter()-t
        oldarr=a.read_npz(base.C/("dev_"+num+"_"+kind+".npz"));a.verify_rows(oldarr,rows);assert base.check(d,oldarr,m)==ref["metrics"][kind]
        mask=s.sha(np.ascontiguousarray(d["alert"]).tobytes()+np.ascontiguousarray(d["hit"]).tobytes());assert mask==old["baseline_masks"][kind];before=d["alert"].copy()
        records=json.loads((C/("trace-"+num+"-"+kind+".json")).read_text())["normal"];assert [x["row"] for x in records]==np.flatnonzero(before).tolist()
        advtrace=[]
        for x in records:
            j=x["row"];lane="parent" if d["parent_alert"][j] else "rescue";path=d["parent_path"][j] if lane=="parent" else d["p"][j];src=str(d["source_tile_sha"][j]);assert lane==x["lane"] and src==x["source_tile_sha256"]
            t=time.perf_counter();tr=s.spiral(rows[j][0],refs[j],path,lane,mm[j],src,fs,pre,bank);out["original_trace_wall"]+=time.perf_counter()-t
            assert tr==x["trace"] and not tr["veto"],"sealed original trace parity"
            tt=tests(tr,lane);t=time.perf_counter();ad=advisory(tr,lane,mm[j],src,fs,pre);out["advisory_wall"]+=time.perf_counter()-t
            assert not ad["veto"] and ad["advanced"]==bool(tt is not None and tt["source"])
            if ad["advanced"]:
                repeat=ad["states"][-1]["snapshot"]["evidence"]["advisory_repeat"]["D4"]
                # Repeated semantic state is explicitly a no-gain stop.
                hh=list(ad["states"]);last=hh[-1]["snapshot"];assert not s.state(last["evidence"],hh,last["unresolved"]) and hh==ad["states"]
            else:repeat={}
            row={"row":j,"kind":kind,"lane":lane,"tests":tt,"original_open":bool(tr["turns"]==2),"advisory_advanced":ad["advanced"],"additional_calls":ad["calls"],"repeat":repeat,"label_evaluation_only":int(d["y"][j]),"localized_evaluation_only":bool(d["hit"][j])}
            out["records"].append(row);advtrace.append({"row":j,"lane":lane,"source_tile_sha256":src,"advisory":ad})
        assert np.array_equal(before,d["alert"]) and all(not v["advisory"]["veto"] for v in advtrace)
        save("trace-"+num+"-"+kind+".json",{"records":advtrace,"baseline_mask":mask,"decision_changes":0})
    out["summary"]=summarize(out["records"]);out["full_wall"]=time.perf_counter()-BOOT;out["cpu"]=time.process_time()-CPU;out["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;save("group-"+num+".json",out)
    cp("group-"+num,{k:v for k,v in out.items() if k not in ("records","summary")}|{"summary":{k:v for k,v in out["summary"].items() if k not in ("cofail","gates")}})
def summary():
    groups=[read("group-"+str(i+1).zfill(2)+".json") for i in range(11)];assert all(g["parity"] for g in groups)
    rows={split:[x for g in groups if g["fit_viewed"]==(split=="fit") for x in g["records"]] for split in ("fit","audit")}
    tables={}
    for split,rr in rows.items():
        tables[split]={"all":summarize(rr),"by_lane":{k:summarize([x for x in rr if x["lane"]==k]) for k in ("parent","rescue")},"by_kind":{k:summarize([x for x in rr if x["kind"]==k]) for k in K},"posthoc_evaluation_only":{"controls":summarize([x for x in rr if x["label_evaluation_only"]==0]),"positives":summarize([x for x in rr if x["label_evaluation_only"]==1])}}
    a0=tables["audit"]["all"];assert (a0["triggered"],a0["pass1"],a0["qualified"],a0["advisory_advanced"],a0["original_open"])==(417,414,409,409,0)
    assert tables["audit"]["posthoc_evaluation_only"]["controls"]["triggered"]==53
    for k,v in {"stationary":206,"native_wide":378,"rival_consensus":2,"rival_stationary_agreement":0,"crossfit_margin":0,"on_ceiling":0,"off_floor":10,"source":409}.items():assert a0["gates"][k]["pass"]==v
    cost={k:(max(g[k] for g in groups) if "rss" in k else sum(g[k] for g in groups)) for k in ("acquisition_wall","baseline_wall","original_trace_wall","advisory_wall","full_wall","cpu","peak_rss_kib","payload_bytes","requests","files")}
    d={"status":"PASS_TRANSITION_DECOUPLED_ADVISORY__HOLD_SCIENTIFIC_SCREEN","groups":11,"rows":3168,"tables":tables,"cost":cost,"original_metrics":{"localized":318,"weak":69,"controls":53,"repeat_stress_cases_original":119},"decision_changes":0,"new_stress_test":False,"source_seal_sha256":L["parent"]["seal_sha256"],"source_result_sha256":L["parent"]["result_sha256"],"protocol_sha256":s.sha(P.read_bytes()),"trace_sha256":{p.name:s.sha(p.read_bytes()) for p in sorted(O.glob("trace-*.json"))},"production_adoption":False,"external_validation":False}
    save("result.json",d);cp("99-summary",{"status":d["status"],"groups":11,"audit":a0,"fit":{k:v for k,v in tables["fit"]["all"].items() if k not in ("gates","cofail")},"cost":cost,"decision_changes":0})
if __name__=="__main__":
    if ARGS[0]=="prepare":prepare()
    elif ARGS[0]=="group":group(int(ARGS[1]))
    elif ARGS[0]=="summary":summary()
    else:raise ValueError(ARGS)
