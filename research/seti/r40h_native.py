#!/usr/bin/env python3
"""R40H native comparison, historical null calibration before reading new spectra."""
import sys
MODE=sys.argv[1];assert MODE in ("bliss","fastdd","turbo")
import r35g_native_compare as native
import json,io,os,pathlib,time
import numpy as np
import r40h_frozen as h
c=h.c;t=h.t
O=pathlib.Path("runtime/r40h_native")/MODE;O.mkdir(parents=True,exist_ok=True)
native.OUT=O/"work";native.OUT.mkdir(exist_ok=True)
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def main():
    fn,setup=native.init();z,n,rec,_=t.parent.restore();cal=[]
    for r in rec["records"]:
        if r["name"].startswith("cal_"):
            q=t.parent.load(z,n,r);s,p,tm=fn(q["V"],pathlib.Path(r["name"]).stem);cal.extend(s[q["y"]==0])
    v=np.array(cal,dtype=float);th=t.parent.cut(v);assert len(v)==256 and int((v>=th).sum())<=12
    lock={"mode":MODE,"threshold":th,"calibration_controls":256,"alerts":int((v>=th).sum()),"source":h.S["historical"],"protocol_sha256":c.sha(h.P.read_bytes())}
    save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes())
    print("R40H_NATIVE_LOCK",json.dumps(lock),flush=True)
    a={"id":int(os.environ["EVAL_ID"]),"sha256":os.environ["EVAL_DIGEST"].removeprefix("sha256:")}
    z,n=c.artifact(a,128);r=json.loads(z.read(n["receipt.json"]))
    assert r["status"]=="PROSPECTIVE_ARCHIVE_INJECTION_ONLY" and len(r["records"])==10
    out=[]
    for rec in r["records"]:
        raw=z.read(n[rec["name"]]);assert c.sha(raw)==rec["sha256"]
        with np.load(io.BytesIO(raw),allow_pickle=False) as q:
            assert c.sha(q["V"].tobytes())==rec["array_sha256"] and q["V"].dtype==np.float32
            st=time.perf_counter();s,p,tm=fn(q["V"],pathlib.Path(rec["name"]).stem);wall=time.perf_counter()-st
            y=q["y"];weak=np.isin(q["amp"],t.S["weak_amplitudes"]);pred=s>=th
            valid=((p>=0)&(p<128)).all(1)
            hit=pred&valid&(((abs(p-q["truth"])<=2)&q["active"]).sum(1)>=.75*q["active"].sum(1))&(y==1)
            b=q["r28_hit"];g=q["r40g_hit"]
            row={"record":rec["name"],"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"control_alerts":int((pred&(y==0)).sum()),"edge_alerts":int((pred&~valid).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum()),"added_vs_r40g":int((hit&~g).sum()),"lost_vs_r40g":int((g&~hit).sum()),"added_vs_r28":int((hit&~b).sum()),"lost_vs_r28":int((b&~hit).sum()),"adapter_wall_seconds":wall,"timing":tm}
            out.append(row);save("partial.json",out);print("R40H_NATIVE_GROUP",json.dumps(row),flush=True)
            np.savez_compressed(O/rec["name"],score=s,path=p,alert=pred,localized=hit)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def tot(p):
        rr=[x for x in out if x["record"].startswith(p)];assert len(rr)==5
        return {k:sum(x[k] for x in rr) for k in ("localized","weak","control_alerts","edge_alerts","positives","controls","added_vs_r40g","lost_vs_r40g","added_vs_r28","lost_vs_r28","adapter_wall_seconds")}
    d={"status":"NATIVE_COMMON_ARRAYS_PROSPECTIVE_REVIEW_ONLY","mode":MODE,"source":a,"lock_sha256":lh,"threshold":th,"eval":tot("eval_"),"shared":tot("shared_"),"setup_seconds":setup,"scope":"Pinned CPU adapters on shared JANUS preprocessing; no whole-native pipeline superiority claim; no external blind qualification","production_adoption":False}
    save("result.json",d);print("R40H_NATIVE_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
