#!/usr/bin/env python3
"""R40H native comparison, historical null calibration before reading new spectra."""
import sys
MODE=sys.argv[1];assert MODE in ("bliss","fastdd","turbo")
import r35g_native_compare as native
import json,io,os,pathlib,time
import numpy as np
import r40h_frozen as h
c=h.c;t=h.t
O=pathlib.Path("runtime/r40ia_matched")/MODE;O.mkdir(parents=True,exist_ok=True)
native.OUT=O/"work";native.OUT.mkdir(exist_ok=True)
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def main():
    fn,setup=native.init();z,n,rec,_=t.parent.restore();cal=[];K=0;L=h.readlocks()
    protocol=json.loads(pathlib.Path('research/seti/r40ia_matched_protocol.json').read_text())
    for r in rec["records"]:
        if r["name"].startswith("cal_"):
            q=t.parent.load(z,n,r);s,p,tm=fn(q["V"],pathlib.Path(r["name"]).stem);cal.extend(s[q["y"]==0]);hh=h.f.historical(q,L["r40c"]);_,_,_,pred=h.replay(hh,L);K+=int((pred&(q["y"]==0)).sum())
    v=np.array(cal,dtype=float);assert len(v)==256 and 0<=K<256;th=float(np.nextafter(np.sort(v)[::-1][K],np.inf));assert int((v>=th).sum())<=K
    lock={"mode":MODE,"threshold":th,"calibration_controls":256,"r40g_calibration_alerts":K,"alerts":int((v>=th).sum()),"source":h.S["historical"],"protocol_sha256":c.sha(pathlib.Path("research/seti/r40ia_matched_protocol.json").read_bytes())}
    save("lock.json",lock);lh=c.sha((O/"lock.json").read_bytes())
    print("R40IA_LOCK",json.dumps(lock),flush=True)
    a={"id":int(os.environ["EVAL_ID"]),"sha256":os.environ["EVAL_DIGEST"].removeprefix("sha256:")}
    assert a=={"id":14476995403,"sha256":"cbb45f9fda24c56b8eddda710e1206c7a1adf1136f90e1d1979458e70a4e31ac"}
    z,n=c.artifact(a,192);r=json.loads(z.read(n["receipt.json"]));nz,nn=c.artifact(protocol["native_artifacts"][MODE],8)
    assert r["status"]=="PROSPECTIVE_ARCHIVE_INJECTION_ONLY" and len(r["records"])==10
    out=[]
    for rec in r["records"]:
        raw=z.read(n[rec["name"]]);assert c.sha(raw)==rec["sha256"]
        with np.load(io.BytesIO(raw),allow_pickle=False) as q:
            assert c.sha(q["V"].tobytes())==rec["array_sha256"] and q["V"].dtype==np.float32
            with np.load(io.BytesIO(nz.read(nn[rec["name"]])),allow_pickle=False) as old:
                s,p=old["score"],old["path"]
            tm={"source":"frozen native artifact"};wall=0.0
            y=q["y"];weak=np.isin(q["amp"],t.S["weak_amplitudes"]);pred=s>=th
            valid=((p>=0)&(p<128)).all(1)
            hit=pred&valid&(((abs(p-q["truth"])<=2)&q["active"]).sum(1)>=.75*q["active"].sum(1))&(y==1)
            b=q["r28_hit"];g=q["r40g_hit"]
            row={"record":rec["name"],"localized":int(hit.sum()),"weak":int((hit&weak).sum()),"control_alerts":int((pred&(y==0)).sum()),"edge_alerts":int((pred&~valid).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum()),"added_vs_r40g":int((hit&~g).sum()),"lost_vs_r40g":int((g&~hit).sum()),"added_vs_r28":int((hit&~b).sum()),"lost_vs_r28":int((b&~hit).sum()),"rescore_wall_seconds":wall,"timing":tm}
            out.append(row);save("partial.json",out);print("R40IA_GROUP",json.dumps(row),flush=True)
            np.savez_compressed(O/rec["name"],score=s,path=p,alert=pred,localized=hit)
        assert c.sha((O/"lock.json").read_bytes())==lh
    def tot(p):
        rr=[x for x in out if x["record"].startswith(p)];assert len(rr)==5
        return {k:sum(x[k] for x in rr) for k in ("localized","weak","control_alerts","edge_alerts","positives","controls","added_vs_r40g","lost_vs_r40g","added_vs_r28","lost_vs_r28","rescore_wall_seconds")}
    d={"status":"MATCHED_HISTORICAL_CALIBRATION_REPLAY_ONLY","mode":MODE,"source":a,"lock_sha256":lh,"threshold":th,"eval":tot("eval_"),"shared":tot("shared_"),"setup_seconds":setup,"calibration_target":K,"actual_calibration_alerts":int((v>=th).sum()),"scope":"Matched historical calibration budget, R40H already viewed, frozen native scores reused; not full raw pipeline or independently matched test FAP","production_adoption":False}
    save("result.json",d);print("R40IA_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
