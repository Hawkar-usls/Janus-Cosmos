#!/usr/bin/env python3
"""R37C physical evidence graph. Retrieval diagnostic only, never a source label or veto."""
import datetime,hashlib,itertools,json,math,pathlib,re
import numpy as np
import r37a_hannah_resolution as a

P=pathlib.Path("research/seti/r37c_physical_identity_protocol.json");S=json.loads(P.read_text());T=S["thresholds"]
O=pathlib.Path("runtime/r37c");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def metadata():
    p={}
    for k,c in a.S["source"].items():
        raw=pathlib.Path(c["sidecar"]).read_bytes();assert sha(raw)==c["sidecar_sha256"]
        d=json.loads(raw)
        for e in d["entries"]:
            ps=[e["raw_evidence"]["passport"]] if k=="stationary" else [e["passport"]]+[v["passport"] for v in e["off_matches"]]
            for v in ps:
                a.mem.validate(v);assert v["id"] not in p;p[v["id"]]=v
    assert len(p)==10
    return p
def epoch(m):
    scan=m["source_scan"];s=re.search(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d+",scan);assert s
    return datetime.datetime.fromisoformat(s[0]).replace(tzinfo=datetime.timezone.utc).timestamp()+m["start_row"]*m["ts"]
def physical(r,v):
    m=v["provenance"];x=r["x"];p=r["p"];df=m["df_hz"];ts=m["ts"]
    assert len(p)==16 and len(x)==16
    fr=m["channel0_mhz"]*1e6+df*p
    assert abs(fr[0]-v["frequency_hz"])<1e-5
    t=np.arange(len(p))*ts;tc=t-t.mean();slope=float(tc@fr/(tc@tc)) if np.all(fr==fr[0]) is False else 0.
    # Center frequencies numerically before slope fit, including exactly stationary case.
    mid=float(fr.mean());slope=float(tc@(fr-mid)/(tc@tc));res=fr-mid-slope*tc
    half=abs(df)/2;error=half+float(np.max(abs(res)))
    slope_u=float(error*np.sum(abs(tc))/(tc@tc))
    if df<0:x=x[:,::-1];p=x.shape[1]-1-p
    df=abs(df);te=np.arange(17)*ts;fe=(np.arange(x.shape[1]+1)-.5)*df
    tt=(te[:-1]+te[1:])/2;ff=p*df
    shape,contrast,cells=a.descriptor(x,te,fe,tt,ff);cells=np.asarray(cells)
    core=np.maximum(cells[:,3]-cells[:,[0,-1]].mean(1),0)
    mean=float(core.mean());assert mean>1e-12
    mod=(core-mean)/mean
    return {"id":r["id"],"group":r["group"],"scope":r["scope"],"hash":r["hash"],"scan":m["source_scan"],"epoch_mid_seconds":epoch(m)+float(t.mean()),"frequency_mid_hz":mid,"frequency_half_width_hz":error,"slope_hz_s":slope,"slope_half_width_hz_s":slope_u,"signed_drift":int(1 if slope>slope_u else -1 if slope<-slope_u else 0),"row_interval_s":ts,"channel_width_hz":df,"modulation_profile":mod.tolist(),"modulation_mean_abs_difference":float(np.mean(abs(np.diff(mod)))),"modulation_rms":float(np.sqrt(np.mean(mod**2))),"shape":shape.tolist(),"contrast_norm":contrast}
def harmonic(x,y):
    delta=abs(x["frequency_mid_hz"]-y["frequency_mid_hz"])/1e6
    ans=[]
    for spacing in (133.33,15.,33.3,128.003,2.000004):
        n=int(round(delta/spacing))
        if 1<=n<=1000:
            residual=abs(delta-n*spacing)
            if residual<=.01:ans.append({"spacing_mhz":spacing,"multiple":n,"residual_mhz":residual,"authority":"DIAGNOSTIC_ONLY"})
    return ans
def pair(x,y):
    shape=a.sim(np.asarray(x["shape"]),np.asarray(y["shape"]))
    mod=float(np.sqrt(np.mean((np.asarray(x["modulation_profile"])-np.asarray(y["modulation_profile"]))**2)))
    dt=abs(x["epoch_mid_seconds"]-y["epoch_mid_seconds"]);t=(x["epoch_mid_seconds"]+y["epoch_mid_seconds"])/2
    def prop(v):
        d=t-v["epoch_mid_seconds"]
        return v["frequency_mid_hz"]+v["slope_hz_s"]*d,v["frequency_half_width_hz"]+v["slope_half_width_hz_s"]*abs(d)
    fx,ux=prop(x);fy,uy=prop(y);dr=abs(x["slope_hz_s"]-y["slope_hz_s"])
    f=dt<=T["max_time_gap_s"] and abs(fx-fy)<=ux+uy+T["frequency_guard_hz"]
    d=dr<=x["slope_half_width_hz_s"]+y["slope_half_width_hz_s"]+T["drift_guard_hz_s"]
    sh=shape>=T["shape"];md=mod<=T["modulation_distance"]
    arms={"shape":sh,"shape_frequency":sh and f,"shape_frequency_drift":sh and f and d,"shape_modulation":sh and md,"joint":sh and f and d and md}
    same=x["group"]==y["group"];known=x["scope"]==y["scope"]=="EXTERNAL_CENTRAL_REPLAY"
    return {"a":x["id"],"b":y["id"],"groups":[x["group"],y["group"]],"same_reported_group":same,"known_central_both":known,"distinct_scan":x["scan"]!=y["scan"],"distinct_tile":x["hash"]!=y["hash"],"shape_cosine":shape,"modulation_distance":mod,"time_separation_s":dt,"conditional_midpoint_delta_hz":abs(fx-fy),"conditional_midpoint_uncertainty_hz":ux+uy,"drift_delta_hz_s":dr,"frequency_compatible":f,"drift_compatible":d,"arms":arms,"counterfactual_frequency_drift_collision":sh and md,"harmonic_diagnostics":harmonic(x,y),"authority":"RETRIEVAL_ONLY"}
def counts(rows):
    return {"pairs":len(rows),"arms":{k:sum(v["arms"][k] for v in rows) for k in S["arms"]},"counterfactual_collision_pass":sum(v["counterfactual_frequency_drift_collision"] for v in rows),"harmonic_matches":sum(bool(v["harmonic_diagnostics"]) for v in rows)}
def main():
    meta=metadata();rr=list(a.records());assert len(rr)==10 and len({v["hash"] for v in rr})==10
    ds=[physical(r,meta[r["id"]]) for r in rr];pp=[pair(x,y) for x,y in itertools.combinations(ds,2)]
    assert len(pp)==45 and all(v["distinct_tile"] for v in pp)
    pos=[v for v in pp if v["known_central_both"] and v["same_reported_group"]];neg=[v for v in pp if v["known_central_both"] and not v["same_reported_group"]]
    review=[v for v in pp if not v["known_central_both"]];cross=[v for v in pp if not v["same_reported_group"]]
    assert len(pos)==6 and len(neg)==9 and len(cross)==37
    assert all(v["distinct_scan"] for v in pos)
    pc=counts(pos);nc=counts(neg);xc=counts(cross);rc=counts(review)
    gates={"central_retention":pc["arms"]["joint"]==6,"central_cross_group_zero":nc["arms"]["joint"]==0}
    d={"status":"PASS_DEVELOPMENT_GROUP_RETRIEVAL" if all(gates.values()) else "HOLD_DEVELOPMENT_GROUP_RETRIEVAL","gates":gates,"measured_tiles":10,"pair_count":45,"known_same":pc,"known_different":nc,"all_different":xc,"review_only":rc,"signed_drift_resolved":sum(v["signed_drift"]!=0 for v in ds),"harmonic_authority":"DIAGNOSTIC_ONLY","counterfactual_scope":"Force frequency/drift compatible; shape/modulation only on distinct measured tiles with artificial metadata collision; not a new measurement.","independent_origin_gate":"HOLD","verified_promoted":0,"production_adoption":False}
    save("descriptors.json",ds);save("pairs.json",pp);save("result.json",d);print("R37C_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
