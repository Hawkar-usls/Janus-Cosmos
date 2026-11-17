#!/usr/bin/env python3
"""R5: frozen R4 spectrum replay, narrowband track tests; descriptive only."""
import hashlib,json,math,pathlib,sys
import numpy as np

def robust(v):
    m=float(np.median(v)); s=float(1.4826*np.median(np.abs(v-m)))
    return m,s

def load(folder,o,p):
    x=np.load(folder/p["artifact"],allow_pickle=False)
    if hashlib.sha256(x.tobytes()).hexdigest()!=p["slice_sha256"]: raise ValueError("slice hash mismatch")
    if list(x.shape)!=p["shape"] or x.ndim!=2 or not np.isfinite(x).all(): raise ValueError("bad array")
    if x.shape[0]!=o["sampled_time_rows"]: raise ValueError("time rows mismatch")
    df0=float(o["channel_hz"]); f0=float(p["channel0_mhz"])
    f=float(p["frequency_mhz"])
    df=np.copysign(df0,f-f0)
    center=(f-f0)*1e6/df
    if abs(center-(x.shape[1]-1)/2)>0.6: raise ValueError("frequency grid ambiguity")
    return x,df,center

def track(x,df,ts,intercept_hz,drift,center):
    t=np.arange(x.shape[0])*ts
    c=np.rint(center+(intercept_hz+drift*t)/df).astype(int)
    if np.any(c<1) or np.any(c>=x.shape[1]-1): return None
    # Row-level robust normalization removes time-varying gain. No calibrated false alarm claim.
    m=np.median(x,axis=1); s=1.4826*np.median(np.abs(x-m[:,None]),axis=1)
    if np.any(s<=0): return None
    y=(x-m[:,None])/s[:,None]
    v=y[np.arange(len(c)),c]
    # Compare with same drift tracks translated in frequency; exclude a guard region.
    null=[]
    for k in range(-x.shape[1],x.shape[1]+1,3):
        if abs(k)<8: continue
        cc=c+k
        if cc.min()<0 or cc.max()>=x.shape[1]: continue
        null.append(float(np.mean(y[np.arange(len(c)),cc])))
    if len(null)<30: return None
    med,scale=robust(np.array(null))
    if scale<=0:return None
    return {"score":float((np.mean(v)-med)/scale),"row_excess_fraction":float(np.mean(v>3)),
            "null_tracks":len(null),"intercept_hz":float(intercept_hz),"drift_hz_s":float(drift),
            "rows":len(c),"channel_hz":float(abs(df))}

def search(x,df,ts,center,expected,drifts,half):
    if len(x)<12:return None
    # Finite fixed grid, common search width in ON and OFF; maximization is descriptive.
    opts=[]
    for delta in np.arange(-half,half+abs(df)*.01,abs(df)):
        for drift in drifts:
            a=track(x,df,ts,expected+float(delta),float(drift),center)
            if a is not None: opts.append(a)
    return max(opts,key=lambda a:a["score"]) if opts else None

def strong(a):
    return bool(a and a["score"]>=8 and a["row_excess_fraction"]>=.5)

def analyze(folder):
    d=json.loads((folder/"result.json").read_text())
    obs=d["observations"]
    result={"target":d["target"],"groups":[],"status":"HOLD_NO_ORIGINAL_OR_CONTROL"}
    for line in d["groups"]:
        g={"group_mhz":line["group_mhz"],"frequency_mhz":line["frequency_mhz"],"catalog_drift_hz_s":line["catalog_drift_hz_s"],"status":"HOLD_MISSING_REQUIRED_DATA","tests":[]}
        result["groups"].append(g)
        def product(o):
            return next((p for p in o.get("products",[]) if p["group_mhz"]==line["group_mhz"] and "artifact" in p),None)
        original=next((o for o in obs if o["role"]=="ORIGINAL_ON"),None)
        if not original or not product(original) or line["catalog_drift_hz_s"] is None:continue
        x,df,c=load(folder,original,product(original));ts=float(original["tsamp_s"]);drift=float(line["catalog_drift_hz_s"])
        if not original["full_time_scan"]:continue
        delta=abs(df)/(max(1,(len(x)-1)*ts))
        drifts=[drift-delta,drift,drift+delta]
        anchor=search(x,df,ts,c,0.,drifts,28.)
        g["anchor"]=anchor
        if not strong(anchor):
            g["status"]="HOLD_ORIGINAL_NARROW_TRACK_NOT_RECOVERED";continue
        for o in obs:
            if o is original:continue
            p=product(o)
            t={"role":o["role"],"url":o["url"],"status":"HOLD_NO_COMPARABLE_SLICE"}
            g["tests"].append(t)
            if not p:continue
            y,df2,c2=load(folder,o,p);ts2=float(o["tsamp_s"])
            if not o["full_time_scan"] or y.shape[0]!=x.shape[0] or abs(ts2-ts)>1e-6 or abs(abs(df2)-abs(df))>1e-5:continue
            elapsed=(float(o["tstart_mjd"])-float(original["tstart_mjd"]))*86400
            t["delta_seconds"]=elapsed
            t["same_fixed_topocentric_window"]=search(y,df2,ts2,c2,anchor["intercept_hz"],drifts,28.)
            if abs(elapsed)<=1800:
                predicted=anchor["intercept_hz"]+anchor["drift_hz_s"]*elapsed
                t["continuous_track"]=search(y,df2,ts2,c2,predicted,drifts,28.)
                t["predicted_intercept_hz"]=predicted
            else:t["continuous_track"]=None
            t["status"]="DESCRIPTIVE_TRACK_COMPARISON"
            t["fixed_window_strong"]=strong(t["same_fixed_topocentric_window"])
            t["continuous_strong"]=strong(t["continuous_track"])
            t["sensitivity_note"]="Equal channel spacing, integration and rows; gain normalized, not flux calibrated."
        off=[t for t in g["tests"] if t["role"] in ("ARCHIVE_LABELED_OFF","NEARBY_OTHER_TARGET_CONTROL")]
        repeat=[t for t in g["tests"] if t["role"]=="SAME_SESSION_ON_REPEAT"]
        if any(t.get("continuous_strong") or t.get("fixed_window_strong") for t in off):
            g["status"]="CONTROL_TRACK_PRESENT_RFI_CONCERN"
        elif any(t.get("continuous_track") is not None and t.get("same_fixed_topocentric_window") is not None for t in off) and any(t.get("continuous_strong") for t in repeat):
            g["status"]="ON_OFF_ON_PRELIMINARY_NEEDS_CALIBRATION_AND_INDEPENDENT_REPEAT"
        else:g["status"]="HOLD_NO_COMPLETE_MATCHED_CADENCE"
    result["status"]="DESCRIPTIVE_ARCHIVE_TRACK_SCREEN_COMPLETE"
    return result

def selftest():
    rng=np.random.default_rng(17); x=rng.normal(0,1,(16,431));df=2.8;ts=18.
    for i in range(16):x[i,215+round(.08*i*ts/df)]+=15
    a=search(x,df,ts,215,0,[.08],5.6)
    assert strong(a),a
    assert not strong(search(rng.normal(0,1,(16,431)),df,ts,215,0,[.08],5.6))
    assert search(x,df,ts,215,10000,[.08],5.6) is None
    print("SYNTHETIC_TRACK_TESTS=PASS")

def main():
    selftest()
    root=pathlib.Path(sys.argv[1]);out=pathlib.Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
    results=[analyze(p.parent) for p in sorted(root.glob("*/result.json"))]
    if len(results)!=7 or sum(len(r["groups"]) for r in results)!=12:raise ValueError("missing frozen groups")
    counts={}
    for r in results:
        for g in r["groups"]:
            counts[g["status"]]=counts.get(g["status"],0)+1
            print("R5_GROUP",json.dumps({"target":r["target"],**g}),flush=True)
    final={"schema":"janus.seti.r5.track_screen.v1","source_run":44041795339,"targets":len(results),"groups":12,
           "counts":counts,"independent_reobservation_confirmed":0,"confirmed_technosignatures":0,
           "caveats":["Robust shifted-track score is descriptive, not calibrated SNR or p-value. Search max selection not calibrated.",
                      "Fixed candidate line selected in R4; limited 1200 Hz window may clip fast drifting signals.",
                      "Nearby different target is not verified beam OFF; no absolute flux sensitivity comparison.",
                      "Different-day comparison at same topocentric frequency is not barycentric re-detection.",
                      "No new observations or independent experimental confirmation."],"results":results}
    (out/"r5-result.json").write_text(json.dumps(final,indent=2)+"\n")
    print("R5_SUMMARY",json.dumps({k:v for k,v in final.items() if k!="results"}))
if __name__=="__main__":main()
