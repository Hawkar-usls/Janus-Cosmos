#!/usr/bin/env python3
"""R6 diagnostic replays, wider windows and frozen start/middle/end conventions.
No calibrated astrophysical claim; other-target controls are not verified beam OFF.
"""
import hashlib,json,pathlib,sys
import numpy as np

def strength(x):
    med=np.median(x,axis=1)
    scale=1.4826*np.median(np.abs(x-med[:,None]),axis=1)
    if not np.isfinite(x).all() or (scale<=0).any(): raise ValueError("nonfinite or zero-noise rows")
    return (x-med[:,None])/scale[:,None]

def load(root,o,p):
    x=np.load(root/p["artifact"],allow_pickle=False)
    if hashlib.sha256(x.tobytes()).hexdigest()!=p["slice_sha256"]:raise ValueError("slice hash mismatch")
    if x.ndim!=2 or list(x.shape)!=p["shape"] or not o["full_time_scan"]:raise ValueError("incomplete spectrum")
    df=np.copysign(float(o["channel_hz"]),float(p["frequency_mhz"])-float(p["channel0_mhz"]))
    center=(p["frequency_mhz"]-p["channel0_mhz"])*1e6/df
    if not 0<=center<x.shape[1]:raise ValueError("frequency grid outside slice")
    return strength(x),df,center,float(o["tsamp_s"])

def score_grid(y,df,center,ts,expected,drifts,offset_bins=10):
    t=np.arange(len(y))*ts
    best=None; n=0
    for drift in drifts:
        for offset in range(-offset_bins,offset_bins+1):
            intercept=expected+offset*abs(df)
            c=np.rint(center+(intercept+drift*t)/df).astype(int)
            if c.min()<0 or c.max()>=y.shape[1]:continue
            vals=y[np.arange(len(y)),c]
            # Same slope at frequency offsets, keeping a guard around tested line.
            ks=np.arange(-y.shape[1],y.shape[1],7)
            ks=ks[(abs(ks)>=14)&(c.min()+ks>=0)&(c.max()+ks<y.shape[1])]
            if len(ks)<30:continue
            null=y[np.arange(len(y))[:,None],c[:,None]+ks[None,:]].mean(axis=0)
            m=float(np.median(null)); mad=float(1.4826*np.median(np.abs(null-m)))
            if mad<=0:continue
            v={"score":float((vals.mean()-m)/mad),"row_excess_fraction":float(np.mean(vals>3)),
               "intercept_hz":float(intercept),"drift_hz_s":float(drift),"null_tracks":len(ks),
               "median_row_excess":float(np.median(vals))}
            n+=1
            if best is None or v["score"]>best["score"]:best=v
    if best is not None:best["tested_tracks"]=n
    return best

def strong(a):return bool(a and a["score"]>=8 and a["row_excess_fraction"]>=.5)

def analyze(root):
    d=json.loads((root/"result.json").read_text());obs=d["observations"]
    out={"target":d["target"],"source_replay_status":d["status"],"hits":[],"groups":{}}
    for line in d["groups"]:
        hit={**line,"status":"HOLD_NO_COMPARABLE_SOURCE","controls":[]}
        out["hits"].append(hit)
        def prod(o):
            return next((p for p in o.get("products",[]) if p["group_mhz"]==line["group_mhz"] and p["catalog_hit_index"]==line["catalog_hit_index"] and "artifact" in p),None)
        orig=next((o for o in obs if o["role"]=="ORIGINAL_ON"),None)
        if not orig or not prod(orig):continue
        try:
            y,df,c,ts=load(root,orig,prod(orig));dr=float(line["catalog_drift_hz_s"])
            span=(len(y)-1)*ts;step=abs(df)/max(span,1)
            drifts=[dr-step,dr,dr+step]
            # Price et al uses observing START time in extrapolation. Mid/end alternatives
            # are explicitly diagnostic; cannot count as a new blind detection.
            variants={"start":0.,"mid":-dr*span/2,"end":-dr*span}
            fits={k:score_grid(y,df,c,ts,v,drifts) for k,v in variants.items()}
            hit["time_convention_fits"]=fits
            anchor=fits["start"]
            hit["status"]="HOLD_START_CONVENTION_TRACK_NOT_RECOVERED"
            if not strong(anchor):
                if any(strong(a) for a in fits.values()):hit["status"]="HOLD_TIME_CONVENTION_AMBIGUOUS"
                continue
            hit["anchor"]=anchor
            hit["status"]="HOLD_NO_MATCHED_CONTROL_OR_REPEAT"
            for o in obs:
                if o is orig:continue
                rec={"role":o["role"],"url":o["url"],"status":"HOLD_NO_SLICE"};hit["controls"].append(rec)
                p=prod(o)
                if not p:continue
                z,df2,c2,ts2=load(root,o,p)
                if len(z)!=len(y) or abs(abs(df2)-abs(df))>1e-5 or abs(ts2-ts)>1e-6:
                    rec["status"]="HOLD_INCOMPARABLE_RESOLUTION";continue
                dt=(float(o["tstart_mjd"])-float(orig["tstart_mjd"]))*86400
                rec["dt_seconds"]=dt
                rec["fixed_window"]=score_grid(z,df2,c2,ts2,anchor["intercept_hz"],drifts)
                if abs(dt)<=1800:
                    pred=anchor["intercept_hz"]+anchor["drift_hz_s"]*dt
                    rec["continuous"]=score_grid(z,df2,c2,ts2,pred,drifts)
                    # Inject in normalized units only: this checks algorithm sensitivity, not flux.
                    injected=[]; amp=max(1.,anchor["median_row_excess"])
                    for k in (-70,-45,45,70):
                        f=pred+k*abs(df2)
                        ch=np.rint(c2+(f+anchor["drift_hz_s"]*np.arange(len(z))*ts2)/df2).astype(int)
                        if ch.min()<0 or ch.max()>=z.shape[1]:continue
                        zi=z.copy();zi[np.arange(len(zi)),ch]+=amp
                        a=score_grid(zi,df2,c2,ts2,f,[anchor["drift_hz_s"]],offset_bins=2)
                        injected.append(strong(a))
                    rec["normalized_injection_recovery"]={"passes":sum(injected),"trials":len(injected),"amplitude_row_noise_units":amp}
                else:rec["continuous"]=None
                rec["status"]="DESCRIPTIVE_MATCHED_TRACK_TEST"
                rec["strong_fixed"]=strong(rec["fixed_window"]);rec["strong_continuous"]=strong(rec["continuous"])
            off=[r for r in hit["controls"] if r["role"] in ("ARCHIVE_LABELED_OFF","NEARBY_OTHER_TARGET_CONTROL")]
            on=[r for r in hit["controls"] if r["role"]=="SAME_SESSION_ON_REPEAT"]
            if any(r.get("strong_fixed") or r.get("strong_continuous") for r in off):
                hit["status"]="CONTROL_TRACK_PRESENT_RFI_CONCERN"
            elif any(r.get("continuous") for r in off) and any(r.get("strong_continuous") for r in on):
                hit["status"]="PRELIMINARY_DIRECTIONAL_PATTERN_HOLD_FLUX_AND_INDEPENDENT_REPEAT"
        except Exception as e:hit["error"]=str(e);hit["status"]="HOLD_ANALYSIS_ERROR"
    for h in out["hits"]:
        key=str(h["group_mhz"]);counts=out["groups"].setdefault(key,{})
        counts[h["status"]]=counts.get(h["status"],0)+1
        print("R6_HIT",json.dumps({"target":out["target"],**h}),flush=True)
    out.update(independent_confirmation=False,confirmed_technosignatures=0,
       limitations=["All catalog hits tested; post-selection data, not independent blind detections.",
                    "Shifted-track score is descriptive, not calibrated SNR or a global p-value.",
                    "Mid/end convention diagnostics do not establish timestamp convention.",
                    "Normalized injection is algorithm sensitivity only, not flux calibration.",
                    "Same-frequency different-day tests are not barycentric comparisons.",
                    "No verified beam-OFF or new telescope measurement."])
    (root/"r6-analysis.json").write_text(json.dumps(out,indent=2)+"\n")
    print("R6_SUMMARY",json.dumps({k:v for k,v in out.items() if k!="hits"}),flush=True)

if __name__=="__main__":
    analyze(pathlib.Path(sys.argv[1]))
