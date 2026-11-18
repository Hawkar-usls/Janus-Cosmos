#!/usr/bin/env python3
"""Conditional exact drift bounds from quantized measured track bins, not source-origin labels."""
import itertools,json,math,pathlib
import numpy as np
import r37c_physical_identity as p
P=pathlib.Path("research/seti/r37d_quantized_drift_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37d");O.mkdir(parents=True,exist_ok=True)
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def bound(t,f,w):
    t,f,w=map(lambda x:np.asarray(x,float),(t,f,w));ix=np.argsort(t,kind="stable");t=t[ix];f=f[ix];w=w[ix]
    lo=-math.inf;hi=math.inf;eps=1e-10
    for i in range(len(t)-1):
        dt=t[i+1:]-t[i];df=f[i+1:]-f[i];ww=w[i+1:]+w[i]
        if np.any(abs(df[dt==0])>ww[dt==0]+1e-8):return {"feasible":False,"sign":None}
        m=dt>0
        if m.any():lo=max(lo,float(np.max((df[m]-ww[m])/dt[m])));hi=min(hi,float(np.min((df[m]+ww[m])/dt[m])))
    if lo>hi+eps:return {"feasible":False,"sign":None}
    assert math.isfinite(lo) and math.isfinite(hi)
    if lo>hi:lo=hi=(lo+hi)/2
    z=1 if lo>eps else -1 if hi<-eps else 0
    # Verify an actual intercept witness at midpoint of feasible slope interval.
    b=(lo+hi)/2;ft=f-f[0]-b*(t-t[0]);al=float(max(ft-w));au=float(min(ft+w))
    assert al<=au+1e-5
    return {"feasible":True,"lower_hz_s":lo,"upper_hz_s":hi,"width_hz_s":hi-lo,"sign":z,"intercept_witness_width_hz":max(0.,au-al)}
def canaries():
    t=np.arange(16)*10.;w=np.ones(16)*.5
    for b in (0.,.04,-.04):
        f=np.round(b*t);r=bound(t,f,w);assert r["feasible"] and r["lower_hz_s"]-1e-10<=b<=r["upper_hz_s"]+1e-10
        assert r["sign"]==int(np.sign(b))
    assert not bound([0,0,10],[0,2,2],[.5,.5,.5])["feasible"]
    assert not bound([0,10,20],[0,10,0],[.5,.5,.5])["feasible"]
    r=bound([0,10,20],[0,1,0],[.5,.5,.5]);assert r["feasible"] and r["sign"]==0
    return "PASS"
def main():
    test=canaries();meta=p.metadata();rr=list(p.a.records());assert len(rr)==10
    tracks=[];out=[]
    for r in rr:
        v=meta[r["id"]];m=v["provenance"];q=np.asarray(v["identity"]["track"]);df=m["df_hz"];t=p.epoch(m)+np.arange(len(q))*m["ts"];f=m["channel0_mhz"]*1e6+q*df;w=np.full(len(q),abs(df)/2)
        full=bound(t,f,w);short=bound(t[:16],f[:16],w[:16]);expected=int(np.sign(v["features"]["drift_hz_s"])) if r["scope"]=="UNKNOWN_REVIEW_ONLY" else 0
        feasible=0;stable=0
        for i in range(len(q)):
            for s in (-1,1):
                z=f.copy();z[i]+=s*abs(df);b=bound(t,z,w);feasible+=b["feasible"];stable+=b["feasible"] and b["sign"]==expected
        out.append({"id":r["id"],"group":r["group"],"scan":m["source_scan"],"scope":r["scope"],"rows":len(q),"expected_track_direction":expected,"full":full,"first16":short,"one_bin_stress":{"trials":2*len(q),"feasible":feasible,"expected_sign_retained":stable}})
        tracks.append((r,m,t,f,w))
    pairs=[]
    for (a,am,t,f,w),(b,bm,u,g,z) in itertools.combinations(tracks,2):
        if a["group"]!=b["group"]:continue
        assert am["source_scan"]!=bm["source_scan"]
        pairs.append({"a":a["id"],"b":b["id"],"group":a["group"],"scope":a["scope"],"joint":bound(np.r_[t,u],np.r_[f,g],np.r_[w,z]),"authority":"CONDITIONAL_AFFINE_COMPATIBILITY_ONLY"})
    assert len(pairs)==8
    known=[r for r in out if r["scope"]=="EXTERNAL_CENTRAL_REPLAY"];unknown=[r for r in out if r["scope"]=="UNKNOWN_REVIEW_ONLY"]
    k=sum(r["full"]["feasible"] and r["full"]["sign"]==0 for r in known);u=sum(r["full"]["feasible"] and r["full"]["sign"]==r["expected_track_direction"] and r["expected_track_direction"]!=0 for r in unknown)
    s=sum(r["first16"]["feasible"] and r["first16"]["sign"]==r["expected_track_direction"] and r["expected_track_direction"]!=0 for r in unknown)
    d={"status":"PASS_CONDITIONAL_DRIFT_DEVELOPMENT" if k==6 and u==4 else "HOLD_CONDITIONAL_DRIFT_DEVELOPMENT","canaries":test,"known_stationary_compatible":k,"known_stationary_total":6,"review_signed_full":u,"review_signed_first16":s,"review_total":4,"joint_pairs":len(pairs),"joint_feasible":sum(r["joint"]["feasible"] for r in pairs),"joint_known_feasible":sum(r["joint"]["feasible"] and r["scope"]=="EXTERNAL_CENTRAL_REPLAY" for r in pairs),"joint_review_feasible":sum(r["joint"]["feasible"] and r["scope"]=="UNKNOWN_REVIEW_ONLY" for r in pairs),"one_bin_stress":{k:sum(r["one_bin_stress"][k] for r in out) for k in ("trials","feasible","expected_sign_retained")},"scope":"Exact conditional bin interval model applied to stored measured paths; no raw source labels, no detector uncertainty or subchannel accuracy claim. Independent origin HOLD.","verified_promoted":0,"production_adoption":False}
    save("tracks.json",out);save("pairs.json",pairs);save("result.json",d);print("R37D_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
