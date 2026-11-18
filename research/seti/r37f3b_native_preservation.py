#!/usr/bin/env python3
"""Native resolution conservation of historical astronomical spectra, isolated read-only audit."""
import hashlib,json,pathlib,math
import numpy as np
P=pathlib.Path("research/seti/r37f3b_native_preservation_protocol.json");S=json.loads(P.read_text())
for n,h in S["pinned_blobs"].items():
    b=(P.parent/n).read_bytes()
    assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==h,n
import r37f2d3_methanol as science
import r37f3a_natural_preservation as guard
from astropy.utils import iers
from astropy.coordinates import solar_system_ephemeris
O=pathlib.Path("runtime/r37f3b");O.mkdir(parents=True,exist_ok=True)
L=S["limits"];assert guard.L==L
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def native(path,v,src,previous):
    with path.open("rb") as f:h,n=science.fil_header(f.read(65536))
    g=science.fil_grid(h,n,path.stat().st_size);assert g["shape"][0]==6
    k=previous["lsrk_factor_start_mid_end"][1];vlo,vhi=src["lsr_window"];margin=science.M["baseline_margin_km_s"];C=science.C;REST=science.REST
    fmin,fmax=sorted([REST*(1-(vlo-margin)/C)/k,REST*(1-(vhi+margin)/C)/k])
    a,b=sorted([(fmin-h["fch1"])/g["df"],(fmax-h["fch1"])/g["df"]]);lo=math.ceil(a);hi=math.floor(b)+1
    B=science.M["block_bins"];hi=lo+(hi-lo)//B*B
    assert 0<=lo<hi<=g["shape"][-1] and hi-lo<=500000
    mm=np.memmap(path,dtype="<f4",mode="r",offset=n,shape=tuple(g["shape"]))
    x=np.array(mm[:,0,lo:hi],dtype=np.float64);del mm
    assert np.isfinite(x).all() and (x>=0).all()
    vc=C*(1-(h["fch1"]+np.arange(lo,hi,dtype=float)*g["df"])*k/REST)
    off=(vc<vlo-margin/2)|(vc>vhi+margin/2)
    scale=np.median(x[:,off],axis=1);assert (scale>0).all();x/=scale[:,None]
    if vc[0]>vc[-1]:vc=vc[::-1].copy();x=x[:,::-1].copy()
    step=abs(C*g["df"]*k/REST);assert step>0 and np.allclose(np.diff(vc),step,rtol=1e-6,atol=1e-10)
    e=vc[0]+step*(np.arange(len(vc)+1)-.5)
    assert np.max(abs((e[:-1]+e[1:])/2-vc))<1e-8
    return vc,x,e,step
def fast_rebin(y,a,b):
    # Exact piecewise-constant integral, evaluated at output boundaries; O(n+m) storage.
    assert a[0]==b[0] and a[-1]==b[-1]
    # Center at the first edge: avoids cancellation from absolute velocities.
    a=a-a[0];b=b-b[0]
    integ=np.r_[0.,np.cumsum(y*np.diff(a),dtype=np.float64)]
    q=np.interp(b,a,integ)
    z=np.diff(q)/np.diff(b)
    err=abs(float(z@np.diff(b))-float(y@np.diff(a)))/max(1.,abs(float(y@np.diff(a))))
    assert err<=L["conservation_rel"] and (z>=0).all()
    return z,err
def grid(e,step,t):
    d=step*t["factor"];q=np.arange(e[0]+d*t["phase"],e[-1],d)
    q=q[(q>e[0]+1e-9)&(q<e[-1]-1e-9)]
    return np.r_[e[0],q,e[-1]]
def selftest():
    old,blob=guard.pinned_rebin();rng=np.random.default_rng(3703);err=0.
    for t in S["transforms"]:
        e=np.linspace(-3,3,129);f=grid(e,e[1]-e[0],t);y=rng.uniform(0,4,len(e)-1)
        ref,_=old(y[None,:],np.array([0.,1.]),e,np.array([0.,1.]),f)
        q,_=fast_rebin(y,e,f);err=max(err,float(np.max(abs(ref[0]-q))))
        assert np.allclose(ref[0],q,atol=1e-10,rtol=1e-10)
    return {"r37a_equivalence_max_abs_error":err,"r37a_blob":blob,"small_grid_cases":8}
def check_group(vc,e,step,y,src,group):
    line=src["lsr_window"];lo,hi=line;bg=float(np.median(y[(vc<lo-10)|(vc>hi+10)]));assert bg>0
    y.setflags(write=False);before=hashlib.sha256(y.tobytes()).hexdigest()
    rows=[]
    for t in S["transforms"]:
        f=grid(e,step,t);z,err=fast_rebin(y,e,f)
        out=guard.compare(e,y,f,z,bg,line)
        rows.append({"source":src["name"],"group":group,"transform":t,"conservation_rel":err,**out})
    # Independent broad damage tests with known expected verdict. No detector training.
    mask=(vc>=lo)&(vc<=hi);pk=vc[np.flatnonzero(mask)[np.argmax(y[mask])]]
    variants={"unchanged":(e,y),"attenuation_50pct":(e,bg+(y-bg)*.5),"peak_notch":(e,np.where(abs(vc-pk)<=.2,bg,y)),"shift_0.5km_s":(e,np.interp(vc-.5,vc,y,left=bg,right=bg)),"nonfinite":(e,np.where(abs(vc-pk)<=.02,np.nan,y)),"incomplete_support":(e[1:],y[1:])}
    can=[]
    for name,(f,z) in variants.items():
        out=guard.compare(e,y,f,z,bg,line);expected="PRESERVED" if name=="unchanged" else "REVIEW"
        can.append({"source":src["name"],"group":group,"name":name,"expected":expected,"correct":out["status"]==expected,**out})
    assert hashlib.sha256(y.tobytes()).hexdigest()==before
    return rows,can
def main():
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
    selftest_result=selftest();rr=science.refs();rows=[];can=[];sources=[]
    for s in science.S["sources"]:
        src={**s,"expected_sha256":S["expected_sha256"][s["name"]]};v=science.choose(src,rr);p=O/"temporary.fil"
        r={"source":s["name"],"day":s["day"],"status":"HOLD"}
        try:
            r["integrity"]=science.download(v,src,p)
            previous=science.inspect(p,v,src);assert previous["status"]=="CONDITIONAL_NATURAL_LINE_SUPPORT"
            vc,x,e,step=native(p,v,src,previous);assert len(x)==6
            r.update(native_df_hz=previous["native_df_hz"],rows=len(x),native_channels_in_science_window=len(vc),science_gate_reproduced=True)
            group_rows=[];group_can=[]
            for name,part in zip(S["groups"],[x,x[:3],x[3:]]):
                a,b=check_group(vc,e,step,part.mean(0),src,name);group_rows.extend(a);group_can.extend(b)
            rows.extend(group_rows);can.extend(group_can)
            r["status"]="PROCESSED";r["comparisons"]=len(group_rows)
        except Exception as e:
            import traceback
            traceback.print_exc();r["error"]=type(e).__name__+": "+str(e)[:400]
        finally:p.unlink(missing_ok=True)
        sources.append(r);print("R37F3B_SOURCE",json.dumps(r),flush=True)
    n=sum(r["status"]=="PRESERVED" for r in rows);k=sum(r["correct"] for r in can)
    d={"status":"PASS_NATIVE_ISOLATED_PRESERVATION" if n==len(rows)==48 and k==len(can)==36 else "HOLD_NATIVE_PRESERVATION","sources":sources,"comparisons":rows,"canaries":can,"preserved":n,"transform_count":len(rows),"canaries_correct":k,"canary_count":len(can),"selftest":selftest_result,"full_integrity_matches":sum("integrity" in r for r in sources),"originals_unchanged":all(r["status"]=="PROCESSED" for r in sources),"native_grid_tested":len(rows)==48,"end_to_end_hannah_tested":False,"exact_bl_frame_verified":False,"spatial_off_verified":False,"production_adoption":False}
    save("result.json",d)
    # Aggregate log; no fine spectral records.
    print("R37F3B_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
