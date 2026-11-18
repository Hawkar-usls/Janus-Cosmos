#!/usr/bin/env python3
"""Mask-preserving historical natural astronomical calibration diagnostic."""
import hashlib,importlib.metadata,json,pathlib,subprocess,sys
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.utils import iers
from astropy.coordinates import solar_system_ephemeris
from dysh.fits.gbtfitsload import GBTFITSLoad
P=pathlib.Path("research/seti/r37f3d1bb1_mask_protocol.json")
S=json.loads(P.read_text());M=S["method"];O=pathlib.Path("runtime/r37f3d1bb1");O.mkdir(parents=True,exist_ok=True)
for n,h in S["pinned_blobs"].items():
    b=(P.parent/n).read_bytes()
    assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==h,n
import r37f3d1bb_average as prev
def rms(x):
    assert np.size(x) and np.isfinite(x).all()
    return float(np.sqrt(np.mean(np.square(x))))
def counts(x):
    x=np.asarray(x,dtype=np.int64)
    return {"total":int(x.sum()),"min_per_row":int(x.min()),"max_per_row":int(x.max()),"rows_affected":int(np.count_nonzero(x))}
def audit_raw(p):
    rows=[];ref={}
    with fits.open(p,memmap=True) as hd:
        d=next(h.data for h in hd if isinstance(h,fits.BinTableHDU) and h.header.get("EXTNAME")=="SINGLE DISH")
        cal=np.char.strip(d["CAL"].astype(str))
        dates=np.char.strip(d["DATE-OBS"].astype(str))
        for s in (152,153):
            for q in M["pols"]:
                for phase in ("T","F"):
                    ix=np.flatnonzero((d["SCAN"]==s)&(d["IFNUM"]==M["ifnum"])&(d["PLNUM"]==q)&(d["FDNUM"]==0)&(cal==phase))
                    assert len(ix)==M["rows"]
                    bad=[];zero=[]
                    for k in ix:
                        v=np.asarray(d["DATA"][k]).ravel()
                        assert len(v)==M["channels"]
                        bad.append(np.count_nonzero(~np.isfinite(v)))
                        zero.append(np.count_nonzero(v==0))
                        if s==152 and phase=="T":
                            ts=str(dates[k]);vf=float(d["VFRAME"][k])
                            assert (q,ts) not in ref and np.isfinite(vf)
                            ref[q,ts]=vf
                    rows.append({"scan":s,"pol":q,"cal":phase,"rows":len(ix),"nonfinite":counts(bad),"zero":counts(zero),"zeros_are_not_assumed_invalid":True})
    return rows,ref
def unpack(spec):
    raw=spec.flux.to_value(u.K)
    y=np.asarray(raw,dtype=float)
    f=np.asarray(spec.spectral_axis.to_value(u.Hz),dtype=float)
    m=np.ma.getmaskarray(raw)|np.broadcast_to(np.asarray(spec.mask,dtype=bool) if spec.mask is not None else False,y.shape)
    assert y.shape==f.shape==m.shape==(M["channels"],)
    assert np.isfinite(f).all()
    assert np.max(np.abs(abs(np.diff(f))-715.2557373046875))<=M["grid_atol_hz"]
    good=np.isfinite(y)&~m
    return y.copy(),f.copy(),good,{"masked":int(m.sum()),"nonfinite":int(np.count_nonzero(~np.isfinite(y))),"invalid_union":int(np.count_nonzero(~good))}
def local_values(target,f,y,good):
    """No extrapolation, no interpolation across a masked neighbour."""
    assert np.all(np.diff(f)>0)
    k=np.searchsorted(f,target,side="left")
    hit=(k<len(f))&(f[np.minimum(k,len(f)-1)]==target)
    result=np.full(target.shape,np.nan)
    exact=np.flatnonzero(hit)
    ix=exact[good[k[exact]]]
    result[ix]=y[k[ix]]
    interior=np.flatnonzero(~hit&(k>0)&(k<len(f)))
    j=interior[good[k[interior]-1]&good[k[interior]]]
    l=k[j]-1;r=k[j]
    t=(target[j]-f[l])/(f[r]-f[l])
    result[j]=(1-t)*y[l]+t*y[r]
    return result,np.isfinite(result)
def selftest():
    f=np.arange(5,dtype=float);y=np.array([1.,2.,np.nan,4.,5.]);g=np.isfinite(y)
    t=np.array([-1.,0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,5.])
    z,ok=local_values(t,f,y,g)
    expected=np.array([False,True,True,True,False,False,False,True,True,True,False])
    assert np.array_equal(ok,expected)
    assert np.allclose(z[ok],[1.,1.5,2.,4.,4.5,5.])
    return True
def measure(sdf,q,ref,row):
    block=sdf.getps(scan=M["scan"],ifnum=M["ifnum"],plnum=q,fdnum=M["fdnum"])
    assert len(block)==1
    ps=block[0];assert int(ps.nint)==M["rows"]
    ys=[];fs=[];gs=[];cs=[];corr=[];used=set()
    for i in range(M["rows"]):
        sp=ps.getspec(i,use_wcs=True);y,f,g,c=unpack(sp)
        ys.append(y);fs.append(f);gs.append(g);cs.append(c)
        ts=prev.text(sp.meta["DATE-OBS"])
        vf=ref[q,ts];assert ts not in used;used.add(ts)
        corr.append(np.sqrt((1+vf/prev.C)/(1-vf/prev.C)))
    y=np.stack(ys);f=np.stack(fs);g=np.stack(gs)
    fully=g.all(axis=0);anyvalid=g.any(axis=0);n=int(fully.sum())
    row.update(calibrated_integrations=len(y),channels=M["channels"],calibrated_invalid={k:counts([c[k] for c in cs]) for k in cs[0]},
               valid_every_integration=n,valid_any_integration=int(anyvalid.sum()),invalid_every_integration=int((~anyvalid).sum()),
               complete_support_fraction=n/M["channels"])
    assert n>=M["min_common_channels"],"no common valid index support"
    avg,axis,agood,ac=unpack(ps.timeaverage(weights=None,use_wcs=True))
    common=fully&agood
    row["dysh_average_invalid"]=ac;row["index_comparison_channels"]=int(common.sum())
    assert np.array_equal(common,fully),"dysh dropped fully valid channels"
    manual=y[:,common].mean(axis=0)
    assert np.all(np.isfinite(manual))
    err=avg[common]-manual
    row.update(index_average_max_abs_K=float(np.max(abs(err))),index_average_rms_K=rms(err),
               index_average_matches=bool(np.allclose(avg[common],manual,atol=M["index_atol_K"],rtol=M["index_rtol"])),
               first_axis_max_error_hz=float(np.max(abs(axis-f[0]))),
               temporal_residual_rms_K=rms(y[:,common]-manual[None,:]),
               temporal_residual_energy_fraction=float(np.mean((y[:,common]-manual[None,:])**2)/max(np.mean(y[:,common]**2),1e-30)))
    assert row["index_average_matches"]
    # Conditional physical-grid derivative; never replace native data.
    hf=f*np.asarray(corr)[:,None]
    if hf[0,0]>hf[0,-1]:
        hf=hf[:,::-1];y=y[:,::-1];g=g[:,::-1];common=common[::-1]
    assert np.all(np.diff(hf,axis=1)>0)
    target=hf[0];zs=[];oks=[]
    for a,b,c in zip(hf,y,g):
        z,ok=local_values(target,a,b,c)
        zs.append(z);oks.append(ok)
    z=np.stack(zs);valid=np.stack(oks).all(axis=0)
    joint=valid&common
    row.update(aligned_common_channels=int(valid.sum()),joint_comparison_channels=int(joint.sum()),
               joint_coverage_fraction=float(joint.mean()),alignment_conditional=True)
    assert int(joint.sum())>=M["min_common_channels"],"no valid common physical/index support"
    aligned=z[:,joint].mean(0);baseline=y[:,joint].mean(0)
    row.update(aligned_vs_index_rms_K=rms(aligned-baseline),
               aligned_vs_index_relative_rms=rms(aligned-baseline)/max(rms(baseline),1e-30),
               aligned_temporal_residual_rms_K=rms(z[:,joint]-aligned[None,:]),
               imputed_samples=0,interpolated_across_gaps=0,extrapolated_samples=0)
def main():
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
    p=O/"original.fits";d={"status":"HOLD_MASK_SUPPORT_DIAGNOSTIC","source":"NGC2415","production_adoption":False,"detector_tested":False,"lossless_average_proven":False,"official_tsys_calibration_qualified":False,"rows":[]}
    try:
        d["synthetic_mask_gap_check"]=selftest()
        assert importlib.metadata.version("dysh")=="1.0.1"
        d["integrity"]=prev.meta.prev.download(p)
        assert d["integrity"]["sha256"]==prev.S["archive_sha256"] and d["integrity"]["bytes"]==prev.S["archive_bytes"]
        p.chmod(0o444);m={};prev.meta.inspect(p,m);assert len(m["gates"])==14 and all(m["gates"].values())
        d["spatial_off_verified"]=True;d["frame"]=prev.framecheck(p);assert all(v["pass"] for v in d["frame"].values())
        d["raw_support"],ref=audit_raw(p)
        sdf=GBTFITSLoad(p)
        for q in M["pols"]:
            row={"pol":q,"status":"HOLD"}
            try:measure(sdf,q,ref,row);row["status"]="DIAGNOSTIC_COMPLETE"
            except Exception as e:
                import traceback;traceback.print_exc();row["error"]=type(e).__name__+": "+str(e)[:500]
            d["rows"].append(row);print("R37F3D1BB1_POL",json.dumps(row),flush=True)
        d["original_unchanged"]=prev.hashfile(p)==prev.S["archive_sha256"]
        d["status"]="PASS_MASK_SUPPORT_DIAGNOSTIC_ONLY" if d["original_unchanged"] and all(r["status"]=="DIAGNOSTIC_COMPLETE" for r in d["rows"]) and len(d["rows"])==2 else "HOLD_MASK_SUPPORT_DIAGNOSTIC"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:p.unlink(missing_ok=True)
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    (O/"pip-freeze.txt").write_text(subprocess.check_output([sys.executable,"-m","pip","freeze"],text=True))
    print("R37F3D1BB1_SUMMARY",json.dumps(d,allow_nan=False),flush=True)
if __name__=="__main__":main()
