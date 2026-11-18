#!/usr/bin/env python3
"""Historical natural astronomy calibration/averaging diagnostic, aggregate results only."""
import hashlib,importlib.metadata,inspect,json,pathlib,subprocess,sys
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord,EarthLocation,solar_system_ephemeris
from astropy.time import Time
from astropy.io import fits
from astropy.utils import iers
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.spectra import scan,spectrum
P=pathlib.Path("research/seti/r37f3d1bb_average_protocol.json")
S=json.loads(P.read_text());M=S["method"];O=pathlib.Path("runtime/r37f3d1bb");O.mkdir(parents=True,exist_ok=True)
for n,h in S["pinned_blobs"].items():
    b=(P.parent/n).read_bytes()
    assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==h,n
import r37f3d1a_metadata as meta
C=299792458.
site=M["physical_site"]
loc=EarthLocation.from_geodetic(site["lon_deg"]*u.deg,site["lat_deg"]*u.deg,site["height_m"]*u.m)
def hashfile(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(4*1024**2),b""):h.update(b)
    return h.hexdigest()
def text(x):return str(x).strip()
def framecheck(p):
    out={}
    with fits.open(p,memmap=True) as hd:
        d=next(h.data for h in hd if isinstance(h,fits.BinTableHDU) and h.header.get("EXTNAME")=="SINGLE DISH")
        for s in (152,153):
            m=(d["SCAN"]==s)&(d["IFNUM"]==0)&(d["PLNUM"]==0)&(np.char.strip(d["CAL"].astype(str))=="T")
            assert int(m.sum())==151
            t=Time(np.char.strip(d["DATE-OBS"][m].astype(str)),format="isot",scale="utc")
            coord=SkyCoord(ra=d["CRVAL2"][m]*u.deg,dec=d["CRVAL3"][m]*u.deg,frame="fk5",equinox="J2000")
            vc=coord.radial_velocity_correction(kind="heliocentric",obstime=t,location=loc).to_value(u.m/u.s)
            vf=np.asarray(d["VFRAME"][m],float)
            assert np.isfinite(vc).all() and np.isfinite(vf).all()
            out[str(s)]={"rows":len(vc),"negation_max_abs_residual_m_s":float(np.max(abs(vf+vc))),"same_sign_max_abs_residual_m_s":float(np.max(abs(vf-vc))),"pass":bool(np.max(abs(vf+vc))<=M["frame_tolerance_m_s"])}
    return out
def code_audit():
    records=[]
    for cls in (scan.PSScan,scan.ScanBlock):
        f=inspect.unwrap(cls.timeaverage);code=inspect.getsource(f)
        records.append({"method":f.__qualname__,"sha256":hashlib.sha256(code.encode()).hexdigest(),"relevant_lines":[x.strip() for x in code.splitlines() if any(k in x.lower() for k in ("getspec(","average_spectra(","average(","align","shift","_calibrated","use_wcs"))]})
    return records
def get_arrays(spec):
    y=np.asarray(spec.flux.to_value(u.K),float)
    f=np.asarray(spec.spectral_axis.to_value(u.Hz),float)
    mask=np.asarray(spec.mask,bool) if spec.mask is not None else np.zeros(y.shape,bool)
    assert y.shape==f.shape==mask.shape==(M["channels"],)
    good=np.isfinite(y)&np.isfinite(f)&~mask
    assert good.all(),"masked or nonfinite channels"
    assert np.max(np.abs(abs(np.diff(f))-715.2557373046875))<=M["grid_tolerance_hz"]
    return y.copy(),f.copy()
def rms(x):return float(np.sqrt(np.mean(np.square(x))))
def measure(sdf,q,frame_ok):
    block=sdf.getps(scan=M["scan"],ifnum=M["ifnum"],plnum=q,fdnum=M["fdnum"])
    assert len(block)==1
    ps=block[0];assert int(ps.nint)==M["rows"]
    ys=[];fs=[];corr=[]
    for i in range(M["rows"]):
        sp=ps.getspec(i,use_wcs=True);y,f=get_arrays(sp);ys.append(y);fs.append(f)
        if frame_ok:
            t=Time(text(sp.meta["DATE-OBS"]),format="isot",scale="utc")
            coord=SkyCoord(ra=float(sp.meta["CRVAL2"])*u.deg,dec=float(sp.meta["CRVAL3"])*u.deg,frame="fk5",equinox="J2000")
            vc=float(coord.radial_velocity_correction(kind="heliocentric",obstime=t,location=loc).to_value(u.m/u.s))
            assert abs(float(sp.meta["VFRAME"])+vc)<=M["frame_tolerance_m_s"]
            corr.append(np.sqrt((1-vc/C)/(1+vc/C)))
    y=np.stack(ys);f=np.stack(fs)
    avg,axis=get_arrays(ps.timeaverage(weights=None,use_wcs=True))
    manual=y.mean(axis=0);err=float(np.max(abs(avg-manual)))
    # Equal weighting is intentionally diagnostic, not the official Tsys output.
    # A floating point allowance is fixed, independent of result.
    index_match=bool(np.allclose(avg,manual,rtol=1e-5,atol=1e-5))
    firsterr=float(np.max(abs(axis-f[0])))
    temporal=rms(y-manual[None,:])
    d={"rows":len(y),"channels":y.shape[1],"dysh_index_average_matches":index_match,"index_average_max_abs_K":err,"dysh_axis_first_max_error_hz":firsterr,"first_axis_matches":firsterr<=M["grid_tolerance_hz"],"temporal_residual_rms_K":temporal,"temporal_variance_fraction":float(np.mean((y-manual[None,:])**2)/max(np.mean(y*y),1e-30)),"frame_diagnostic_performed":False}
    assert index_match
    if frame_ok:
        hf=f*np.asarray(corr)[:,None]
        # Restrict to common interior support. Interpolation is a labelled
        # diagnostic derivative, not a replacement for any native observation.
        if hf[0,0]>hf[0,-1]:hf=hf[:,::-1];y=y[:,::-1];manual=manual[::-1]
        assert np.all(np.diff(hf,axis=1)>0)
        lo=float(np.max(hf[:,0]));hi=float(np.min(hf[:,-1]));keep=(hf[0]>=lo)&(hf[0]<=hi)
        target=hf[0,keep];assert len(target)>=M["channels"]-4
        z=np.stack([np.interp(target,a,b) for a,b in zip(hf,y)])
        aligned=z.mean(0);baseline=manual[keep]
        d.update(frame_diagnostic_performed=True,common_channels=len(target),discarded_edge_channels=M["channels"]-len(target),aligned_vs_index_rms_K=rms(aligned-baseline),aligned_vs_index_relative_rms=rms(aligned-baseline)/max(rms(baseline),1e-30),frame_aligned_derivative_only=True,temporal_residual_aligned_rms_K=rms(z-aligned[None,:]))
    return d
def main():
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
    p=O/"original.fits";d={"status":"HOLD_AVERAGE_DIAGNOSTIC","source":"NGC2415","production_adoption":False,"detector_tested":False,"lossless_average_proven":False,"official_tsys_calibration_qualified":False,"rows":[]}
    try:
        assert importlib.metadata.version("dysh")==S["dysh"]
        d["source_audit"]=code_audit()
        d["integrity"]=meta.prev.download(p)
        assert d["integrity"]["sha256"]==S["archive_sha256"] and d["integrity"]["bytes"]==S["archive_bytes"]
        p.chmod(0o444);m={};meta.inspect(p,m);assert len(m["gates"])==14 and all(m["gates"].values())
        d["spatial_off_verified"]=True;d["frame"]=framecheck(p);ok=all(x["pass"] for x in d["frame"].values())
        d["frame_sign_consistent_with_independent_heliocentric_projection"]=ok
        sdf=GBTFITSLoad(p)
        for q in M["pols"]:
            row={"pol":q,"status":"HOLD"}
            try:
                row.update(measure(sdf,q,ok));row["status"]="DIAGNOSTIC_COMPLETE"
            except Exception as e:
                import traceback;traceback.print_exc();row["error"]=type(e).__name__+": "+str(e)[:500]
            d["rows"].append(row);print("R37F3D1BB_POL",json.dumps(row),flush=True)
        d["original_unchanged"]=hashfile(p)==S["archive_sha256"]
        d["status"]="PASS_ISOLATED_AVERAGE_DIAGNOSTIC" if d["original_unchanged"] and ok and len(d["rows"])==2 and all(r["status"]=="DIAGNOSTIC_COMPLETE" for r in d["rows"]) else "HOLD_AVERAGE_DIAGNOSTIC"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:p.unlink(missing_ok=True)
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    (O/"pip-freeze.txt").write_text(subprocess.check_output([sys.executable,"-m","pip","freeze"],text=True))
    print("R37F3D1BB_SUMMARY",json.dumps(d,allow_nan=False),flush=True)
if __name__=="__main__":main()
