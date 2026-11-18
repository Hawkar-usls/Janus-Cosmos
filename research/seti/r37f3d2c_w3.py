#!/usr/bin/env python3
"""W3_1 isolated arithmetic and native WCS; no grid harmonization."""
import json,pathlib,importlib.metadata
import numpy as np
from astropy.io import fits
from astropy.utils import iers
from astropy.coordinates import solar_system_ephemeris
from dysh.fits.gbtfitsload import GBTFITSLoad
import r37f3d2b_exposure as prior
import r37f3d2_archive_inventory as inv
P=pathlib.Path("research/seti/r37f3d2c_w3_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f3d2c");O.mkdir(parents=True,exist_ok=True)
def close(a,b,atol,rtol):
    a=np.asarray(a,float);b=np.asarray(b,float);assert a.shape==b.shape and a.size and np.isfinite(a).all() and np.isfinite(b).all()
    err=float(np.max(abs(a-b)));assert np.allclose(a,b,atol=atol,rtol=rtol),str(err);return err
def measure(p,d):
    meta=prior.metadata(p);ixs=np.asarray(meta["kept_row_indices"]);assert ixs.shape==(4,30)
    sdf=GBTFITSLoad(p);bl=sdf.getps(scan=60,fdnum=0,ifnum=0,plnum=0,apply_flags=True);assert len(bl)==1
    ps=bl[0];assert ps.nint==30 and ps.nchan==32768 and not ps._nocal and ps._smoothref==1 and ps.tscale.lower()=="ta"
    arrays=[];masks=[];exp=[];df=[];tc=[];grids=[]
    with fits.open(p,memmap=True) as hd:
        h=next(h for h in hd if isinstance(h,fits.BinTableHDU) and "DATA" in h.columns.names);tab=h.data
        for j,(a,r) in enumerate((("sigcalon","sigonrows"),("sigcaloff","sigoffrows"),("refcalon","refonrows"),("refcaloff","refoffrows"))):
            ix=np.asarray(getattr(ps,"_"+r),int);assert np.array_equal(ix,ixs[j]),a+" row mismatch"
            raw=np.stack([np.asarray(tab["DATA"][i],float).ravel() for i in ix]);lib=np.ma.asarray(getattr(ps,"_"+a));m=np.ma.getmaskarray(lib)
            assert raw.shape==m.shape==(30,32768)
            close(raw[~m],lib.data[~m],0,0)
            arrays.append(raw);masks.append(m);exp.append(np.asarray(tab["EXPOSURE"][ix],float));df.append(np.asarray(tab["CDELT1"][ix],float));tc.append(np.asarray(tab["TCAL"][ix],float))
            grids.append(np.asarray(tab["CRVAL1"][ix],float)[:,None]+(np.arange(32768)[None,:]+1-np.asarray(tab["CRPIX1"][ix],float)[:,None])*df[-1][:,None])
    d["row_ids_match"]=True;d["raw_unmasked_match"]=True;d["input_masks"]=[int(m.sum()) for m in masks]
    edge=int(32768*.1);sl=slice(edge,-(edge-1));ts=[]
    assert np.array_equal(tc[2],tc[3]) and (tc[2]>0).all()
    for i in range(30):
        off=arrays[3][i,sl];on=arrays[2][i,sl];mo=~masks[3][i,sl];md=mo&~masks[2][i,sl]
        assert mo.any() and md.any();diff=np.mean((on-off)[md]);assert diff!=0
        ts.append(float(tc[3][i]*(np.mean(off[mo])/diff+.5)))
    ts=np.asarray(ts);d["tsys_error_K"]=close(ps.tsys,ts,1e-9,1e-9)
    e1=exp[0]+exp[1];e2=exp[2]+exp[3];e=e1*e2/(e1+e2);width=abs(sum(df)/4)
    d["exposure_error_s"]=close(ps.exposure,e,1e-10,1e-10);d["width_error_hz"]=close(abs(ps.delta_freq),width,1e-10,1e-10)
    w=e*width/ts**2;d["weight_error"]=close(ps.tsys_weight,w,1e-10,1e-9)
    sig=(arrays[0]+arrays[1])/2;ref=(arrays[2]+arrays[3])/2
    data=np.ma.asarray(ps._calibrated);g=~np.ma.getmaskarray(data)&np.isfinite(data.data)
    expect=~np.logical_or.reduce(masks)&np.isfinite(ref)&(ref!=0)&np.isfinite(sig)
    d["mask_exact"]=bool(np.array_equal(g,expect));assert d["mask_exact"]
    d["calibrated_error_K"]=close(data.data[g],np.broadcast_to(ts[:,None],sig.shape)[g]*(sig[g]-ref[g])/ref[g],.0001,.000001)
    d["valid_calibrated_samples"]=int(g.sum())
    axes=[]
    for i in range(30):
        sp=ps.getspec(i,use_wcs=True);axes.append(sp.spectral_axis.to_value("Hz"))
    d["native_wcs_error_hz"]=close(np.stack(axes),grids[1],.01,0)
    delta=(grids[0][:,0]-grids[2][:,0])/df[0]
    d["fractional_grid"]={"offset_channels_min":float(delta.min()),"offset_channels_max":float(delta.max()),"max_residual_channels":float(np.max(abs(delta-np.rint(delta)))),"integer_common_grid":bool((abs(delta-np.rint(delta))<.001).all())}
    d["imputed_samples"]=0;d["interpolated_samples"]=0
def main():
    d={"status":"HOLD","qualified":False,"production_adoption":False,"detector_tested":False};p=O/"original.fits"
    try:
        assert importlib.metadata.version("dysh")=="1.0.1";iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
        src=next(x for x in inv.S["sources"] if x["id"]=="w3_1");d["integrity"]=inv.download(src,p);assert d["integrity"]["sha256"]==S["source_sha256"];p.chmod(0o444)
        measure(p,d);d["status"]="PASS_ISOLATED_ARITHMETIC__HOLD_PHYSICAL_GRID"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:
        if p.exists():
            d["original_unchanged"]=prior.sha(p)==S["source_sha256"];p.unlink()
            if not d["original_unchanged"]:d["status"]="HOLD_INTEGRITY"
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");print("R37F3D2C_SUMMARY",json.dumps(d,allow_nan=False),flush=True)
if __name__=="__main__":main()
