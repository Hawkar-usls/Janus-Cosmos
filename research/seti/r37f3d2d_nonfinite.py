#!/usr/bin/env python3
"""R37F3D2D: distinguish nonfinite archive pixels, library mask and value mismatch."""
import json,pathlib,importlib.metadata
import numpy as np
from astropy.io import fits
from astropy.utils import iers
from astropy.coordinates import solar_system_ephemeris
from dysh.fits.gbtfitsload import GBTFITSLoad
import r37f3d2b_exposure as prior
import r37f3d2_archive_inventory as inv
O=pathlib.Path("runtime/r37f3d2d");O.mkdir(parents=True,exist_ok=True)
def main():
    p=O/"original.fits";d={"status":"HOLD","rows":[],"qualified":False,"production_adoption":False,"imputed_samples":0,"interpolated_samples":0}
    try:
        assert importlib.metadata.version("dysh")=="1.0.1";iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
        src=next(x for x in inv.S["sources"] if x["id"]=="w3_1");d["integrity"]=inv.download(src,p);assert d["integrity"]["sha256"]==prior.S["sha256"];p.chmod(0o444)
        m=prior.metadata(p);ixs=np.asarray(m["kept_row_indices"])
        sdf=GBTFITSLoad(p);ps=sdf.getps(scan=60,fdnum=0,ifnum=0,plnum=0,apply_flags=True)[0];assert ps.nint==30 and ps.nchan==32768
        masks=[];arrays=[]
        with fits.open(p,memmap=True) as hd:
            h=next(h for h in hd if isinstance(h,fits.BinTableHDU) and "DATA" in h.columns.names);tab=h.data
            for j,(a,r) in enumerate((("sigcalon","sigonrows"),("sigcaloff","sigoffrows"),("refcalon","refonrows"),("refcaloff","refoffrows"))):
                ix=np.asarray(getattr(ps,"_"+r),int);assert np.array_equal(ix,ixs[j])
                raw=np.stack([np.asarray(tab["DATA"][i],float).ravel() for i in ix]);lib=np.ma.asarray(getattr(ps,"_"+a));mask=np.ma.getmaskarray(lib);v=np.asarray(lib.data,float)
                assert raw.shape==v.shape==mask.shape==(30,32768)
                bad=~np.isfinite(raw);lb=~np.isfinite(v);flags=np.stack([np.asarray(tab["FLAGS"][i],bool).ravel() for i in ix]) if "FLAGS" in tab.names else np.zeros_like(mask)
                safe=~mask&~bad&~lb;diff=int(np.sum(raw[safe]!=v[safe]))
                mismatch=int(np.sum((bad!=lb)&~mask));unmasked=bad&~mask
                row={"phase":a,"raw_nan":int(np.isnan(raw).sum()),"raw_inf":int(np.isinf(raw).sum()),"raw_bad_unmasked":int(unmasked.sum()),"library_nonfinite_unmasked":int((lb&~mask).sum()),"nonfinite_support_mismatch":mismatch,"mask_count":int(mask.sum()),"archive_flags":int(flags.sum()),"raw_bad_flagged":int((bad&flags).sum()),"finite_compared":int(safe.sum()),"unequal_finite":diff,"unmasked_bad_by_raw_row":[{"row":int(i),"count":int(k)} for i,k in zip(ix,unmasked.sum(1)) if k]}
                d["rows"].append(row);masks.append(mask|bad|lb);arrays.append(raw)
                print("R37F3D2D_PHASE",json.dumps(row),flush=True)
        sig=(arrays[0]+arrays[1])/2;ref=(arrays[2]+arrays[3])/2
        expected=~np.logical_or.reduce(masks)&np.isfinite(sig)&np.isfinite(ref)&(ref!=0)
        data=np.ma.asarray(ps._calibrated);finite=~np.ma.getmaskarray(data)&np.isfinite(data.data)
        d["calibrated_support"]={"expected":int(expected.sum()),"official_finite":int(finite.sum()),"official_extra":int((finite&~expected).sum()),"official_missing":int((expected&~finite).sum()),"exact_match":bool(np.array_equal(expected,finite))}
        d["status"]="PASS_NONFINITE_DIAGNOSIS__HOLD_ARITHMETIC_GRID"
        assert all(x["unequal_finite"]==0 and x["nonfinite_support_mismatch"]==0 for x in d["rows"])
    except Exception as e:
        import traceback;traceback.print_exc();d["status"]="HOLD_DIAGNOSTIC";d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:
        if p.exists():
            d["original_unchanged"]=prior.sha(p)==prior.S["sha256"];p.unlink()
            if not d["original_unchanged"]:d["status"]="HOLD_INTEGRITY"
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");print("R37F3D2D_SUMMARY",json.dumps(d,allow_nan=False),flush=True)
if __name__=="__main__":main()
