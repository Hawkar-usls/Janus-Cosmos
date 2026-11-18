#!/usr/bin/env python3
"""W3 physical frequency metadata audit; never shift or interpolate power."""
import json,pathlib
import numpy as np
from astropy.io import fits
import r37f3d2b_exposure as prior
import r37f3d2_archive_inventory as inv
O=pathlib.Path("runtime/r37f3d2f");O.mkdir(parents=True,exist_ok=True)
def plain(x):
    if isinstance(x,(bytes,np.bytes_)):return x.decode(errors="replace").strip()
    if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
    if isinstance(x,(int,np.integer)):return int(x)
    return str(x).strip()
def main():
    p=O/"original.fits";d={"status":"HOLD","qualified":False,"production_adoption":False,"power_read":False,"interpolated_samples":0}
    try:
        src=next(x for x in inv.S["sources"] if x["id"]=="w3_1");d["integrity"]=inv.download(src,p);assert d["integrity"]["sha256"]==prior.S["sha256"]
        p.chmod(0o444);m=prior.metadata(p);ix=np.asarray(m["kept_row_indices"]);assert ix.shape==(4,30)
        with fits.open(p,memmap=True) as hd:
            h=next(h for h in hd if isinstance(h,fits.BinTableHDU) and "DATA" in h.columns.names);tab=h.data
            keys=["CTYPE1","CUNIT1","SPECSYS","SSYSOBS","VELDEF","VELFRAME","VFRAME","RVSYS","VELOCITY","RESTFREQ","OBSFREQ","CRVAL1","CRPIX1","CDELT1","FREQRES","BANDWID","LST","DATE-OBS","CRVAL2","CRVAL3","CTYPE2","CTYPE3","EQUINOX","RADESYS","OBJECT","OBSMODE","SIG","CAL","FDNUM","SAMPLER","FEED","IFNUM","PLNUM","LO1","LO2","LO1FREQ","LO2FREQ","SIDEBAND"]
            d["column_names"]=list(tab.names)
            d["header"]= {k:plain(h.header[k]) for k in keys if k in h.header}
            d["primary_header"]={k:plain(hd[0].header[k]) for k in keys if k in hd[0].header}
            d["phase_metadata"]={k:[[plain(x) for x in tab[k][ids]] for ids in ix] for k in keys if k in tab.names}
            v=np.asarray(tab["CRVAL1"][ix],float);pix=np.asarray(tab["CRPIX1"][ix],float);df=np.asarray(tab["CDELT1"][ix],float)
            f=v+(1-pix)*df;off=f[0]-f[2];frac=off/df[0]-np.rint(off/df[0])
            assert (np.abs(frac)>.001).all()
            d["frequency_offsets_hz"]={"min":float(off.min()),"max":float(off.max()),"values":off.tolist(),"phase_pair_max":float(max(np.max(abs(f[0]-f[1])),np.max(abs(f[2]-f[3]))))}
            d["fractional_residual_channels"]={"min":float(frac.min()),"max":float(frac.max())}
            d["candidate_explanations"]={}
            # Diagnostic differences only: these are not proof of a reference frame transformation.
            for k in ("OBSFREQ","RESTFREQ","VFRAME","RVSYS","VELOCITY","LO1","LO2","LO1FREQ","LO2FREQ"):
                if k not in tab.names:continue
                a=np.asarray(tab[k][ix],float);diff=a[0]-a[2]
                if np.isfinite(diff).all():
                    d["candidate_explanations"][k]={"unit":str(h.columns[k].unit),"min_difference":float(diff.min()),"max_difference":float(diff.max())}
                    if k in ("OBSFREQ","RESTFREQ","LO1","LO2","LO1FREQ","LO2FREQ"):
                        d["candidate_explanations"][k]["max_abs_residual_vs_native_hz"]=float(np.max(abs(off-diff)))
            d["status"]="PASS_PHYSICAL_METADATA_INVENTORY__HOLD_FRAME_CAUSAL_ATTRIBUTION"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:
        if p.exists():
            d["original_unchanged"]=prior.sha(p)==prior.S["sha256"];p.unlink()
            if not d["original_unchanged"]:d["status"]="HOLD_INTEGRITY"
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");print("R37F3D2F_SUMMARY",json.dumps({k:v for k,v in d.items() if k!="phase_metadata"},allow_nan=False),flush=True)
if __name__=="__main__":main()
