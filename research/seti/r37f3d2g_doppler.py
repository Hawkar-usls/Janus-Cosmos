#!/usr/bin/env python3
"""W3 metadata-only Doppler differential, no power or grid resampling."""
import json,pathlib
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord,EarthLocation,ICRS,LSRK,CartesianDifferential,solar_system_ephemeris
from astropy.time import Time
from astropy.utils import iers
import astropy.units as u
import r37f3d2b_exposure as prior
import r37f3d2_archive_inventory as inv
P=pathlib.Path("research/seti/r37f3d2g_doppler_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f3d2g");O.mkdir(parents=True,exist_ok=True)
C=299792458.
def stats(x):
    x=np.asarray(x,float);assert np.isfinite(x).all()
    return {"min":float(x.min()),"max":float(x.max()),"max_abs":float(abs(x).max()),"rms":float(np.sqrt(np.mean(x*x)))}
def main():
    p=O/"original.fits";d={"status":"HOLD","qualified":False,"production_adoption":False,"power_read":False,"interpolated_samples":0}
    try:
        src=next(x for x in inv.S["sources"] if x["id"]=="w3_1");d["integrity"]=inv.download(src,p);assert d["integrity"]["sha256"]==S["source_sha256"];p.chmod(0o444)
        m=prior.metadata(p);ix=np.asarray(m["kept_row_indices"])
        with fits.open(p,memmap=True) as hd:
            h=next(h for h in hd if isinstance(h,fits.BinTableHDU) and "DATA" in h.columns.names);tab=h.data
            def num(k):return np.asarray(tab[k][ix],float)
            def st(k):return np.char.upper(np.char.strip(np.asarray(tab[k][ix]).astype(str)))
            d["definitions"]={k:sorted(set(st(k).ravel().tolist())) for k in ("VELDEF","CTYPE1","CTYPE2","CTYPE3","RADESYS","SIDEBAND","OBSMODE","PROCTYPE","OBJECT")}
            d["numeric_ranges"]={k:stats(num(k)) for k in ("DOPFREQ","RESTFREQ","VFRAME","RVSYS","VELOCITY","OBSFREQ","CRVAL1","CRPIX1","CDELT1","FOFFREF1","SITELONG","SITELAT","SITEELEV","TRGTLONG","TRGTLAT","CRVAL2","CRVAL3")}
            observed=num("OBSFREQ");v=num("VELOCITY");vf=num("VFRAME");f0=num("DOPFREQ");rest=num("RESTFREQ");delta=lambda a:a[0]-a[2]
            assert np.isfinite(f0).all() and (f0>0).all()
            pred=f0*(1-(v+vf)/C);p2=rest*(1-(v+vf)/C)
            beta=vf/C;rel=f0*(1-v/C)*np.sqrt((1-beta)/(1+beta))
            d["observed_delta_hz"]=stats(delta(observed))
            d["radio_dopfreq_residual_hz"]=stats(delta(observed)-delta(pred));d["radio_restfreq_residual_hz"]=stats(delta(observed)-delta(p2));d["relativistic_frame_residual_hz"]=stats(delta(observed)-delta(rel))
            d["absolute_radio_residual_hz"]=stats(observed-pred)
            d["radio_definition"]=all(k.startswith(("RADI","VRAD")) for k in d["definitions"]["VELDEF"])
            d["radio_differential_pass"]=d["radio_definition"] and d["radio_dopfreq_residual_hz"]["max_abs"]<=5
            # Independent projected-motion diagnostic; only use explicitly equatorial coordinates.
            iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
            d["frame_motion"]={"status":"HOLD_UNSUPPORTED_OR_AMBIGUOUS_COORDINATES"}
            defs=d["definitions"]
            if all(k.startswith("RA") for k in defs["CTYPE2"]) and all(k.startswith("DEC") for k in defs["CTYPE3"]) and all(k in ("FK5","ICRS") for k in defs["RADESYS"]):
                # Site longitude in east-positive convention per GBT SDFITS; report raw values.
                loc=EarthLocation.from_geodetic(num("SITELONG")[0]*u.deg,num("SITELAT")[0]*u.deg,num("SITEELEV")[0]*u.m)
                for n,ra,dec in (("pointing",num("CRVAL2"),num("CRVAL3")),("target",num("TRGTLONG"),num("TRGTLAT"))):
                    if not (np.isfinite(ra).all() and np.isfinite(dec).all() and (abs(dec)<=90).all()):continue
                    tt=[];bb=[]
                    for j in (0,2):
                        t=Time(st("DATE-OBS")[j].tolist(),scale="utc");site=EarthLocation.from_geodetic(num("SITELONG")[j]*u.deg,num("SITELAT")[j]*u.deg,num("SITEELEV")[j]*u.m)
                        coord=SkyCoord(ra=ra[j]*u.deg,dec=dec[j]*u.deg,frame="fk5" if defs["RADESYS"]==["FK5"] else "icrs",**({"equinox":"J2000"} if defs["RADESYS"]==["FK5"] else {}))
                        bary=coord.radial_velocity_correction(obstime=t,location=site,kind="barycentric").to_value(u.m/u.s)
                        # LSRK additive projection from static fictitious zero-velocity ICRS source.
                        ic=coord.icrs
                        zero=SkyCoord(ra=ic.ra,dec=ic.dec,distance=np.ones(len(ic))*u.pc,pm_ra_cosdec=np.zeros(len(ic))*u.mas/u.yr,pm_dec=np.zeros(len(ic))*u.mas/u.yr,radial_velocity=np.zeros(len(ic))*u.km/u.s,frame="icrs")
                        lsr=zero.transform_to(LSRK()).radial_velocity.to_value(u.m/u.s)
                        bb.append(bary);tt.append(bary+lsr)
                    d["frame_motion"][n]={"barycentric_delta_residual_vs_negative_vframe_m_s":stats((bb[0]-bb[1])+delta(vf)),"lsrk_delta_residual_vs_negative_vframe_m_s":stats((tt[0]-tt[1])+delta(vf)),"barycentric_delta_m_s":stats(bb[0]-bb[1]),"lsrk_delta_m_s":stats(tt[0]-tt[1])}
                if len(d["frame_motion"])>1:d["frame_motion"]["status"]="DIAGNOSTIC_ONLY"
            d["status"]="PASS_METADATA_DOPPLER_DIFFERENTIAL__HOLD_PHYSICAL_QUALIFICATION" if d["radio_differential_pass"] else "HOLD_DOPPLER_DIFFERENTIAL"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:
        if p.exists():
            d["original_unchanged"]=prior.sha(p)==S["source_sha256"];p.unlink()
            if not d["original_unchanged"]:d["status"]="HOLD_INTEGRITY"
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");print("R37F3D2G_SUMMARY",json.dumps(d,allow_nan=False),flush=True)
if __name__=="__main__":main()
