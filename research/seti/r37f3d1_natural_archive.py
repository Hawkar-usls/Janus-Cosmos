#!/usr/bin/env python3
"""Qualification of documented historical astronomy archive metadata only."""
import hashlib,json,pathlib,urllib.parse
import numpy as np
import requests
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
P=pathlib.Path("research/seti/r37f3d1_natural_archive_protocol.json")
S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f3d1");O.mkdir(parents=True,exist_ok=True)
def download(p):
    src=S["source"];h=hashlib.sha256();n=0
    with requests.get(src["archive"],stream=True,timeout=(20,90),allow_redirects=False) as r:
        r.raise_for_status();assert r.status_code==200
        with p.open("wb") as f:
            for b in r.iter_content(4*1024**2):
                n+=len(b);assert n<=src["max_bytes"];h.update(b);f.write(b)
        if r.headers.get("Content-Length"):assert n==int(r.headers["Content-Length"])
    assert n>0
    return {"bytes":n,"sha256":h.hexdigest(),"published_checksum_verified":False,"https_archive":True}
def inspect(p):
    src=S["source"]
    with fits.open(p,memmap=True) as hd:
        hd.verify("exception")
        assert "NRAO" in str(hd[0].header.get("TELESCOP","")).upper()
        tables=[x for x in hd if isinstance(x,fits.BinTableHDU) and x.header.get("EXTNAME")=="SINGLE DISH"]
        assert len(tables)==1
        t=tables[0];d=t.data
        def col(n):return np.asarray(d[n])
        def st(n):return np.char.strip(col(n).astype(str))
        checks={}
        for i,x in enumerate(hd):
            checks[str(i)]={k:bool(getattr(x,"verify_"+k.lower())()==1) if k in x.header else None for k in ("CHECKSUM","DATASUM")}
        assert all(v is not False for m in checks.values() for v in m.values())
        on,off=src["scans"];scan=col("SCAN");obj=st("OBJECT");proc=st("OBSMODE")
        assert set(np.unique(scan))=={on,off}
        assert set(obj)=={src["name"]}
        assert all("ONOFF" in x.upper() for x in np.unique(proc))
        assert set(st("SIG")[scan==on])=={"T"} and set(st("SIG")[scan==off])=={"F"}
        assert set(col("PROCSEQN")[scan==on])=={1} and set(col("PROCSEQN")[scan==off])=={2}
        assert len(set(col("IFNUM")))==src["expected_windows"]
        assert all(np.isfinite(col(n)).all() for n in ("CDELT1","DURATION","EXPOSURE","CRVAL2","CRVAL3"))
        assert (col("DURATION")>0).all() and (col("EXPOSURE")>0).all() and (col("CDELT1")!=0).all()
        # Fixed IF and pol groups; do not count calibration phases or duplicate time stamps as extra integrations.
        counts=[];stamps=st("DATE-OBS")
        for s in (on,off):
            for i in sorted(set(col("IFNUM"))):
                for q in sorted(set(col("PLNUM"))):
                    mask=(scan==s)&(col("IFNUM")==i)&(col("PLNUM")==q)
                    for cal in sorted(set(st("CAL")[mask])):
                        m=mask&(st("CAL")==cal)
                        count=len(set(stamps[m]));assert count==int(m.sum())
                        counts.append(count)
        assert counts and min(counts)>=src["min_unique_integrations"]
        # Spatial reference: coordinate types must be identical equatorial frames,
        # actual pointing separation, not the mere existence of a SIG flag.
        assert all(x.upper().startswith("RA") for x in set(st("CTYPE2")))
        assert all(x.upper().startswith("DEC") for x in set(st("CTYPE3")))
        assert len(set(st("RADESYS")))==1
        coords=[]
        for s in (on,off):
            m=scan==s
            a=SkyCoord(ra=col("CRVAL2")[m]*u.deg,dec=col("CRVAL3")[m]*u.deg)
            coords.append(a)
        sep=float(coords[0][0].separation(coords[1][0]).deg)
        spread=max(float(np.max(a.separation(a[0]).deg)) for a in coords)
        assert sep>0.1 and spread<sep/4
        return {"rows":len(d),"unique_time_integrations_min":min(counts),"unique_time_integrations_max":max(counts),"groups":len(counts),"windows":len(set(col("IFNUM"))),"polarizations":len(set(col("PLNUM"))),"native_channel_spacing_hz_values":sorted(set(np.abs(col("CDELT1")).tolist())),"spatial_off_separation_deg":sep,"max_pointing_spread_deg":spread,"scan_procedure_verified":True,"fits_checksums":checks,"metadata_qualified":True,"calibrated_line_tested":False}
def main():
    p=O/"temporary.fits";d={"status":"HOLD_ARCHIVE_QUALIFICATION","source":S["source"]["name"],"production_adoption":False,"detector_tested":False,"independent_algorithm_validation":False}
    try:
        d["integrity"]=download(p);d["metadata"]=inspect(p)
        d["status"]="PASS_NEW_NATURAL_ARCHIVE_METADATA"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:p.unlink(missing_ok=True)
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    print("R37F3D1_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
