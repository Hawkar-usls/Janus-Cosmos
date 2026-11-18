#!/usr/bin/env python3
"""Independent metadata checks for historical natural galaxy SDFITS."""
import hashlib,json,pathlib
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
P=pathlib.Path("research/seti/r37f3d1a_metadata_protocol.json")
S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f3d1a");O.mkdir(parents=True,exist_ok=True)
b=(P.parent/"r37f3d1_natural_archive.py").read_bytes()
assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==S["pinned_blob"]
import r37f3d1_natural_archive as prev
def inspect(p,out):
    gates={}
    out["gates"]=gates
    def gate(n,v):gates[n]=bool(v)
    with fits.open(p,memmap=True) as hd:
        hd.verify("exception")
        t=[h for h in hd if isinstance(h,fits.BinTableHDU) and h.header.get("EXTNAME")=="SINGLE DISH"]
        assert len(t)==1
        t=t[0];d=t.data
        def col(n):return np.asarray(d[n])
        def st(n):return np.char.upper(np.char.strip(col(n).astype(str)))
        def vals(n,m):return sorted(set(st(n)[m].tolist()))
        scan=col("SCAN");on,off=prev.S["source"]["scans"]
        gate("telescope","NRAO" in str(hd[0].header.get("TELESCOP","")).upper())
        gate("source_scan",set(st("OBJECT"))=={"NGC2415"} and set(scan)=={on,off})
        out["rows"]=len(d)
        names=("SIG","CAL","PROCSCAN","PROCSEQN","PROCSIZE","PROCTYPE","OBSMODE","CTYPE2","CTYPE3","RADESYS","LASTON","LASTOFF")
        out["scan_fields"]={str(s):{n:vals(n,scan==s) for n in names if n in d.names} for s in (on,off)}
        out["sig_format"]=str(t.columns["SIG"].format)
        out["cal_format"]=str(t.columns["CAL"].format)
        gate("procedure",all("ONOFF" in x for x in set(st("OBSMODE"))) and set(col("PROCSEQN")[scan==on])=={1} and set(col("PROCSEQN")[scan==off])=={2} and set(col("PROCSIZE"))=={2})
        # SIG is a switching phase; support the documented total-power
        # position switch only if its flags are consistent within each scan.
        sig=[set(st("SIG")[scan==s]) for s in (on,off)]
        gate("phase_consistent",sig in [[{"T"},{"F"}],[{"T"},{"T"}]] and set(st("CAL"))=={"T","F"})
        ifnum,plnum,fdnum=[col(n) for n in ("IFNUM","PLNUM","FDNUM")]
        gate("dimensions",len(set(ifnum))==5 and len(set(plnum))==2 and len(set(fdnum))==1 and len(d)==6040)
        gate("duration",all(np.isfinite(col(n)).all() and (col(n)>0).all() for n in ("DURATION","EXPOSURE")))
        gate("native_grid",np.isfinite(col("CDELT1")).all() and (col("CDELT1")!=0).all())
        out["native_channel_spacing_hz"]=sorted(set(abs(col("CDELT1")).tolist()))
        # Time groups: no counting simultaneous polarizations, IFs or cal phases twice.
        groups=[];alltimes={};pair_ok=True;no_overlap=True
        stamps=np.char.strip(col("DATE-OBS").astype(str))
        for s in (on,off):
            for i in sorted(set(ifnum)):
                for q in sorted(set(plnum)):
                    for f in sorted(set(fdnum)):
                        base=(scan==s)&(ifnum==i)&(plnum==q)&(fdnum==f)
                        phase={}
                        for cal in ("T","F"):
                            m=base&(st("CAL")==cal);raw=stamps[m];unique=sorted(set(raw.tolist()))
                            count=len(unique);phase[cal]=unique
                            if count:
                                time=Time(unique,format="isot",scale="utc")
                                dif=np.diff(time.tai.unix_tai)
                                duration=float(np.median(col("DURATION")[m]))
                                no_overlap &= bool(np.all(dif>=duration-.01))
                                span=float((time[-1]-time[0]).sec)+duration
                            else:duration=span=0.
                            groups.append({"scan":int(s),"if":int(i),"pol":int(q),"feed":int(f),"cal":cal,"rows":int(m.sum()),"unique":count,"duration_sec":duration,"span_sec":span})
                        pair_ok &= phase["T"]==phase["F"]
                        k=str(s);alltimes.setdefault(k,phase["T"])
                        pair_ok &= alltimes[k]==phase["T"]
        out["time_groups"]=groups
        out["unique_integrations_by_scan"]={k:len(v) for k,v in alltimes.items()}
        gate("time_counts",len(groups)==40 and all(g["unique"]==g["rows"]==151 for g in groups))
        gate("cal_and_window_time_alignment",pair_ok)
        gate("nonoverlap_within_phase",no_overlap)
        gate("disjoint_scan_timestamps",not set(alltimes[str(on)])&set(alltimes[str(off)]))
        # Same declared equatorial frame, actual pointing coordinates.
        frame=set(st("RADESYS"));out["coordinate_frame"]=sorted(frame)
        gate("coordinate_frame",len(frame)==1 and all(x.startswith("RA") for x in set(st("CTYPE2"))) and all(x.startswith("DEC") for x in set(st("CTYPE3"))))
        coords=[]
        for s in (on,off):
            m=scan==s;ra=col("CRVAL2")[m];dec=col("CRVAL3")[m]
            assert np.isfinite(ra).all() and np.isfinite(dec).all()
            coords.append(SkyCoord(ra=ra*u.deg,dec=dec*u.deg))
        sep=float(coords[0][0].separation(coords[1][0]).deg)
        spread=max(float(np.max(a.separation(a[0]).deg)) for a in coords)
        out["spatial_off_separation_deg"]=sep;out["max_pointing_spread_deg"]=spread
        gate("spatial_off",gates["coordinate_frame"] and sep>0.1 and spread<sep/4)
        checks={}
        for i,h in enumerate(hd):
            checks[str(i)]={k:bool(getattr(h,"verify_"+k.lower())()==1) if k in h.header else None for k in ("CHECKSUM","DATASUM")}
        out["fits_checksums"]=checks
        gate("internal_checksums",all(v is not False for m in checks.values() for v in m.values()))
def main():
    p=O/"original.fits";d={"status":"HOLD_METADATA_DIAGNOSTIC","source":"NGC2415","production_adoption":False,"detector_tested":False,"calibrated_line_tested":False,"independent_algorithm_validation":False}
    try:
        d["integrity"]=prev.download(p)
        assert d["integrity"]["bytes"]==S["archive_bytes"] and d["integrity"]["sha256"]==S["archive_sha256"]
        inspect(p,d)
        d["status"]="PASS_NATURAL_SPATIAL_OFF_METADATA" if all(d["gates"].values()) else "HOLD_METADATA_DIAGNOSTIC"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:p.unlink(missing_ok=True)
    d["qualified_sessions"]=int(d["status"]=="PASS_NATURAL_SPATIAL_OFF_METADATA")
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    print("R37F3D1A_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
