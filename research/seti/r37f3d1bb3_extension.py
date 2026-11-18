#!/usr/bin/env python3
"""Repeat the frozen scientific calibration audit on four historical galaxy windows."""
import hashlib,json,pathlib,subprocess,sys,importlib.metadata
from astropy.utils import iers
from astropy.coordinates import solar_system_ephemeris
from dysh.fits.gbtfitsload import GBTFITSLoad
P=pathlib.Path("research/seti/r37f3d1bb3_extension_protocol.json")
S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f3d1bb3");O.mkdir(parents=True,exist_ok=True)
def pins():
    for n,h in S["pinned_blobs"].items():
        b=(P.parent/n).read_bytes()
        assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==h,n
pins()
import r37f3d1bb2b_header as audit
def main():
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
    p=O/"original.fits";d={"status":"HOLD_NATURAL_WINDOW_EXTENSION","source":"NGC2415","rows":[],"production_adoption":False,"detector_tested":False,"lossless_average_proven":False,"absolute_astrophysical_accuracy_proven":False}
    before=dict(audit.M)
    try:
        assert importlib.metadata.version("dysh")=="1.0.1"
        d["integrity"]=audit.prev.meta.prev.download(p)
        assert d["integrity"]["sha256"]==audit.prev.S["archive_sha256"] and d["integrity"]["bytes"]==audit.prev.S["archive_bytes"]
        p.chmod(0o444);m={};audit.prev.meta.inspect(p,m);assert len(m["gates"])==14 and all(m["gates"].values())
        d["spatial_off_verified"]=True
        for i in S["ifnums"]:
            audit.M["ifnum"]=i
            sdf=GBTFITSLoad(p)
            for q in S["pols"]:
                row={"if":i,"pol":q,"status":"HOLD"}
                try:
                    audit.measure(sdf,p,q,row)
                    row["status"]="PASS_TSYS_ARITHMETIC_SUPPORT"
                except Exception as e:
                    import traceback;traceback.print_exc();row["error"]=type(e).__name__+": "+str(e)[:500]
                d["rows"].append(row);print("R37F3D1BB3_WINDOW",json.dumps(row),flush=True)
            del sdf
        d["original_unchanged"]=audit.prev.hashfile(p)==audit.prev.S["archive_sha256"]
        pins()
        passed=[r for r in d["rows"] if r["status"]=="PASS_TSYS_ARITHMETIC_SUPPORT"]
        d["passed_spectra"]=len(passed);d["planned_spectra"]=8
        d["valid_calibrated_samples_passed"]=sum(r["calibration_valid_samples"] for r in passed)
        d["status"]="PASS_ISOLATED_NATURAL_WINDOW_EXTENSION" if d["original_unchanged"] and len(passed)==8 else "HOLD_NATURAL_WINDOW_EXTENSION"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:
        audit.M.clear();audit.M.update(before);p.unlink(missing_ok=True)
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    (O/"pip-freeze.txt").write_text(subprocess.check_output([sys.executable,"-m","pip","freeze"],text=True))
    print("R37F3D1BB3_SUMMARY",json.dumps(d,allow_nan=False),flush=True)
if __name__=="__main__":main()
