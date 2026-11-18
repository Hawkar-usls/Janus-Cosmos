#!/usr/bin/env python3
"""Calibrate a historical natural galaxy's position-switch data with dysh."""
import hashlib,importlib.metadata,json,pathlib,subprocess,sys
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.utils import iers
from dysh.fits.gbtfitsload import GBTFITSLoad
P=pathlib.Path("research/seti/r37f3d1b_natural_dysh_protocol.json")
S=json.loads(P.read_text());M=S["method"];O=pathlib.Path("runtime/r37f3d1b");O.mkdir(parents=True,exist_ok=True)
for n,h in S["pinned_blobs"].items():
    b=(P.parent/n).read_bytes()
    assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==h,n
import r37f3d1a_metadata as meta
def hashfile(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(4*1024**2),b""):h.update(b)
    return h.hexdigest()
def raw_grids(p):
    grids={}
    with fits.open(p,memmap=True) as hd:
        d=next(h.data for h in hd if isinstance(h,fits.BinTableHDU) and h.header.get("EXTNAME")=="SINGLE DISH")
        for i in M["ifnums"]:
            for q in M["pols"]:
                sel=(d["SCAN"]==M["scan"])&(d["IFNUM"]==i)&(d["PLNUM"]==q)&(d["FDNUM"]==M["fdnum"])
                assert int(sel.sum())==302
                w={n:np.asarray(d[n][sel],dtype=float) for n in ("CRVAL1","CRPIX1","CDELT1")}
                assert all(np.isfinite(v).all() for v in w.values())
                assert np.allclose(w["CDELT1"],w["CDELT1"][0],rtol=0,atol=M["frequency_spacing_atol_hz"])
                f0=w["CRVAL1"]+(1-w["CRPIX1"])*w["CDELT1"]
                assert np.ptp(f0)<=M["frequency_spacing_atol_hz"],"ON reference grid not fixed"
                grids[i,q]=f0[0]+np.arange(M["channels"])*w["CDELT1"][0]
    return grids
def measure(sdf,i,q,grid):
    ps=sdf.getps(scan=M["scan"],ifnum=i,plnum=q,fdnum=M["fdnum"])
    assert len(ps)==1
    ta=ps.timeaverage(weights=M["weights"])
    x=np.asarray(ta.flux.to_value(u.K),dtype=float)
    f=np.asarray(ta.spectral_axis.to_value(u.Hz),dtype=float)
    mask=np.asarray(ta.mask,dtype=bool) if ta.mask is not None else np.zeros(x.shape,dtype=bool)
    assert x.shape==f.shape==mask.shape==(M["channels"],)
    assert np.isfinite(f).all() and np.isfinite(x).all() and not mask.any()
    spacing=np.abs(np.diff(f));err=float(np.max(np.abs(spacing-M["native_df_hz"])))
    assert err<=M["frequency_spacing_atol_hz"],"native channel spacing changed"
    griderr=float(np.max(np.abs(f-grid)))
    assert griderr<=M["frequency_spacing_atol_hz"],"native observer-grid changed"
    return x.copy(),f.copy(),{"channels":len(x),"median_K":float(np.median(x)),"robust_spread_K":float(1.4826*np.median(np.abs(x-np.median(x)))),"spacing_max_error_hz":err,"grid_max_error_hz":griderr,"masked_channels":int(mask.sum()),"native_grid_preserved":True}
def main():
    iers.conf.auto_download=False
    p=O/"original.fits";d={"status":"HOLD_NATURAL_DYSH_CALIBRATION","source":"NGC2415","production_adoption":False,"detector_tested":False,"independent_algorithm_validation":False,"baseline_subtracted":False,"flux_density_calibrated":False,"rows":[]}
    try:
        assert importlib.metadata.version("dysh")==M["dysh"]
        d["dysh"]=M["dysh"];d["astropy"]=importlib.metadata.version("astropy")
        d["integrity"]=meta.prev.download(p)
        assert d["integrity"]["bytes"]==S["archive_bytes"] and d["integrity"]["sha256"]==S["archive_sha256"]
        m={};meta.inspect(p,m);assert len(m["gates"])==14 and all(m["gates"].values())
        d["spatial_off_verified"]=True;d["unique_integrations_by_scan"]=m["unique_integrations_by_scan"]
        grids=raw_grids(p);p.chmod(0o444);sdf=GBTFITSLoad(p)
        for i in M["ifnums"]:
            for q in M["pols"]:
                row={"if":i,"pol":q,"status":"HOLD"}
                try:
                    x,f,r=measure(sdf,i,q,grids[i,q]);row.update(r)
                    if i==q==0:
                        x2,f2,_=measure(sdf,i,q,grids[i,q])
                        d["tutorial_replay_exact"]=bool(np.array_equal(x,x2) and np.array_equal(f,f2))
                        assert d["tutorial_replay_exact"]
                    row["status"]="CALIBRATED_NATIVE"
                except Exception as e:
                    import traceback;traceback.print_exc();row["error"]=type(e).__name__+": "+str(e)[:400]
                d["rows"].append(row);print("R37F3D1B_WINDOW",json.dumps(row),flush=True)
        d["original_unchanged"]=hashfile(p)==S["archive_sha256"]
        d["calibrated"]=sum(r["status"]=="CALIBRATED_NATIVE" for r in d["rows"])
        if d["calibrated"]==10 and d.get("tutorial_replay_exact") and d["original_unchanged"]:
            d["status"]="PASS_NATURAL_DYSH_NATIVE_CALIBRATION"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:p.unlink(missing_ok=True)
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    (O/"pip-freeze.txt").write_text(subprocess.check_output([sys.executable,"-m","pip","freeze"],text=True))
    print("R37F3D1B_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
