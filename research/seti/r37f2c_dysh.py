#!/usr/bin/env python3
"""Official dysh frequency-switch natural W3OH control; no detector or source classification."""
import hashlib,importlib.metadata,json,pathlib,subprocess,sys
import numpy as np, requests
from astropy import units as u
from astropy.io import fits
from astropy.utils import iers
from dysh.fits.gbtfitsload import GBTFITSLoad
P=pathlib.Path("research/seti/r37f2c_dysh_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f2c");O.mkdir(parents=True,exist_ok=True)
C=S["source"];M=S["method"]
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def main():
    assert importlib.metadata.version("dysh")==M["dysh"]
    iers.conf.auto_download=False
    f=O/"W3OH_input.fits";b=bytearray()
    with requests.get(C["url"],stream=True,timeout=(15,60),allow_redirects=False) as r:
        r.raise_for_status();assert r.status_code==200
        for x in r.iter_content(65536):
            b.extend(x);assert len(b)<=C["bytes"]
    assert len(b)==C["bytes"] and sha(b)==C["sha256"];f.write_bytes(b)
    with fits.open(f,memmap=False) as h:
        a=next(t.data for t in h if isinstance(t,fits.BinTableHDU) and "DATA" in t.columns.names)
        assert len(a)==240 and set(str(v).strip() for v in a["OBJECT"])=={"W3OH"}
        assert set(int(v) for v in a["SCAN"])==set(C["scans"])
        assert set(float(v) for v in a["RESTFREQ"])=={float(C["rest_hz"])}
        assert set(str(v).strip() for v in a["CTYPE1"])=={"FREQ-OBS"}
        assert set(str(v).strip() for v in a["VELDEF"])=={"OPTI-LSR"}
    sdf=GBTFITSLoad(f);rows=[];arrays={}
    def band(v,intervals):
        z=np.zeros(v.shape,dtype=bool)
        for lo,hi in intervals:z|=(v>=lo)&(v<=hi)
        return z
    for scan in C["scans"]:
        for pol in C["polarizations"]:
            row={"scan":scan,"pol":pol,"status":"HOLD"}
            try:
                ta=sdf.getfs(scan=scan,ifnum=C["ifnum"],plnum=pol,fdnum=C["fdnum"],fold=True).timeaverage(weights=M["weights"])
                x=np.asarray(ta.flux.to_value(u.K),dtype=float).copy()
                v=np.asarray(ta.velocity_axis_to("km/s",toframe=M["velocity_frame"],doppler_convention=M["velocity_convention"]).to_value(u.km/u.s))
                assert x.shape==v.shape==(32768,)
                ta.baseline(model=M["baseline"],degree=M["degree"],include=[tuple(k) for k in M["include_channels"]],remove=True)
                y=np.asarray(ta.flux.to_value(u.K),dtype=float)
                mask=np.asarray(ta.mask,dtype=bool) if ta.mask is not None else np.zeros(y.shape,bool)
                assert mask.shape==y.shape
                ok=np.isfinite(v)&np.isfinite(y)&~mask
                line=ok&band(v,M["line_intervals_km_s"]);noise=ok&band(v,M["noise_intervals_km_s"])
                assert line.sum()>=10 and noise.sum()>=100
                med=float(np.median(y[noise]));sigma=1.4826*float(np.median(np.abs(y[noise]-med)));assert sigma>0
                j=int(np.flatnonzero(line)[np.argmax(y[line])]);dv=float(np.median(np.abs(np.diff(v))))
                z=float((y[j]-med)/sigma)
                integral=float(np.sum(y[line]-med)*dv)
                row.update(status="CALIBRATED",peak_radio_lsrk_km_s=float(v[j]),peak_K=float(y[j]),noise_median_K=med,noise_robust_K=sigma,descriptive_z=z,integrated_K_km_s=integral,velocity_bin_km_s=dv,line_channels=int(line.sum()),noise_channels=int(noise.sum()),supported=bool(z>=M["descriptive_threshold"] and integral>0))
                arrays[f"scan{scan}_pol{pol}"]=np.stack([v,x,y])
            except Exception as e:
                import traceback
                traceback.print_exc()
                row["error"]=type(e).__name__+": "+str(e)[:400]
            rows.append(row);print("R37F2C_SCAN",json.dumps(row),flush=True)
    if arrays:np.savez_compressed(O/"calibrated_spectra.npz",**arrays)
    f.unlink()
    ok=[r for r in rows if r["status"]=="CALIBRATED"]
    supported=[r for r in ok if r["supported"]]
    scans=[s for s in C["scans"] if all(any(r["scan"]==s and r["pol"]==p for r in supported) for p in C["polarizations"])]
    d={"status":"INDEPENDENT_EPOCH_NATURAL_LINE_SUPPORT" if len(ok)==10 and len(scans)>=M["minimum_scans_supported"] else "HOLD_CALIBRATION_OR_LINE_CRITERIA","input_sha256":C["sha256"],"dysh":importlib.metadata.version("dysh"),"astropy":importlib.metadata.version("astropy"),"date":C["date"],"rows":rows,"calibrated":len(ok),"supported":len(supported),"scans_both_polarizations_supported":scans,"distinct_observation_dates":1,"spatial_off_verified":False,"fine_grid_equivalence":False,"bl_frame_verified":False,"production_adoption":False}
    save("result.json",d)
    (O/"pip-freeze.txt").write_text(subprocess.check_output([sys.executable,"-m","pip","freeze"],text=True))
    print("R37F2C_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
