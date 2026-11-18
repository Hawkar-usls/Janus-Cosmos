#!/usr/bin/env python3
"""R37F3D2B: keep original bytes; evaluate exposure masks and official calibration."""
import hashlib,inspect,importlib.metadata,json,pathlib
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.utils import iers
from astropy.coordinates import solar_system_ephemeris
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.spectra.scan import PSScan
import r37f3d2_archive_inventory as inv
P=pathlib.Path("research/seti/r37f3d2b_exposure_protocol.json");S=json.loads(P.read_text())["source"]
O=pathlib.Path("runtime/r37f3d2b");O.mkdir(parents=True,exist_ok=True)
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(4*1024**2),b""):h.update(b)
    return h.hexdigest()
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def vec(x):
    x=np.asarray(x,dtype=float);return {"shape":list(x.shape),"finite":int(np.isfinite(x).sum()),"positive":int((x>0).sum()),"zero":int((x==0).sum()),"values":[float(v) if np.isfinite(v) else None for v in x.ravel()] if x.size<=128 else None}
def metadata(p):
    with fits.open(p,memmap=True) as hd:
        hd.verify("exception");tab=[h for h in hd if isinstance(h,fits.BinTableHDU) and "DATA" in h.columns.names];assert len(tab)==1;d=tab[0].data
        mask=(d["IFNUM"]==S["ifnum"])&(d["PLNUM"]==S["plnum"])&(d["FDNUM"]==S["fdnum"])
        cal=np.char.upper(np.char.strip(np.asarray(d["CAL"]).astype(str)));dates=np.char.strip(np.asarray(d["DATE-OBS"]).astype(str))
        indices=[];times=[];phase=[]
        for scan in S["scans"]:
            rr=[]
            for c in ("T","F"):
                i=np.flatnonzero(mask&(d["SCAN"]==scan)&(cal==c));i=i[np.argsort(dates[i])]
                assert len(i)==31 and len(set(dates[i]))==31;rr.append(i)
                v=np.asarray(d["EXPOSURE"][i],float);g=np.asarray(d["DURATION"][i],float)
                assert np.isfinite(g).all() and (g>0).all()
                phase.append({"scan":scan,"cal":c,"indices":i.tolist(),"exposure":vec(v),"duration":vec(g)})
            assert np.array_equal(dates[rr[0]],dates[rr[1]])
            times.append(Time(dates[rr[0]],format="isot",scale="utc").tai.unix_tai);indices.extend(rr)
        ix=np.stack(indices);exp=np.asarray(d["EXPOSURE"][ix],float);valid=np.isfinite(exp)&(exp>0)
        keep=valid.all(0);assert int(keep.sum())==30 and sorted(ix[~valid].tolist())==[61,123]
        dt=[np.diff(t) for t in times];assert np.allclose(dt[0],dt[1],atol=.01,rtol=0) and times[0][-1]<times[1][0]
        f=np.asarray(d["CRVAL1"][ix],float);pix=np.asarray(d["CRPIX1"][ix],float);df=np.asarray(d["CDELT1"][ix],float)
        assert np.isfinite(f).all() and np.isfinite(pix).all() and np.isfinite(df).all() and (df!=0).all()
        assert np.allclose(df,df[0,0],atol=0,rtol=0)
        # Frequency of first channel under FITS 1-based CRPIX convention.
        f0=f+(1-pix)*df;delta=(f0[0]-f0[2])/df[0]
        dist=np.abs(delta-np.rint(delta))
        return {"paired_valid":int(keep.sum()),"paired_total":len(keep),"excluded_ordinals":np.flatnonzero(~keep).tolist(),"excluded_all_four_raw_rows":ix[:,~keep].tolist(),"invalid_exposure_rows":ix[~valid].tolist(),"kept_row_indices":ix[:,keep].tolist(),"original_phase_metadata":phase,"cadence_s":[float(np.median(x)) for x in dt],"channel_hz":float(df[0,0]),"on_off_offset_channels_min":float(delta[keep].min()),"on_off_offset_channels_max":float(delta[keep].max()),"max_fractional_channel_residual":float(dist[keep].max()),"integer_common_grid":bool((dist[keep]<.001).all()),"power_read":False}
def official(p):
    assert importlib.metadata.version("dysh")=="1.0.1"
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
    sources={}
    for name in ("_calc_exposure","calibrate","timeaverage"):
        f=getattr(PSScan,name,None)
        if f is not None:
            src=inspect.getsource(f);sources[name]={"sha256":hashlib.sha256(src.encode()).hexdigest(),"relevant_lines":[ln for ln in src.splitlines() if any(v in ln.lower() for v in ("exposure","weight","mask","nan","blank"))][:80]}
    save("official_source_audit.json",sources)
    sdf=GBTFITSLoad(p);block=sdf.getps(scan=60,fdnum=0,ifnum=0,plnum=0,apply_flags=True)
    rows=[]
    for scan in block:
        e=np.asarray(scan.exposure,dtype=float);ts=np.asarray(scan.tsys,dtype=float);w=np.asarray(scan.weights,dtype=float)
        rows.append({"class":type(scan).__name__,"nint":int(scan.nint),"nchan":int(scan.nchan),"exposure":vec(e),"tsys":vec(ts),"weights":vec(w)})
    return {"version":importlib.metadata.version("dysh"),"scan_count":len(rows),"scans":rows,"official_success":True,"parity_proven":False,"notes":"Official outputs inspected without changing raw rows. Zero raw exposure handling must be reconciled with the paired mask; no automatic qualification."}
def main():
    src=next(x for x in inv.S["sources"] if x["id"]==S["id"]);p=O/"original.fits"
    d={"status":"HOLD","source":S,"qualified":False,"detector_tested":False,"production_adoption":False}
    try:
        d["integrity"]=inv.download(src,p);assert d["integrity"]["sha256"]==S["sha256"] and d["integrity"]["bytes"]==S["bytes"]
        p.chmod(0o444);d["metadata"]=metadata(p);d["status"]="PASS_EXPOSURE_PAIR_MASK__HOLD_CALIBRATION"
        try:d["official"]=official(p)
        except Exception as e:
            import traceback;traceback.print_exc();d["official"]={"official_success":False,"error":type(e).__name__+": "+str(e)[:500]}
        assert sha(p)==S["sha256"];d["original_unchanged"]=True
    except Exception as e:
        import traceback;traceback.print_exc();d["status"]="HOLD_AUDIT_FAILURE";d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:p.unlink(missing_ok=True)
    save("result.json",d)
    meta=d.get("metadata",{});summary={**d,"metadata":{k:v for k,v in meta.items() if k not in ("kept_row_indices","original_phase_metadata")}}
    print("R37F3D2B_SUMMARY",json.dumps(summary,allow_nan=False),flush=True)
if __name__=="__main__":main()
