#!/usr/bin/env python3
"""Historical natural astronomy WCS diagnostics; never read spectral DATA samples."""
import hashlib,importlib.metadata,inspect,json,pathlib,subprocess,sys
import numpy as np
from astropy.io import fits
from dysh.spectra import scan,spectrum
P=pathlib.Path("research/seti/r37f3d1ba_wcs_protocol.json")
S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f3d1ba");O.mkdir(parents=True,exist_ok=True)
for n,h in S["pinned_blobs"].items():
    b=(P.parent/n).read_bytes()
    assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==h,n
import r37f3d1a_metadata as meta
def summary(v):
    v=np.asarray(v,dtype=float);assert v.size and np.isfinite(v).all()
    return {"min":float(np.min(v)),"max":float(np.max(v)),"span":float(np.ptp(v)),"unique":int(len(np.unique(v)))}
def sha_file(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(4*1024**2),b""):h.update(b)
    return h.hexdigest()
def audit(p):
    groups=[];total=0
    with fits.open(p,memmap=True) as hd:
        t=next(h for h in hd if isinstance(h,fits.BinTableHDU) and h.header.get("EXTNAME")=="SINGLE DISH")
        d=t.data
        def col(n):return np.asarray(d[n])
        def st(n):return np.char.strip(col(n).astype(str))
        textfields=("CTYPE1","CUNIT1","VELDEF","SPECSYS","SSYSOBS")
        fields={n:sorted(set(st(n).tolist())) if n in d.names else t.header.get(n) for n in textfields}
        for s in (152,153):
            for i in range(5):
                for q in (0,1):
                    for cal in ("T","F"):
                        m=(col("SCAN")==s)&(col("IFNUM")==i)&(col("PLNUM")==q)&(col("FDNUM")==0)&(st("CAL")==cal)
                        assert int(m.sum())==151
                        total+=int(m.sum())
                        v={n:np.asarray(col(n)[m],float) for n in ("CRVAL1","CRPIX1","CDELT1","VFRAME","RESTFREQ","VELOCITY","OBSFREQ") if n in d.names}
                        assert all(n in v for n in ("CRVAL1","CRPIX1","CDELT1","VFRAME"))
                        f0=v["CRVAL1"]+(1-v["CRPIX1"])*v["CDELT1"]
                        assert np.all(v["CDELT1"]!=0)
                        beta=v["VFRAME"]/299792458.0;assert np.max(abs(beta))<1
                        fac=np.sqrt((1+beta)/(1-beta))
                        # Only conditional algebra; not a validated frame transformation.
                        corr={}
                        for name,k in (("multiply",fac),("divide",1/fac)):
                            corr[name]={"first_channel_span_hz":float(np.ptp(f0*k)),"channel_spacing_span_hz":float(np.ptp(abs(v["CDELT1"]*k)))}
                        groups.append({"scan":s,"if":i,"pol":q,"cal":cal,"rows":151,"fields":{n:summary(a) for n,a in v.items()},"observer_first_channel_span_hz":float(np.ptp(f0)),"observer_spacing_span_hz":float(np.ptp(abs(v["CDELT1"]))),"observer_fixed_0p01_hz":bool(np.ptp(f0)<=.01),"conditional_frame_spans":corr})
        assert total==len(d)==6040 and len(groups)==40
    return {"rows_audited":total,"declared_frames":fields,"groups":groups,"group_count":len(groups),
            "observer_max_first_channel_span_hz":max(g["observer_first_channel_span_hz"] for g in groups),
            "observer_max_spacing_span_hz":max(g["observer_spacing_span_hz"] for g in groups),
            "observer_fixed_groups":sum(g["observer_fixed_0p01_hz"] for g in groups),
            "conditional_max_first_channel_span_hz":{k:max(g["conditional_frame_spans"][k]["first_channel_span_hz"] for g in groups) for k in ("multiply","divide")}}
def code_audit():
    rows=[]
    for clsname in ("PSScan","ScanBlock","SpectralAverageMixin","Scan"):
        cls=getattr(scan,clsname,None)
        if cls is None:continue
        for name in ("timeaverage","_timeaverage","_average","getspec","calibrate"):
            f=getattr(cls,name,None)
            if f is None:continue
            try:
                code=inspect.getsource(inspect.unwrap(f))
                rows.append({"class":clsname,"method":name,"signature":str(inspect.signature(f)),"defined_in":getattr(f,"__qualname__",name),"source_sha256":hashlib.sha256(code.encode()).hexdigest(),"references_wcs":"wcs" in code.lower(),"references_align":"align" in code.lower(),"references_shift":"shift" in code.lower(),"references_average":"average" in code.lower(),"references_first_index":"[0]" in code})
            except (TypeError,OSError,ValueError):continue
    assert any(r["method"]=="timeaverage" for r in rows)
    return {"methods":rows,"scan_module_sha256":sha_file(pathlib.Path(scan.__file__)),"spectrum_module_sha256":sha_file(pathlib.Path(spectrum.__file__))}
def main():
    p=O/"original.fits";d={"status":"HOLD_WCS_DIAGNOSTIC","source":"NGC2415","production_adoption":False,"detector_tested":False,"calibration_tested":False,"averaging_lossless_proven":False,"frame_convention_validated":False,"prior_calibration_status":"HOLD_NATURAL_DYSH_CALIBRATION"}
    try:
        assert importlib.metadata.version("dysh")==S["dysh"]
        d["dysh"]=S["dysh"]
        d["integrity"]=meta.prev.download(p);assert d["integrity"]["sha256"]==S["archive_sha256"] and d["integrity"]["bytes"]==S["archive_bytes"]
        p.chmod(0o444);m={};meta.inspect(p,m);assert len(m["gates"])==14 and all(m["gates"].values())
        d["spatial_off_verified"]=True;d["wcs"]=audit(p);d["dysh_source"]=code_audit()
        d["original_unchanged"]=sha_file(p)==S["archive_sha256"];assert d["original_unchanged"]
        d["status"]="PASS_WCS_DIAGNOSTIC_ONLY"
    except Exception as e:
        import traceback;traceback.print_exc();d["error"]=type(e).__name__+": "+str(e)[:500]
    finally:p.unlink(missing_ok=True)
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    (O/"pip-freeze.txt").write_text(subprocess.check_output([sys.executable,"-m","pip","freeze"],text=True))
    print("R37F3D1BA_SUMMARY",json.dumps(d,allow_nan=False),flush=True)
if __name__=="__main__":main()
