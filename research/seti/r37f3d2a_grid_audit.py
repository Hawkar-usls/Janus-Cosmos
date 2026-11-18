#!/usr/bin/env python3
"""R37F3D2A: isolate native grid / exposure failure without replacing bad metadata."""
import hashlib,json,pathlib
import numpy as np
from astropy.io import fits
import r37f3d2_archive_inventory as inv
OUT=pathlib.Path("runtime/r37f3d2a");OUT.mkdir(parents=True,exist_ok=True)
EXPECTED="f48450007a9d4a58c4ea6338be0d4a92175e53853bf875c05d0101affa9ce9a3"
def audit(src,p):
    with fits.open(p,memmap=True) as hd:
        hd.verify("exception")
        tables=[h for h in hd if isinstance(h,fits.BinTableHDU) and "DATA" in h.columns.names];assert len(tables)==1
        t=tables[0];d=t.data
        scan=np.asarray(d["SCAN"])
        cal=np.char.upper(np.char.strip(np.asarray(d["CAL"]).astype(str)))
        keys=["IFNUM","PLNUM","FDNUM"]
        cols=[np.asarray(d[k]) for k in keys]
        pairs=[]
        for g in sorted(set(zip(*[c[np.isin(scan,src["pair"])] for c in cols]))):
            mask=np.logical_and.reduce([c==v for c,v in zip(cols,g)])
            groups=[]
            for sc in src["pair"]:
                for c in ("T","F"):
                    m=mask&(scan==sc)&(cal==c);idx=np.flatnonzero(m)
                    fields={}
                    for key in ("CRVAL1","CRPIX1","CDELT1","EXPOSURE","DURATION"):
                        v=np.asarray(d[key][m],dtype=float);ok=np.isfinite(v);finite=v[ok]
                        bad=~ok
                        if key=="CDELT1":bad|=v==0
                        if key in ("EXPOSURE","DURATION"):bad|=v<=0
                        fields[key]={"finite":int(ok.sum()),"rows":len(v),"invalid":int(bad.sum()),"invalid_row_indices":idx[bad].tolist(),"min":float(finite.min()) if len(finite) else None,"max":float(finite.max()) if len(finite) else None}
                    flags={k:sorted(set(np.char.strip(np.asarray(d[k][m]).astype(str)).tolist())) for k in ("CTYPE1","CUNIT1","SPECSYS","VELDEF","BACKEND") if k in d.names}
                    groups.append({"scan":sc,"cal":c,"rows":len(idx),"date_count":len(set(np.asarray(d["DATE-OBS"][m]).astype(str))),"fields":fields,"frame_fields":flags,"native_grid_fields_valid":len(idx)>0 and all(fields[k]["invalid"]==0 for k in ("CRVAL1","CRPIX1","CDELT1","EXPOSURE"))})
            pairs.append({"ifnum":int(g[0]),"plnum":int(g[1]),"fdnum":int(g[2]),"groups":groups,"native_grid_fields_valid":all(v["native_grid_fields_valid"] for v in groups)})
        q=next(x for x in pairs if [x["ifnum"],x["plnum"],x["fdnum"]]==[src["ifnum"],src["plnum"],src["fdnum"]])
        return {"original_selected":q,"all_groups":pairs,"alternatives_valid":sum(x["native_grid_fields_valid"] for x in pairs),"spectral_column_format":str(t.columns["DATA"].format),"spectral_column_dim":str(t.columns["DATA"].dim),"power_read":False}
def main():
    src=next(s for s in inv.S["sources"] if s["id"]=="w3_1");p=OUT/"w3.fits"
    try:
        integ=inv.download(src,p);assert integ["sha256"]==EXPECTED and integ["bytes"]==130800960
        d=audit(src,p);out={"status":"DIAGNOSTIC_METADATA_ONLY","source":src,"integrity":integ,**d,"automatically_qualified":False,"detector_tested":False,"production_adoption":False}
        (OUT/"result.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
        print("R37F3D2A_SUMMARY",json.dumps({k:v for k,v in out.items() if k!="all_groups"},allow_nan=False),flush=True)
    finally:p.unlink(missing_ok=True)
if __name__=="__main__":main()
