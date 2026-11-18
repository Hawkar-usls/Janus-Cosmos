#!/usr/bin/env python3
"""Metadata-only qualification of historical natural astronomy controls."""
import hashlib,json,pathlib,sys
import numpy as np
import requests
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
P=pathlib.Path("research/seti/r37f3d2_archive_inventory_protocol.json")
S=json.loads(P.read_text());L=S["limits"]
O=pathlib.Path("runtime/r37f3d2");O.mkdir(parents=True,exist_ok=True)
def download(src,p):
    h=hashlib.sha256();n=0
    with requests.get(src["url"],stream=True,timeout=(20,75),allow_redirects=False) as r:
        r.raise_for_status();assert r.status_code==200
        if r.headers.get("Content-Length"):assert int(r.headers["Content-Length"])<=L["max_bytes_per_archive"],"archive exceeds byte cap"
        with p.open("wb") as f:
            for b in r.iter_content(4*1024**2):
                n+=len(b);assert n<=L["max_bytes_per_archive"],"archive exceeds byte cap"
                h.update(b);f.write(b)
        if r.headers.get("Content-Length"):assert n==int(r.headers["Content-Length"])
    assert n>0
    return {"bytes":n,"sha256":h.hexdigest(),"published_checksum_verified":False}
def inspect(src,p,out):
    with fits.open(p,memmap=True) as hd:
        hd.verify("exception")
        ts=[h for h in hd if isinstance(h,fits.BinTableHDU) and "DATA" in h.columns.names]
        assert ts,"no spectral table"
        out["spectral_tables"]=len(ts)
        checks=[]
        for h in hd:
            checks.append({k:(bool(getattr(h,"verify_"+k.lower())()==1) if k in h.header else None) for k in ("CHECKSUM","DATASUM")})
        out["internal_checksums"]=checks
        assert all(v is not False for m in checks for v in m.values())
        inventory=[];pairs=[]
        for t in ts:
            d=t.data
            def col(n):return np.asarray(d[n])
            def st(n):return np.char.upper(np.char.strip(col(n).astype(str)))
            required={"OBJECT","SCAN","DATE-OBS","IFNUM","PLNUM","FDNUM","CAL"}
            assert required.issubset(d.names),str(sorted(required-set(d.names)))
            scan=col("SCAN");obj=st("OBJECT");date=np.char.strip(col("DATE-OBS").astype(str))
            ifn=col("IFNUM");pol=col("PLNUM");fd=col("FDNUM");cal=st("CAL")
            modes=st("OBSMODE") if "OBSMODE" in d.names else np.full(len(d),"UNKNOWN")
            out.setdefault("dates",[]);out["dates"]=sorted(set(out["dates"])|set(x[:10] for x in date))
            for sc in sorted(set(scan)):
                m=scan==sc
                group_counts=[]
                for key in set(zip(ifn[m],pol[m],fd[m],cal[m])):
                    g=m&(ifn==key[0])&(pol==key[1])&(fd==key[2])&(cal==key[3])
                    group_counts.append(len(set(date[g])))
                inventory.append({"scan":int(sc),"source":sorted(set(obj[m])),"modes":sorted(set(modes[m])),"raw_rows":int(m.sum()),"max_unique_times_fixed_group":max(group_counts),"min_unique_times_fixed_group":min(group_counts)})
            if "pair" not in src:continue
            on,off=src["pair"]
            gated={"table_has_pair":set(src["pair"]).issubset(set(scan))}
            if not gated["table_has_pair"]:continue
            needed={"PROCSEQN","PROCSCAN","CRVAL2","CRVAL3","CTYPE2","CTYPE3","RADESYS","DURATION","EXPOSURE","CDELT1","CRPIX1","CRVAL1"}
            if not needed.issubset(d.names):
                pairs.append({"gates":{"metadata_columns":False},"missing":sorted(needed-set(d.names))});continue
            select=(ifn==src["ifnum"])&(pol==src["plnum"])&(fd==src["fdnum"])
            masks=[select&(scan==s) for s in (on,off)]
            gated["selected_both"]=all(m.any() for m in masks)
            if not gated["selected_both"]:pairs.append({"gates":gated});continue
            gated["same_source"]=all(set(obj[m])=={src["name"].upper()} for m in masks)
            gated["procedure"]=all("ONOFF" in v for m in masks for v in set(modes[m])) and all(set(col("PROCSEQN")[m])=={i+1} for i,m in enumerate(masks)) and all(set(st("PROCSCAN")[m])=={role} for m,role in zip(masks,("ON","OFF")))
            gated["equatorial_frame"]=len(set(st("RADESYS")[masks[0]|masks[1]]))==1 and all(v.startswith("RA") for m in masks for v in set(st("CTYPE2")[m])) and all(v.startswith("DEC") for m in masks for v in set(st("CTYPE3")[m]))
            counts=[];timeline=[];align=True;nonoverlap=True
            for m in masks:
                phases=[]
                for c in ("T","F"):
                    g=m&(cal==c);raw=date[g];uniq=sorted(set(raw));counts.append(len(uniq));phases.append(uniq)
                    if not uniq:nonoverlap=False;continue
                    ti=Time(uniq,format="isot",scale="utc");dur=col("DURATION")[g]
                    nonoverlap &= len(raw)==len(uniq) and bool(np.isfinite(dur).all() and (dur>0).all() and np.all(np.diff(ti.tai.unix_tai)>=float(np.max(dur))-.01))
                align &= phases[0]==phases[1]
                timeline.append(phases[0])
            gated["time_counts"]=min(counts)>=L["min_unique_per_scan"]
            gated["cal_times"]=align and all(set(cal[m])=={"T","F"} for m in masks)
            gated["nonoverlap"]=nonoverlap and not set(timeline[0])&set(timeline[1])
            gated["native_grid"]=all(np.isfinite(col(n)[m]).all() for m in masks for n in ("CRVAL1","CRPIX1","CDELT1","EXPOSURE")) and all((col("CDELT1")[m]!=0).all() and (col("EXPOSURE")[m]>0).all() for m in masks)
            coords=[SkyCoord(ra=col("CRVAL2")[m]*u.deg,dec=col("CRVAL3")[m]*u.deg) for m in masks]
            sep=float(coords[0][0].separation(coords[1][0]).deg)
            spread=max(float(np.max(a.separation(a[0]).deg)) for a in coords)
            gated["spatial_off"]=gated["equatorial_frame"] and sep>L["min_separation_deg"] and spread<sep/4
            pairs.append({"gates":{k:bool(v) for k,v in gated.items()},"time_counts":counts,"separation_deg":sep,"pointing_spread_deg":spread})
        out["scan_inventory"]=inventory
        out["pair_checks"]=pairs
        out["position_switch_metadata_qualified"]=not src.get("inventory_only",False) and len(pairs)==1 and all(pairs[0]["gates"].values()) and len(pairs[0]["gates"])>=9
def main():
    wanted=sys.argv[1]
    src=next(s for s in S["sources"] if s["id"]==wanted)
    p=O/(wanted+".fits")
    out={"source_id":wanted,"project":src["project"],"status":"HOLD_ARCHIVE_INVENTORY","position_switch_metadata_qualified":False,"calibrated_line_tested":False,"detector_tested":False,"production_adoption":False}
    try:
        out["integrity"]=download(src,p)
        inspect(src,p,out)
        out["status"]="PASS_POSITION_SWITCH_METADATA_ONLY" if out["position_switch_metadata_qualified"] else "INVENTORY_ONLY"
    except Exception as e:
        import traceback
        traceback.print_exc();out["error"]=type(e).__name__+": "+str(e)[:400]
    finally:p.unlink(missing_ok=True)
    (O/(wanted+".json")).write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print("R37F3D2_SUMMARY",json.dumps({k:v for k,v in out.items() if k!="scan_inventory"},allow_nan=False),flush=True)
if __name__=="__main__":main()
