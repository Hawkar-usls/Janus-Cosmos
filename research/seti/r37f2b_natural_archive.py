#!/usr/bin/env python3
"""Independent astronomical control metadata and frequency-frame provenance."""
import collections,hashlib,io,json,pathlib,re,urllib.parse,concurrent.futures
import numpy as np
import requests
from astropy.io import fits
P=pathlib.Path("research/seti/r37f2b_natural_archive_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f2b");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,x):(O/n).write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def download(u,cap,hosts,params=None):
    assert urllib.parse.urlparse(u).scheme=="https" and urllib.parse.urlparse(u).hostname in hosts
    with requests.get(u,params=params,stream=True,timeout=(15,60),allow_redirects=False) as r:
        r.raise_for_status();assert r.status_code==200
        if r.headers.get("Content-Length","").isdigit():assert int(r.headers["Content-Length"])<=cap
        b=bytearray()
        for x in r.iter_content(65536):
            b.extend(x);assert len(b)<=cap
        return bytes(b)
def scalar(v):
    if isinstance(v,bytes):return v.decode("ascii",errors="replace").strip()
    if isinstance(v,(str,bool,int,np.integer,np.bool_)):return v.item() if hasattr(v,"item") else v
    if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
    return str(v)
def values(a,k):
    if k not in (a.dtype.names or ()):return []
    return sorted({str(scalar(v)) for v in a[k]})
def sdfits(spec):
    d={"id":spec["id"],"url":spec["url"],"status":"HOLD"}
    try:
        b=download(spec["url"],S["limits"]["download_bytes"],{"www.gb.nrao.edu","safe.nrao.edu"})
        d.update(file_bytes=len(b),sha256=sha(b));tables=[]
        with fits.open(io.BytesIO(b),memmap=False,checksum=True) as hd:
            for h in hd:
                if not isinstance(h,fits.BinTableHDU):continue
                a=h.data;names=set(a.dtype.names or ());assert len(a)<=100000
                if not {"OBJECT","SCAN","DATA"}<=names:continue
                sel=np.array([str(scalar(v)).strip().upper()=="W3OH" for v in a["OBJECT"]])
                q=a[sel];assert len(q)>0
                scans=sorted({int(v) for v in q["SCAN"]})
                frames={k:values(q,k) for k in ("CTYPE1","VELDEF","SPECSYS","SSYSOBS","CUNIT1","RADESYS")}
                head={k:scalar(h.header[k]) for k in ("SPECSYS","SSYSOBS","CTYPE1","VELDEF","CUNIT1","RADESYS") if k in h.header}
                fields={k:values(q,k) for k in ("DATE-OBS","SIG","CAL","RESTFREQ","OBSFREQ","CDELT1","CRVAL1","CRPIX1","EXPOSURE","IFNUM","PLNUM","FDNUM","PROJID","OBSERVER") if k in names and k!="OBSERVER"}
                phases=[]
                for scan in scans:
                    s=q[q["SCAN"]==scan]
                    sig=values(s,"SIG");cal=values(s,"CAL")
                    phases.append({"scan":scan,"rows":len(s),"sig":sig,"cal":cal,"both_sig_ref":len(sig)==2,"frequency_setting_count":len(values(s,"CRVAL1"))})
                tables.append({"rows":len(q),"scans":scans,"expected_scans_present":set(spec["expected_scans"])<=set(scans),"frame_fields":frames,"header_frame":head,"fields":fields,"phase_summary":phases,"rows_with_both_phases":sum(v["both_sig_ref"] for v in phases),"data_shape":list(q["DATA"].shape),"data_columns_decoded_for_analysis":False})
        assert tables
        d.update(status="SDFITS_METADATA_VERIFIED",tables=tables,frequency_switch_control=any(v["rows_with_both_phases"]>0 and v["expected_scans_present"] for v in tables),spatial_off_verified=False,natural_line_verified=False,archive_frame_scope="SDFITS only")
    except Exception as e:d["error"]=type(e).__name__+": "+str(e)[:300]
    return d
def catalogue():
    cfg=S["berkeley"];seen={};queries=[]
    for alias in cfg["queries"]:
        row={"alias":alias}
        try:
            b=download(cfg["api"],S["limits"]["metadata_bytes"],{"seti.berkeley.edu"},{"target":alias,"telescope":"GBT","file-type":"HDF5"})
            j=json.loads(b);rr=j if isinstance(j,list) else next(j[k] for k in ("data","files","results") if isinstance(j.get(k),list))
            n=0
            for r in rr:
                if not isinstance(r,dict):continue
                u=r.get("url");t=str(r.get("target","")).strip().upper()
                if not isinstance(u,str) or not re.fullmatch(r"(?:DIAG_)?W3(?:OH)?(?:_OFF[A-Z]?)?",t):continue
                try:m=float(r["mjd"]);assert np.isfinite(m)
                except (ValueError,TypeError,AssertionError,KeyError):continue
                role="OFF" if "_OFF" in t else "ON";seen[u]={"target":t,"day":int(m),"role":role,"format":pathlib.PurePosixPath(urllib.parse.urlparse(u).path).suffix.lower(),"mjd":m,"url":u}
                n+=1
            row.update(status="PASS",rows=len(rr),exact_source_rows=n,sha256=sha(b))
        except Exception as e:row.update(status="HOLD",error=str(e)[:200])
        queries.append(row)
    dates=[]
    for day in cfg["days"]:
        q=[v for v in seen.values() if v["day"]==day]
        dates.append({"day":day,"unique":len(q),"on":sum(v["role"]=="ON" for v in q),"off":sum(v["role"]=="OFF" for v in q),"off_targets":sorted({v["target"] for v in q if v["role"]=="OFF"}),"formats":dict(collections.Counter(v["format"] for v in q))})
    offdays=sorted({v["day"] for v in seen.values() if v["role"]=="OFF"})
    return {"queries":queries,"unique_urls":len(seen),"off_catalogue_days":offdays,"requested_dates":dates,"additional_off_metadata_count":sum(v["off"] for v in dates),"qualified_spatial_pairs":0}
def main():
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        a=pool.submit(catalogue);f=[pool.submit(sdfits,v) for v in S["sdfits"]]
        r={"status":"NATURAL_CONTROL_METADATA_AUDIT_ONLY","berkeley":a.result(),"sdfits":[x.result() for x in f],"breakthrough_frame":{"status":"DOCUMENTED_TOPocentric_CONVENTION_NOT_EXACT_OBSERVATION_LINEAGE","evidence":S["frame_evidence"],"exact_input_frame_verified":False},"production_adoption":False,"independent_natural_line_verified":False}
    save("result.json",r)
    # Show duration/group summaries, retain full headers in artifact.
    def compact(t):
        return {**t,"tables":[{**x,"fields":{k:(v[:8] if k=="DATE-OBS" else v[:12]) for k,v in x["fields"].items()}} for x in t.get("tables",[])]}
    print("R37F2B_SUMMARY",json.dumps({**r,"sdfits":[compact(t) for t in r["sdfits"]]}),flush=True)
if __name__=="__main__":main()
