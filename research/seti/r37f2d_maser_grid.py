#!/usr/bin/env python3
"""R37F2D: published natural maser archive metadata only; no spectral power read."""
import collections,hashlib,json,math,pathlib,re,urllib.parse,concurrent.futures
import requests,h5py,hdf5plugin
import r4_chunk_replay as base
import r37b_source_header as geom
P=pathlib.Path("research/seti/r37f2d_maser_grid_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f2d");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def norm(n):return re.sub(r"[\s_]","",str(n).upper())
def role(n):
    n=re.sub(r"^(?:DIAG|CAL[0-9]+)_","",str(n).strip().upper())
    if n.endswith("_OFF"):return n[:-4],"OFF"
    return n,"ON"
def identity(t,s):
    n,r=role(t)
    opts={norm(v) for v in [s["name"]]+s["aliases"] if not v.upper().endswith("_OFF")}
    opts|={x.removeprefix("G") for x in list(opts)}
    return r if norm(n) in opts or norm(n).removeprefix("G") in opts else None
def fmt(u):return pathlib.PurePosixPath(urllib.parse.urlparse(u).path).suffix.lower()
def source(s):
    queries=[];files={}
    for alias in s["aliases"]:
        q={"alias":alias}
        try:
            b,u=base.fetch("https://seti.berkeley.edu/opendata/api/query-files",{"target":alias,"telescope":"GBT","file-type":"HDF5"},limit=S["selection"]["max_metadata_bytes"])
            assert urllib.parse.urlparse(u).hostname=="seti.berkeley.edu"
            j=json.loads(b);rr=j if isinstance(j,list) else next(j[k] for k in ("data","files","results") if isinstance(j.get(k),list))
            count=0
            for v in rr:
                if not isinstance(v,dict):continue
                t=str(v.get("target","")).strip().upper();r=identity(t,s);url=v.get("url")
                if not r or not isinstance(url,str):continue
                p=urllib.parse.urlparse(url)
                if p.hostname!="bldata.berkeley.edu" and not re.fullmatch(r"blpd\d+\.ssl\.berkeley\.edu",p.hostname or ""):continue
                try:m=float(v["mjd"]);assert math.isfinite(m) and 50000<=m<=65000
                except (KeyError,TypeError,ValueError,AssertionError):continue
                if int(m) in S["selection"]["exclude_mjd"]:continue
                ff=fmt(url)
                if ff not in (".h5",".fil",".raw"):continue
                f={"source":s["name"],"catalogue_target":t,"role":r,"mjd":m,"day":int(m),"format":ff,"url":url,"catalogue":{k:v.get(k) for k in ("ra","decl","center_freq","size","md5sum","id")}}
                if url in files and files[url]!=f:files[url]["conflict"]=True
                else:files[url]=f
                count+=1
            q.update(status="PASS",rows=len(rr),exact_rows=count,sha256=sha(b))
        except Exception as e:q.update(status="HOLD",error=type(e).__name__+": "+str(e)[:240])
        queries.append(q)
    vals=[v for v in files.values() if not v.get("conflict")]
    def rank(v):
        try:c=float(v["catalogue"]["center_freq"]);dist=abs(c-S["reference_rest_mhz"])
        except (ValueError,TypeError):dist=1e9
        return (".gpuspec.0000." not in v["url"],dist,v["mjd"],v["url"])
    fs=[v for v in vals if v["format"]==".h5" and base.archive_url(v["url"])]
    chosen=[];seen=set()
    for day in sorted({v["day"] for v in fs}):
        for r in ("ON","OFF"):
            q=sorted([v for v in fs if v["day"]==day and v["role"]==r],key=rank)
            if q and q[0]["url"] not in seen:chosen.append(q[0]);seen.add(q[0]["url"])
    headers=[]
    for v in chosen[:S["selection"]["max_headers_per_source"]]:
        row={**v,"status":"HOLD"}
        try:
            with base.RangeFile(v["url"],8*base.MIB) as remote:
                with h5py.File(remote,"r") as h:
                    d=h["data"];a=d.attrs
                    name=a["source_name"];name=name.decode() if isinstance(name,bytes) else str(name)
                    assert name.strip().upper()==v["catalogue_target"],"header source mismatch"
                    t=float(a["tstart"]);df=float(a["foff"]);ts=float(a["tsamp"]);f0=float(a["fch1"])
                    assert all(math.isfinite(z) for z in (t,df,ts,f0)) and abs(t-v["mjd"])*86400<=3
                    assert d.ndim==3 and df!=0 and ts>0
                    hdr={"src_raj":str(a.get("src_raj","")),"src_dej":str(a.get("src_dej",""))}
                    try:point=geom.pointing(hdr,v["catalogue"])
                    except Exception:point=None
                    bounds=sorted([f0,f0+(d.shape[-1]-1)*df]);g=S["grid"]
                    grid=abs(df)*1e6<=g["max_channel_hz"] and ts<=g["max_integration_s"] and d.shape[0]>=g["min_rows"]
                    lo,hi=S["coverage_window_mhz"]
                    flags={k:str(a[k]) for k in ("barycentric","pulsarcentric","SPECSYS","specsys","SSYSOBS","ssysobs","VELDEF","veldef","CTYPE1","ctype1") if k in a}
                    row.update(status="HEADER_VERIFIED",channel_hz=abs(df)*1e6,integration_s=ts,rows=int(d.shape[0]),nominal_duration_s=ts*d.shape[0],pointing_deg=point,comparable_grid=grid,exact_reference_grid=abs(abs(df)*1e6-g["reference_df_hz"])<1e-7 and abs(ts-g["reference_ts_s"])<1e-7,covers_methanol=bool(bounds[0]<=lo and bounds[1]>=hi),frequency_bounds_mhz=bounds,frame_fields=flags,bytes_read=remote.total)
        except Exception as e:row["error"]=type(e).__name__+": "+str(e)[:200]
        headers.append(row)
    # Legacy products are catalogue candidates only; no fabricated HTTPS mirror.
    pairs=[]
    for a in headers:
        if a["role"]!="ON" or not (a.get("comparable_grid") and a.get("covers_methanol") and a.get("pointing_deg")):continue
        for b in headers:
            if b["role"]!="OFF" or b["day"]!=a["day"] or not (b.get("comparable_grid") and b.get("covers_methanol") and b.get("pointing_deg")):continue
            sep=geom.sep(a["pointing_deg"],b["pointing_deg"]);dt=abs(a["mjd"]-b["mjd"])*86400
            if sep>=.15 and 0<dt<=600 and abs(a["channel_hz"]-b["channel_hz"])<1e-7 and abs(a["integration_s"]-b["integration_s"])<1e-7:
                pairs.append({"day":a["day"],"on":a["url"],"off":b["url"],"separation_deg":sep,"time_gap_s":dt,"status":"METADATA_PAIR_CANDIDATE_ONLY"})
    good=[r for r in headers if r.get("comparable_grid") and r.get("covers_methanol") and r["role"]=="ON"]
    return {"source":s["name"],"queries":queries,"unique_urls":len(vals),"formats":dict(collections.Counter(v["format"] for v in vals)),"dates_by_format":{k:sorted({v["day"] for v in vals if v["format"]==k}) for k in (".h5",".fil",".raw")},"on_count":sum(v["role"]=="ON" for v in vals),"off_count":sum(v["role"]=="OFF" for v in vals),"headers":headers,"fine_on_covering_methanol":len(good),"spatial_metadata_pairs":pairs,"legacy_candidates":[{k:v[k] for k in ("catalogue_target","role","day","format","url","catalogue")} for v in sorted([v for v in vals if v["format"]==".fil"],key=rank)[:6]]}
def main():
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(source,S["known_sources"]))
    d={"status":"NATURAL_FINE_GRID_METADATA_AUDIT_ONLY","sources":rows,"unique_urls_sum":sum(v["unique_urls"] for v in rows),"fine_on_covering_methanol":sum(v["fine_on_covering_methanol"] for v in rows),"pair_candidates":sum(len(v["spatial_metadata_pairs"]) for v in rows),"natural_line_verified":False,"exact_bl_frame_verified":False,"raw_power_read":False,"production_adoption":False}
    save("result.json",d)
    print("R37F2D_SUMMARY",json.dumps({**d,"sources":[{k:v for k,v in r.items() if k!="legacy_candidates"} for r in rows]}),flush=True)
if __name__=="__main__":main()
