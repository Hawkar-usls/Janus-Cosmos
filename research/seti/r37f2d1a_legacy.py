#!/usr/bin/env python3
"""R37F2D1 natural astronomical catalogue and SIGPROC metadata only."""
import collections,hashlib,json,math,pathlib,re,urllib.parse,concurrent.futures
import requests
import astropy.units as u
from astropy.coordinates import SkyCoord
import r4_chunk_replay as base
from r37f1de_archive import fil_header,fil_grid
P=pathlib.Path("research/seti/r37f2d1a_legacy_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f2d1a");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def official(url):
    p=urllib.parse.urlparse(url)
    return p.scheme in ("http","https") and (p.hostname=="bldata.berkeley.edu" or bool(re.fullmatch(r"blpd\d+\.ssl\.berkeley\.edu",p.hostname or ""))) and not p.username and not p.password
def fetchcat(alias):
    b,url=base.fetch("https://seti.berkeley.edu/opendata/api/query-files",{"target":alias,"telescope":"GBT","file-type":"HDF5"},limit=S["limits"]["catalogue_bytes"])
    assert urllib.parse.urlparse(url).hostname=="seti.berkeley.edu"
    j=json.loads(b);rr=j if isinstance(j,list) else next(j[k] for k in ("data","files","results") if isinstance(j.get(k),list))
    return b,rr,url
def exact(t,f):
    t=re.sub(r"^(?:DIAG|CAL[0-9]+)_","",t)
    role="OFF" if t.endswith("_OFF") else "ON"
    if role=="OFF":t=t[:-4]
    if "galactic" not in f:return (t.replace("_","")==f["name"].replace("_",""),role)
    if t in f.get("verified_catalogue_aliases",[]):return True,role
    # Preserve sign and both numeric coordinates. No prefix-only fuzzy matching.
    m=re.fullmatch(r"G?(\d{1,3})[._](\d{2,3})[+_](\d{1,2})[._](\d{2,3})",t)
    if not m:return False,role
    l=float(m[1]+"."+m[2]);b=float(m[3]+"."+m[4])
    return abs(l-f["galactic"][0])<.005 and abs(b-f["galactic"][1])<.005,role
def coord(v,f):
    if "galactic" not in f:return None
    try:
        ra=float(v["ra"]);dec=float(v["decl"]);assert 0<=ra<360 and abs(dec)<=90
        g=SkyCoord(ra=ra*u.deg,dec=dec*u.deg,frame="icrs").galactic
        ref=SkyCoord(l=f["galactic"][0]*u.deg,b=f["galactic"][1]*u.deg,frame="galactic")
        return float(g.separation(ref).deg)
    except Exception:return None
def probe(v):
    url=v["url"];p=urllib.parse.urlparse(url);urls=[urllib.parse.urlunparse(p._replace(scheme="https"))]
    if p.scheme=="http":urls.append(url)
    attempts=[]
    for q in urls:
        a={"scheme":urllib.parse.urlparse(q).scheme}
        try:
            cap=S["limits"]["header_prefix_bytes"]
            with requests.get(q,headers={"Range":f"bytes=0-{cap-1}","Accept-Encoding":"identity"},stream=True,timeout=(10,S["limits"]["timeout_s"]),allow_redirects=False) as r:
                a["http_status"]=r.status_code
                if r.status_code!=206:raise IOError("server did not return exact range")
                m=re.fullmatch(r"bytes 0-(\d+)/(\d+)",r.headers.get("Content-Range",""));assert m and int(m[1])<cap
                b=bytearray()
                for c in r.iter_content(8192):b.extend(c);assert len(b)<=cap
                assert len(b)==int(m[1])+1
                h,n=fil_header(bytes(b));g=fil_grid(h,n,int(m[2]))
                assert g["source"].strip().upper()==v["target"],"exact header source mismatch"
                assert abs(g["mjd"]-v["mjd"])*86400<=3,"header time mismatch"
                df=abs(g["df"])*1e6;ts=g["ts"];gg=S["grid"];assert math.isfinite(df) and math.isfinite(ts) and df>0 and ts>0
                lo,hi=S["coverage_mhz"]
                a.update(status="HEADER_PARSED",header_sha256=sha(bytes(b[:n])),channel_hz=df,integration_s=ts,rows=g["shape"][0],bounds_mhz=g["bounds"],covers_methanol=bool(g["bounds"][0]<=lo and g["bounds"][1]>=hi),fine_grid=bool(df<=gg["df_hz_max"] and ts<=gg["ts_s_max"] and g["shape"][0]>=gg["min_rows"]),exact_reference_grid=bool(abs(df-gg["df_reference"])<1e-7 and abs(ts-gg["ts_reference"])<1e-7),authenticated_transport=a["scheme"]=="https",frame_flags={k:h.get(k) for k in ("barycentric","pulsarcentric")})
        except Exception as e:a.update(status="HOLD",error=type(e).__name__+": "+str(e)[:220])
        attempts.append(a)
        if a.get("authenticated_transport"):break
    return {**v,"attempts":attempts}
def source(f):
    b,rr,final=fetchcat(f["alias"]);names=collections.Counter();files={};reject=collections.Counter();coord_rows=[]
    for r in rr:
        if not isinstance(r,dict):reject["not_object"]+=1;continue
        t=str(r.get("target","")).strip().upper();names[t]+=1;ok,role=exact(t,f);d=coord(r,f)
        if "galactic" in f:coord_rows.append({"target":t,"separation_deg":d,"exact_name":ok})
        if not ok:reject["name"]+=1;continue
        if "galactic" in f and (d is None or d>.05):reject["pointing"]+=1;continue
        url=r.get("url")
        if not isinstance(url,str) or not official(url):reject["archive_host"]+=1;continue
        try:m=float(r["mjd"]);assert math.isfinite(m)
        except Exception:reject["epoch"]+=1;continue
        fmt=pathlib.PurePosixPath(urllib.parse.urlparse(url).path).suffix.lower()
        if fmt not in (".fil",".h5",".raw"):reject["format"]+=1;continue
        v={"url":url,"target":t,"role":role,"mjd":m,"day":int(m),"format":fmt,"metadata_center_freq":r.get("center_freq"),"catalogue_md5_present":bool(re.fullmatch(r"[a-fA-F0-9]{32}",str(r.get("md5sum","")))),"catalogue_pointing_separation_deg":d}
        if url in files and files[url]!=v:reject["conflict"]+=1;files[url]["conflict"]=True
        else:files[url]=v
    fs=[v for v in files.values() if not v.get("conflict")]
    def rank(v):
        try:
            c=float(v["metadata_center_freq"]);dist=abs(c-S["rest_mhz"])
            if not math.isfinite(dist):dist=1e12
        except Exception:dist=1e12
        return (".gpuspec.0000." not in v["url"],dist,v["mjd"],v["url"])
    fil=sorted([v for v in fs if v["format"]==".fil" and v["role"]=="ON"],key=rank)
    picks=[];seen=set()
    for v in fil:
        if v["day"] in seen:continue
        picks.append(v);seen.add(v["day"])
    for v in fil:
        if len(picks)>=S["limits"]["per_source_headers"]:break
        if v not in picks:picks.append(v)
    probes=[probe(v) for v in picks[:S["limits"]["per_source_headers"]]]
    good=[(v,a) for v in probes for a in v["attempts"] if a.get("fine_grid") and a.get("covers_methanol")]
    sep=[v["separation_deg"] for v in coord_rows if v["separation_deg"] is not None]
    return {"source":f["name"],"query_rows":len(rr),"catalog_sha256":sha(b),"catalog_final_scheme":urllib.parse.urlparse(final).scheme,"target_names":dict(names),"coordinate_range_deg":[min(sep),max(sep)] if sep else None,"rejected_rows":dict(reject),"unique_qualified_catalogue_urls":len(fs),"formats":dict(collections.Counter(v["format"] for v in fs)),"on":sum(v["role"]=="ON" for v in fs),"off":sum(v["role"]=="OFF" for v in fs),"md5_present":sum(v["catalogue_md5_present"] for v in fs),"header_probes":probes,"fine_covering_http":sum(not a["authenticated_transport"] for _,a in good),"fine_covering_https":sum(a["authenticated_transport"] for _,a in good),"fine_covering_days":sorted({v["day"] for v,a in good}),"full_file_checksum_verified":False}
def main():
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(source,S["fields"]))
    d={"status":"NATURAL_LEGACY_METADATA_ONLY","sources":rows,"additional_rows":sum(r["query_rows"] for r in rows[:2]),"additional_qualified_urls":sum(r["unique_qualified_catalogue_urls"] for r in rows[:2]),"fine_covering_http":sum(r["fine_covering_http"] for r in rows),"fine_covering_https":sum(r["fine_covering_https"] for r in rows),"independent_spatial_pairs":0,"full_file_checksum_verified":False,"natural_line_verified":False,"production_adoption":False}
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");print("R37F2D1A_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
