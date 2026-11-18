#!/usr/bin/env python3
"""R37F1D/E: archive chronology metadata and natural sky header qualification."""
import collections,datetime as dt,hashlib,io,json,math,os,pathlib,re,struct,sys,urllib.parse,zipfile
import requests
import r37f1bc_metadata as prev
P=pathlib.Path("research/seti/r37f1de_archive_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r37f1de");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(name,x):(O/name).write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def parent(key,name):
    a=S["parent"][key]
    r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['artifact']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=(15,90))
    r.raise_for_status();b=r.content;assert len(b)<16*1024**2 and sha(b)==a["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        n=[i for i in z.infolist() if pathlib.PurePosixPath(i.filename).name==name]
        assert len(n)==1 and n[0].file_size<8*1024**2
        return json.loads(z.read(n[0]))
def cassini():
    p=parent("cassini","cassini.json");rows=[]
    for day in p["days"]:
        for r in day["labels"]:
            if r["status"]!="LABEL_VERIFIED":continue
            a={"day_of_year":day["day_of_year"],"status":"HOLD"}
            try:
                b=prev.fetch(r["url"]);assert sha(b)==r["sha256"]
                q=prev.detached_label(b);s=b.decode("ascii",errors="strict")
                assert q["start"]==r["start"] and q["stop"]==r["stop"]
                f=q["fields"];n=str(f.get("PRODUCT_ID") or "")+" "+pathlib.PurePosixPath(urllib.parse.urlparse(r["url"]).path).name
                upper=s.upper()
                is_rsr=bool(re.search(r"\b(?:RADIO SCIENCE RECEIVER|RSR)\b",upper)) and bool(re.search(r"(?im)^\s*(?:\^TABLE|OBJECT\s*=\s*TABLE|RECORD_BYTES)\s*=",s))
                # Never infer a record cadence from a sample rate or from endpoint division.
                intervals=[]
                for m in re.finditer(r"(?im)^\s*(?:RECORD_INTERVAL|RECORD_DURATION|RECORD_TIME_INTERVAL)\s*=\s*([0-9.eE+-]+)\s*<(?:S|SECOND|SECONDS)>\s*$",s):
                    v=float(m[1])
                    if math.isfinite(v) and v>0:intervals.append(v)
                nominal=None
                if len(set(intervals))==1 and f.get("FILE_RECORDS") and is_rsr:
                    try:
                        count=int(f["FILE_RECORDS"]);v=intervals[0];assert count>0
                        nominal={"declared_interval_s":v,"count":count,"count_span_s":count*v,"endpoint_span_s":q["nominal_duration_s"],"difference_s":count*v-q["nominal_duration_s"]}
                    except (ValueError,AssertionError):pass
                a.update(status="LABEL_REAUTHENTICATED",label_sha256=sha(b),nominal_duration_s=q["nominal_duration_s"],rsr_candidate=is_rsr,has_time_field_definition=bool(re.search(r"(?i)(?:TIME_TAG|TIME_STAMP|SECOND_OF_DAY|TIME OF FIRST SAMPLE)",s)),has_record_count=bool(f.get("FILE_RECORDS")),has_record_bytes=bool(f.get("RECORD_BYTES")),declared_interval_audit=nominal)
            except Exception as e:a["error"]=type(e).__name__+": "+str(e)[:180]
            rows.append(a)
    ok=[v for v in rows if v["status"]=="LABEL_REAUTHENTICATED"]
    rsr=[v for v in ok if v["rsr_candidate"]]
    d={"status":"LABEL_CHRONOLOGY_AUDIT_ONLY","parent_count":p["verified_labels"],"reauthenticated":len(ok),"rsr_candidates":len(rsr),"rsr_nominal_duration_min_s":min((v["nominal_duration_s"] for v in rsr),default=None),"rsr_nominal_duration_max_s":max((v["nominal_duration_s"] for v in rsr),default=None),"rsr_with_record_count":sum(v["has_record_count"] for v in rsr),"rsr_with_documented_time_field":sum(v["has_time_field_definition"] for v in rsr),"explicit_interval_comparisons":sum(v["declared_interval_audit"] is not None for v in rsr),"actual_continuity_verified":False,"raw_receiver_data_read":False,"rows":rows}
    save("cassini.json",d);print("R37F1D_SUMMARY",json.dumps({k:v for k,v in d.items() if k!="rows"}),flush=True)
STR={"source_name","rawdatafile"}
I32={"machine_id","telescope_id","data_type","barycentric","pulsarcentric","ibeam","nbeams","nbits","nchans","nifs","nsamples"}
F64={"az_start","za_start","src_raj","src_dej","tstart","tsamp","fch1","foff","refdm","period"}
def fil_header(b):
    f=io.BytesIO(b)
    def read(n):
        q=f.read(n)
        if len(q)!=n:raise ValueError("truncated header")
        return q
    def token():
        n=struct.unpack("<i",read(4))[0]
        if not 1<=n<=4096:raise ValueError("invalid string length")
        return read(n).decode("ascii")
    if token()!="HEADER_START":raise ValueError("missing HEADER_START")
    d={}
    for _ in range(128):
        k=token()
        if k=="HEADER_END":return d,f.tell()
        if k in d:raise ValueError("duplicate header key")
        if k in STR:v=token()
        elif k in I32:v=struct.unpack("<I",read(4))[0]
        elif k in F64:v=struct.unpack("<d",read(8))[0]
        elif k=="signed":v=struct.unpack("<b",read(1))[0]
        else:raise ValueError("unsupported SIGPROC key "+k)
        if isinstance(v,float) and not math.isfinite(v):raise ValueError("nonfinite field")
        d[k]=v
    raise ValueError("unterminated header")
def fil_grid(d,offset,total):
    if not (d["nchans"]>0 and d["nifs"]>0 and d["nbits"] in (1,2,4,8,16,32,64) and d["tsamp"]>0 and d["foff"]!=0 and total>offset):raise ValueError("invalid grid")
    bits=(total-offset)*8;rowbits=d["nchans"]*d["nifs"]*d["nbits"]
    if bits%rowbits:raise ValueError("nonintegral row count")
    rows=bits//rowbits
    if d.get("nsamples",0) not in (0,rows):raise ValueError("nsamples disagreement")
    return {"source":d["source_name"],"mjd":d["tstart"],"ts":d["tsamp"],"df":d["foff"],"shape":[rows,d["nifs"],d["nchans"]],"bounds":sorted([d["fch1"],d["fch1"]+(d["nchans"]-1)*d["foff"]]),"header_bytes":offset,"size_bytes":total,"row_count_from_file_size":True}
def canaries():
    def token(v):v=v.encode();return struct.pack("<i",len(v))+v
    d={"source_name":"W3","nchans":8,"nifs":1,"nbits":32,"tstart":57424.92,"tsamp":20.,"fch1":1668.,"foff":-0.000003}
    b=token("HEADER_START")+b"".join(token(k)+(token(v) if k in STR else struct.pack("<I",v) if k in I32 else struct.pack("<d",v)) for k,v in d.items())+token("HEADER_END")
    q,n=fil_header(b);assert q==d and n==len(b) and fil_grid(q,n,n+8*4*5)["shape"]==[5,1,8]
    for bad in (b[:-1],token("HEADER_START")+token("unknown")):
        try:fil_header(bad)
        except ValueError:pass
        else:raise AssertionError("invalid header accepted")
    try:fil_grid(q,n,n+1)
    except ValueError:pass
    else:raise AssertionError("invalid row size accepted")
def source_role(s):
    s=str(s).strip().upper();s=re.sub(r"^(?:DIAG|CAL[0-9]+)_","",s)
    role="OFF" if s.endswith("_OFF") else "ON"
    if role=="OFF":s=s[:-4]
    return s,role
def get_fil(u,eng):
    cap=S["natural"]["max_fil_prefix_bytes"]
    with requests.get(u,headers={"Range":f"bytes=0-{cap-1}","Accept-Encoding":"identity"},stream=True,timeout=(15,45)) as r:
        r.raise_for_status();assert eng.allowed(r.url) and r.status_code==206
        m=re.fullmatch(r"bytes 0-(\d+)/(\d+)",r.headers.get("Content-Range",""));assert m and int(m[1])<cap
        b=bytearray()
        for c in r.iter_content(8192):
            b.extend(c);assert len(b)<=cap
        assert len(b)==int(m[1])+1
    d,n=fil_header(bytes(b));return fil_grid(d,n,int(m[2]))
def natural():
    import r9_expanded_search as eng
    canaries();c=S["natural"];p=parent("natural","natural.json")
    queries=[];files={}
    for alias in c["aliases"]:
        try:
            b=prev.fetch("https://seti.berkeley.edu/opendata/api/query-files",{"target":alias,"telescope":"GBT","file-type":"HDF5"})
            j=json.loads(b);rr=j if isinstance(j,list) else next(j[k] for k in ("data","files","results") if isinstance(j.get(k),list))
            n=0
            for r in rr:
                if not isinstance(r,dict):continue
                u=r.get("url");t=str(r.get("target","")).strip().upper()
                if not isinstance(u,str) or not eng.allowed(u):continue
                field,role=source_role(t)
                if field not in c["aliases"]:continue
                ext=pathlib.PurePosixPath(urllib.parse.urlparse(u).path).suffix.lower()
                if ext not in (".fil",".h5",".hdf5"):continue
                try:m=float(r["mjd"]);assert 50000<=m<=65000
                except (ValueError,TypeError,AssertionError,KeyError):continue
                v={"url":u,"catalogue_target":t,"field":field,"role":role,"mjd":m,"format":"fil" if ext==".fil" else "h5"}
                if u in files and files[u]!=v:files[u]["conflict"]=True
                else:files[u]=v
                n+=1
            queries.append({"alias":alias,"rows":len(rr),"qualified_catalogue_rows":n,"sha256":sha(b)})
        except Exception as e:queries.append({"alias":alias,"error":str(e)[:180]})
    vals=[v for v in files.values() if not v.get("conflict")];chosen=[]
    def order(v):return (".gpuspec.0000." not in v["url"],v["mjd"],v["url"])
    for fmt in ("h5","fil"):
        fs=[v for v in vals if v["format"]==fmt];sel={}
        exact=sorted([v for v in fs if v["role"]=="ON" and v["field"]=="W3" and abs(v["mjd"]-c["epoch"])<=.02],key=order)
        for v in exact[:3]:sel[v["url"]]=v
        for field in c["aliases"]:
            on=[v for v in fs if v["field"]==field and v["role"]=="ON"]
            days=sorted({int(v["mjd"]) for v in on})
            for day in sorted(set(days[:1]+days[-1:])):
                q=sorted([v for v in on if int(v["mjd"])==day],key=order)
                if q:sel[q[0]["url"]]=q[0]
                off=sorted([v for v in fs if v["field"]==field and v["role"]=="OFF" and int(v["mjd"])==day],key=order)
                if off:sel[off[0]["url"]]=off[0]
        chosen.extend(list(sel.values())[:c["max_headers_per_format"]])
    results=[]
    for v in chosen:
        r={**v,"status":"HOLD_HEADER"}
        try:
            u=v["url"];p=urllib.parse.urlparse(u)
            if p.scheme=="http":u=urllib.parse.urlunparse(p._replace(scheme="https"))
            h=get_fil(u,eng) if v["format"]=="fil" else eng.hdr({"url":u})
            assert str(h["source"]).strip().upper()==v["catalogue_target"],"exact catalogue/header source disagreement"
            assert abs(h["mjd"]-v["mjd"])<=3/86400 and len(h["shape"])==3,"header epoch/shape disagreement"
            df=abs(h["df"])*1e6;ts=h["ts"];assert math.isfinite(df) and math.isfinite(ts) and df>0 and ts>0
            g=c["grid"];grid=df<=g["max_channel_hz"] and ts<=g["max_integration_s"] and h["shape"][0]>=g["min_rows"]
            r.update(status="HEADER_VERIFIED",channel_hz=df,integration_s=ts,rows=h["shape"][0],nominal_duration_s=ts*h["shape"][0],comparable_grid=grid,oh_coverage=[f for f in c["lines_mhz"] if h["bounds"][0]<=f<=h["bounds"][1]],exact_publication_epoch=v["field"]=="W3" and v["role"]=="ON" and abs(v["mjd"]-c["epoch"])<=.02,line_presence_verified=False)
        except Exception as e:r["error"]=type(e).__name__+": "+str(e)[:200]
        results.append(r)
    ok=[v for v in results if v["status"]=="HEADER_VERIFIED"]
    good=[v for v in ok if v["role"]=="ON" and v["comparable_grid"] and v["oh_coverage"]]
    d={"status":"NATURAL_HEADER_ONLY","parser_canaries":"PASS","queries":queries,"unique_candidate_urls":len(vals),"selected":len(results),"verified_headers":len(ok),"on_oh_grid_headers":len(good),"on_oh_grid_dates":sorted({int(v["mjd"]) for v in good}),"publication_epoch_matches":sum(v["exact_publication_epoch"] for v in ok),"headers":results,"power_decoded":False,"natural_line_verified":False,"production_adoption":False}
    save("natural.json",d);print("R37F1E_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":
    {"cassini":cassini,"natural":natural}[sys.argv[1]]()
