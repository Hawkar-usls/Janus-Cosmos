#!/usr/bin/env python3
"""R37F1BC: independent historical science metadata and natural calibration header audits."""
import collections,csv,datetime as dt,hashlib,io,json,pathlib,re,sys,urllib.parse
import requests
P=pathlib.Path("research/seti/r37f1bc_metadata_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r37f1bc");O.mkdir(parents=True,exist_ok=True)
ALLOWED={"atmos.nmsu.edu","pds-atmospheres.nmsu.edu","seti.berkeley.edu"}
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def fetch(u,params=None):
    assert urllib.parse.urlparse(u).hostname in ALLOWED
    with requests.get(u,params=params,stream=True,timeout=(15,40)) as r:
        r.raise_for_status();assert urllib.parse.urlparse(r.url).hostname in ALLOWED
        b=bytearray()
        for c in r.iter_content(65536):
            b.extend(c);assert len(b)<=16*1024**2
        return bytes(b)
def hrefs(b):return sorted(set(re.findall(r'href=["\x27]([^"\x27]+)["\x27]',b.decode("utf8",errors="replace"),re.I)))
def parse_pds_time(v):
    v=v.strip().strip('"').rstrip("Z")
    try:return dt.datetime.fromisoformat(v)
    except ValueError:
        for fmt in ("%Y-%jT%H:%M:%S.%f","%Y-%jT%H:%M:%S"):
            try:return dt.datetime.strptime(v,fmt)
            except ValueError:pass
        raise
def detached_label(b):
    s=b.decode("ascii",errors="strict")
    def field(k):
        m=re.search(r"(?m)^\s*"+re.escape(k)+r"\s*=\s*([^\r\n]+)",s)
        return m[1].strip().strip('"') if m else None
    a=field("START_TIME");z=field("STOP_TIME")
    assert a and z and field("DATA_SET_ID")
    start=parse_pds_time(a);stop=parse_pds_time(z);span=(stop-start).total_seconds()
    assert span>0
    fields={k:field(k) for k in ("DATA_SET_ID","PRODUCT_ID","RECORD_TYPE","RECORD_BYTES","FILE_RECORDS","INSTRUMENT_HOST_NAME","INSTRUMENT_ID","TARGET_NAME")}
    assert "CO-S-RSS" in fields["DATA_SET_ID"].upper()
    expected=None
    if fields["RECORD_BYTES"] and fields["FILE_RECORDS"]:
        try:expected=int(fields["RECORD_BYTES"])*int(fields["FILE_RECORDS"])
        except ValueError:pass
    return {"start":start.isoformat(),"stop":stop.isoformat(),"nominal_duration_s":span,"fields":fields,"expected_bytes":expected,"continuity_verified":False}
def cassini():
    c=S["cassini"];b=fetch(c["pds4"]);rr=list(csv.reader(io.StringIO(b.decode("utf8-sig") if b.startswith(b"\xef\xbb\xbf") else b.decode("utf8"))))
    # Index inventory only, no inferred path from logical identifier.
    catalog={"sha256":sha(b),"rows":len(rr),"nonempty_rows":sum(bool(r) for r in rr),"column_counts":dict(collections.Counter(str(len(r)) for r in rr)),"unique_rows":len({tuple(r) for r in rr}),"examples":rr[:5]}
    days=[]
    for day in c["days"]:
        root=c["base"]+f"sroc1_{day}/";queue=[(root,0)];visited=set();seen=[];errors=[];label_urls=[]
        while queue and len(visited)<c["max_dirs_per_day"]:
            u,depth=queue.pop(0)
            if u in visited:continue
            visited.add(u)
            try:
                raw=fetch(u);ls=hrefs(raw);seen.append({"url":u,"sha256":sha(raw),"links":len(ls)})
                for v in ls:
                    if re.fullmatch(r"[A-Za-z0-9_.-]+\.lbl",v,re.I):label_urls.append(u+v)
                    elif depth<c["max_depth"] and re.fullmatch(r"[A-Za-z0-9_-]+/",v):queue.append((u+v,depth+1))
            except Exception as e:errors.append({"url":u,"error":str(e)[:220]})
        labels=[]
        for u in sorted(set(label_urls))[:c["max_labels_per_day"]]:
            try:
                raw=fetch(u);d=detached_label(raw)
                # Archive date controls grouping; do not treat individual labels as independent sessions.
                assert dt.datetime.fromisoformat(d["start"]).timetuple().tm_yday==day
                labels.append({"url":u,"status":"LABEL_VERIFIED","sha256":sha(raw),**d})
            except Exception as e:labels.append({"url":u,"status":"HOLD","error":type(e).__name__+": "+str(e)[:220]})
        days.append({"day_of_year":day,"directories":len(visited),"directory_receipts":seen,"label_candidates":len(set(label_urls)),"labels":labels,"errors":errors})
    valid=[v for d in days for v in d["labels"] if v["status"]=="LABEL_VERIFIED"]
    d={"status":"CATALOG_AND_DETACHED_LABEL_AUDIT_ONLY","catalog":catalog,"days":days,"verified_labels":len(valid),"days_with_verified_labels":sum(any(x["status"]=="LABEL_VERIFIED" for x in d["labels"]) for d in days),"nominal_durations_s":sorted(set(v["nominal_duration_s"] for v in valid)),"continuity_verified":False,"raw_data_read":False}
    save("cassini.json",d);print("R37F1B_SUMMARY",json.dumps({**d,"days":[{k:v for k,v in a.items() if k not in ("directory_receipts","labels","errors")} for a in days]}),flush=True)
def natural():
    import r9_expanded_search as eng
    c=S["berkeley"];queries=[];files={}
    def norm(x):return re.sub(r"[^A-Z0-9]","",str(x).upper())
    for alias in c["aliases"]:
        q={"alias":alias}
        try:
            b=fetch(c["api"],{"target":alias,"telescope":"GBT","file-type":"HDF5"})
            j=json.loads(b);rows=j if isinstance(j,list) else next(j[k] for k in ("data","files","results") if isinstance(j.get(k),list))
            suffixes=collections.Counter();keys=collections.Counter();url_count=0;ok=0;examples=[]
            for r in rows:
                if not isinstance(r,dict):continue
                keys.update(r.keys())
                u=r.get("url")
                if not isinstance(u,str):continue
                url_count+=1;p=urllib.parse.urlparse(u);name=urllib.parse.unquote(p.path.rsplit("/",1)[-1])
                suffixes.update([pathlib.PurePosixPath(name).suffix.lower()])
                if len(examples)<5:examples.append({"name":name,"mjd":r.get("mjd"),"target":r.get("target"),"source_name":r.get("source_name"),"metadata_keys":sorted(r)})
                if not eng.allowed(u) or not name.lower().endswith((".h5",".hdf5")):continue
                try:m=float(r["mjd"]);assert 50000<=m<=65000
                except (TypeError,ValueError,KeyError,AssertionError):continue
                # Do not require a guessed filename grammar. Header source name is the qualification.
                f=files.setdefault(u,{"url":u,"name":name,"mjd":m,"aliases":[]})
                if f["mjd"]!=m:f["mjd_conflict"]=True
                if alias not in f["aliases"]:f["aliases"].append(alias)
                ok+=1
            q.update({"status":"PASS","sha256":sha(b),"rows":len(rows),"url_count":url_count,"hdf5_with_mjd":ok,"suffix_counts":dict(suffixes),"keys":dict(keys),"examples":examples})
        except Exception as e:q.update({"status":"HOLD","error":str(e)[:240]})
        queries.append(q)
    vals=[v for v in files.values() if not v.get("mjd_conflict")]
    exact=sorted([v for v in vals if any(a in ("W3","W3OH") for a in v["aliases"]) and abs(v["mjd"]-c["epoch"])<=.02],key=lambda v:(abs(v["mjd"]-c["epoch"]),v["url"]))
    chosen=exact[:8];used={v["url"] for v in chosen}
    for alias in c["aliases"]:
        q=sorted([v for v in vals if alias in v["aliases"]],key=lambda v:(v["mjd"],v["url"]))
        if q:
            for v in (q[0],q[-1]):
                if v["url"] not in used:chosen.append(v);used.add(v["url"])
    results=[]
    for v in chosen[:c["max_headers"]]:
        r={**v,"status":"HOLD_HEADER"}
        try:
            u=v["url"];p=urllib.parse.urlparse(u)
            if p.scheme=="http":u=urllib.parse.urlunparse(p._replace(scheme="https"))
            h=eng.hdr({"url":u});assert len(h["shape"])==3
            assert any(norm(h["source"])==norm(a) for a in v["aliases"]),"catalogue and actual source disagree"
            assert abs(h["mjd"]-v["mjd"])<=3/86400,"catalogue and header epoch disagree"
            df=abs(h["df"])*1e6;ts=h["ts"];g=c["grid"]
            assert df>0 and ts>0
            r.update({"status":"HEADER_VERIFIED","header":h,"channel_hz":df,"integration_s":ts,"duration_s":ts*h["shape"][0],"comparable_grid":df<=g["max_df_hz"] and ts<=g["max_ts_s"] and h["shape"][0]>=g["min_rows"],"oh_coverage":[f for f in c["line_mhz"] if h["bounds"][0]<=f<=h["bounds"][1]],"line_presence_verified":False})
        except Exception as e:r["error"]=type(e).__name__+": "+str(e)[:240]
        results.append(r)
    d={"status":"NATURAL_HEADER_AUDIT_ONLY","queries":queries,"unique_hdf5_urls":len(vals),"exact_epoch_candidates":len(exact),"headers":results,"verified_headers":sum(r["status"]=="HEADER_VERIFIED" for r in results),"comparable_grid_headers":sum(bool(r.get("comparable_grid")) for r in results),"oh_grid_headers":sum(bool(r.get("comparable_grid") and r.get("oh_coverage")) for r in results),"raw_power_read":False,"natural_line_verified":False}
    save("natural.json",d);print("R37F1C_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":
    assert len(sys.argv)==2
    {"cassini":cassini,"natural":natural}[sys.argv[1]]()
