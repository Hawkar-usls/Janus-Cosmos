#!/usr/bin/env python3
"""Historical science archive metadata checks; no receiver voltage or power processing."""
import datetime as dt,hashlib,html,json,pathlib,re,urllib.parse,concurrent.futures
import requests
P=pathlib.Path("research/seti/r37f1_science_archive_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f1");O.mkdir(parents=True,exist_ok=True)
def save(n,x):(O/n).write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def sha(b):return hashlib.sha256(b).hexdigest()
ALLOWED={"pds-ppi.igpp.ucla.edu","pds.nasa.gov","seti.berkeley.edu"}
def get(u,params=None):
    assert urllib.parse.urlparse(u).hostname in ALLOWED
    with requests.get(u,params=params,stream=True,timeout=(15,45)) as r:
        r.raise_for_status();assert urllib.parse.urlparse(r.url).hostname in ALLOWED
        b=bytearray()
        for x in r.iter_content(65536):
            b.extend(x);assert len(b)<=S["selection"]["max_metadata_bytes"]
        raw=bytes(b);return raw,r.url
def label(raw):
    keys=["RECORD_TYPE","RECORD_BYTES","FILE_RECORDS","INSTRUMENT_HOST_NAME","TARGET_NAME","INSTRUMENT_ID","DATA_SET_ID","PRODUCT_ID","FILE_NAME","START_TIME","STOP_TIME","DSN_STATION_NUMBER","^DESCRIPTION"]
    s=raw.decode("ascii",errors="strict");d={}
    for k in keys:
        m=re.search(r"(?m)^\s*"+re.escape(k)+r"\s*=\s*(?:\"([^\"]*)\"|([^\r\n]+))",s)
        if not m:raise ValueError("missing field "+k)
        d[k]=(m[1] if m[1] is not None else m[2]).strip()
    return d
def head(u):
    try:
        r=requests.head(u,allow_redirects=True,timeout=(15,30));r.raise_for_status()
        assert urllib.parse.urlparse(r.url).hostname in ALLOWED
        return {"status":r.status_code,"content_length":int(r.headers["Content-Length"]) if r.headers.get("Content-Length","").isdigit() else None}
    except Exception as e:return {"error":str(e)[:250]}
def archive(cfg):
    root=cfg["url"];durl=root+"DATA/";raw,_=get(durl)
    hrefs=set(html.unescape(v) for v in re.findall(r'href=["\']([^"\']+)["\']',raw.decode("utf-8",errors="replace"),re.I))
    labels=sorted(v for v in hrefs if re.fullmatch(r"[A-Z0-9]+\.LBL",v,re.I) and v[:-4]+".ODR" in hrefs)
    assert labels
    count=min(len(labels),S["selection"]["max_labels_per_archive"])
    ids=sorted(set(round(i*(len(labels)-1)/max(1,count-1)) for i in range(count)))
    try:
        checksums,_=get(root+"CHECKSUMS.TXT");(O/(cfg["id"]+"-checksums.txt")).write_bytes(checksums)
        cs=checksums.decode("ascii",errors="replace")
    except Exception:checksums=None;cs=""
    def audit(n):
        row={"name":n,"status":"HOLD"}
        try:
            b,_=get(durl+n);v=label(b)
            assert v["INSTRUMENT_HOST_NAME"]==cfg["host"] and v["DATA_SET_ID"]==cfg["id"]
            assert v["FILE_NAME"]==n[:-4]+".ODR" and v["PRODUCT_ID"]==v["FILE_NAME"]
            assert v["RECORD_TYPE"]=="FIXED_LENGTH"
            rb=int(v["RECORD_BYTES"]);nr=int(v["FILE_RECORDS"]);assert rb==5056 and nr>0
            t=dt.datetime.fromisoformat(v["START_TIME"]);u=dt.datetime.fromisoformat(v["STOP_TIME"])
            span=(u-t).total_seconds();assert span>0
            h=head(durl+v["FILE_NAME"])
            matching=[line for line in cs.splitlines() if re.search(r"\b"+re.escape(v["FILE_NAME"])+r"\b",line,re.I)]
            row.update({"status":"LABEL_VERIFIED","label_sha256":sha(b),"label":v,"nominal_span_s":span,"expected_bytes":rb*nr,"head":h,"head_size_matches":h.get("content_length")==rb*nr,"published_checksum_mentions":len(matching),"continuity":"UNVERIFIED_RECORD_COUNTERS","independent_session":"UNVERIFIED_GROUPED_BY_CAMPAIGN"})
        except Exception as e:row["error"]=type(e).__name__+": "+str(e)[:300]
        return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as p:rows=list(p.map(audit,[labels[i] for i in ids]))
    valid=[r for r in rows if r["status"]=="LABEL_VERIFIED"]
    return {"archive":cfg["id"],"listing_sha256":sha(raw),"paired_labels":len(labels),"sampled_labels":len(rows),"verified_labels":len(valid),"head_sizes_match":sum(r["head_size_matches"] for r in valid),"checksum_mentions":sum(r["published_checksum_mentions"]>0 for r in valid),"checksum_manifest_sha256":sha(checksums) if checksums else None,"nominal_spans_s":[r["nominal_span_s"] for r in valid],"declared_days":sorted({r["label"]["START_TIME"][:10] for r in valid}),"rows":rows}
def natural():
    c=S["natural_control"];out=[]
    for name in c["query_aliases"]:
        row={"alias":name}
        try:
            raw,u=get(c["api"],{"target":name,"telescope":"GBT","file-type":"HDF5"})
            d=json.loads(raw);rr=d if isinstance(d,list) else next(d[k] for k in ("data","files","results") if isinstance(d.get(k),list))
            rec=[]
            for r in rr:
                if not isinstance(r,dict):continue
                try:mjd=float(r.get("mjd"));assert abs(mjd-c["mjd"])<=.02
                except (ValueError,TypeError,AssertionError):continue
                url=r.get("url","");p=urllib.parse.urlparse(url)
                if p.scheme!="https" or not url.lower().endswith(".h5"):continue
                rec.append({"mjd":mjd,"url":url,"metadata":r})
            row.update({"status":"METADATA_QUERY_OK","response_sha256":sha(raw),"count":len(rr),"publication_epoch_candidates":rec[:40]})
        except Exception as e:row.update({"status":"HOLD_QUERY","error":str(e)[:250]})
        out.append(row)
    return out
def main():
    archives=[]
    for c in S["archives"]:
        try:archives.append(archive(c))
        except Exception as e:archives.append({"archive":c["id"],"status":"HOLD","error":type(e).__name__+": "+str(e)[:300]})
    w=natural()
    d={"status":"METADATA_AUDIT_COMPLETE__CONTINUITY_AND_MATCHED_GRID_HOLD","archives":archives,"natural_queries":w,"natural_epoch_candidate_count":len({v["url"] for r in w for v in r.get("publication_epoch_candidates",[])}),"cassini":S["cassini"],"natural_control_publication":S["natural_control"]["publication"],"continuity_verified":False,"matched_frequency_verified":False,"independent_same_source_sessions_verified":False,"natural_raw_signal_verified":False,"raw_power_read":False,"production_adoption":False}
    save("result.json",d)
    print("R37F1_SUMMARY",json.dumps({**d,"archives":[{k:v for k,v in r.items() if k!="rows"} for r in archives],"natural_queries":[{**r,"publication_epoch_candidates":len(r.get("publication_epoch_candidates",[]))} for r in w]}),flush=True)
if __name__=="__main__":main()
