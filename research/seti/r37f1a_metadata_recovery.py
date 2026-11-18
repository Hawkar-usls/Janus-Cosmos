#!/usr/bin/env python3
"""R37F1A: bounded archival metadata recovery and natural sky control headers."""
import datetime as dt,hashlib,json,pathlib,re,urllib.parse,xml.etree.ElementTree as ET,concurrent.futures
import requests
import r9_expanded_search as eng
import r37f1_science_archive as prev
P=pathlib.Path("research/seti/r37f1a_metadata_recovery_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f1a");O.mkdir(parents=True,exist_ok=True)
HOSTS={"pds-ppi.igpp.ucla.edu","atmos.nmsu.edu","pds-atmospheres.nmsu.edu","seti.berkeley.edu"}
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,x):(O/n).write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def get(url,params=None):
    assert urllib.parse.urlparse(url).hostname in HOSTS
    with requests.get(url,params=params,stream=True,timeout=(15,40)) as r:
        r.raise_for_status();assert urllib.parse.urlparse(r.url).hostname in HOSTS
        b=bytearray()
        for x in r.iter_content(65536):
            b.extend(x);assert len(b)<=S["limits"]["metadata_bytes"]
        return bytes(b),r.url
def probe(u,parse=lambda x:{}):
    try:
        b,v=get(u);return {"status":"PASS","url":v,"sha256":sha(b),"bytes":len(b),**parse(b)}
    except Exception as e:return {"status":"HOLD","url":u,"error":type(e).__name__+": "+str(e)[:240]}
def vg(c):
    root="https://pds-ppi.igpp.ucla.edu/data/"+c["dataset"]+"/"
    def parse(b):
        d=prev.label(b);assert d["DATA_SET_ID"]==c["dataset"];assert int(d["RECORD_BYTES"])==5056
        span=(dt.datetime.fromisoformat(d["STOP_TIME"])-dt.datetime.fromisoformat(d["START_TIME"])).total_seconds()
        assert span==c["nominal_span_s"]
        return {"dataset":d["DATA_SET_ID"],"nominal_span_s":span,"record_count":int(d["FILE_RECORDS"]),"expected_bytes":int(d["FILE_RECORDS"])*5056,"continuity_verified":False}
    return {"dataset":c["dataset"],"listing":probe(root+"DATA/"),"direct_label":probe(root+"DATA/"+c["file"],parse),"checksum_manifest":probe(root+"CHECKSUMS.TXT")}
def hrefs(b):return sorted(set(re.findall(r'href=["\']([^"\']+)["\']',b.decode(errors="replace"),re.I)))
def xml_info(b):
    assert b"<!DOCTYPE" not in b.upper() and b"<!ENTITY" not in b.upper()
    root=ET.fromstring(b)
    def values(tag):return [e.text.strip() for e in root.iter() if e.tag.rsplit("}",1)[-1]==tag and e.text and e.text.strip()]
    return {"root":root.tag.rsplit("}",1)[-1],"logical_identifier":values("logical_identifier")[:3],"start":values("start_date_time")[:4],"stop":values("stop_date_time")[:4],"files_declared":len(values("file_name"))}
def cassini():
    c=S["cassini"];p=c["pds4"];q=c["pds3"]
    a={"pds4_collection":probe(p+c["collection"],xml_info),"pds4_index":probe(p+c["index"],lambda b:{"nonempty_lines":len([s for s in b.splitlines() if s.strip()])}),"pds3_listing":probe(q,lambda b:{"day_directories":[x for x in hrefs(b) if re.fullmatch(r"sroc1_\d{3}/",x)]}),"products":[]}
    # Follow only official, bounded year/day directory names; one detached label per chosen day.
    for day in (123,141,159,177):
        r={"day_of_year":day,"attempts":[]}
        try:
            b,_=get(p+"2005/");dirs=[x for x in hrefs(b) if re.fullmatch(r"(?:\d{3}|2005[-_]?\d{3})/",x) and int(x.rstrip("/")[-3:])==day]
            for d in dirs[:1]:
                raw,_=get(p+"2005/"+d);xs=[x for x in hrefs(raw) if re.fullmatch(r"[A-Za-z0-9_.-]+\.xml",x,re.I)]
                if xs:r["attempts"].append(probe(p+"2005/"+d+xs[0],xml_info))
        except Exception as e:r["attempts"].append({"status":"HOLD","error":str(e)[:240]})
        r["legacy_day_directory"]=probe(q+f"sroc1_{day}/",lambda b:{"child_directory_count":sum(bool(re.fullmatch(r"[A-Za-z0-9_-]+/",x)) for x in hrefs(b))})
        a["products"].append(r)
    return a
def natural():
    cfg=S["natural"];queries=[];files={}
    def norm(s):return re.sub(r"[^A-Z0-9]","",str(s).upper())
    def query(alias):
        b,_=get("https://seti.berkeley.edu/opendata/api/query-files",{"target":alias,"telescope":"GBT","file-type":"HDF5"})
        d=json.loads(b);rows=d if isinstance(d,list) else next(d[k] for k in ("data","files","results") if isinstance(d.get(k),list))
        out=[]
        for r in rows:
            if not isinstance(r,dict):continue
            u=r.get("url")
            if not isinstance(u,str) or not eng.allowed(u):continue
            n=urllib.parse.unquote(urllib.parse.urlparse(u).path.rsplit("/",1)[-1])
            if not n.lower().endswith(".h5"):continue
            m=re.fullmatch(r"GBT_(\d{5})_(\d{5})_(.+)_fine\.h5",n,re.I)
            if not m:
                m=re.search(r"guppi_(\d{5})_(\d{5})_(.+?)(?:_fine|\.fine)\.h5$",n,re.I)
            if not m or norm(m[3])!=norm(alias):continue
            fm=int(m[1])+int(m[2])/86400
            try:
                v=float(r.get("mjd"))
                if abs(v-fm)>3/86400:continue
            except (TypeError,ValueError):continue
            out.append({"url":u,"name":n,"mjd":fm,"target":m[3],"alias":alias})
        return {"alias":alias,"status":"PASS","rows":len(rows),"fine_exact_alias":len(out),"epoch_candidates":sum(abs(v["mjd"]-cfg["epoch"])<=.02 for v in out),"sha256":sha(b)},out
    for alias in cfg["aliases"]:
        try:r,fs=query(alias);queries.append(r);files.update({v["url"]:v for v in fs})
        except Exception as e:queries.append({"alias":alias,"status":"HOLD","error":str(e)[:240]})
    vals=list(files.values());oh={norm(x) for x in cfg["oh_aliases"]}
    exact=sorted([v for v in vals if norm(v["target"]) in oh and abs(v["mjd"]-cfg["epoch"])<=.02],key=lambda x:(abs(x["mjd"]-cfg["epoch"]),x["url"]))
    chosen=exact[:8];used={v["url"] for v in chosen}
    for alias in ["W3","W3OH","DR21","NGC7538"]:
        q=sorted([v for v in vals if norm(v["target"])==norm(alias)],key=lambda x:(x["mjd"],x["url"]))
        if q:
            for v in (q[0],q[-1]):
                if v["url"] not in used:chosen.append(v);used.add(v["url"])
    headers=[]
    for v in chosen[:cfg["max_headers"]]:
        r={**v,"status":"HOLD_HEADER"}
        try:
            u=v["url"];p=urllib.parse.urlparse(u)
            if p.scheme=="http":u=urllib.parse.urlunparse(p._replace(scheme="https"))
            h=eng.hdr({"url":u});assert norm(h["source"])==norm(v["target"])
            assert abs(h["mjd"]-v["mjd"])<=3/86400 and len(h["shape"])==3
            df=abs(h["df"])*1e6;ts=h["ts"];assert df>0 and ts>0
            r.update({"status":"HEADER_VERIFIED","header":h,"channel_hz":df,"integration_s":ts,"nominal_duration_s":ts*h["shape"][0],"comparable_grid":df<=10 and ts<=60 and h["shape"][0]>=2,"covers_oh":[f for f in cfg["frequency_mhz"] if h["bounds"][0]<=f<=h["bounds"][1]],"exact_literature_epoch":norm(v["target"]) in oh and abs(v["mjd"]-cfg["epoch"])<=.02,"line_presence_verified":False})
        except Exception as e:r["error"]=type(e).__name__+": "+str(e)[:250]
        headers.append(r)
    return {"queries":queries,"unique_exact_alias_fine_files":len(vals),"exact_epoch_candidates":len(exact),"headers":headers,"qualified_grid_headers":sum(bool(r.get("comparable_grid")) for r in headers),"oh_coverage_headers":sum(bool(r.get("comparable_grid") and r.get("covers_oh") and norm(r["target"]) in oh) for r in headers)}
def main():
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as p:
        a=p.submit(lambda:[vg(c) for c in S["voyager"]]);b=p.submit(cassini);c=p.submit(natural)
        result={"status":"METADATA_REPAIR_AUDIT_ONLY","voyager":a.result(),"cassini":b.result(),"natural":c.result(),"raw_power_read":False,"continuity_verified":False,"source_identity_verified":False,"production_adoption":False}
    save("result.json",result);print("R37F1A_SUMMARY",json.dumps(result),flush=True)
if __name__=="__main__":main()
