#!/usr/bin/env python3
"""Metadata-only alternative archive candidates. Not a new benchmark or label source."""
import concurrent.futures,hashlib,io,json,os,pathlib,re,urllib.parse,zipfile
import requests
P=pathlib.Path("research/seti/r37a_hannah_resolution_protocol.json");S=json.loads(P.read_text())["archive_scout"]
F=json.loads(pathlib.Path("research/seti/r36f2_pipeline_protocol.json").read_text())
OUT=pathlib.Path("runtime/r37_sources");OUT.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def body(url,params=None):
    with requests.get(url,params=params,stream=True,timeout=(15,45)) as r:
        r.raise_for_status();assert urllib.parse.urlparse(r.url).hostname=="seti.berkeley.edu"
        data=bytearray()
        for b in r.iter_content(65536):
            data.extend(b);assert len(data)<=32*1024**2
        return bytes(data),r.url
def walk(o):
    if isinstance(o,str):yield o
    elif isinstance(o,list):
        for x in o:yield from walk(x)
    elif isinstance(o,dict):
        for x in o.values():yield from walk(x)
def exclusion():
    ts=set(F["exclude_extra"])|{"HIP93069"};days=set(F["shield"]["exclude_days"])|{58000};names=set();urls=set()
    for a in F["r36f2"]["prior_manifests"]:
        r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=60);r.raise_for_status();b=r.content
        assert len(b)<16*1024**2 and sha(b)==a["digest"]
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            ns=[v for v in z.infolist() if pathlib.PurePosixPath(v.filename).name=="manifest.json"];assert len(ns)==1 and ns[0].file_size<8*1024**2
            d=json.loads(z.read(ns[0]))
        for q in d["selected"]:
            ts.add(q["target"]);days.add(int(q.get("day",q["files"][0]["day"])))
            for f in q["files"]:names.add(f["name"]);urls.add(f["url"])
    return ts,days,names,urls
def main():
    ts,days,names,urls=exclusion();raw,index=body(S["index"]);(OUT/"target_index.json").write_bytes(raw)
    strings=set(walk(json.loads(raw)));all_names=sorted({x[:-2].upper() for x in strings if re.fullmatch(r"[A-Za-z0-9_-]{2,48}_S",x,re.I)})
    aliases=[q for q in all_names if not re.fullmatch(r"HIP\d+",q) and q not in ts]
    aliases=sorted(aliases,key=lambda x:sha(("R37-SOURCE|"+x).encode()))[:32]
    gbt=[x for x in S["gbt_targets"] if x not in ts]
    jobs=[("Parkes",x,x) for x in aliases]+[("GBT",x,y) for x in gbt for y in (x,x+"_OFF")]
    def run(job):
        telescope,target,query=job;row={"telescope":telescope,"target":target,"query":query}
        try:
            b,u=body("https://seti.berkeley.edu/opendata/api/query-files",{"target":query,"telescope":telescope,"file-type":"HDF5"})
            d=json.loads(b);rr=d if isinstance(d,list) else next(d[k] for k in ("data","files","results") if isinstance(d.get(k),list))
            nm="index-"+sha(u.encode())[:16]+".json";(OUT/nm).write_bytes(b)
            matches=[]
            for f in rr:
                if not isinstance(f,dict):continue
                a=f.get("url","")
                if not isinstance(a,str):continue
                p=urllib.parse.urlparse(a)
                if not (p.hostname=="bldata.berkeley.edu" or re.fullmatch(r"blpd\d+\.ssl\.berkeley\.edu",p.hostname or "")):continue
                n=urllib.parse.unquote(p.path.rsplit("/",1)[-1])
                try:mjd=float(f.get("mjd"));day=int(mjd)
                except (ValueError,TypeError):continue
                if day in days or n in names or a in urls or not n.lower().endswith(".h5"):continue
                if telescope=="Parkes" and not re.fullmatch(r"Parkes_\d{5}_\d{5}_"+re.escape(target)+r"_[SR]_fine\.h5",n,re.I):continue
                matches.append({"name":n,"url":a,"mjd":mjd,"day":day,"metadata":f})
            row.update({"url":u,"index_sha256":sha(b),"rows":len(rr),"nonexcluded_files":len(matches),"eligible_dates":sorted({m["day"] for m in matches}),"files":matches[:200]})
        except Exception as e:row["error"]=str(e)[:300]
        return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:res=list(pool.map(run,jobs))
    pp=[r for r in res if r["telescope"]=="Parkes"];gg=[r for r in res if r["telescope"]=="GBT"]
    # Count same-date possible ON/OFF groups, explicitly not a validated cadence.
    groups=[]
    for target in gbt:
        a=[r for r in gg if r["target"]==target and r["query"]==target];b=[r for r in gg if r["target"]==target and r["query"]==target+"_OFF"]
        if a and b:
            dd=set(a[0].get("eligible_dates",[]))&set(b[0].get("eligible_dates",[]))
            if dd:groups.append({"target":target,"days":sorted(dd)})
    result={"status":"METADATA_ONLY_REQUIRES_HEADER_AND_POINTING_QUALIFICATION","prior_excluded_targets":len(ts),"prior_excluded_days":len(days),"archive_index_sha256":sha(raw),"nonhip_aliases_total":len([q for q in all_names if not re.fullmatch(r"HIP\d+",q) and q not in ts]),"parkes_queries":len(pp),"parkes_nonexcluded_files":sum(r.get("nonexcluded_files",0) for r in pp),"gbt_queries":len(gg),"gbt_nonexcluded_files":sum(r.get("nonexcluded_files",0) for r in gg),"gbt_possible_on_off_dates":groups,"successful_queries":sum("error" not in r for r in res),"queries":len(res),"no_power_access":True,"exam_started":False,"production_adoption":False}
    save("queries.json",res);save("result.json",result);print("R37_SOURCE_SUMMARY",json.dumps(result),flush=True)
if __name__=="__main__":main()
