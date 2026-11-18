#!/usr/bin/env python3
"""R33B4 audit: source identity and bounded archive discovery, no spectral read."""
import hashlib,json,pathlib,urllib.request,urllib.parse,concurrent.futures,re
P=pathlib.Path("research/seti/r33b4_source_identity_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r33b4");OUT.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    raw=pathlib.Path(S["parent"]["path"]).read_bytes();assert sha(raw)==S["parent"]["sha256"]
    d=json.loads(raw);assert len(d["entries"])==2
    audit=[]
    for e in d["entries"]:
        t=e["on_track"];p=e["passport"]["provenance"];nom=e["published"]["rounded_mhz"]*1e6
        residual=t["scan_start_frequency_hz"]-nom;resolution=abs(p["df_hz"])/((t["rows"]-1)*p["ts"])
        audit.append({"id":e["published"]["id"],"frequency_residual_hz":residual,"within_half_khz_at_scan_start":abs(residual)<=500,"scan_start_hz":t["scan_start_frequency_hz"],"measured_drift_hz_s":t["drift_hz_s"],"endpoint_resolution_hz_s":resolution,"off_residuals_hz":[o["frequency_residual_hz"] for o in e["off_matches"]],"status":"HOLD_ORIGINAL_HIT_AND_MEMBERSHIP_MAPPING_UNAVAILABLE","verified_promoted":0})
    jobs=[]
    for group,targets,types in [("catalog",S["catalog_scout"]["targets"],S["catalog_scout"]["file_types"]),("figure6",S["figure6"]["index_targets"],S["figure6"]["types"]),("r34",S["r34_scout"]["targets"],S["r34_scout"]["types"])]:
        for target in targets:
            for typ in types:jobs.append((group,target,typ))
    def query(item):
        group,target,typ=item;url="https://seti.berkeley.edu/opendata/api/query-files?"+urllib.parse.urlencode({"target":target,"telescope":"Parkes","file-type":typ})
        row={"group":group,"target":target,"type":typ,"url":url}
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"JANUS COSMOS R33B4 academic source metadata"})
            with urllib.request.urlopen(req,timeout=40) as resp:
                assert urllib.parse.urlparse(resp.url).hostname=="seti.berkeley.edu"
                raw=resp.read(32*1024*1024+1);assert len(raw)<=32*1024*1024
            d=json.loads(raw);rr=d if isinstance(d,list) else next(d[k] for k in ("data","files","results") if isinstance(d.get(k),list))
            path=OUT/("index-"+sha(url.encode())[:16]+".json");path.write_bytes(raw)
            matches=[]
            for v in rr:
                if not isinstance(v,dict):continue
                u=v.get("url","")
                if not isinstance(u,str):continue
                parsed=urllib.parse.urlparse(u)
                if not (parsed.hostname=="bldata.berkeley.edu" or re.fullmatch(r"blpd\d+\.ssl\.berkeley\.edu",parsed.hostname or "")):continue
                name=urllib.parse.unquote(parsed.path.rsplit("/",1)[-1])
                try:day=int(float(v.get("mjd",0)))
                except (TypeError,ValueError):day=0
                if group=="catalog":
                    if day!=58602 and "2019-04-29" not in name:continue
                    if not any(s in name for s in S["catalog_scout"]["exact_scans"]):continue
                elif group=="figure6":
                    if day!=58598 and "2019-04-25" not in name:continue
                matches.append({"name":name,"url":u,"mjd":day,"file_type":v.get("file_type"),"size":v.get("size")})
            row.update(count=len(rr),matching_count=len(matches),index_sha256=sha(raw),matches=matches[:300],catalog_like_count=sum(m["name"].lower().endswith((".dat",".csv",".json")) for m in matches))
        except Exception as e:row["error"]=str(e)[:300]
        return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(query,jobs))
    # Exact filename only; do not infer the title/filename discrepancy away.
    f6=[m for r in results if r["group"]=="figure6" for m in r.get("matches",[]) if "blc03_HIP98959_S_2019-04-25T19:24:30.000" in m["name"]]
    f6=list({m["url"]:m for m in f6}.values())
    summary={"status":"HOLD_ORIGINAL_HIT_MAPPING","audits":audit,"query_count":len(results),"queries_ok":sum("error" not in r for r in results),"catalog_like_count":sum(r.get("catalog_like_count",0) for r in results if r["group"]=="catalog"),"figure6_exact_file_matches":f6,"figure6_caption_filename_discrepancy":True,"r34_metadata_matches":sum(r.get("matching_count",0) for r in results if r["group"]=="r34"),"queries":[{k:v for k,v in r.items() if k!="matches"} for r in results],"source_labels_promoted":0,"no_spectral_read":True,"exam_started":False,"production_adoption":False}
    (OUT/"queries.json").write_text(json.dumps(results,indent=2)+"\n")
    (OUT/"result.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("R33B4_SUMMARY",json.dumps(summary),flush=True)
if __name__=="__main__":main()
