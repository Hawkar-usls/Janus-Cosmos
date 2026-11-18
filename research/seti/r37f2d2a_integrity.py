#!/usr/bin/env python3
"""Natural archive integrity only: match exact HTTPS reference and stream one complete file."""
import hashlib,io,json,os,pathlib,urllib.parse,zipfile,time
import requests
import r37f2d2_integrity as base
P=pathlib.Path("research/seti/r37f2d2a_integrity_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f2d2a");O.mkdir(parents=True,exist_ok=True)
def save(name,d):
    (O/name).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def load():
    p=S["parent"];r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['artifact']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90)
    r.raise_for_status();b=r.content;assert len(b)<16*1024**2 and hashlib.sha256(b).hexdigest()==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        n=[v for v in z.infolist() if pathlib.PurePosixPath(v.filename).name=="result.json"];assert len(n)==1 and n[0].file_size<8*1024**2
        d=json.loads(z.read(n[0]))
    assert len(d["files"])==11
    return d["files"]
def main():
    rows=load();pages=[];errors=[]
    wanted={base.locator(v["url"]):v for v in rows};assert len(wanted)==11 and None not in wanted
    refs={}
    for k,v in wanted.items():
        if v["status"]=="HTTPS_PUBLISHED_MD5_JOIN":refs[k]=[v["reference"]]
    epochs=sorted({v["mjd"] for v in rows if v["status"]!="HTTPS_PUBLISHED_MD5_JOIN"})
    for epoch in epochs:
        for page in range(1,S["max_pages_per_epoch"]+1):
            q=urllib.parse.urlencode({"project":"GBT","file_type":"filterbank","mjd":f"{epoch:.8f}","mjd_range":S["narrow_mjd_range_days"],"search":"Search","perPage":100,"page":page})
            try:
                soup,rr,receipt=base.portal_page(S["portal"]+"?"+q);pages.append(receipt)
                for r in rr:
                    k=base.locator(r["url"])
                    if k in wanted and r["target"].strip().upper()==wanted[k]["target"] and abs(r["mjd"]-wanted[k]["mjd"])*86400<=S["match_mjd_tolerance_s"]:
                        refs.setdefault(k,[]).append(r)
                if len(rr)<100:break
            except Exception as e:errors.append({"epoch":epoch,"page":page,"error":type(e).__name__+": "+str(e)[:250]});break
    result=[]
    for k,v in wanted.items():
        r={k2:v[k2] for k2 in ("url","target","mjd","day","source")}
        matches=refs.get(k,[])
        if matches and len({(x["size"],x["md5"]) for x in matches})==1:
            r.update(status="HTTPS_PUBLISHED_MD5_JOIN",reference=matches[0])
        elif matches:r["status"]="HOLD_CONFLICTING_REFERENCE"
        else:r["status"]="HOLD_NO_EXACT_HTTPS_REFERENCE"
        result.append(r)
    eligible=sorted([r for r in result if r["status"]=="HTTPS_PUBLISHED_MD5_JOIN" and 0<r["reference"]["size"]<=S["full_stream_bytes"]],key=lambda r:(r["reference"]["size"],r["url"]))
    selected=eligible[0] if eligible else None
    save("frozen_selection.json",{"files":result,"selected":selected,"selection_before_full_read":True,"pages":pages,"errors":errors})
    receipt={"status":"HOLD_NO_ELIGIBLE_FILE"}
    if selected:
        ref=selected["reference"];n=0;md5=hashlib.md5();sha=hashlib.sha256();start=time.monotonic()
        receipt={"status":"HOLD","source":selected["source"],"target":selected["target"],"expected_size":ref["size"],"expected_md5":ref["md5"],"url":selected["url"]}
        try:
            assert base.locator(selected["url"])==base.locator(ref["url"])
            with requests.get(selected["url"],stream=True,timeout=(15,90),allow_redirects=False,headers={"Accept-Encoding":"identity"}) as f:
                f.raise_for_status();assert f.status_code==200
                if f.headers.get("Content-Length","").isdigit():assert int(f.headers["Content-Length"])==ref["size"]
                assert f.headers.get("Content-Encoding","identity")=="identity"
                for chunk in f.iter_content(4*1024**2):
                    n+=len(chunk);assert n<=ref["size"] and n<=S["full_stream_bytes"]
                    md5.update(chunk);sha.update(chunk)
                    if n%(256*1024**2)==0:print("R37F2D2A_BYTES",n,flush=True)
            assert n==ref["size"],"size mismatch"
            receipt.update(bytes_read=n,computed_md5=md5.hexdigest(),local_sha256=sha.hexdigest(),seconds=time.monotonic()-start)
            receipt["status"]="FULL_SIZE_AND_PUBLISHED_MD5_MATCH" if md5.hexdigest()==ref["md5"] else "HOLD_MD5_MISMATCH"
        except Exception as e:receipt.update(error=type(e).__name__+": "+str(e)[:300],bytes_read=n,seconds=time.monotonic()-start)
    d={"status":"FULL_FILE_INTEGRITY_MATCH" if receipt["status"]=="FULL_SIZE_AND_PUBLISHED_MD5_MATCH" else "ARCHIVE_INTEGRITY_HOLD","reference_matches":sum(r["status"]=="HTTPS_PUBLISHED_MD5_JOIN" for r in result),"full_file_matches":int(receipt["status"]=="FULL_SIZE_AND_PUBLISHED_MD5_MATCH"),"files":result,"stream":receipt,"pages":pages,"errors":errors,"natural_line_verified":False,"independent_spatial_pairs":0,"production_adoption":False}
    save("result.json",d);print("R37F2D2A_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
