#!/usr/bin/env python3
"""R37F2D2: official natural archive metadata and bounded complete checksum audit."""
import hashlib,io,json,math,os,pathlib,re,urllib.parse,zipfile
import requests
from bs4 import BeautifulSoup
P=pathlib.Path("research/seti/r37f2d2_integrity_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f2d2");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def locator(u):
    p=urllib.parse.urlparse(u)
    if p.scheme not in ("http","https") or not (p.hostname=="bldata.berkeley.edu" or re.fullmatch(r"blpd\d+\.ssl\.berkeley\.edu",p.hostname or "")):return None
    if p.username or p.password or p.query or p.fragment:return None
    return (p.hostname,urllib.parse.unquote(p.path))
def parent(p):
    r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['artifact']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90)
    r.raise_for_status();b=r.content;assert len(b)<16*1024**2 and sha(b)==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        n=[v for v in z.infolist() if pathlib.PurePosixPath(v.filename).name=="result.json"];assert len(n)==1 and n[0].file_size<8*1024**2
        d=json.loads(z.read(n[0]))
    return [dict(url=h["url"],target=h["target"],mjd=h["mjd"],day=h["day"],source=s["source"]) for s in d["sources"] for h in s["header_probes"] if any(a.get("fine_grid") and a.get("covers_methanol") for a in h["attempts"])]
def get(u,hosts):
    for _ in range(5):
        p=urllib.parse.urlparse(u);assert p.scheme=="https" and p.hostname in hosts
        with requests.get(u,timeout=(15,45),stream=True,allow_redirects=False) as r:
            if r.status_code in (301,302,303,307,308):u=urllib.parse.urljoin(u,r.headers["Location"]);continue
            r.raise_for_status();assert r.status_code==200
            b=bytearray()
            for x in r.iter_content(65536):b.extend(x);assert len(b)<=S["limits"]["metadata_bytes"]
            return bytes(b),u
    raise ValueError("redirect limit")
def portal_page(url):
    b,final=get(url,S["allowed_portal_hosts"]);s=BeautifulSoup(b,"html.parser");rows=[]
    for tr in s.find_all("tr"):
        cells=[x.get_text(" ",strip=True) for x in tr.find_all("td")]
        if len(cells)<10:continue
        md5=[x for x in cells if re.fullmatch(r"(?i)[0-9a-f]{32}",x)]
        if len(md5)!=1:continue
        try:mjd=float(cells[1]);size=int(cells[-3].replace(",","").replace(" ",""))
        except (ValueError,IndexError):continue
        for a in tr.find_all("a",href=True):
            u=urllib.parse.urljoin(final,a["href"]);k=locator(u)
            if k:rows.append({"url":u,"target":cells[3],"mjd":mjd,"size":size,"md5":md5[0].lower(),"portal_page":final,"page_sha256":sha(b)})
    return s,rows,{"url":final,"sha256":sha(b),"rows_parsed":len(rows)}
def main():
    vals=[v for p in S["parents"] for v in parent(p)]
    wanted={locator(v["url"]):v for v in vals};assert None not in wanted and len(wanted)==len(vals)==11
    release,ru=get(S["official_release"],["seti.berkeley.edu"])
    soup=BeautifulSoup(release,"html.parser");links=[]
    for a in soup.find_all("a",href=True):
        url=urllib.parse.urljoin(ru,a["href"]);p=urllib.parse.urlparse(url)
        if p.hostname in S["allowed_portal_hosts"] and p.path=="/opendatasearch" and p.scheme=="https":links.append(url)
        elif p.hostname in S["allowed_portal_hosts"] and p.path=="/opendatasearch" and p.scheme=="http":
            links.append(urllib.parse.urlunparse(p._replace(scheme="https")))
    links=sorted(set(links));assert len(links)>=4
    hits={};pages=[];errors=[]
    for url in links:
        p=urllib.parse.urlparse(url);params=urllib.parse.parse_qs(p.query)
        params["perPage"]=["100"];params["file_type"]=["filterbank"]
        for i in range(1,S["limits"]["maximum_pages_per_query"]+1):
            params["page"]=[str(i)];q=urllib.parse.urlunparse(p._replace(query=urllib.parse.urlencode(params,doseq=True)))
            try:
                s,rows,receipt=portal_page(q);pages.append(receipt)
                for r in rows:
                    k=locator(r["url"])
                    if k not in wanted:continue
                    v=wanted[k]
                    if r["target"].strip().upper()!=v["target"] or abs(r["mjd"]-v["mjd"])*86400>3:continue
                    hits.setdefault(k,[]).append(r)
                if not rows or len(rows)<100:break
            except Exception as e:errors.append({"url":q,"error":type(e).__name__+": "+str(e)[:200]});break
    result=[];remaining=S["limits"]["full_file_bytes"]
    for k,v in wanted.items():
        matches=hits.get(k,[]);r={**v,"status":"HOLD_NO_EXACT_HTTPS_REFERENCE"}
        if matches and len({(x["size"],x["md5"]) for x in matches})==1:
            ref=matches[0];r.update(status="HTTPS_PUBLISHED_MD5_JOIN",reference=ref,file_integrity_verified=False)
            if ref["size"]>remaining:r["download_status"]="HOLD_FIXED_FULL_FILE_BUDGET"
            else:
                remaining-=ref["size"]
                try:
                    md5=hashlib.md5();h=hashlib.sha256();n=0
                    with requests.get(v["url"],stream=True,timeout=(15,60),allow_redirects=False,headers={"Accept-Encoding":"identity"}) as f:
                        f.raise_for_status();assert f.status_code==200
                        for chunk in f.iter_content(1024**2):
                            n+=len(chunk);assert n<=ref["size"];md5.update(chunk);h.update(chunk)
                    assert n==ref["size"] and md5.hexdigest()==ref["md5"]
                    r.update(download_status="FULL_SIZE_AND_PUBLISHED_MD5_MATCH",file_integrity_verified=True,local_sha256=h.hexdigest(),bytes_read=n)
                except Exception as e:r.update(download_status="HOLD",download_error=type(e).__name__+": "+str(e)[:200])
        elif matches:r["status"]="HOLD_CONFLICTING_HTTPS_REFERENCES"
        result.append(r)
    d={"status":"ARCHIVE_INTEGRITY_AUDIT_ONLY","release_sha256":sha(release),"portal_links":links,"pages":pages,"errors":errors,"files":result,"reference_matches":sum(v["status"]=="HTTPS_PUBLISHED_MD5_JOIN" for v in result),"full_file_matches":sum(bool(v.get("file_integrity_verified")) for v in result),"scientific_line_verified":False,"independent_spatial_pairs":0,"production_adoption":False}
    (O/"result.json").write_text(json.dumps(d,indent=2)+"\n")
    print("R37F2D2_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
