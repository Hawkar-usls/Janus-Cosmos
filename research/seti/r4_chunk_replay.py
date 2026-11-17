#!/usr/bin/env python3
"""Budgeted, chunk-aware replay of frozen Breakthrough Listen event groups.
No autonomous promotions or claims. Archive observations only, not new telescope time.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import pathlib
import re
import tempfile
from urllib.parse import urljoin, urlparse

import bs4
import h5py
import hdf5plugin  # registers the archive's bitshuffle decoder
import numpy as np
import requests

MIB = 1024 * 1024
ARCHIVE = "https://breakthroughinitiatives.org/opendatasearch"
UA = "JANUS COSMOS scientific archive replay/4.0"
NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")
SCAN = re.compile(r"_guppi_(\d{5})_(\d{5})_(?:\d+_)?(.+?)_(\d{4})\.gpuspec\.0000\.h5$", re.I)

def numbers(s):
    return [float(v) for v in NUMBER.findall(str(s))]

def fetch(url, params=None, limit=8*MIB):
    with requests.get(url, params=params, headers={"User-Agent":UA},
                      timeout=(15,40), stream=True) as r:
        r.raise_for_status()
        data = bytearray()
        for b in r.iter_content(65536):
            data.extend(b)
            if len(data)>limit:
                raise ValueError("metadata resource exceeds budget")
        return bytes(data), r.url

def basename(url):
    return pathlib.PurePosixPath(urlparse(url).path).name

def archive_url(url):
    return urlparse(url).scheme=="https" and urlparse(url).hostname=="bldata.berkeley.edu"

def scan(url):
    name=basename(url)
    m=SCAN.search(name)
    if not m:
        return None
    return dict(day=int(m[1]), seconds=int(m[2]), target=m[3].upper(), scan=int(m[4]),
                prefix=name.split("_guppi_")[0], session=urlparse(url).path.rsplit("/",1)[0])

def norm(name):
    return name.upper().removesuffix("_S").removesuffix("_OFF")

def query_archive(params):
    raw, url=fetch(ARCHIVE, params=dict({"file_type":"HDF5","perPage":100,
                                           "search":"Search","order":"asc","sort":"utc"},**params))
    soup=bs4.BeautifulSoup(raw,"html.parser")
    links=sorted({urljoin(url,a["href"]) for a in soup.select("a[href]")
                  if archive_url(urljoin(url,a["href"])) and basename(urljoin(url,a["href"])).endswith(".0000.h5")})
    return links,{"query":url,"sha256":hashlib.sha256(raw).hexdigest(),"matched_links":len(links)}

class RangeFile(io.RawIOBase):
    """Exact, bounded GET ranges with a disk-backed per-file exact-request cache."""
    def __init__(self,url,budget):
        self.url=url; self.budget=budget; self.pos=0; self.total=0
        self.requests_count=0; self.cache={}; self.cache_hits=0
        self.cache_file=tempfile.TemporaryFile()
        self.session=requests.Session()
        self.session.headers.update({"User-Agent":UA,"Accept-Encoding":"identity"})
        try:
            with self.session.get(url,headers={"Range":"bytes=0-0"},stream=True,timeout=(15,45)) as r:
                self.requests_count+=1
                if r.status_code!=206 or not r.headers.get("Content-Range","").startswith("bytes 0-0/"):
                    raise IOError("archive refused exact range")
                self.size=int(r.headers["Content-Range"].split("/")[-1])
        except Exception:
            self.close()
            raise
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self,offset,whence=0):
        pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset
        if pos<0: raise ValueError("negative seek")
        self.pos=pos
        return pos
    def readinto(self,b):
        n=min(len(b),self.size-self.pos)
        if n<=0: return 0
        start=self.pos; key=(start,n)
        if key in self.cache:
            self.cache_file.seek(self.cache[key])
            memoryview(b)[:n]=self.cache_file.read(n)
            self.cache_hits+=1; self.pos+=n
            return n
        if self.total+n>self.budget:
            raise IOError(f"range budget exceeded: fetched={self.total}, requested={n}, budget={self.budget}")
        dest=memoryview(b).cast("B")
        # One contiguous HTTP request per HDF5 range, streamed in small blocks.
        with self.session.get(self.url,headers={"Range":f"bytes={start}-{start+n-1}"},stream=True,timeout=(15,60)) as r:
            self.requests_count+=1
            if r.status_code!=206 or r.headers.get("Content-Range")!=f"bytes {start}-{start+n-1}/{self.size}":
                raise IOError("range integrity check failed")
            offset=0; self.cache_file.seek(0,2); at=self.cache_file.tell()
            for chunk in r.iter_content(1024*1024):
                end=offset+len(chunk)
                if end>n: raise IOError("overlong range")
                dest[offset:end]=chunk
                self.cache_file.write(chunk); offset=end
            if offset!=n: raise IOError("short range")
        self.cache[key]=at
        self.pos+=n; self.total+=n
        return n
    def close(self):
        if not self.closed:
            self.session.close()
            self.cache_file.close()
        super().close()

def window(d,f0,df,f,half_hz):
    center=int(round((f-f0)/df))
    rad=max(3,int(math.ceil(half_hz/(abs(df)*1e6))))
    lo=max(0,center-rad); hi=min(d.shape[-1],center+rad+1)
    if not lo<=center<hi: return None
    return (lo,hi,center)

def preflight(d,windows,nt):
    if d.chunks is None: return None
    ct,ci,cf=d.chunks
    coords={(t,0,k*cf) for t in range(0,nt,ct) for lo,hi,_ in windows
            for k in range(lo//cf,(hi-1)//cf+1)}
    size=0
    for coord in coords:
        size+=int(d.id.get_chunk_info_by_coord(coord).size)
    return {"unique_chunks":len(coords),"compressed_bytes":size}

def inspect(url, lines, role, out, budget):
    assert archive_url(url)
    result={"url":url,"role":role,"status":"HOLD_UNREAD","budget_bytes":budget}
    remote=None
    try:
        remote=RangeFile(url,budget)
        with h5py.File(remote,"r",rdcc_nbytes=128*MIB,rdcc_nslots=10007) as h:
            d=h["data"]; a=d.attrs
            f0=float(a["fch1"]); df=float(a["foff"]); ts=float(a.get("tsamp",0))
            result.update(shape=list(d.shape),chunks=list(d.chunks) if d.chunks else None,
                          file_bytes=remote.size, channel_hz=abs(df)*1e6,tsamp_s=ts,
                          tstart_mjd=float(a.get("tstart",0)),source_name=str(a.get("source_name","")),
                          frequency_range_mhz=sorted([f0,f0+(d.shape[-1]-1)*df]))
            print("SPECTRUM_HEADER",json.dumps(result),flush=True)
            if d.ndim!=3 or df==0 or abs(df)*1e6>10 or ts<=0:
                raise ValueError("not supported high resolution time spectrum")
            nw=[(line,window(d,f0,df,line["frequency_mhz"],600)) for line in lines]
            nw=[(line,w) for line,w in nw if w is not None]
            if not nw: raise ValueError("no target frequency in receiver coverage")
            # Plan unique compressed chunks before reading. A partial scan is explicitly labeled.
            attempts=sorted({int(d.shape[0]),min(32,d.shape[0]),min(16,d.shape[0]),
                             min(8,d.shape[0]),min(4,d.shape[0]),1},reverse=True)
            plans=[]
            for nt in attempts:
                p=preflight(d,[w for _,w in nw],nt)
                plans.append({"rows":nt,"plan":p})
                if p is None or p["compressed_bytes"]+remote.total+8*MIB<=budget:
                    break
            else:
                raise IOError("even one time-row tile exceeds budget")
            result["chunk_plan"]=plans; result["sampled_time_rows"]=nt
            result["total_time_rows"]=int(d.shape[0]); result["full_time_scan"]=(nt==d.shape[0])
            print("CHUNK_PLAN",json.dumps({k:result[k] for k in
                ("role","sampled_time_rows","total_time_rows","chunk_plan")}),flush=True)
            products=[]
            for line,(lo,hi,center) in nw:
                try:
                    slab=np.asarray(d[:nt,0,lo:hi],dtype=np.float32)
                    vals=slab[np.isfinite(slab)]
                    if vals.size==0: raise ValueError("nonfinite data")
                    med=float(np.median(vals)); mad=float(np.median(np.abs(vals-med)))
                    peak=float(np.max(vals))
                    # Descriptive strength only; not a calibrated SNR or discovery statistic.
                    z=(peak-med)/(1.4826*mad) if mad>0 else None
                    name=hashlib.sha256((url+"|"+str(line["frequency_mhz"])).encode()).hexdigest()[:24]
                    np.save(out/f"{name}.npy",slab,allow_pickle=False)
                    products.append({**line,"shape":list(slab.shape),"center_channel":center,
                                     "channel0_mhz":f0+lo*df,"median":med,"mad":mad,
                                     "peak":peak,"peak_robust_z_descriptive":z,
                                     "slice_sha256":hashlib.sha256(slab.tobytes()).hexdigest(),
                                     "artifact":f"{name}.npy"})
                except Exception as e:
                    products.append({**line,"error":str(e)})
            result["products"]=products
            result["status"]="SLICES_READ" if any("artifact" in p for p in products) else "HOLD_NO_SLICES"
    except Exception as e:
        result["error"]=str(e)
    finally:
        if remote is not None:
            result.update(bytes_read=remote.total,range_requests=remote.requests_count,
                          cache_hits=remote.cache_hits)
            remote.close()
    print("OBSERVATION_RESULT",json.dumps(result,ensure_ascii=False),flush=True)
    return result

def run(target, out):
    spec=json.loads(pathlib.Path("research/seti/frozen_r4_targets.json").read_text())
    obj=spec["targets"][target]; catalog=spec["catalogs"][obj["band"]]
    raw,_=fetch(catalog["url"])
    if hashlib.sha256(raw).hexdigest()!=catalog["sha256"]: raise ValueError("catalog hash changed")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    lines=[]
    for mid in obj["frequencies"]:
        r=[r for r in rows if r["Source"]==target and abs(float(r["FreqMid"])-mid)<1e-7]
        if len(r)!=1 or r[0]["FileID"]!=obj["filename"]: raise ValueError("frozen identity mismatch")
        r=r[0]
        fs=numbers(r["Freqs"]); sn=numbers(r["SNR"]); ds=numbers(r["DriftRates"])
        if not fs or not all(abs(f-mid)<2 for f in fs):
            raise ValueError("unexpected event frequencies")
        # A fixed representative line per group; no result-based post-selection.
        k=max(range(len(fs)),key=lambda i:sn[i]) if len(sn)==len(fs) else len(fs)//2
        lines.append({"group_mhz":mid,"frequency_mhz":fs[k],"catalog_snr":sn[k] if len(sn)==len(fs) else None,
                      "catalog_drift_hz_s":ds[k] if len(ds)==len(fs) else None,
                      "selection":"max_catalog_snr" if len(sn)==len(fs) else "middle_catalog_index",
                      "group_event_count":int(float(r["Nevent"]))})
    out.mkdir(parents=True,exist_ok=True)
    receipt={"schema":"janus.seti.r4.chunk_replay.v1","target":target,"groups":lines,
             "catalog_sha256":catalog["sha256"],"archive_queries":[],"observations":[],
             "status":"HOLD","independent_confirmation":False,"technosignature_confirmed":False,
             "limitations":["A fixed representative event line per group, not a full blind search.",
                            "A peak in a slice is not a detection; no calibrated false alarm p-value.",
                            "OFF/other target controls and repeats must be interpreted with frequency drift and receiver sensitivity.",
                            "No new observations or claim of extraterrestrial origin."]}
    def save():
        (out/"result.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
    def lookup(params):
        try:
            links,q=query_archive(params); receipt["archive_queries"].append(q); return links
        except Exception as e:
            receipt["archive_queries"].append({"params":params,"error":str(e)}); return []
    print("FROZEN_LINES",json.dumps(lines),flush=True)
    url=obj["url"]; s0=scan(obj["filename"])
    if not s0: raise ValueError("unrecognized frozen file name")
    project="GBT" if obj["band"].startswith("GBT") else "all"
    epoch=s0["day"]+s0["seconds"]/86400
    target_links=lookup({"project":project,"target":target})
    # The exact time query also permits recovering the Parkes record and nearby scans.
    session_links=lookup({"project":project,"mjd":f"{epoch:.7f}","mjd_range":"0.025"})
    all_links=sorted(set(target_links+session_links))
    if url is None:
        matched=[u for u in all_links if basename(u)==obj["filename"]]
        url=matched[0] if len(matched)==1 else None
    if not url:
        receipt["status"]="HOLD_NO_EXACT_ARCHIVE_URL"; save(); return
    if basename(url)!=obj["filename"]: raise ValueError("mismatched archived source")
    s0=scan(url)
    receipt["source_file"]=url; save()
    original=inspect(url,lines,"ORIGINAL_ON",out,384*MIB)
    receipt["observations"].append(original); save()
    if original["status"]!="SLICES_READ":
        receipt["status"]="HOLD_ORIGINAL_READ_FAILED"; save(); return
    candidates=[]
    for u in all_links:
        v=scan(u)
        if u==url or not v or v["prefix"]!=s0["prefix"]: continue
        delta=(v["day"]-s0["day"])*86400+v["seconds"]-s0["seconds"]
        same_session=v["session"]==s0["session"] and abs(delta)<=1800
        if same_session and v["target"]!=s0["target"] and norm(v["target"])==norm(s0["target"]) and "_OFF" in v["target"]:
            role="ARCHIVE_LABELED_OFF"
        elif same_session and norm(v["target"])!=norm(s0["target"]):
            role="NEARBY_OTHER_TARGET_CONTROL"
        elif same_session and v["target"]==s0["target"] and v["scan"]!=s0["scan"]:
            role="SAME_SESSION_ON_REPEAT"
        elif v["target"]==s0["target"] and v["day"]!=s0["day"]:
            role="DIFFERENT_DAY_SAME_TARGET_REPEAT"
        else: continue
        candidates.append({"role":role,"url":u,"delta_seconds":delta,"scan":v["scan"],
                           "target":v["target"],"day":v["day"]})
    receipt["control_candidates"]=candidates
    # Deterministic nearest control, then nearest different-day observation, no signal-based choice.
    controls=[c for c in candidates if c["role"] in ("ARCHIVE_LABELED_OFF","NEARBY_OTHER_TARGET_CONTROL")]
    controls.sort(key=lambda c:(c["role"]!="ARCHIVE_LABELED_OFF",abs(c["delta_seconds"]),c["url"]))
    repeats=[c for c in candidates if c["role"] in ("DIFFERENT_DAY_SAME_TARGET_REPEAT","SAME_SESSION_ON_REPEAT")]
    repeats.sort(key=lambda c:(c["role"]!="DIFFERENT_DAY_SAME_TARGET_REPEAT",abs(c["delta_seconds"]),c["url"]))
    for c in controls[:1]+repeats[:1]:
        receipt["observations"].append(inspect(c["url"],lines,c["role"],out,192*MIB)); save()
    receipt["status"]="ARCHIVE_SLICES_AVAILABLE_REVIEW_REQUIRED"
    save()

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--target",required=True)
    parser.add_argument("--out",required=True,type=pathlib.Path)
    args=parser.parse_args()
    run(args.target,args.out)
