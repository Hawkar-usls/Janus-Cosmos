"""SIT2E Parkes archived three-date PHYSICAL HEADER probe; NO radio power reads.

Real filenames are discovered anew from the published 2025 archive listing.
July's official products are SIGPROC .fil, Sept/Oct use .h5. A bounded
HTTP RangeFile fetches physical HDF5 metadata only; provenance is
FIRST-SEEN header-byte evidence, NOT a whole-file SHA256 certificate.
"""
from __future__ import annotations
import hashlib
import io
import json
import pathlib
import re
import struct
import sys
from urllib.parse import urljoin,unquote

import h5py
import hdf5plugin
import requests

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"research"/"seti"))
from r4_chunk_replay import RangeFile
import janus_sit2e_parkes_preflight as q

OUT=ROOT/"runtime"/"sit2e_parkes"
OUT.mkdir(parents=True,exist_ok=True)
INT={"nchans","nifs","nbits","ibeam","nbeams","telescope_id","machine_id","data_type",
     "barycentric","pulsarcentric","signed","nbins","nsamples"}
DOUBLE={"tsamp","tstart","fch1","foff","src_raj","src_dej","az_start","za_start","refdm","period","fchannel"}
STRING={"source_name","rawdatafile","project_id","observer","hostname","date_obs"}
EXPECTED=("fch1","foff","tstart","tsamp")

def parse_sigproc_header(raw:bytes):
    """Header only; disallow unknown offsets instead of inventing layout."""
    buf=io.BytesIO(raw)
    def str1():
        raw=buf.read(4)
        if len(raw)!=4:raise ValueError("SIGPROC_TRUNCATED_LENGTH")
        k=struct.unpack("<i",raw)[0]
        if not 1<=k<=512:raise ValueError("SIGPROC_INVALID_STRING_LENGTH")
        v=buf.read(k)
        if len(v)!=k:raise ValueError("SIGPROC_TRUNCATED_STRING")
        return v.decode("utf-8",errors="replace")
    if str1()!="HEADER_START":raise ValueError("NO_SIGPROC_HEADER_START")
    h={}
    while True:
        n=str1()
        if n=="HEADER_END":break
        if n in INT:
            x=buf.read(4)
            if len(x)!=4:raise ValueError("TRUNCATED_INT")
            h[n]=struct.unpack("<i",x)[0]
        elif n in DOUBLE:
            x=buf.read(8)
            if len(x)!=8:raise ValueError("TRUNCATED_DOUBLE")
            h[n]=struct.unpack("<d",x)[0]
        elif n in STRING:h[n]=str1()
        else:raise ValueError("UNKNOWN_SIGPROC_PARAMETER:"+n)
        if buf.tell()>8192:raise ValueError("SIGPROC_HEADER_OVER_8KB")
    absent=[n for n in EXPECTED if n not in h]
    if absent:raise ValueError("SIGPROC_MISSING_AXES:"+",".join(absent))
    if h.get("nchans",0)<128 or h["tsamp"]<=0 or h["foff"]==0:
        raise ValueError("INVALID_SIGPROC_SPECTRAL_GRID")
    h["frequency_range_mhz"]=sorted((h["fch1"],h["fch1"]+(h["nchans"]-1)*h["foff"]))
    h["channel_width_hz"]=abs(h["foff"])*1e6
    h["header_bytes"]=buf.tell()
    h["no_radio_power_decoded"]=True
    return h

def attr_scalar(v):
    if isinstance(v,bytes):return v.decode("utf-8",errors="replace")
    if hasattr(v,"tolist"):v=v.tolist()
    if isinstance(v,list):return [attr_scalar(x) for x in v[:5]]
    if isinstance(v,(int,float,str,bool)) or v is None:return v
    return str(v)

def hdf5_remote_header(url):
    source=None
    try:
        source=RangeFile(url,budget=12*1024*1024)
        with h5py.File(source,"r") as h:
            if "data" not in h:raise ValueError("HDF5_NO_DATA")
            d=h["data"]
            attrs={str(k).lower():attr_scalar(v) for k,v in h.attrs.items()}
            attrs.update({str(k).lower():attr_scalar(v) for k,v in d.attrs.items()})
            missing=[k for k in EXPECTED if k not in attrs]
            if missing:raise ValueError("HDF5_MISSING_PHYSICAL_AXES:"+",".join(missing))
            n=int(d.shape[-1]);f0=float(attrs["fch1"]);df=float(attrs["foff"])
            lo=min(f0,f0+(n-1)*df);hi=max(f0,f0+(n-1)*df)
            return {"source_name":attrs.get("source_name"),"src_raj":attrs.get("src_raj"),
                    "src_dej":attrs.get("src_dej"),"tstart":float(attrs["tstart"]),
                    "tsamp":float(attrs["tsamp"]),"fch1":f0,"foff":df,"nchans":n,
                    "nifs":attrs.get("nifs"),"data_shape":list(d.shape),
                    "channel_width_hz":abs(df)*1e6,
                    "frequency_range_mhz":[lo,hi],
                    "header_mode":"HTTP_RANGE_HDF5_METADATA_ONLY__NOT_FULL_SHA256",
                    "read_radio_power":False}
    finally:
        if source is not None:source.close()

def bytes_range_header(url):
    if not q.safe_url(url,q.BASE):raise ValueError("NON_OFFICIAL_HOST")
    r=requests.get(url,headers={"Range":"bytes=0-8191","Accept-Encoding":"identity"},
                   timeout=(15,45))
    if r.status_code!=206:raise ValueError("NO_EXACT_HTTP_RANGE")
    if not r.headers.get("Content-Range","").startswith("bytes 0-8191/"):
        raise ValueError("BAD_RANGE_CONTENT_RANGE")
    if len(r.content)!=8192:raise ValueError("SHORT_HTTP_RANGE")
    header=parse_sigproc_header(r.content)
    return {"physical":header,"first_8192_byte_range_sha256":hashlib.sha256(r.content).hexdigest(),
            "total_source_file_sha256":None,"full_file_downloaded":False}

def files_from_real_index(root):
    with requests.get(root,timeout=(12,35)) as r:
        r.raise_for_status()
        if len(r.content)>7_000_000:raise ValueError("INDEX_TOO_LARGE")
        parsed=q.LinkParser();parsed.feed(r.text)
        out=[]
        for h in parsed.hrefs:
            h=unquote(h).split("?",1)[0]
            if not h or h.endswith("/") or h.startswith("../") or h.startswith("/"):continue
            url=urljoin(root,h)
            if q.safe_url(url,root) and (h.endswith(".fil") or h.endswith(".h5")):
                out.append({"filename":h,"url":url,"role":q.role_hint(h),
                            "kind":"SIGPROC_FIL" if h.endswith(".fil") else "HDF5"})
        return out

def probe_date(spec):
    root=urljoin(q.BASE,spec["archive_path"])
    # Fixed backend sample: same backend blc00 in every date, then
    # frequency grid must be authenticated rather than assumed.
    backend=urljoin(root,"blc00/")
    entries=files_from_real_index(backend)
    mode="SIGPROC_FIL" if spec["date"]=="2025-07-31" else "HDF5"
    candidates=[]
    for role in ("SOURCE_NAME_ONLY","REFERENCE_NAME_ONLY"):
        rr=[z for z in entries if z["role"]==role and z["kind"]==mode
            and re.search(r"\.0002\.(?:fil|h5)$",z["filename"])]
        # One actual file per unique scan time, no duplicates from different
        # backend banks, no inference based on source-name alone.
        candidates.extend(rr[:10])
    report={"date":spec["date"],"archive_root":root,"actual_probe_backend":backend,
            "files_available_in_blc00":len(entries),
            "source":None,"reference":None,
            "physical_headers_read":0,"full_source_sha256_certified":False,
            "antenna_pointing_independently_certified":False,
            "no_HDF5_radio_power_read":True,
            "role_labels_from_archive_filename_only":True,
            "range_probe_errors":[]}
    with requests.Session() as sess:
        sizeprobe=[]
        for x in candidates:
            try:
                response=sess.head(x["url"],timeout=(10,25),allow_redirects=True)
                response.raise_for_status()
                if not q.safe_url(response.url,root):raise ValueError("UNEXPECTED_HOST")
                size=int(response.headers.get("Content-Length","-1"))
                sizeprobe.append({**x,"bytes":size})
            except Exception as e:report["range_probe_errors"].append({"name":x["filename"],"stage":"HEAD","error":str(e)[:150]})
        report["sampled_head_sizes"]=sizeprobe
        # 0-length and 7624-byte placeholders were observed in Sep;
        # don't count them as genuine spectra. Prefer valid files >1MB.
        for role,key in (("SOURCE_NAME_ONLY","source"),("REFERENCE_NAME_ONLY","reference")):
            possible=[x for x in sizeprobe if x["role"]==role and x["bytes"]>1_000_000]
            if not possible:continue
            candidate=min(possible,key=lambda x:x["bytes"])
            try:
                if mode=="SIGPROC_FIL":
                    meta=bytes_range_header(candidate["url"])
                    physical=meta["physical"]
                    del meta["physical"]
                else:
                    physical=hdf5_remote_header(candidate["url"])
                    meta={"total_source_file_sha256":None,"full_file_downloaded":False}
                report[key]={**candidate,"physical":physical,**meta,
                             "role_measured_telemetry":False}
                report["physical_headers_read"]+=1
            except Exception as e:
                report["range_probe_errors"].append({"name":candidate["filename"],
                                                     "stage":"HEADER","error_type":type(e).__name__,
                                                     "error":str(e)[:260]})
    if report["source"] and report["reference"]:
        a=report["source"]["physical"]["frequency_range_mhz"]
        b=report["reference"]["physical"]["frequency_range_mhz"]
        overlap=max(a[0],b[0])<min(a[1],b[1])
        report["ON_OFF_same_backend_physical_frequency_overlap"]=overlap
        report["status"]="PAIR_HEADER_RANGE_PASS__FULL_BYTES_AND_INDEPENDENT_POINTING_HOLD" if overlap else "HOLD_SAME_BACKEND_PHYSICAL_FREQUENCY_OVERLAP"
    else:
        report["status"]="HOLD_NO_VALID_TWO_ROLE_PHYSICAL_HEADER_PAIR"
    return report

def run():
    result={"schema":"janus.cosmos.sit2e.parkes_mixed_format_range_header_pilot.v1",
            "status":"THREE_DATE_METADATA_PREFLIGHT_ONLY",
            "source":"https://bldata.berkeley.edu/ATLAS/",
            "all_radio_power_unread":True,"october_algorithmic_scoring_unopened":True,
            "thresholds_fitted":False,"whole_source_SHA256_computed":False,
            "actual_antenna_pointing_authenticated":False,"science_gain_claim":None,
            "dates":[]}
    for spec in q.P["archive_sessions"]:
        try:r=probe_date(spec)
        except Exception as e:r={"date":spec["date"],"status":"HOLD_INDEX_OR_FETCH_ERROR","error":str(e)[:300]}
        result["dates"].append(r)
        print("SIT2E_RANGE_DATE",json.dumps(r,sort_keys=True,ensure_ascii=False),flush=True)
    ranges=[x for x in result["dates"] if x.get("source") and x.get("reference")]
    if len(ranges)==3:
        f=[x["source"]["physical"]["frequency_range_mhz"] for x in ranges]
        all_overlap=max(y[0] for y in f)<min(y[1] for y in f)
        result["three_date_source_band_overlap"]=all_overlap
    else:
        result["three_date_source_band_overlap"]=False
    result["status"]=("THREE_DATE_HEADER_PARTIAL_PASS__REMAINING_SHA_AND_POINTING_HOLD"
        if all(x.get("status","").startswith("PAIR_HEADER") for x in result["dates"])
        else "THREE_DATE_SOURCE_HEADER_QUALIFICATION_HOLD")
    return result

def main():
    result=run()
    (OUT/"sit2e_parkes_range_headers.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("SIT2E_RANGE_FINAL",json.dumps({"status":result["status"],
        "three_date_source_band_overlap":result["three_date_source_band_overlap"],
        "october_algorithmic_scoring_unopened":result["october_algorithmic_scoring_unopened"]
    },sort_keys=True))

if __name__=="__main__":main()
