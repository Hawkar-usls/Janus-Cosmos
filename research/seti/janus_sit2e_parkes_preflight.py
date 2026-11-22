"""JANUS SIT2E — exact public Parkes/Murriyang 3-date HDF5 qualification.

Research-only. This qualifies original filenames and physical headers; reads
NO HDF5 radio-power elements, certainly no heldout October test statistic.
Separate temporal dates are necessary but not sufficient for independent
physically-qualified pointing / technosignature evidence.
"""
from __future__ import annotations
import hashlib
import html.parser
import json
import pathlib
import re
import sys
import time
from urllib.parse import urljoin, urlparse, unquote
import requests
import hdf5plugin  # register hdf5 bitshuffle filters
import h5py

ROOT=pathlib.Path(__file__).resolve().parents[2]
P=json.loads((ROOT/"research/seti/janus_sit2e_parkes_preflight_protocol.json").read_text())
BASE=P["source"].rstrip("/")+"/"
OUT=ROOT/"runtime"/"sit2e_parkes"
OUT.mkdir(parents=True,exist_ok=True)
MAXFILE=int(P["resources"]["maximum_file_bytes"])

class LinkParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__();self.hrefs=[]
    def handle_starttag(self,tag,attrs):
        if tag.lower()=="a":
            u=dict(attrs).get("href")
            if u is not None:self.hrefs.append(u)

def safe_url(url,base):
    p=urlparse(url);b=urlparse(base)
    return p.scheme=="https" and p.netloc==b.netloc and p.path.startswith(b.path)

def role_hint(name):
    s=name.lower()
    # Absence of an explicit label is NOT evidence of source/reference.
    if re.search(r"(?:^|[_\-.])(reference|ref|off|r)(?:[_\-.]|$)",s):
        return "REFERENCE_NAME_ONLY"
    if re.search(r"(?:^|[_\-.])(source|src|on|s)(?:[_\-.]|$)",s):
        return "SOURCE_NAME_ONLY"
    return "UNKNOWN"

def product_hint(name):
    s=name.lower()
    if not (s.endswith(".h5") or s.endswith(".hdf5")):return "not_hdf5"
    for x in ("fine","mid","time"):
        if re.search(r"(?:^|[_\-.])"+x+r"(?:[_\-.]|$)",s):return x
    return "other_hdf5"

def session_index(spec,session:requests.Session):
    root=urljoin(BASE,spec["archive_path"])
    if not safe_url(root,BASE):raise ValueError("BAD_DATE_INDEX_PATH")
    pending=[(root,0)]
    visited=set();files=[];issues=[];indexed=[]
    while pending and len(visited)<int(P["resources"]["directory_max"]):
        u,depth=pending.pop(0)
        if u in visited:continue
        visited.add(u)
        try:
            r=session.get(u,timeout=(12,35))
            r.raise_for_status()
            if len(r.content)>int(P["resources"]["index_max_bytes"]):
                raise ValueError("DIRECTORY_HTML_TOO_LARGE")
            content_type=r.headers.get("Content-Type","")
            if "html" not in content_type.lower() and "<a " not in r.text.lower():
                raise ValueError("NOT_DIRECTORY_HTML:"+content_type)
            parsed=LinkParser();parsed.feed(r.text)
            indexed.append({"url":u,"status_code":r.status_code,"index_sha256":hashlib.sha256(r.content).hexdigest(),
                            "href_count":len(parsed.hrefs)})
            for h in parsed.hrefs:
                h=unquote(h).split("?",1)[0].split("#",1)[0]
                if not h or h in ("../","./") or h.startswith("/") or h.startswith("../"):continue
                dst=urljoin(u,h)
                if not safe_url(dst,root):continue
                name=dst.rstrip("/").rsplit("/",1)[-1]
                if h.endswith("/"):
                    if depth<int(P["resources"]["recursive_depth"]) and dst not in visited:
                        pending.append((dst,depth+1))
                elif product_hint(name)!="not_hdf5":
                    files.append({"filename":name,"url":dst,"product":product_hint(name),"role_name_only":role_hint(name)})
        except Exception as e:
            issues.append({"url":u,"error_type":type(e).__name__,"message":str(e)[:220]})
    files=sorted({v["url"]:v for v in files}.values(),key=lambda v:v["url"])
    return {"session_date":spec["date"],"role":spec["role"],"index_root":root,"directories":indexed,
            "index_errors":issues,"filenames_found":len(files),"filenames":files[:1000],
            "truncated_file_listing":len(files)>1000}

def safe_head(c,session):
    result={"filename":c["filename"],"url":c["url"],"name_role":c["role_name_only"],"product":c["product"]}
    try:
        r=session.head(c["url"],allow_redirects=True,timeout=(10,28))
        if not safe_url(r.url,BASE):raise ValueError("UNEXPECTED_REDIRECT")
        r.raise_for_status()
        n=r.headers.get("Content-Length")
        result.update({"status":r.status_code,"bytes_header":int(n) if n and n.isdigit() else None,
                       "accept_ranges":r.headers.get("Accept-Ranges"),
                       "content_type":r.headers.get("Content-Type")})
    except Exception as e:
        result.update({"head_error":str(e)[:240],"head_error_type":type(e).__name__})
    return result

def source_header(path):
    with h5py.File(path,"r") as f:
        if "data" not in f:raise ValueError("MISSING_DATA_DATASET")
        d=f["data"]
        attrs={str(k).lower():v for k,v in f.attrs.items()}
        for k,v in d.attrs.items():attrs.setdefault(str(k).lower(),v)
        def norm(v):
            if isinstance(v,bytes):return v.decode("utf-8",errors="replace")
            if hasattr(v,"tolist"):v=v.tolist()
            if isinstance(v,list):return [norm(z) for z in v[:8]]
            if isinstance(v,(int,float,str,bool)) or v is None:return v
            return str(v)[:100]
        names=("tstart","tsamp","fch1","foff","nchans","nifs","nbits","source_name",
               "src_raj","src_dej","telescope_id","machine_id","rawdatafile")
        out={k:norm(attrs[k]) for k in names if k in attrs}
        required=("tstart","tsamp","fch1","foff")
        missing=[k for k in required if k not in out]
        if missing:raise ValueError("MISSING_PHYSICAL_AXES:"+",".join(missing))
        nfreq=int(d.shape[-1]);start=float(out["fch1"]);step=float(out["foff"])
        lo=min(start,start+(nfreq-1)*step)
        hi=max(start,start+(nfreq-1)*step)
        out["frequency_start_end_mhz"]=[lo,hi]
        out["channel_width_hz"]=abs(step)*1e6
        out["data_shape"]=list(d.shape)
        out["data_dtype"]=str(d.dtype)
        out["compression"]=str(d.compression)
        out["header_only_radio_power_not_read"]=True
        return out

def qualify_candidates(report,session):
    # Probe a bounded number of source/reference candidates from official HTML
    # only; never synthesize archive paths or infer role from target name.
    files=report["filenames"]
    candidate=[x for x in files if x["product"] in ("mid","time","fine","other_hdf5")]
    order={"mid":0,"time":1,"fine":2,"other_hdf5":3}
    candidate.sort(key=lambda x:(order[x["product"]],x["role_name_only"]=="UNKNOWN",x["filename"]))
    for role in ("SOURCE_NAME_ONLY","REFERENCE_NAME_ONLY"):
        selected=[x for x in candidate if x["role_name_only"]==role]
        # Early HEAD evidence on up to 4 files per role. No arbitrary blind download.
        for x in selected[:4]:
            report.setdefault("metadata_probe",[]).append(safe_head(x,session))
    # If no explicit role labels exist, inspect up to 4 files as a diagnostic
    # but do NOT treat them as a qualified Source/Reference pair.
    if not any(x["role_name_only"]=="SOURCE_NAME_ONLY" for x in candidate) or not any(x["role_name_only"]=="REFERENCE_NAME_ONLY" for x in candidate):
        for x in candidate[:4]:
            if x["url"] not in {y["url"] for y in report.get("metadata_probe",[])}:
                report.setdefault("metadata_probe",[]).append(safe_head(x,session))
    selected=[]
    for role in ("SOURCE_NAME_ONLY","REFERENCE_NAME_ONLY"):
        options=[x for x in report.get("metadata_probe",[]) if x["name_role"]==role and
                 x.get("bytes_header") and x["bytes_header"]<=MAXFILE]
        if options:
            selected.append(min(options,key=lambda x:(x["bytes_header"],order[x["product"]])))
    if len(selected)!=2:
        report["qualification"]="HOLD_NO_EXPLICIT_SMALL_HDF5_SOURCE_AND_REFERENCE_PAIR"
        return
    report["downloaded_headers"]=[]
    for x in selected:
        name=x["filename"]
        path=OUT/name
        digest=hashlib.sha256();size=0
        try:
            with session.get(x["url"],stream=True,timeout=(15,60)) as r:
                r.raise_for_status()
                if r.url!=x["url"]:raise ValueError("UNEXPECTED_REDIRECT")
                with path.open("wb") as f:
                    for block in r.iter_content(1024*1024):
                        if not block:continue
                        size+=len(block)
                        if size>MAXFILE:raise ValueError("DOWNLOAD_EXCEEDS_BUDGET")
                        digest.update(block);f.write(block)
            if x["bytes_header"]!=size:raise ValueError("HEAD_SIZE_MISMATCH")
            h=source_header(path)
            report["downloaded_headers"].append({"filename":name,"role_name_only":x["name_role"],
                "bytes":size,"first_seen_sha256":digest.hexdigest(),"physical":h,
                "independent_antenna_pointing":False})
        except Exception as e:
            report["downloaded_headers"].append({"filename":name,"error":str(e)[:400],
                                                   "error_type":type(e).__name__})
        finally:path.unlink(missing_ok=True)
    if len(report["downloaded_headers"])==2 and all("physical" in x for x in report["downloaded_headers"]):
        a,b=[x["physical"] for x in report["downloaded_headers"]]
        overlap=max(a["frequency_start_end_mhz"][0],b["frequency_start_end_mhz"][0])<min(a["frequency_start_end_mhz"][1],b["frequency_start_end_mhz"][1])
        report["source_reference_frequency_overlap"]=overlap
        report["qualification"]="PAIRED_ARCHIVE_HEADERS_AND_SHA256_OBSERVED__POINTING_HOLD" if overlap else "HOLD_FREQUENCY_GRIDS_DO_NOT_OVERLAP"
    else:
        report["qualification"]="HOLD_HDF5_DOWNLOAD_OR_HEADER_FAILURE"

def run():
    session=requests.Session()
    session.headers.update({"User-Agent":"JANUS-SIT2E-Parkes-Academic-Archive-Metadata-Qualification/1.0"})
    result={"schema":"janus.cosmos.sit2e.parkes_real_source_preflight_result.v1",
            "status":"ARCHIVE_PREFLIGHT_RUNNING","official":BASE,
            "sessions":[],"radio_power_read":False,
            "October_labels_or_scores_read":False,"thresholds_fitted":False,
            "external_independent_pointing_verified":False,
            "scientific_sensitivity_gain":None,"model_or_production_changed":False}
    for spec in P["archive_sessions"]:
        report=session_index(spec,session)
        qualify_candidates(report,session)
        result["sessions"].append(report)
        print("PARKES_SIT2E_DATE",json.dumps({
            "date":spec["date"],"directory_count":len(report["directories"]),
            "candidate_hdf5_files":report["filenames_found"],
            "candidate_role_counts":{r:sum(x["role_name_only"]==r for x in report["filenames"]) for r in ("SOURCE_NAME_ONLY","REFERENCE_NAME_ONLY","UNKNOWN")},
            "qualification":report.get("qualification"),
            "headers":report.get("downloaded_headers",[])},ensure_ascii=False,allow_nan=False),flush=True)
    count=sum(x.get("qualification")=="PAIRED_ARCHIVE_HEADERS_AND_SHA256_OBSERVED__POINTING_HOLD" for x in result["sessions"])
    result["paired_dates_qualified_at_header_level"]=count
    result["status"]=("ARCHIVE_HDF5_PAIR_METADATA_3_OF_3_PASS__ACTUAL_ANTENNA_POINTING_AND_SCIENCE_HOLD" if count==3
        else "ARCHIVE_HDF5_QUALIFICATION_INCOMPLETE__SIT2E_EXAM_HOLD")
    return result

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    try:result=run()
    except Exception as e:
        result={"status":"PREFLIGHT_SCRIPT_ERROR__EXAM_HOLD","exception":repr(e),
                "thresholds_fitted":False,"radio_power_read":False,"October_labels_or_scores_read":False}
    (OUT/"sit2e_parkes_preflight.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("PARKES_SIT2E_FINAL",json.dumps({
        "status":result["status"],"paired_dates":result.get("paired_dates_qualified_at_header_level",0),
        "radio_power_read":result["radio_power_read"],
        "thresholds_fitted":result["thresholds_fitted"]
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()
