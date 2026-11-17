#!/usr/bin/env python3
"""Find exact Parkes source via public BL API, verify renamed records by headers."""
import hashlib,json,pathlib,re,sys
from urllib.parse import urlparse,unquote
import requests
import h5py,hdf5plugin
import r4_chunk_replay as base
import r6_track_diagnostics as diag

OUT=pathlib.Path("runtime/HIP47425_S")
OUT.mkdir(parents=True,exist_ok=True)
OLD="spliced_blc000102030405_guppi_57854_41329_56305256_HIP47425_S_0001.gpuspec.0000.h5"
MJD=57854+41329/86400
FREQS=[3426.021167,3426.129478,3426.208251]
API="https://seti.berkeley.edu/opendata/api/query-files"
# Allow only public Berkeley archive hosts, no credentials.
def archive_url(u):
    p=urlparse(u)
    return p.scheme in ("http","https") and (p.hostname=="bldata.berkeley.edu" or bool(re.fullmatch(r"blpd\d+\.ssl\.berkeley\.edu",p.hostname or "")))
base.archive_url=archive_url
receipt={"schema":"janus.seti.r7.parkes_api_recovery.v1","target":"HIP47425_S","old_file":OLD,"expected_mjd":MJD,"queries":[],"files":[],"headers":[],"observations":[],"status":"HOLD","independent_confirmation":False,"confirmed_technosignatures":0}
def save():
    (OUT/"result.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
def query(target):
    q={"target":target,"telescope":"Parkes"}
    try:
        raw,url=base.fetch(API,q,limit=16*base.MIB)
        data=json.loads(raw)
        rows=data if isinstance(data,list) else next((data[k] for k in ("data","files","results") if isinstance(data.get(k),list)),[])
        (OUT/("api-"+target+".json")).write_bytes(raw)
        receipt["queries"].append({"url":url,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),"rows":len(rows)})
        return rows
    except Exception as e:
        receipt["queries"].append({"params":q,"error":str(e)})
        return []
def parse(u):
    n=unquote(base.basename(u))
    a=re.fullmatch(r"Parkes_(\d{5})_(\d{5})_(HIP47425)(?:_([SR]))?_(fine|mid|time)\.h5",n,re.I)
    if a:return {"name":n,"day":int(a[1]),"seconds":int(a[2]),"position":(a[4] or "?").upper(),"resolution":a[5].lower(),"url":u}
    b=re.search(r"_guppi_(\d{5})_(\d{5})_(?:\d+_)?HIP47425(?:_([SR]))?_\d{4}\.gpuspec\.0000\.h5$",n,re.I)
    if b:return {"name":n,"day":int(b[1]),"seconds":int(b[2]),"position":(b[3] or "?").upper(),"resolution":"fine","url":u}
    return None
def read_header(f):
    r={"file":f,"status":"HOLD"}
    try:
        with base.RangeFile(f["url"],8*base.MIB) as remote:
            with h5py.File(remote,"r") as h:
                d=h["data"]; a=d.attrs
                f0=float(a["fch1"]);df=float(a["foff"]);ts=float(a["tsamp"])
                mjd=float(a["tstart"])
                source=a.get("source_name","")
                if isinstance(source,bytes):source=source.decode(errors="replace")
                bounds=sorted([f0,f0+(d.shape[-1]-1)*df])
                r.update(status="HEADER_READ",tstart_mjd=mjd,delta_from_original_seconds=(mjd-MJD)*86400,
                         source_name=str(source),src_raj=str(a.get("src_raj","")),src_dej=str(a.get("src_dej","")),
                         telescope_id=str(a.get("telescope_id","")),channel_hz=abs(df)*1e6,
                         frequency_range_mhz=bounds,tsamp_s=ts,shape=list(d.shape),
                         all_lines_covered=bool(all(bounds[0]<v<bounds[1] for v in FREQS)),
                         high_resolution=bool(0<abs(df)*1e6<=10),bytes_read=remote.total)
    except Exception as e:r["error"]=str(e)
    print("R7_HEADER",json.dumps(r),flush=True)
    return r
def eligible(h):
    return h["status"]=="HEADER_READ" and h["all_lines_covered"] and h["high_resolution"]
# Frozen catalog rows and all three original event frequencies.
spec=json.loads(pathlib.Path("research/seti/frozen_r4_targets.json").read_text())
raw,_=base.fetch(spec["catalogs"]["PARKES_10CM"]["url"])
assert hashlib.sha256(raw).hexdigest()==spec["catalogs"]["PARKES_10CM"]["sha256"]
import csv,io
rs=[r for r in csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))) if r["FileID"]==OLD and abs(float(r["FreqMid"])-3426.114709)<1e-7]
assert len(rs)==1
r=rs[0];fs=base.numbers(r["Freqs"]);ds=base.numbers(r["DriftRates"]);sn=base.numbers(r["SNR"])
assert len(fs)==len(ds)==len(sn)==3 and sorted(fs)==sorted(FREQS)
lines=[{"group_mhz":3426.114709,"frequency_mhz":fs[i],"catalog_drift_hz_s":ds[i],"catalog_snr":sn[i],"catalog_hit_index":i,"selection":"all_frozen_catalog_hits","group_event_count":3} for i in range(3)]
receipt["groups"]=lines
files={}
for target in ("HIP47425","HIP47425_S","HIP47425_R"):
    for row in query(target):
        if not isinstance(row,dict):continue
        u=row.get("url","")
        if not isinstance(u,str) or not archive_url(u):continue
        f=parse(u)
        if f:
            f["archive_metadata"]=row
            files[u]=f
receipt["files"]=list(files.values())
fine=[f for f in files.values() if f["resolution"]=="fine"]
print("R7_INDEX",json.dumps({"queries":receipt["queries"],"file_count":len(files),"fine_count":len(fine),"days":sorted({f["day"] for f in fine}),"original_day_files":[f for f in fine if f["day"]==57854]},ensure_ascii=False),flush=True)
save()
# Exact filename preferred, then same day and same start second. A renamed
# product is a time/target-matched proxy, not a proven byte-identical old file.
candidates=sorted((f for f in fine if f["day"]==57854 and abs(f["seconds"]-41329)<=2 and f["position"]!="R"),key=lambda f:(f["name"]!=OLD,abs(f["seconds"]-41329),f["url"]))
for f in candidates[:4]:
    h=read_header(f);receipt["headers"].append(h)
save()
valid=[h for h in receipt["headers"] if eligible(h) and abs(h["delta_from_original_seconds"])<=2 and h["source_name"].upper()=="HIP47425_S"]
if not valid:
    receipt["status"]="HOLD_EXACT_OR_TIME_MATCHED_SOURCE_NOT_VERIFIED";save()
else:
    h=valid[0];original=h["file"]
    receipt["identity_grade"]="EXACT_CATALOG_FILENAME" if original["name"]==OLD else "HEADER_MATCHED_RENAMED_PRODUCT_NOT_BYTE_IDENTITY"
    # Original plus same-day ON and OFF controls. Select chronologically, blind to intensities.
    same_day=[f for f in fine if f["day"]==57854 and f["url"]!=original["url"] and abs(f["seconds"]-41329)<=1800]
    controls=sorted((f for f in same_day if f["position"]=="R"),key=lambda f:(abs(f["seconds"]-41329),f["url"]))
    repeats=sorted((f for f in same_day if f["position"]=="S"),key=lambda f:(abs(f["seconds"]-41329),f["url"]))
    epochs=sorted((f for f in fine if f["day"]!=57854 and f["position"]=="S"),key=lambda f:(abs(f["day"]-57854),f["seconds"],f["url"]))
    selected=[(original,"ORIGINAL_ON")]+[(f,"ARCHIVE_LABELED_OFF") for f in controls[:2]]+[(f,"SAME_SESSION_ON_REPEAT") for f in repeats[:2]]+[(f,"DIFFERENT_DAY_SAME_TARGET_REPEAT") for f in epochs[:1]]
    for f,role in selected:
        h=read_header(f);receipt["headers"].append(h)
        if not eligible(h):continue
        if role=="ARCHIVE_LABELED_OFF" and h["source_name"].upper()!="HIP47425_R":continue
        if role!="ARCHIVE_LABELED_OFF" and h["source_name"].upper()!="HIP47425_S":continue
        sp=base.inspect(f["url"],lines,role,OUT,384*base.MIB if role=="ORIGINAL_ON" else 192*base.MIB)
        receipt["observations"].append(sp);save()
    receipt["status"]="ARCHIVE_REPLAY_COMPLETE_REVIEW_REQUIRED" if any(o["role"]=="ORIGINAL_ON" and o["status"]=="SLICES_READ" for o in receipt["observations"]) else "HOLD_ORIGINAL_READ_FAILED"
    save()
    diag.analyze(OUT)
print("R7_SUMMARY",json.dumps({k:v for k,v in receipt.items() if k not in ("files","headers","observations","groups")}),flush=True)
