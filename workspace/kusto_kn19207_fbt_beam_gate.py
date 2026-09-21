#!/usr/bin/env python3
import gzip, hashlib, json, math, struct
from pathlib import Path
import requests

TARGET={"id":"CD169_TOBI_LOCUS","lat":-3.865418,"lon":-12.14244}
BASENAME="sb20080107081434.xse.mb94"
RAW_BASE="https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/"
GEN_BASE=RAW_BASE+"generated/"
OUT=Path("workspace/kusto_open_seafloor_out")
OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session()
S.headers.update({"User-Agent":"JANUS-KUSTO-fbt-beam-gate/1.0"})

def fetch(url):
    r=S.get(url,timeout=180)
    r.raise_for_status()
    return r.content

def sha(b): return hashlib.sha256(b).hexdigest()

def local_nav_xy(lon,lat):
    east=(lon-TARGET["lon"])*111320.0*math.cos(math.radians(TARGET["lat"]))
    north=(lat-TARGET["lat"])*111320.0
    return east,north

def parse_fbt(data):
    off=0
    records=[]
    ids={b"V4":(4,90), b"V5":(5,98), b"cc":("comment",2), b"##":("oldcomment",30)}
    while off+2 <= len(data):
        tag=data[off:off+2]
        if tag not in ids:
            raise RuntimeError(f"Unknown FBT tag {tag!r} at byte {off}")
        version,hsize=ids[tag]
        if version in ("comment","oldcomment"):
            if off+hsize+128 > len(data): break
            off += hsize+128
            continue
        if off+hsize > len(data): break
        h=data[off:off+hsize]
        # MBF_MBLDEOIH uses network/big-endian byte order in generated FBT.
        time_d,lon,lat,sensordepth,altitude=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,beam_xwidth,beam_lwidth=struct.unpack_from(">7f",h,42)
        if version==4:
            nbath,namp,nss,sensorhead=struct.unpack_from(">4h",h,70)
            depth_scale,distance_scale=struct.unpack_from(">2f",h,78)
            tail=86
        else:
            nbath,namp,nss,sensorhead=struct.unpack_from(">4i",h,70)
            depth_scale,distance_scale=struct.unpack_from(">2f",h,86)
            tail=94
        if min(nbath,namp,nss) < 0 or nbath>100000 or namp>100000 or nss>1000000:
            raise RuntimeError(f"Implausible dimensions at {off}: {nbath},{namp},{nss}")
        pos=off+hsize
        need=nbath + 2*nbath*3 + 2*namp + 2*nss*3
        if pos+need>len(data):
            raise RuntimeError(f"Truncated data at {off}, need {need}")
        flags=list(data[pos:pos+nbath]); pos+=nbath
        def shorts(n):
            nonlocal pos
            if n==0:return []
            vals=list(struct.unpack_from(f">{n}h",data,pos)); pos+=2*n
            return vals
        bath=shorts(nbath)
        across=shorts(nbath)
        along=shorts(nbath)
        _amp=shorts(namp)
        _ss=shorts(nss); _ssx=shorts(nss); _ssy=shorts(nss)
        navx,navy=local_nav_xy(lon,lat)
        sh=math.sin(math.radians(heading)); ch=math.cos(math.radians(heading))
        beams=[]
        for i in range(nbath):
            if flags[i] == 1: # MB_FLAG_NULL
                continue
            xtrack=distance_scale*across[i]
            ltrack=distance_scale*along[i]
            east=navx + ltrack*sh + xtrack*ch
            north=navy + ltrack*ch - xtrack*sh
            dist=math.hypot(east,north)
            depth=depth_scale*bath[i]+sensordepth
            beams.append({
                "beam":i,"flag":flags[i],"good":flags[i]==0,
                "depth_m":depth,"acrosstrack_m":xtrack,"alongtrack_m":ltrack,
                "target_distance_m":dist,"target_east_m":east,"target_north_m":north
            })
        records.append({
            "version":version,"time_d":time_d,"lon":lon,"lat":lat,
            "sensordepth_m":sensordepth,"altitude_m":altitude,
            "heading_deg":heading,"speed":speed,"roll":roll,"pitch":pitch,"heave":heave,
            "beam_xwidth_deg":beam_xwidth,"beam_lwidth_deg":beam_lwidth,
            "nbath":nbath,"namp":namp,"nss":nss,"sensorhead":sensorhead,
            "depth_scale":depth_scale,"distance_scale":distance_scale,
            "beams":beams
        })
        off=pos
    return records

fbt_url=GEN_BASE+BASENAME+".fbt"
inf_url=GEN_BASE+BASENAME+".inf"
fnv_url=GEN_BASE+BASENAME+".fnv"
raw_url=RAW_BASE+BASENAME+".gz"

fbt=fetch(fbt_url)
inf=fetch(inf_url).decode("utf-8","replace")
fnv=fetch(fnv_url)
raw_gz=fetch(raw_url)
records=parse_fbt(fbt)

# exact FNV gate bracket found in prior frozen sweep
t0=1199695282.441006
t1=1199695295.051087
pad=20.0
window=[r for r in records if t0-pad <= r["time_d"] <= t1+pad]
allbeams=[(r,b) for r in window for b in r["beams"]]
good=[(r,b) for r,b in allbeams if b["good"]]
nonnull=allbeams
good.sort(key=lambda rb:rb[1]["target_distance_m"])
nonnull.sort(key=lambda rb:rb[1]["target_distance_m"])

def compact(rb):
    if not rb:return None
    r,b=rb
    return {
      "time_d":r["time_d"],"nav_lon":r["lon"],"nav_lat":r["lat"],
      "heading_deg":r["heading_deg"],"beam_count":r["nbath"],
      "beam":b
    }

nearest_good=compact(good[0] if good else None)
nearest_nonnull=compact(nonnull[0] if nonnull else None)

# local depth statistics using valid/good beams near target from pings in the bracket
neighborhoods={}
for rad in [50,100,250,500,1000]:
    vals=[b["depth_m"] for r,b in good if b["target_distance_m"]<=rad]
    neighborhoods[str(rad)]={
      "n_good_beams":len(vals),
      "depth_min_m":min(vals) if vals else None,
      "depth_max_m":max(vals) if vals else None,
      "depth_mean_m":sum(vals)/len(vals) if vals else None,
      "relief_m":max(vals)-min(vals) if vals else None
    }

receipt={
 "artifact_id":"JANUS-KUSTO-KN19207-FBT-REAL-BEAM-GATE-2026-09-21-v1.0",
 "target":TARGET,
 "source":{
   "survey":"KN192-07","platform":"Knorr","instrument":"SeaBeam 3012",
   "raw_file":BASENAME+".gz","raw_url":raw_url,"raw_gz_bytes":len(raw_gz),"raw_gz_sha256":sha(raw_gz),
   "fbt_file":BASENAME+".fbt","fbt_url":fbt_url,"fbt_bytes":len(fbt),"fbt_sha256":sha(fbt),
   "fnv_file":BASENAME+".fnv","fnv_sha256":sha(fnv),
   "inf_file":BASENAME+".inf","inf_sha256":sha(inf.encode()),"inf_text":inf[:30000]
 },
 "parser_authority":{
   "format":"MBF_MBLDEOIH / MB-System fast bathymetry (.fbt)",
   "implementation_basis":[
      "src/mbio/mbr_mbldeoih.c: V4/V5 headers and binary arrays",
      "src/mbio/mbsys_ldeoih.c: bath = depth_scale*stored_bath + sensordepth; track offsets = distance_scale*stored offsets"
   ],
   "beam_position_transform":"ship nav + MB-System across/along offsets rotated by heading in local tangent plane",
   "flag_rule":"MB_FLAG_NULL(1) excluded; nearest_good additionally requires flag==0"
 },
 "record_count":len(records),
 "window_record_count":len(window),
 "window_epoch":[t0-pad,t1+pad],
 "nearest_good_beam":nearest_good,
 "nearest_nonnull_beam":nearest_nonnull,
 "depth_neighborhoods":neighborhoods,
 "direct_beam_gate_pass":bool(nearest_good and nearest_good["beam"]["target_distance_m"]<=100.0),
 "claim_ceiling":"INDEPENDENT_PROCESSED_FAST_BATHYMETRY_BEAM_NEAR_FROZEN_EARTH_CELL__NOT_RAW_REPROCESSING__NOT_OBJECT_IDENTITY",
 "next_if_pass":"FREEZE_LOCAL_BEAM_CLOUD__COMPARE_DEPTH_WITH_GEBCO_AND_CD169__THEN_PREREGISTER_BLIND_LOCAL_MORPHOLOGY_METRICS",
 "raw_note":"Raw XSE bytes are preserved only by hash/size in this run; FBT is an MB-System generated derivative and is not equivalent to independent raw reprocessing."
}
p=OUT/"JANUS-KUSTO-KN19207-FBT-REAL-BEAM-GATE-2026-09-21-v1.0.json"
p.write_text(json.dumps(receipt,indent=2),encoding="utf-8")
print(json.dumps({
 "record_count":receipt["record_count"],
 "window_record_count":receipt["window_record_count"],
 "nearest_good_beam":receipt["nearest_good_beam"],
 "nearest_nonnull_beam":receipt["nearest_nonnull_beam"],
 "depth_neighborhoods":receipt["depth_neighborhoods"],
 "direct_beam_gate_pass":receipt["direct_beam_gate_pass"],
 "raw_gz_sha256":receipt["source"]["raw_gz_sha256"],
 "fbt_sha256":receipt["source"]["fbt_sha256"]
},indent=2))
