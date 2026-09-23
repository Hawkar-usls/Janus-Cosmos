#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures,hashlib,json,math,re,struct,time
from pathlib import Path
from urllib.parse import urljoin
import numpy as np
import requests
from scipy.stats import wasserstein_distance
from shapely.geometry import shape,LineString
from shapely.ops import unary_union,transform as shp_transform
from shapely import intersects_xy
from pyproj import Transformer,Geod

OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
ARC="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
ROOTS={
 "PF0501":"https://data.ngdc.noaa.gov/platforms/ocean/ships/pathfinder/PF0501/multibeam/data/version1/MB/generated/",
 "RB1202":"https://data.ngdc.noaa.gov/platforms/ocean/ships/ronald_h._brown/RB1202/multibeam/data/version1/MB/generated/"
}
EXPECTED={"PF0501":{"files":18,"rows":9056},"RB1202":{"files":26,"rows":9941}}
CELL=100.0
BINS=[("INNER",0.0,3000.0),("MID",3000.0,4500.0),("OUTER",4500.0,5500.0),("EXTREME",5500.0,float("inf"))]
UA={"User-Agent":"JANUS-KUSTO-first-unseen-transfer/1.0"}
GEOD=Geod(ellps="WGS84")
TO_UTM=Transformer.from_crs("EPSG:4326","EPSG:32618",always_xy=True)
CHANNELS=[
 ("LINE18_ERA","KN192-07_2008"),
 ("LINE18_ERA","KNOX15RR_2008"),
 ("LINE32_ERA","KN192-07_2008"),
 ("LINE32_ERA","KNOX15RR_2008")
]

def get(u,params=None,timeout=180):
    last=None
    for i in range(8):
        try:
            r=requests.get(u,params=params,headers=UA,timeout=timeout)
            if r.status_code==429:
                time.sleep(min(30,2**i));continue
            r.raise_for_status();return r
        except Exception as e:
            last=e
            if i==7:raise
            time.sleep(min(30,2**i))
    raise last

def footprint(sid):
    d=get(ARC,{
      "where":f"SURVEY_ID='{sid}'","outFields":"SURVEY_ID",
      "returnGeometry":"true","outSR":"4326","f":"geojson"
    }).json()
    gs=[shape(f["geometry"]) for f in d.get("features",[]) if f.get("geometry")]
    if not gs:raise RuntimeError("missing footprint "+sid)
    return unary_union(gs)

INTER_WGS=footprint("PF0501").intersection(footprint("RB1202"))
if INTER_WGS.is_empty:raise RuntimeError("empty footprint intersection")
INTER_UTM=shp_transform(TO_UTM.transform,INTER_WGS)

def list_fnv(root):
    text=get(root,timeout=90).text
    return sorted(set(urljoin(root,h) for h in re.findall(r'href=["\']([^"\']+\.fnv)["\']',text,re.I)))

def fnv_hits(u):
    rows=[]
    for line in get(u,timeout=90).text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:
            epoch=float(p[6]);plon=float(p[15]);plat=float(p[16]);slon=float(p[17]);slat=float(p[18])
        except:continue
        if LineString([(plon,plat),(slon,slat)]).intersects(INTER_WGS):
            rows.append(epoch)
    return u,rows

def freeze_file_set(sid):
    fnvs=list_fnv(ROOTS[sid]);hit_files=[];row_count=0
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        futs=[ex.submit(fnv_hits,u) for u in fnvs]
        for f in concurrent.futures.as_completed(futs):
            u,rows=f.result()
            if rows:
                row_count+=len(rows)
                hit_files.append({"fnv_url":u,"fnv_file":u.rsplit("/",1)[-1],"hit_rows":len(rows),
                                  "fbt_url":u[:-4]+".fbt"})
    hit_files.sort(key=lambda x:x["fnv_file"])
    exp=EXPECTED[sid]
    if len(hit_files)!=exp["files"] or row_count!=exp["rows"]:
        raise RuntimeError(f"FROZEN_FNV_SET_MISMATCH {sid}: files={len(hit_files)} rows={row_count} expected={exp}")
    return hit_files

FILESETS={sid:freeze_file_set(sid) for sid in ROOTS}
print("frozen file sets", {k:(len(v),sum(x["hit_rows"] for x in v)) for k,v in FILESETS.items()},flush=True)

def parse_fbt_to_grid(data,sid,grid):
    off=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("comment",30)}
    records=0;good_inside=0
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids:
            raise RuntimeError(f"FBT_BAD_TAG {sid} {tag!r} @{off}")
        ver,hs=ids[tag]
        if isinstance(ver,str):
            off+=hs+128;continue
        if off+hs>len(data):raise RuntimeError("truncated header")
        h=data[off:off+hs]
        t,lon,lat,sensordepth,alt=struct.unpack_from(">5d",h,2)
        heading=struct.unpack_from(">f",h,42)[0]
        if ver==4:
            nb,na,ns,head=struct.unpack_from(">4h",h,70);ds,xs=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,head=struct.unpack_from(">4i",h,70);ds,xs=struct.unpack_from(">2f",h,86)
        p=off+hs
        if min(nb,na,ns)<0 or nb>10000 or na>10000 or ns>100000:raise RuntimeError("bad FBT dimensions")
        flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p).copy();p+=nb
        def arr(n):
            nonlocal p
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(np.float64);p+=2*n;return a
        bath=arr(nb);across=arr(nb);along=arr(nb)
        p+=2*na+6*ns
        good=np.where(flags==0)[0]
        if len(good):
            z=ds*bath[good]+sensordepth
            xt=xs*across[good];lt=xs*along[good]
            hr=math.radians(heading)
            east=lt*math.sin(hr)+xt*math.cos(hr)
            north=lt*math.cos(hr)-xt*math.sin(hr)
            dist=np.hypot(east,north)
            az=np.degrees(np.arctan2(east,north))
            lons,lats,_=GEOD.fwd(np.full(len(good),lon),np.full(len(good),lat),az,dist)
            xx,yy=TO_UTM.transform(lons,lats)
            xx=np.asarray(xx);yy=np.asarray(yy)
            inside=np.asarray(intersects_xy(INTER_UTM,xx,yy),dtype=bool)
            for x0,y0,z0,a0 in zip(xx[inside],yy[inside],z[inside],np.abs(xt[inside])):
                key=(math.floor(x0/CELL),math.floor(y0/CELL))
                cell=grid.setdefault(key,{"depth":[],"abs_across":[]})
                cell["depth"].append(float(z0))
                if sid=="PF0501":cell["abs_across"].append(float(a0))
            good_inside+=int(np.sum(inside))
        records+=1;off=p
    return records,good_inside

def load_grid(sid):
    grid={};manifest=[]
    def one(rec):
        b=get(rec["fbt_url"],timeout=240).content
        return rec,hashlib.sha256(b).hexdigest(),len(b),b
    # modest concurrency for download; parsing sequentially keeps memory bounded
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        futs=[ex.submit(one,r) for r in FILESETS[sid]]
        for k,f in enumerate(concurrent.futures.as_completed(futs),1):
            rec,sha,n,b=f.result()
            nr,ni=parse_fbt_to_grid(b,sid,grid)
            manifest.append({"fbt_file":rec["fbt_url"].rsplit("/",1)[-1],"sha256":sha,"bytes":n,
                             "records":nr,"good_beams_inside_polygon":ni})
            print(sid,"parsed",k,"/",len(FILESETS[sid]),"cells",len(grid),flush=True)
    final={}
    for key,v in grid.items():
        final[key]={
          "depth_m":float(np.median(v["depth"])),
          "beam_count":len(v["depth"]),
          "abs_across_m":float(np.median(v["abs_across"])) if sid=="PF0501" and v["abs_across"] else None
        }
    manifest.sort(key=lambda x:x["fbt_file"])
    return final,manifest

PF,PFMAN=load_grid("PF0501")
RB,RBMAN=load_grid("RB1202")
common=sorted(set(PF)&set(RB))
if not common:raise RuntimeError("NO_COMMON_DEPTH_CELLS")
unseen=[]
for k in common:
    ac=PF[k]["abs_across_m"]
    resid=PF[k]["depth_m"]-RB[k]["depth_m"]
    label=None
    for name,lo,hi in BINS:
        if lo<=ac<hi:label=name;break
    unseen.append({"key":[k[0],k[1]],"abs_across_m":ac,"bin":label,
                   "PF0501_depth_m":PF[k]["depth_m"],"RB1202_depth_m":RB[k]["depth_m"],
                   "signed_residual_m":resid,"absolute_residual_m":abs(resid),
                   "PF_beams":PF[k]["beam_count"],"RB_beams":RB[k]["beam_count"]})

# Load frozen CD169 library.
lib_path=Path("data/cousteau/JANUS-KUSTO-CD169-EM12-EMPIRICAL-RESIDUAL-BOOTSTRAP-LIBRARY-2026-09-23-v1.0.json")
lib_bytes=lib_path.read_bytes();LIB=json.loads(lib_bytes)
libsha=hashlib.sha256(lib_bytes).hexdigest()

unseen_by={name:[] for name,_,_ in BINS}
for r in unseen:unseen_by[r["bin"]].append(r["signed_residual_m"])

def arr_abs(raw,surv,binname):
    x=LIB["library"][raw][surv][binname]
    vals=x.get("absolute_residual_m")
    if vals is None:
        vals=[abs(v) for v in x.get("signed_residual_m",[])]
    return np.asarray(vals,dtype=float)

# Shared bin eligibility exactly as preregistered.
elig=[]
elig_detail={}
for name,lo,hi in BINS:
    un=np.asarray(unseen_by[name],dtype=float)
    calib_ns={}
    ok=len(un)>=50
    for raw,surv in CHANNELS:
        n=len(arr_abs(raw,surv,name));calib_ns[f"{raw} x {surv}"]=n
        ok=ok and n>=30
    elig_detail[name]={"unseen_n":len(un),"calibration_n":calib_ns,"eligible":bool(ok)}
    if ok:elig.append(name)

support_pass=(len(elig)>=3 and "EXTREME" in elig)
channel_results={}
passes=0
if support_pass:
    total_unseen=sum(len(unseen_by[b]) for b in elig)
    weights={b:len(unseen_by[b])/total_unseen for b in elig}
    for raw,surv in CHANNELS:
        pooled=np.concatenate([arr_abs(raw,surv,b) for b,_,_ in BINS])
        pool_tail=float(np.mean(pooled>=100))
        v2w=basew=v2t=baset=0.0
        per={}
        for b in elig:
            u=np.abs(np.asarray(unseen_by[b],dtype=float))
            c=arr_abs(raw,surv,b)
            w=weights[b]
            wv=float(wasserstein_distance(u,c))
            wb=float(wasserstein_distance(u,pooled))
            up=float(np.mean(u>=100));cp=float(np.mean(c>=100))
            tv=abs(up-cp);tb=abs(up-pool_tail)
            v2w+=w*wv;basew+=w*wb;v2t+=w*tv;baset+=w*tb
            per[b]={"weight":w,"unseen_n":len(u),"unseen_tail_ge100":up,
                    "calib_bin_tail_ge100":cp,"calib_pool_tail_ge100":pool_tail,
                    "V2_wasserstein":wv,"BASELINE_wasserstein":wb,
                    "V2_tail_error":tv,"BASELINE_tail_error":tb}
        cp=bool(v2w<basew and v2t<baset)
        passes+=int(cp)
        channel_results[f"{raw} x {surv}"]={
          "eligible_bins":elig,"per_bin":per,
          "V2_weighted_wasserstein":v2w,"BASELINE_weighted_wasserstein":basew,
          "V2_weighted_tail_error":v2t,"BASELINE_weighted_tail_error":baset,
          "wasserstein_ratio_V2_over_baseline":None if basew==0 else v2w/basew,
          "tail_error_ratio_V2_over_baseline":None if baset==0 else v2t/baset,
          "channel_pass":cp
        }

def stats(vals):
    a=np.asarray(vals,dtype=float);aa=np.abs(a)
    if len(a)==0:return {"n":0}
    return {"n":len(a),"median_signed_m":float(np.median(a)),"median_abs_m":float(np.median(aa)),
            "p90_abs_m":float(np.percentile(aa,90)),"p99_abs_m":float(np.percentile(aa,99)),
            "fraction_ge50":float(np.mean(aa>=50)),"fraction_ge100":float(np.mean(aa>=100)),
            "fraction_ge150":float(np.mean(aa>=150))}
unseen_stats={b:stats(unseen_by[b]) for b,_,_ in BINS}
ext=np.abs(np.asarray(unseen_by["EXTREME"],float));non=np.abs(np.asarray([v for b in ["INNER","MID","OUTER"] for v in unseen_by[b]],float))
tailratio=None
if len(ext) and len(non):
    ep=float(np.mean(ext>=100));np0=float(np.mean(non>=100))
    tailratio=ep/np0 if np0>0 else (float("inf") if ep>0 else None)

if not support_pass: verdict="INSUFFICIENT_UNSEEN_BIN_SUPPORT"
elif passes>=3: verdict="PASS_UNSEEN_TRANSFER"
elif passes==2: verdict="MIXED_UNSEEN_TRANSFER"
else: verdict="FAIL_UNSEEN_TRANSFER"

out={
 "artifact_id":"JANUS-KUSTO-FIRST-UNSEEN-ACROSS-TRACK-RESIDUAL-TRANSFER-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-FIRST-UNSEEN-ACROSS-TRACK-RESIDUAL-TRANSFER-PREREG-2026-09-23-v1.0.json",
 "depth_read_after_prereg":True,
 "frozen_intersection":{"area_deg2":float(INTER_WGS.area),"bounds":list(INTER_WGS.bounds)},
 "frozen_file_sets":{"PF0501":FILESETS["PF0501"],"RB1202":FILESETS["RB1202"]},
 "fbt_manifests":{"PF0501":PFMAN,"RB1202":RBMAN},
 "grid":{"crs":"EPSG:32618","cell_m":CELL,"PF_cells":len(PF),"RB_cells":len(RB),"common_cells":len(common)},
 "calibration_library":{"path":str(lib_path),"sha256":libsha},
 "bin_eligibility":elig_detail,
 "shared_eligible_bins":elig,
 "support_pass":support_pass,
 "unseen_bin_statistics":unseen_stats,
 "unseen_extreme_to_nonextreme_tail_ge100_ratio":tailratio,
 "channel_results":channel_results,
 "channel_pass_count":passes,
 "authoritative_verdict":verdict,
 "signed_direction_used_for_transfer_verdict":False,
 "claim_ceiling":"FIRST_UNSEEN_TRANSFER_TEST_OF_CD169_ACROSS_TRACK_CONDITIONED_RESIDUAL_MAGNITUDE_MODEL__NOT_UNIVERSAL_PREDICTIVE_SKILL"
}
p=OUT/"JANUS-KUSTO-FIRST-UNSEEN-ACROSS-TRACK-RESIDUAL-TRANSFER-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "grid":out["grid"],"bin_eligibility":elig_detail,"shared_eligible_bins":elig,
 "support_pass":support_pass,"unseen_bin_statistics":unseen_stats,
 "unseen_extreme_to_nonextreme_tail_ge100_ratio":tailratio,
 "channel_results":channel_results,"channel_pass_count":passes,
 "authoritative_verdict":verdict
},indent=2))
