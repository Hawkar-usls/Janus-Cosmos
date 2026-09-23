#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures, json, math, re, struct, time, os
from pathlib import Path
from urllib.parse import urljoin, urlparse

import numpy as np
import requests
from scipy.spatial import cKDTree
from shapely.geometry import shape, LineString, Point
from shapely.ops import unary_union
from shapely.prepared import prep

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace/kusto_global_anomaly_out"; OUT.mkdir(parents=True,exist_ok=True)
STAGE0=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-GLOBAL-NCEI-MULTISURVEY-OPPORTUNITY-MAP-RECEIPT-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-global-stage1-blind-morphology/1.0"}
FOOT="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
GRID=250.0
RADII=[250,500,1000,2000]
MIN_N=100
TAIL=0.005
EXTREME_ACROSS=5500.0
PATCH_CAP=5000
M_PER_DEG=111320.0

FAMILIES={
 "RELIEF_ROUGHNESS":["relief_p05_p95_m","plane_rmse_m","quadratic_rmse_m"],
 "CURVATURE":["curvature_anisotropy","curvature_strength"],
 "RADIAL":["radial_quadratic_r2","annular_contrast_abs"],
 "ANGULAR":["m2","m3","m4","m5","m6"],
 "LINEARITY":["structure_tensor_anisotropy","dominant_lineament_strength"],
 "TERRACE_STEP":["radial_step_strength"],
 "TOPOLOGIC":["radial_profile_extrema_count","isolated_extremum_persistence"]
}
STRUCTURAL=set(FAMILIES)-{"RELIEF_ROUGHNESS"}

def get(url,stream=False):
    last=None
    for a in range(7):
        try:
            r=requests.get(url,headers=UA,timeout=120,stream=stream)
            if r.status_code==429:
                time.sleep(min(30,2**a)); continue
            r.raise_for_status(); return r
        except Exception as e:
            last=e
            if a==6: raise
            time.sleep(min(30,2**a))
    raise last

def hrefs(u):
    r=get(u)
    return sorted(set(urljoin(u,h) for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I)))

def derive_base(download_url):
    p=urlparse(download_url); parts=[x for x in p.path.split("/") if x]
    ship=parts[-2]; survey=parts[-1][:-8]
    return f"https://data.ngdc.noaa.gov/platforms/ocean/ships/{ship}/{survey}/multibeam/data/"

def list_fnv(base,maxdepth=4):
    seen=set(); out=[]
    def walk(u,d):
        if u in seen or d>maxdepth:return
        seen.add(u)
        for x in hrefs(u):
            if not x.startswith(base) or x.rstrip("/")==u.rstrip("/"):continue
            if x.endswith("/"): walk(x,d+1)
            elif x.lower().endswith(".fnv"): out.append(x)
    walk(base,0)
    return sorted(set(out))

def footprint(sid):
    p={"where":f"SURVEY_ID='{sid}'","outFields":"SURVEY_ID","returnGeometry":"true","outSR":"4326","f":"geojson"}
    d=get(FOOT+"?"+requests.compat.urlencode(p)).json()
    gs=[shape(f["geometry"]) for f in d.get("features",[]) if f.get("geometry")]
    g=unary_union(gs)
    if not g.is_valid:g=g.buffer(0)
    return g

def parse_fnv_text(text):
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:
            epoch=float(p[6]); plon=float(p[15]); plat=float(p[16]); slon=float(p[17]); slat=float(p[18])
        except:continue
        rows.append((epoch,plon,plat,slon,slat))
    return rows

def make_xy(lat0):
    cos0=math.cos(math.radians(lat0))
    def xy(lon,lat): return lon*M_PER_DEG*cos0,lat*M_PER_DEG
    def ll(x,y): return x/(M_PER_DEG*cos0),y/M_PER_DEG
    return xy,ll

def sample_line_cells(plon,plat,slon,slat,xy):
    ax,ay=xy(plon,plat); bx,by=xy(slon,slat)
    d=math.hypot(bx-ax,by-ay)
    n=max(1,int(math.ceil(d/(GRID/2))))
    out=set()
    for k in range(n+1):
        t=k/n
        x=ax+(bx-ax)*t; y=ay+(by-ay)*t
        out.add((math.floor(x/GRID),math.floor(y/GRID)))
    return out

def fnv_scan(url,inter,xy):
    text=get(url).text
    rows=parse_fnv_text(text)
    hit=False; cells=set()
    for _,plon,plat,slon,slat in rows:
        seg=LineString([(plon,plat),(slon,slat)])
        if seg.intersects(inter):
            hit=True
            for cell in sample_line_cells(plon,plat,slon,slat,xy):
                cx=(cell[0]+0.5)*GRID; cy=(cell[1]+0.5)*GRID
                lon,lat=ll_global(cx,cy) if False else (None,None)
                cells.add(cell)
    return url,hit,cells,len(rows)

def fbt_url_from_fnv(u):
    return u[:-4]+".fbt"

def parse_fbt(data,name,xy,inter):
    off=0; parts=[]; recs=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    pint=prep(inter)
    minx,miny,maxx,maxy=inter.bounds
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids: raise RuntimeError(f"{name}: bad tag {tag!r} @{off}")
        ver,hs=ids[tag]
        if isinstance(ver,str):
            off+=hs+128; continue
        h=data[off:off+hs]
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,bx,bl=struct.unpack_from(">7f",h,42)
        if ver==4:
            nb,na,ns,shd=struct.unpack_from(">4h",h,70); dscale,xscale=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,shd=struct.unpack_from(">4i",h,70); dscale,xscale=struct.unpack_from(">2f",h,86)
        p=off+hs
        flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p).copy();p+=nb
        def arr(n):
            nonlocal p
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(float);p+=2*n;return a
        bath=arr(nb); across=arr(nb); along=arr(nb);p+=2*na+6*ns
        navx,navy=xy(lon,lat); hr=math.radians(heading); sh=math.sin(hr); ch=math.cos(hr)
        xt=xscale*across; lt=xscale*along
        east=navx+lt*sh+xt*ch; north=navy+lt*ch-xt*sh; depth=dscale*bath+sd
        good=np.where(flags==0)[0]
        if len(good):
            rows=[]
            for i in good:
                # lon/lat only for polygon clipping after cheap bbox
                lonb,latb=ll_current(float(east[i]),float(north[i]))
                if lonb<minx or lonb>maxx or latb<miny or latb>maxy: continue
                if not pint.intersects(Point(lonb,latb)): continue
                rows.append((float(east[i]),float(north[i]),float(depth[i]),abs(float(xt[i]))))
            if rows:parts.append(np.asarray(rows,float))
        off=p;recs+=1
    return (np.vstack(parts) if parts else np.empty((0,4),float)),recs

def deterministic_cap(p,cx,cy):
    if len(p)<=PATCH_CAP:return p
    dx=p[:,0]-cx;dy=p[:,1]-cy
    ang=np.arctan2(dy,dx); rr=np.hypot(dx,dy)
    order=np.lexsort((rr,ang))
    idx=np.linspace(0,len(order)-1,PATCH_CAP,dtype=int)
    return p[order[idx]]

def plane_residual(x,y,z):
    A=np.column_stack([np.ones(len(x)),x,y])
    coef=np.linalg.lstsq(A,z,rcond=None)[0]; res=z-A@coef
    return res,float(np.sqrt(np.mean(res*res)))

def safe_ratio(a,b):
    return None if b is None or abs(b)<1e-12 else float(a/b)

def patch_metrics(p,cx,cy,rad):
    if len(p)<MIN_N:return None
    p=deterministic_cap(p,cx,cy)
    x=p[:,0]-cx;y=p[:,1]-cy;z=p[:,2]; across=p[:,3];rr=np.hypot(x,y)
    p05,p95=np.percentile(z,[5,95])
    resid,prmse=plane_residual(x,y,z)
    mad=np.median(np.abs(resid-np.median(resid))); scale=max(1e-9,1.4826*mad)
    Q=np.column_stack([np.ones(len(x)),x,y,x*x,x*y,y*y])
    qc=np.linalg.lstsq(Q,z,rcond=None)[0]; qres=z-Q@qc
    H=np.array([[2*qc[3],qc[4]],[qc[4],2*qc[5]]]); eig=np.linalg.eigvalsh(H); ab=np.abs(eig)
    anis=float(ab.max()/max(1e-12,ab.min()))
    curv=float(np.linalg.norm(H))
    Rm=np.column_stack([np.ones(len(rr)),rr,rr*rr])
    rc=np.linalg.lstsq(Rm,resid,rcond=None)[0]; pred=Rm@rc
    sst=float(np.sum((resid-resid.mean())**2)); ssr=float(np.sum((resid-pred)**2))
    radial=0.0 if sst<=0 else float(1-ssr/sst)
    # inner-vs-outer annular contrast
    q1=np.quantile(rr,.33);q2=np.quantile(rr,.67)
    inn=np.median(resid[rr<=q1]) if np.any(rr<=q1) else 0
    out=np.median(resid[rr>=q2]) if np.any(rr>=q2) else 0
    ann=float(abs(out-inn)/scale)
    theta=np.arctan2(y,x)
    amps={}
    rms=float(np.sqrt(np.mean(resid*resid)))
    for m in range(2,7):
        amp=2*abs(np.mean(resid*np.exp(-1j*m*theta)))
        amps[f"m{m}"]=0.0 if rms==0 else float(amp/rms)
    # structure tensor proxy on high residual geometry
    high=np.abs(resid)>=np.median(np.abs(resid))
    if np.sum(high)>=10:
        C=np.cov(np.column_stack([x[high],y[high]]).T)
        ce=np.linalg.eigvalsh(C); lin=float(ce.max()/max(1e-9,ce.min()))
        lstrength=float((ce.max()-ce.min())/max(1e-9,ce.max()+ce.min()))
    else: lin=1.0;lstrength=0.0
    # radial bins for steps/topologic proxy
    edges=np.arange(0,rad+50,50.0); meds=[]; centers=[]
    for a,b in zip(edges[:-1],edges[1:]):
        m=(rr>=a)&(rr<b)
        if np.sum(m)>=5:
            meds.append(float(np.median(z[m])));centers.append((a+b)/2)
    jumps=[abs(meds[i+1]-meds[i]) for i in range(len(meds)-1)]
    step=(max(jumps)/scale) if jumps else 0.0
    ext=0
    persistence=0.0
    for i in range(1,len(meds)-1):
        if (meds[i]>meds[i-1] and meds[i]>meds[i+1]) or (meds[i]<meds[i-1] and meds[i]<meds[i+1]):
            ext+=1
            persistence=max(persistence,min(abs(meds[i]-meds[i-1]),abs(meds[i]-meds[i+1]))/scale)
    center_mask=rr<=max(50,rad*0.15)
    outer_mask=(rr>=rad*0.6)&(rr<=rad)
    if np.sum(center_mask)>=5 and np.sum(outer_mask)>=5:
        persistence=max(persistence,abs(float(np.median(z[center_mask])-np.median(z[outer_mask])))/scale)
    return {
      "n_good_beams":int(len(p)),
      "extreme_beam_fraction":float(np.mean(across>=EXTREME_ACROSS)),
      "relief_p05_p95_m":float(p95-p05),
      "plane_rmse_m":prmse,
      "quadratic_rmse_m":float(np.sqrt(np.mean(qres*qres))),
      "curvature_anisotropy":anis,
      "curvature_strength":curv,
      "radial_quadratic_r2":radial,
      "annular_contrast_abs":ann,
      **amps,
      "structure_tensor_anisotropy":lin,
      "dominant_lineament_strength":lstrength,
      "radial_step_strength":float(step),
      "radial_profile_extrema_count":int(ext),
      "isolated_extremum_persistence":float(persistence)
    }

def empirical_rank(vals):
    order=np.argsort(vals,kind="mergesort"); ranks=np.empty(len(vals),float)
    ranks[order]=(np.arange(len(vals))+0.5)/len(vals)
    return ranks

regions_out=[]
rank_filter=os.environ.get("KUSTO_REGION_RANK")
selected_regions=STAGE0["stage1_regions"]
if rank_filter:
    selected_regions=[r for r in selected_regions if str(r["rank"])==str(rank_filter)]
    if len(selected_regions)!=1:
        raise RuntimeError(f"Frozen region rank {rank_filter} not found exactly once")
for reg in selected_regions:
    rank=reg["rank"]; A=reg["survey_a"]; B=reg["survey_b"]
    sidA=A["survey_id"];sidB=B["survey_id"]
    print("REGION",rank,sidA,sidB,flush=True)
    inter=footprint(sidA).intersection(footprint(sidB))
    if inter.is_empty: raise RuntimeError("frozen pair intersection vanished")
    lat0=float(inter.representative_point().y)
    xy,ll=make_xy(lat0)
    globals()["ll_current"]=ll
    baseA=derive_base(A["download_url"]);baseB=derive_base(B["download_url"])
    fnvA=list_fnv(baseA);fnvB=list_fnv(baseB)
    def scan_set(urls):
        hits=[];cells=set();meta=[]
        def one(u):
            text=get(u).text;rows=parse_fnv_text(text);hit=False;cc=set()
            for _,plon,plat,slon,slat in rows:
                seg=LineString([(plon,plat),(slon,slat)])
                if seg.intersects(inter):
                    hit=True;cc.update(sample_line_cells(plon,plat,slon,slat,xy))
            return u,hit,cc,len(rows)
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
            futs=[ex.submit(one,u) for u in urls]
            for k,f in enumerate(concurrent.futures.as_completed(futs),1):
                u,hit,cc,nr=f.result()
                if hit:
                    hits.append(u);cells.update(cc);meta.append({"fnv":u,"rows":nr,"cells":len(cc)})
                if k%100==0:print(" fnv",k,"/",len(urls),flush=True)
        return sorted(hits),cells,meta
    hitA,cellsA,metaA=scan_set(fnvA);hitB,cellsB,metaB=scan_set(fnvB)
    common=sorted(cellsA & cellsB)
    centers=[]
    pint=prep(inter)
    for ix,iy in common:
        cx=(ix+0.5)*GRID;cy=(iy+0.5)*GRID;lon,lat=ll(cx,cy)
        if pint.intersects(Point(lon,lat)):
            centers.append({"ix":ix,"iy":iy,"x":cx,"y":cy,"lon":lon,"lat":lat})
    print(" common centers",len(centers),"hit files",len(hitA),len(hitB),flush=True)
    # Discovery depth only from survey A
    def fetch_fbt(u):
        fu=fbt_url_from_fnv(u);b=get(fu).content
        pts,recs=parse_fbt(b,fu.rsplit("/",1)[-1],xy,inter)
        return fu,len(b),pts,recs
    parts=[];manifest=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
        futs=[ex.submit(fetch_fbt,u) for u in hitA]
        for k,f in enumerate(concurrent.futures.as_completed(futs),1):
            fu,nb,pts,recs=f.result();parts.append(pts);manifest.append({"fbt":fu,"bytes":nb,"records":recs,"good_points_in_overlap":int(len(pts))})
            if k%20==0:print(" fbt",k,"/",len(hitA),flush=True)
    points=np.vstack(parts) if parts else np.empty((0,4),float)
    tree=cKDTree(points[:,:2]) if len(points) else None
    records=[]
    for ci,c in enumerate(centers):
        rec={**c,"metrics":{}}
        for rad in RADII:
            idx=tree.query_ball_point([c["x"],c["y"]],rad) if tree is not None else []
            pp=points[np.asarray(idx,dtype=int)] if len(idx) else np.empty((0,4))
            rec["metrics"][str(rad)]=patch_metrics(pp,c["x"],c["y"],rad)
        records.append(rec)
        if (ci+1)%1000==0:print(" metrics",ci+1,"/",len(centers),flush=True)
    # Empirical extremes and robust generic score per radius
    thresholds={}
    for rad in RADII:
        rr=str(rad);thresholds[rr]={}
        keys=sorted({x for xs in FAMILIES.values() for x in xs})
        robust_z={i:[] for i in range(len(records))}
        for key in keys:
            vals=[];refs=[]
            for i,r in enumerate(records):
                m=r["metrics"][rr]
                if m is not None and m.get(key) is not None and np.isfinite(m[key]):
                    vals.append(float(m[key]));refs.append(i)
            if not vals:continue
            arr=np.asarray(vals,float);med=float(np.median(arr));mad=float(np.median(np.abs(arr-med)));scale=max(1e-12,1.4826*mad)
            lo=float(np.quantile(arr,TAIL));hi=float(np.quantile(arr,1-TAIL));pcts=empirical_rank(arr)
            thresholds[rr][key]={"n":len(arr),"q005":lo,"q995":hi,"median":med,"mad":mad}
            for j,i in enumerate(refs):
                m=records[i]["metrics"][rr]
                m.setdefault("_percentiles",{})[key]=float(pcts[j])
                m.setdefault("_extremes",{})[key]=bool(arr[j]<=lo or arr[j]>=hi)
                robust_z[i].append(abs((arr[j]-med)/scale))
        scores=[];refs=[]
        for i,zs in robust_z.items():
            m=records[i]["metrics"][rr]
            if m is None or not zs:continue
            score=0.5*max(zs)+0.5*float(np.mean(zs))
            m["generic_robust_score"]=score;scores.append(score);refs.append(i)
        if scores:
            arr=np.asarray(scores,float);q=float(np.quantile(arr,1-TAIL));thresholds[rr]["generic_robust_score"]={"n":len(arr),"q995":q}
            for val,i in zip(arr,refs):
                records[i]["metrics"][rr]["generic_extreme"]=bool(val>=q)
    flags=[]
    for r in records:
        for rad in RADII:
            rr=str(rad);m=r["metrics"][rr]
            if m is None:continue
            votes=[];details={}
            for fam,keys in FAMILIES.items():
                ex=[k for k in keys if m.get("_extremes",{}).get(k,False)]
                if ex:votes.append(fam);details[fam]=ex
            generic=bool(m.get("generic_extreme",False))
            follow=(len(votes)>=2 and any(v in STRUCTURAL for v in votes)) or generic
            dominated=m["extreme_beam_fraction"]>0.50
            flags.append({
              "lon":r["lon"],"lat":r["lat"],"radius_m":rad,
              "family_votes":votes,"extreme_metrics_by_family":details,
              "generic_extreme":generic,"generic_robust_score":m.get("generic_robust_score"),
              "nuisance_dominated":dominated,
              "promotion_suppressed":bool(follow and dominated),
              "followup_flag":bool(follow and not dominated),
              "metrics":{k:v for k,v in m.items() if not k.startswith("_")}
            })
    candidates=[x for x in flags if x["followup_flag"]]
    candidates.sort(key=lambda x:(-len(x["family_votes"]),-float(x.get("generic_robust_score") or 0),x["lat"],x["lon"],x["radius_m"]))
    kept=[]
    for z in candidates:
        zx,zy=xy(z["lon"],z["lat"])
        if any(math.hypot(zx-xy(q["lon"],q["lat"])[0],zy-xy(q["lon"],q["lat"])[1])<2000 for q in kept):continue
        kept.append(z)
        if len(kept)>=20:break
    suppressed=[x for x in flags if x["promotion_suppressed"]]
    regions_out.append({
      "rank":rank,"discovery_survey":sidA,"heldout_survey":sidB,
      "intersection_area_deg2":float(inter.area),
      "fnv_discovery_hits":len(hitA),"fnv_heldout_hits":len(hitB),
      "common_support_center_count":len(centers),
      "discovery_good_beam_points_in_overlap":int(len(points)),
      "thresholds":thresholds,
      "pre_dedup_followup_count":len(candidates),
      "spatially_deduplicated_followup_count":len(kept),
      "followup_candidates":kept,
      "nuisance_suppressed_count":len(suppressed),
      "nuisance_suppressed_examples":suppressed[:20],
      "discovery_manifest":manifest,
      "heldout_depth_read":False
    })

out={
 "artifact_id":"JANUS-KUSTO-GLOBAL-STAGE1-BLIND-MORPHOLOGY-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-GLOBAL-STAGE1-BLIND-MORPHOLOGY-EXECUTOR-PREREG-2026-09-23-v1.0.json",
 "implementation_addendum":"data/cousteau/JANUS-KUSTO-GLOBAL-STAGE1-IMPLEMENTATION-ADDENDUM-2026-09-23-v1.0.json",
 "stage0_receipt":"data/cousteau/JANUS-KUSTO-GLOBAL-NCEI-MULTISURVEY-OPPORTUNITY-MAP-RECEIPT-2026-09-23-v1.0.json",
 "stage0p5_receipt":"data/cousteau/JANUS-KUSTO-GLOBAL-STAGE0P5-SOURCE-RESOLVER-RECEIPT-2026-09-23-v1.0.json",
 "regions":regions_out,
 "heldout_survey_depth_read":False,
 "claim_ceiling":"GLOBAL_STAGE1_BLIND_MORPHOLOGY_CANDIDATE_GENERATION_ONLY__NO_CROSS_SURVEY_PERSISTENCE_YET"
}
suffix=f"-RANK{rank_filter}" if rank_filter else ""
p=OUT/f"JANUS-KUSTO-GLOBAL-STAGE1-BLIND-MORPHOLOGY-RUN-2026-09-23-v1.0{suffix}.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "regions":[{
   "rank":r["rank"],"discovery_survey":r["discovery_survey"],"heldout_survey":r["heldout_survey"],
   "common_support_center_count":r["common_support_center_count"],
   "discovery_good_beam_points_in_overlap":r["discovery_good_beam_points_in_overlap"],
   "pre_dedup_followup_count":r["pre_dedup_followup_count"],
   "spatially_deduplicated_followup_count":r["spatially_deduplicated_followup_count"],
   "nuisance_suppressed_count":r["nuisance_suppressed_count"],
   "top_candidates":r["followup_candidates"][:10]
 } for r in regions_out],
 "heldout_survey_depth_read":False
},indent=2))
