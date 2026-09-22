#!/usr/bin/env python3
import concurrent.futures, hashlib, json, math, re, struct
from pathlib import Path
from urllib.parse import urljoin
import numpy as np
import requests
from scipy.stats import spearmanr

CANDS=[
 {"id":"KN19207_CAND_002","lat":-3.9727527956056825,"lon":-12.272824298723462},
 {"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432},
]
SURVEYS={
 "KN192-07":{
   "platform":"Knorr","instrument":"SeaBeam 3012",
   "base":"https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/generated/"
 },
 "KNOX15RR":{
   "platform":"Roger Revelle","instrument":"Kongsberg EM120",
   "base":"https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/data/version1/MB/generated/"
 }
}
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-cross-survey-beam-replication/1.0"}
M=111320.0
LAT0=-4.0
COS0=math.cos(math.radians(LAT0))
SELECT_M=1600.0
GRID_R=1000.0
CELL=100.0

def gxy(lon,lat):
    return lon*M*COS0, lat*M

def dist_point_seg(px,py,ax,ay,bx,by):
    vx,vy=bx-ax,by-ay; wx,wy=px-ax,py-ay
    den=vx*vx+vy*vy
    t=0.0 if den==0 else max(0.0,min(1.0,(wx*vx+wy*vy)/den))
    qx,qy=ax+t*vx,ay+t*vy
    return math.hypot(px-qx,py-qy)

def parse_fnv(text):
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19: continue
        try:
            rows.append({
              "t":float(p[6]),"navlon":float(p[7]),"navlat":float(p[8]),
              "portlon":float(p[15]),"portlat":float(p[16]),
              "stbdlon":float(p[17]),"stbdlat":float(p[18])
            })
        except: pass
    return rows

def fnv_min_cross_distance(rows,c):
    px,py=gxy(c["lon"],c["lat"])
    best=float("inf")
    for r in rows:
        ax,ay=gxy(r["portlon"],r["portlat"]); bx,by=gxy(r["stbdlon"],r["stbdlat"])
        d=dist_point_seg(px,py,ax,ay,bx,by)
        if d<best: best=d
    return best

def fetch(url,timeout=120):
    r=requests.get(url,headers=UA,timeout=timeout)
    r.raise_for_status()
    return r.content

def list_files(base,ext):
    html=fetch(base,60).decode("latin1","replace")
    return sorted(set(re.findall(r'href="([^"]+\'+ext+r')"',html,re.I)))

def fetch_fnv(base,name):
    b=fetch(urljoin(base,name),90)
    return name,hashlib.sha256(b).hexdigest(),parse_fnv(b.decode("utf-8","replace"))

def parse_fbt(data,name):
    off=0; chunks=[]; nrec=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids: raise RuntimeError(f"{name}: bad tag {tag!r} at {off}")
        ver,hs=ids[tag]
        if isinstance(ver,str):
            off += hs+128
            continue
        if off+hs>len(data): break
        h=data[off:off+hs]
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,bxw,blw=struct.unpack_from(">7f",h,42)
        if ver==4:
            nb,na,ns,head=struct.unpack_from(">4h",h,70)
            dscale,xscale=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,head=struct.unpack_from(">4i",h,70)
            dscale,xscale=struct.unpack_from(">2f",h,86)
        p=off+hs
        if min(nb,na,ns)<0 or nb>100000 or na>100000 or ns>1000000:
            raise RuntimeError(f"{name}: implausible dims {nb},{na},{ns}")
        flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p).copy(); p+=nb
        def arr(n):
            nonlocal p
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(float); p+=2*n; return a
        bath=arr(nb); across=arr(nb); along=arr(nb)
        p += 2*na + 6*ns
        gx,gy=gxy(lon,lat)
        hr=math.radians(heading); sh=math.sin(hr); ch=math.cos(hr)
        xt=xscale*across; lt=xscale*along
        east=gx + lt*sh + xt*ch
        north=gy + lt*ch - xt*sh
        depth=dscale*bath + sd
        good=np.where(flags==0)[0]
        if len(good):
            chunks.append(np.column_stack([east[good],north[good],depth[good]]))
        off=p; nrec+=1
    return np.vstack(chunks) if chunks else np.empty((0,3)), nrec

def select_files(survey):
    base=SURVEYS[survey]["base"]
    names=list_files(base,"fnv")
    selected={c["id"]:[] for c in CANDS}
    manifest=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as ex:
        futs=[ex.submit(fetch_fnv,base,n) for n in names]
        for fut in concurrent.futures.as_completed(futs):
            name,sha,rows=fut.result()
            rec={"fnv":name,"sha256":sha,"rows":len(rows),"candidate_min_cross_m":{}}
            for c in CANDS:
                d=fnv_min_cross_distance(rows,c)
                rec["candidate_min_cross_m"][c["id"]]=d
                if d<=SELECT_M:
                    selected[c["id"]].append(name[:-4]+".fbt")
            manifest.append(rec)
    for k in selected: selected[k]=sorted(set(selected[k]))
    manifest.sort(key=lambda x:x["fnv"])
    return selected,manifest

def load_selected_fbt(survey,selected):
    base=SURVEYS[survey]["base"]
    names=sorted(set(n for xs in selected.values() for n in xs))
    manifest=[]; parts=[]
    def one(n):
        b=fetch(urljoin(base,n),180)
        pts,nrec=parse_fbt(b,n)
        return n,hashlib.sha256(b).hexdigest(),len(b),pts,nrec
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
        futs=[ex.submit(one,n) for n in names]
        for fut in concurrent.futures.as_completed(futs):
            n,sha,size,pts,nrec=fut.result()
            parts.append((n,pts))
            manifest.append({"fbt":n,"sha256":sha,"bytes":size,"records":nrec,"good_beams":int(len(pts))})
    manifest.sort(key=lambda x:x["fbt"])
    return parts,manifest

def candidate_points(c, selected_names, parts, radius=1600):
    cx,cy=gxy(c["lon"],c["lat"])
    wanted=set(selected_names)
    out=[]
    for n,pts in parts:
        if n not in wanted or len(pts)==0: continue
        rr=np.hypot(pts[:,0]-cx,pts[:,1]-cy)
        m=rr<=radius
        if m.any(): out.append(pts[m])
    return np.vstack(out) if out else np.empty((0,3))

def basic_support(c,pts):
    cx,cy=gxy(c["lon"],c["lat"])
    if len(pts)==0:return {"nearest_good_beam_m":None,"n_1000m":0}
    rr=np.hypot(pts[:,0]-cx,pts[:,1]-cy)
    return {"nearest_good_beam_m":float(rr.min()),"n_1000m":int(np.sum(rr<=1000))}

def grid_medians(c,pts):
    cx,cy=gxy(c["lon"],c["lat"])
    x=pts[:,0]-cx; y=pts[:,1]-cy; z=pts[:,2]
    rr=np.hypot(x,y)
    m=rr<=GRID_R
    x=x[m];y=y[m];z=z[m]
    cells={}
    for xx,yy,zz in zip(x,y,z):
        ix=math.floor((xx+GRID_R)/CELL); iy=math.floor((yy+GRID_R)/CELL)
        cells.setdefault((ix,iy),[]).append(float(zz))
    out={}
    for (ix,iy),vals in cells.items():
        xc=-GRID_R+(ix+0.5)*CELL; yc=-GRID_R+(iy+0.5)*CELL
        if math.hypot(xc,yc)<=GRID_R:
            out[(ix,iy)]={"x":xc,"y":yc,"z":float(np.median(vals)),"n":len(vals)}
    return out

def fit_stats(x,y,z):
    A=np.column_stack([np.ones(len(x)),x,y])
    coef=np.linalg.lstsq(A,z,rcond=None)[0]
    pred=A@coef; res=z-pred
    rmse=float(np.sqrt(np.mean(res*res)))
    grad=np.array([coef[1],coef[2]],float)
    gm=float(np.hypot(*grad))
    az=(math.degrees(math.atan2(grad[0],grad[1]))+360)%360
    q=np.column_stack([np.ones(len(x)),x,y,x*x,x*y,y*y])
    qc=np.linalg.lstsq(q,z,rcond=None)[0]
    H=np.array([[2*qc[3],qc[4]],[qc[4],2*qc[5]]])
    eig=np.linalg.eigvalsh(H)
    r=np.hypot(x,y)
    R=np.column_stack([np.ones(len(r)),r,r*r])
    rc=np.linalg.lstsq(R,res,rcond=None)[0]
    rp=R@rc
    sst=float(np.sum((res-res.mean())**2)); ssr=float(np.sum((res-rp)**2))
    radial=None if sst<=0 else float(1-ssr/sst)
    p05,p95=np.percentile(z,[5,95])
    return {
      "coef":[float(v) for v in coef],
      "residuals":res,
      "plane_rmse_m":rmse,
      "plane_gradient_magnitude":gm,
      "plane_gradient_azimuth_deg":az,
      "p05_p95_relief_m":float(p95-p05),
      "quadratic_hessian_eigenvalues_per_m":[float(v) for v in eig],
      "radial_quadratic_r2_after_plane":radial
    }

def angle_diff(a,b):
    d=abs(a-b)%360
    return min(d,360-d)

def compare(c,pa,pb):
    sa=basic_support(c,pa); sb=basic_support(c,pb)
    support=(sa["nearest_good_beam_m"] is not None and sb["nearest_good_beam_m"] is not None and
             sa["nearest_good_beam_m"]<=100 and sb["nearest_good_beam_m"]<=100 and
             sa["n_1000m"]>=100 and sb["n_1000m"]>=100)
    ga=grid_medians(c,pa); gb=grid_medians(c,pb)
    common=sorted(set(ga)&set(gb))
    out={"support_A":sa,"support_B":sb,"support_gate_pass":support,
         "grid_cells_A":len(ga),"grid_cells_B":len(gb),"common_cell_count":len(common),
         "common_cell_gate_pass":len(common)>=30}
    if not support or len(common)<30:
        out["primary_reproduction_gate_pass"]=False
        out["verdict"]="NO_CROSS_SURVEY_MORPHOLOGY_TEST" if not support else "INSUFFICIENT_COMMON_GRID_SUPPORT"
        return out
    x=np.array([ga[k]["x"] for k in common],float)
    y=np.array([ga[k]["y"] for k in common],float)
    za=np.array([ga[k]["z"] for k in common],float)
    zb=np.array([gb[k]["z"] for k in common],float)
    fa=fit_stats(x,y,za); fb=fit_stats(x,y,zb)
    ra=fa.pop("residuals"); rb=fb.pop("residuals")
    pear=float(np.corrcoef(ra,rb)[0,1]) if np.std(ra)>0 and np.std(rb)>0 else None
    spear=float(spearmanr(ra,rb).statistic) if len(ra)>=3 else None
    sign=float(np.mean(np.sign(ra)==np.sign(rb)))
    primary=bool(pear is not None and spear is not None and pear>=0.60 and spear>=0.60 and sign>=0.65)
    medoff=float(np.median(zb-za))
    relief_ratio=None if fa["p05_p95_relief_m"]==0 else fb["p05_p95_relief_m"]/fa["p05_p95_relief_m"]
    grad_ratio=None if fa["plane_gradient_magnitude"]==0 else fb["plane_gradient_magnitude"]/fa["plane_gradient_magnitude"]
    rmse_ratio=None if fa["plane_rmse_m"]==0 else fb["plane_rmse_m"]/fa["plane_rmse_m"]
    out.update({
      "lineage_A_common_grid_stats":fa,
      "lineage_B_common_grid_stats":fb,
      "residual_pearson":pear,
      "residual_spearman":spear,
      "residual_sign_concordance":sign,
      "median_depth_offset_B_minus_A_m":medoff,
      "relief_ratio_B_over_A":relief_ratio,
      "plane_gradient_magnitude_ratio_B_over_A":grad_ratio,
      "plane_gradient_azimuth_difference_deg":angle_diff(fa["plane_gradient_azimuth_deg"],fb["plane_gradient_azimuth_deg"]),
      "plane_rmse_ratio_B_over_A":rmse_ratio,
      "primary_reproduction_gate_pass":primary,
      "verdict":"LOCAL_BATHYMETRIC_MORPHOLOGY_REPRODUCED_ACROSS_INDEPENDENT_SURVEYS" if primary else "NO_CROSS_SURVEY_MORPHOLOGY_REPLICATION_AT_FROZEN_CENTER"
    })
    return out

inventory={}
parts={}
for survey in SURVEYS:
    print("SELECT",survey)
    selected,fnvman=select_files(survey)
    print(survey,"selected",json.dumps({k:len(v) for k,v in selected.items()}))
    p,fbtman=load_selected_fbt(survey,selected)
    inventory[survey]={"selected_fbt_by_candidate":selected,"fnv_manifest":fnvman,"fbt_manifest":fbtman}
    parts[survey]=p
    print(survey,"loaded FBT",len(p))

results=[]
for c in CANDS:
    pa=candidate_points(c,inventory["KN192-07"]["selected_fbt_by_candidate"][c["id"]],parts["KN192-07"])
    pb=candidate_points(c,inventory["KNOX15RR"]["selected_fbt_by_candidate"][c["id"]],parts["KNOX15RR"])
    comp=compare(c,pa,pb)
    comp.update(c)
    comp["point_count_within_1600m_A"]=int(len(pa)); comp["point_count_within_1600m_B"]=int(len(pb))
    results.append(comp)
    print("\nRESULT",c["id"])
    print(json.dumps(comp,indent=2))

out={
 "artifact_id":"JANUS-KUSTO-KN19207-VS-KNOX15RR-CROSS-SURVEY-BEAM-REPLICATION-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-KN19207-VS-KNOX15RR-OPERATIONAL-PREREG-2026-09-22-v1.0.json",
 "lineages":SURVEYS,
 "inventory":inventory,
 "results":results,
 "pass_candidates":[r["id"] for r in results if r.get("primary_reproduction_gate_pass")],
 "claim_ceiling":"CROSS_SURVEY_SEAFLOOR_REPRODUCIBILITY_ONLY__NO_NOVELTY__NO_IDENTITY__NO_ARTIFICIALITY"
}
p=OUT/"JANUS-KUSTO-KN19207-VS-KNOX15RR-CROSS-SURVEY-BEAM-REPLICATION-2026-09-22-v1.0.json"
p.write_text(json.dumps(out,indent=2),encoding="utf-8")
print("\nSUMMARY",json.dumps({"pass_candidates":out["pass_candidates"]},indent=2))
