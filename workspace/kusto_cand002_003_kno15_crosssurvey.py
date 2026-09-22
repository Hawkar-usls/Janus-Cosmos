#!/usr/bin/env python3
import concurrent.futures, hashlib, json, math, re, struct
from pathlib import Path
from urllib.parse import urljoin
import numpy as np
import requests
from scipy.stats import spearmanr

OUT=Path("workspace/kusto_crosssurvey_out"); OUT.mkdir(parents=True,exist_ok=True)
SURVEYS={
 "KN192-07":{
   "generated":"https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/generated/",
   "raw":"https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/"
 },
 "KNOX15RR":{
   "generated":"https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/data/version1/MB/generated/",
   "raw":"https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/data/version1/MB/"
 }
}
CANDS=[
 {"id":"KN19207_CAND_002","lat":-3.9727527956056825,"lon":-12.272824298723462},
 {"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432},
]
RADS=[500,1000]
GRID=100.0
SELECT_PAD=2500.0
UA={"User-Agent":"JANUS-KUSTO-crosssurvey-replication/1.0"}
S=requests.Session(); S.headers.update(UA)

def fetch(url,timeout=180):
    r=S.get(url,timeout=timeout); r.raise_for_status(); return r.content

def listing(url,suffix):
    txt=fetch(url,90).decode("utf-8","replace")
    return sorted(set(re.findall(r'href="([^"]+'+re.escape(suffix)+r')"',txt,re.I)))

def localxy(lat0,lon0,lat,lon):
    return ((lon-lon0)*111320.0*math.cos(math.radians(lat0)),(lat-lat0)*111320.0)

def segdist0(a,b):
    ax,ay=a; bx,by=b; vx,vy=bx-ax,by-ay
    den=vx*vx+vy*vy
    t=0.0 if den==0 else max(0.0,min(1.0,-(ax*vx+ay*vy)/den))
    return math.hypot(ax+t*vx,ay+t*vy)

def parse_fnv(txt,c):
    rows=[]
    for line in txt.splitlines():
        p=line.split()
        if len(p)<19: continue
        try:
            rows.append({
              "epoch":float(p[6]),"navlon":float(p[7]),"navlat":float(p[8]),"heading":float(p[9]),
              "portlon":float(p[15]),"portlat":float(p[16]),"stbdlon":float(p[17]),"stbdlat":float(p[18])
            })
        except: pass
    best=1e99
    for r in rows:
        a=localxy(c["lat"],c["lon"],r["portlat"],r["portlon"])
        b=localxy(c["lat"],c["lon"],r["stbdlat"],r["stbdlon"])
        best=min(best,segdist0(a,b))
    return rows,best

def choose_files(sid,c):
    g=SURVEYS[sid]["generated"]
    names=listing(g,".fnv")
    selected=[]
    best=[]
    def one(n):
        b=fetch(urljoin(g,n),90)
        rows,d=parse_fnv(b.decode("utf-8","replace"),c)
        return n,d,hashlib.sha256(b).hexdigest(),len(rows)
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as ex:
        for n,d,h,nrows in ex.map(one,names):
            best.append((d,n))
            if d<=SELECT_PAD:
                selected.append({"fnv":n,"min_cross_swath_distance_m":d,"fnv_sha256":h,"rows":nrows})
    best.sort()
    selected.sort(key=lambda z:z["fnv"])
    return selected,best[:10],len(names)

def parse_fbt(data,c):
    off=0; pts=[]; records=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids: raise RuntimeError(f"bad fbt tag {tag!r} at {off}")
        ver,hs=ids[tag]
        if isinstance(ver,str):
            off += hs+128; continue
        if off+hs>len(data): break
        h=data[off:off+hs]
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,bx,bl=struct.unpack_from(">7f",h,42)
        if ver==4:
            nb,na,ns,shd=struct.unpack_from(">4h",h,70); dscale,xscale=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,shd=struct.unpack_from(">4i",h,70); dscale,xscale=struct.unpack_from(">2f",h,86)
        p=off+hs
        flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p).copy(); p+=nb
        def arr(n):
            nonlocal p
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(float); p+=2*n; return a
        bath=arr(nb); across=arr(nb); along=arr(nb)
        p += 2*na + 6*ns
        navx,navy=localxy(c["lat"],c["lon"],lat,lon)
        hr=math.radians(heading); sh=math.sin(hr); ch=math.cos(hr)
        xtrack=xscale*across; ltrack=xscale*along
        east=navx+ltrack*sh+xtrack*ch
        north=navy+ltrack*ch-xtrack*sh
        depth=dscale*bath+sd
        good=np.where(flags==0)[0]
        if len(good):
            rr=np.hypot(east[good],north[good])
            keep=good[rr<=3000]
            if len(keep):
                pts.append(np.column_stack([east[keep],north[keep],depth[keep]]))
        off=p; records+=1
    return np.vstack(pts) if pts else np.empty((0,3)),records

def fbt_name_from_fnv(n):
    return n[:-4]+".fbt"

def raw_name_from_fbt(n):
    return n[:-4]+".gz"

def load_points(sid,c,selected):
    g=SURVEYS[sid]["generated"]; raw=SURVEYS[sid]["raw"]
    parts=[]; manifest=[]
    for item in selected:
        fn=fbt_name_from_fnv(item["fnv"])
        u=urljoin(g,fn)
        b=fetch(u,180)
        pts,nrec=parse_fbt(b,c)
        parts.append(pts)
        rawname=raw_name_from_fbt(fn)
        rawurl=urljoin(raw,rawname)
        rawmeta={"raw_file":rawname,"raw_url":rawurl}
        try:
            rb=fetch(rawurl,180)
            rawmeta.update({"raw_bytes":len(rb),"raw_sha256":hashlib.sha256(rb).hexdigest()})
        except Exception as e:
            rawmeta.update({"raw_error":repr(e)})
        manifest.append({
          **item,"fbt":fn,"fbt_url":u,"fbt_bytes":len(b),"fbt_sha256":hashlib.sha256(b).hexdigest(),
          "fbt_records":nrec,"points_within_3km":int(len(pts)),**rawmeta
        })
    return (np.vstack(parts) if parts else np.empty((0,3))),manifest

def make_grid(points,rad):
    grid={}
    if len(points)==0:return grid
    for x,y,z in points:
        if x*x+y*y>rad*rad: continue
        ix=int(math.floor((x+rad)/GRID)); iy=int(math.floor((y+rad)/GRID))
        cx=-rad+(ix+0.5)*GRID; cy=-rad+(iy+0.5)*GRID
        if cx*cx+cy*cy>rad*rad: continue
        grid.setdefault((ix,iy),[]).append(float(z))
    return {k:float(np.median(v)) for k,v in grid.items()}

def plane_residual(keys,grid):
    xy=[]; z=[]
    for ix,iy in keys:
        x=ix*GRID; y=iy*GRID
        xy.append((x,y)); z.append(grid[(ix,iy)])
    xy=np.asarray(xy,float); z=np.asarray(z,float)
    A=np.column_stack([np.ones(len(z)),xy[:,0],xy[:,1]])
    coef=np.linalg.lstsq(A,z,rcond=None)[0]
    return z-A@coef

def centered_compare(g1,g2,min_common):
    keys=sorted(set(g1)&set(g2))
    if len(keys)<min_common:return {"common_cell_count":len(keys),"sufficient":False}
    r1=plane_residual(keys,g1); r2=plane_residual(keys,g2)
    pear=float(np.corrcoef(r1,r2)[0,1]) if np.std(r1)>0 and np.std(r2)>0 else None
    spr=float(spearmanr(r1,r2).statistic) if len(r1)>=3 else None
    return {"common_cell_count":len(keys),"sufficient":True,"pearson_r":pear,"spearman_rho":spr}

def shifted_compare(g1,g2,min_common):
    best=None
    for sx in range(-3,4):
      for sy in range(-3,4):
        # compare g1[k] with g2[k shifted back]
        pairs=[]
        for k in g1:
            k2=(k[0]+sx,k[1]+sy)
            if k2 in g2:pairs.append((k,k2))
        if len(pairs)<min_common:continue
        # independent plane fits on matched coordinates in their own grids
        k1=[a for a,b in pairs]; k2=[b for a,b in pairs]
        r1=plane_residual(k1,g1); r2=plane_residual(k2,g2)
        if np.std(r1)==0 or np.std(r2)==0:continue
        pear=float(np.corrcoef(r1,r2)[0,1])
        rec={"sx_cells":sx,"sy_cells":sy,"shift_x_m":sx*GRID,"shift_y_m":sy*GRID,
             "shift_distance_m":math.hypot(sx*GRID,sy*GRID),"pearson_r":pear,"common_cell_count":len(pairs)}
        if best is None or pear>best["pearson_r"]:best=rec
    return best

def grid_metrics(g,rad):
    if not g:return None
    keys=list(g); vals=np.array([g[k] for k in keys],float)
    xy=np.array([[(-rad+(k[0]+.5)*GRID),(-rad+(k[1]+.5)*GRID)] for k in keys],float)
    A=np.column_stack([np.ones(len(vals)),xy[:,0],xy[:,1]])
    coef=np.linalg.lstsq(A,vals,rcond=None)[0]; res=vals-A@coef
    rr=np.hypot(xy[:,0],xy[:,1])
    R=np.column_stack([np.ones(len(vals)),rr,rr*rr])
    rc=np.linalg.lstsq(R,res,rcond=None)[0]; pred=R@rc
    sst=float(np.sum((res-res.mean())**2)); ssr=float(np.sum((res-pred)**2))
    radial=None if sst<=0 else 1-ssr/sst
    ann=(rr>=250)&(rr<=1000)
    m4=None
    if ann.sum()>=20:
        theta=np.arctan2(xy[ann,1],xy[ann,0]); a=2*abs(np.mean(res[ann]*np.exp(-4j*theta))); rms=np.sqrt(np.mean(res[ann]**2))
        if rms>0:m4=float(a/rms)
    p05,p95=np.percentile(vals,[5,95])
    return {"cell_count":len(vals),"p05_p95_relief_m":float(p95-p05),"plane_rmse_m":float(np.sqrt(np.mean(res*res))),
            "radial_quadratic_r2":float(radial) if radial is not None else None,"angular_m4":m4}

out={"artifact_id":"JANUS-KUSTO-CAND002-003-KN19207-X-KNOX15RR-CROSSSURVEY-RUN-2026-09-22-v1.0",
     "prereg":"data/cousteau/JANUS-KUSTO-CAND002-003-KN19207-X-KNOX15RR-CROSSSURVEY-PREREG-2026-09-22-v1.0.json",
     "results":[]}
for c in CANDS:
    print("\nCANDIDATE",c["id"])
    survey_data={}
    for sid in SURVEYS:
        selected,nearest,nfnv=choose_files(sid,c)
        print(sid,"FNV",nfnv,"selected",len(selected),"nearest",nearest[:3])
        pts,manifest=load_points(sid,c,selected)
        nearest_beam=float(np.min(np.hypot(pts[:,0],pts[:,1]))) if len(pts) else None
        survey_data[sid]={"selected_files":manifest,"point_count_within_3km":int(len(pts)),"nearest_good_beam_m":nearest_beam,"points":pts}
        print(sid,"points",len(pts),"nearest beam",nearest_beam)
    cres={"candidate":c,"surveys":{},"radii":{}}
    for sid,d in survey_data.items():
        cres["surveys"][sid]={k:v for k,v in d.items() if k!="points"}
    for rad in RADS:
        g1=make_grid(survey_data["KN192-07"]["points"],rad)
        g2=make_grid(survey_data["KNOX15RR"]["points"],rad)
        minc=20 if rad==500 else 80
        cen=centered_compare(g1,g2,minc)
        sh=shifted_compare(g1,g2,minc)
        gm1=grid_metrics(g1,rad); gm2=grid_metrics(g2,rad)
        cres["radii"][str(rad)]={"KN192-07_grid_metrics":gm1,"KNOX15RR_grid_metrics":gm2,
                                 "centered":cen,"best_shift":sh}
    p=cres["radii"]["1000"]; cen=p["centered"]; sh=p["best_shift"]
    n1=cres["surveys"]["KN192-07"]["nearest_good_beam_m"]; n2=cres["surveys"]["KNOX15RR"]["nearest_good_beam_m"]
    if not cen.get("sufficient",False) or n1 is None or n2 is None:
        verdict="INSUFFICIENT_COMMON_COVERAGE"
    else:
        req=[
          n1<=150,n2<=150,cen["common_cell_count"]>=80,
          cen["pearson_r"] is not None and cen["pearson_r"]>=0.60,
          cen["spearman_rho"] is not None and cen["spearman_rho"]>=0.55,
          sh is not None and sh["shift_distance_m"]<=200,
          sh is not None and cen["pearson_r"]>=sh["pearson_r"]-0.10
        ]
        verdict="PASS_CROSS_SURVEY_MORPHOLOGY_REPLICATION" if all(req) else "FAIL_CROSS_SURVEY_MORPHOLOGY_REPLICATION"
        cres["gate_requirements"]=req
    cres["verdict"]=verdict
    out["results"].append(cres)
    print(json.dumps({"id":c["id"],"verdict":verdict,"surveys":cres["surveys"],"radii":cres["radii"]},indent=2)[:30000])

out["summary"]={r["candidate"]["id"]:r["verdict"] for r in out["results"]}
out["claim_ceiling"]="CROSS_SURVEY_SEAFLOOR_MORPHOLOGY_REPLICATION_ONLY__NO_OBJECT_IDENTITY__NO_ARTIFICIALITY"
p=OUT/"JANUS-KUSTO-CAND002-003-KN19207-X-KNOX15RR-CROSSSURVEY-RUN-2026-09-22-v1.0.json"
p.write_text(json.dumps(out,indent=2),encoding="utf-8")
print("\nSUMMARY",json.dumps(out["summary"],indent=2))
