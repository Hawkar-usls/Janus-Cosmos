#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures, ftplib, hashlib, json, math, re, struct
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import requests
from scipy.stats import spearmanr

C={"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432}
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
M=111320.0; LAT0=-4.0; COS0=math.cos(math.radians(LAT0))
GRID_R=1000.0; SELECT_R=1600.0; CELL=100.0
UA={"User-Agent":"JANUS-KUSTO-CAND003-three-year-replication/1.0"}

CD_HOST='livftp.noc.ac.uk'
CD_ROOT='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/EM12'
SURVEYS={
 "KN192-07_2008":"https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/generated/",
 "KNOX15RR_2008":"https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/data/version1/MB/generated/"
}

def gxy(lon,lat): return lon*M*COS0,lat*M
def cxy(): return gxy(C["lon"],C["lat"])

def ftp():
    f=ftplib.FTP(timeout=120);f.connect(CD_HOST,21);f.login('anonymous','janus-kusto@example.invalid');f.voidcmd('TYPE I');return f

def parse_xyz(raw):
    s=raw.decode('ascii','ignore').strip()
    if not s or s[0] in '#;!':return None
    p=s.replace(',',' ').split()
    if len(p)<3:return None
    try:lon=float(p[0]);lat=float(p[1]);z=float(p[2])
    except:return None
    if not all(map(math.isfinite,(lon,lat,z))):return None
    return lon,lat,z

def list_cd():
    f=ftp()
    try:
        out=[]
        for p in f.nlst(CD_ROOT):
            name=p.rstrip('/').split('/')[-1]
            if name.lower().endswith('.xyz.ascii'):
                full=p if p.startswith('/') else CD_ROOT+'/'+name
                try:size=f.size(full)
                except:size=None
                out.append((name,full,size))
        return sorted(out)
    finally:
        try:f.quit()
        except:f.close()

def load_cd169():
    cx,cy=cxy(); selected=[]; pts=[]; manifest=[]
    # Geometry selection is performed using all files, independent of depth values.
    file_rows=[]
    for name,path,size in list_cd():
        f=ftp();h=hashlib.sha256();nearest=float('inf');xy=[];zrows=[]
        try:
            sock=f.transfercmd('RETR '+path);stream=sock.makefile('rb')
            try:
                for raw in stream:
                    h.update(raw);q=parse_xyz(raw)
                    if q is None:continue
                    lon,lat,z=q;x,y=gxy(lon,lat);d=math.hypot(x-cx,y-cy)
                    xy.append(d);zrows.append((x,y,z,d))
                    nearest=min(nearest,d)
            finally:stream.close();sock.close()
        finally:f.close()
        rec={"file":name,"path":path,"bytes":size,"sha256":h.hexdigest(),"nearest_m":nearest}
        file_rows.append(rec)
        if nearest<=SELECT_R:
            selected.append(name)
            pts.extend((x,y,z) for x,y,z,d in zrows if d<=SELECT_R)
    file_rows.sort(key=lambda r:r["nearest_m"])
    return np.asarray(pts,float),{"files":file_rows,"selected_files":selected}

def fetch(url,timeout=120):
    r=requests.get(url,headers=UA,timeout=timeout);r.raise_for_status();return r.content

def list_files(base,ext):
    html=fetch(base,60).decode("latin1","replace")
    hrefs=re.findall(r'href="([^"]+)"',html,re.I)
    suffix='.'+ext.lower()
    return sorted(set(h for h in hrefs if h.lower().endswith(suffix)))

def parse_fnv(text):
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try:rows.append({"portlon":float(p[15]),"portlat":float(p[16]),"stbdlon":float(p[17]),"stbdlat":float(p[18])})
        except:pass
    return rows

def segdist(px,py,ax,ay,bx,by):
    vx,vy=bx-ax,by-ay; wx,wy=px-ax,py-ay; den=vx*vx+vy*vy
    t=0 if den==0 else max(0,min(1,(wx*vx+wy*vy)/den))
    return math.hypot(px-(ax+t*vx),py-(ay+t*vy))

def min_fnv(rows):
    px,py=cxy();best=float('inf')
    for r in rows:
        ax,ay=gxy(r["portlon"],r["portlat"]);bx,by=gxy(r["stbdlon"],r["stbdlat"])
        best=min(best,segdist(px,py,ax,ay,bx,by))
    return best

def parse_fbt(data,name):
    off=0; chunks=[]; nrec=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids:raise RuntimeError(f"{name} bad tag {tag!r} @{off}")
        ver,hs=ids[tag]
        if isinstance(ver,str):off+=hs+128;continue
        h=data[off:off+hs]
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,bxw,blw=struct.unpack_from(">7f",h,42)
        if ver==4:
            nb,na,ns,head=struct.unpack_from(">4h",h,70);dscale,xscale=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,head=struct.unpack_from(">4i",h,70);dscale,xscale=struct.unpack_from(">2f",h,86)
        p=off+hs
        flags=np.frombuffer(data,dtype=np.uint8,count=nb,offset=p).copy();p+=nb
        def arr(n):
            nonlocal p
            a=np.frombuffer(data,dtype=">i2",count=n,offset=p).astype(float);p+=2*n;return a
        bath=arr(nb);across=arr(nb);along=arr(nb);p+=2*na+6*ns
        gx,gy=gxy(lon,lat);hr=math.radians(heading);sh=math.sin(hr);ch=math.cos(hr)
        xt=xscale*across;lt=xscale*along
        x=gx+lt*sh+xt*ch;y=gy+lt*ch-xt*sh;z=dscale*bath+sd
        good=np.where(flags==0)[0]
        if len(good):chunks.append(np.column_stack([x[good],y[good],z[good]]))
        off=p;nrec+=1
    return (np.vstack(chunks) if chunks else np.empty((0,3))),nrec

def load_2008(label,base):
    names=list_files(base,"fnv"); selected=[]; fnv_manifest=[]
    def getfnv(n):
        b=fetch(urljoin(base,n),90);return n,hashlib.sha256(b).hexdigest(),parse_fnv(b.decode("utf-8","replace"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as ex:
        for n,sha,rows in (f.result() for f in concurrent.futures.as_completed([ex.submit(getfnv,n) for n in names])):
            d=min_fnv(rows);fnv_manifest.append({"file":n,"sha256":sha,"min_cross_m":d})
            if d<=SELECT_R:selected.append(n[:-4]+".fbt")
    selected=sorted(set(selected));fnv_manifest.sort(key=lambda r:r["file"])
    parts=[];fbt_manifest=[]
    def getfbt(n):
        b=fetch(urljoin(base,n),180);pts,nrec=parse_fbt(b,n);return n,hashlib.sha256(b).hexdigest(),len(b),pts,nrec
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
        for n,sha,size,pts,nrec in (f.result() for f in concurrent.futures.as_completed([ex.submit(getfbt,n) for n in selected])):
            parts.append(pts);fbt_manifest.append({"file":n,"sha256":sha,"bytes":size,"records":nrec,"good_beams":int(len(pts))})
    allpts=np.vstack(parts) if parts else np.empty((0,3))
    cx,cy=cxy()
    rr=np.hypot(allpts[:,0]-cx,allpts[:,1]-cy) if len(allpts) else np.array([])
    allpts=allpts[rr<=SELECT_R] if len(allpts) else allpts
    return allpts,{"selected_fbt":selected,"fbt_manifest":sorted(fbt_manifest,key=lambda r:r["file"]),"fnv_manifest":fnv_manifest}

def support(pts):
    cx,cy=cxy()
    if not len(pts):return {"nearest_m":None,"n1000":0}
    rr=np.hypot(pts[:,0]-cx,pts[:,1]-cy)
    return {"nearest_m":float(rr.min()),"n1000":int(np.sum(rr<=1000))}

def grid(pts):
    cx,cy=cxy(); x=pts[:,0]-cx;y=pts[:,1]-cy;z=pts[:,2];r=np.hypot(x,y);m=r<=GRID_R;x=x[m];y=y[m];z=z[m]
    d={}
    for xx,yy,zz in zip(x,y,z):
        ix=math.floor((xx+GRID_R)/CELL);iy=math.floor((yy+GRID_R)/CELL)
        d.setdefault((ix,iy),[]).append(float(zz))
    out={}
    for k,v in d.items():
        ix,iy=k;xc=-GRID_R+(ix+.5)*CELL;yc=-GRID_R+(iy+.5)*CELL
        if math.hypot(xc,yc)<=GRID_R:out[k]={"x":xc,"y":yc,"z":float(np.median(v)),"n":len(v)}
    return out

def fit(x,y,z):
    A=np.column_stack([np.ones(len(x)),x,y]);coef=np.linalg.lstsq(A,z,rcond=None)[0];res=z-A@coef
    return coef,res

def pair(labelA,ptsA,labelB,ptsB):
    sa,sb=support(ptsA),support(ptsB);ga,gb=grid(ptsA),grid(ptsB);common=sorted(set(ga)&set(gb))
    out={"A":labelA,"B":labelB,"support_A":sa,"support_B":sb,"common_cells":len(common)}
    if sa["nearest_m"] is None or sb["nearest_m"] is None or sa["nearest_m"]>100 or sb["nearest_m"]>100 or sa["n1000"]<100 or sb["n1000"]<100 or len(common)<30:
        out.update({"gate_pass":False,"verdict":"INSUFFICIENT_SUPPORT_OR_COMMON_CELLS"});return out
    x=np.array([ga[k]["x"] for k in common]);y=np.array([ga[k]["y"] for k in common]);za=np.array([ga[k]["z"] for k in common]);zb=np.array([gb[k]["z"] for k in common])
    ca,ra=fit(x,y,za);cb,rb=fit(x,y,zb)
    pear=float(np.corrcoef(ra,rb)[0,1]);spear=float(spearmanr(ra,rb).statistic);sign=float(np.mean(np.sign(ra)==np.sign(rb)))
    p05a,p95a=np.percentile(za,[5,95]);p05b,p95b=np.percentile(zb,[5,95])
    aza=(math.degrees(math.atan2(ca[1],ca[2]))+360)%360;azb=(math.degrees(math.atan2(cb[1],cb[2]))+360)%360;ad=abs(aza-azb)%360;ad=min(ad,360-ad)
    gate=pear>=.60 and spear>=.60 and sign>=.65
    out.update({
      "gate_pass":bool(gate),"residual_pearson":pear,"residual_spearman":spear,"residual_sign_concordance":sign,
      "median_depth_offset_B_minus_A_m":float(np.median(zb-za)),
      "relief_A_m":float(p95a-p05a),"relief_B_m":float(p95b-p05b),
      "plane_gradient_A":[float(ca[1]),float(ca[2])],"plane_gradient_B":[float(cb[1]),float(cb[2])],
      "plane_gradient_azimuth_A_deg":aza,"plane_gradient_azimuth_B_deg":azb,"plane_gradient_azimuth_difference_deg":ad,
      "verdict":"PASS_PERSISTENT_LOCAL_MORPHOLOGY" if gate else "FAIL_PERSISTENT_LOCAL_MORPHOLOGY"
    })
    return out

print("Loading CD169 EM12...")
cd,cdman=load_cd169();print("CD points",len(cd),"support",support(cd))
loaded={"CD169_EM12_2005":cd};manifest={"CD169_EM12_2005":cdman}
for lab,base in SURVEYS.items():
    print("Loading",lab);pts,man=load_2008(lab,base);loaded[lab]=pts;manifest[lab]=man;print(lab,len(pts),support(pts))

pairs=[
 pair("CD169_EM12_2005",loaded["CD169_EM12_2005"],"KN192-07_2008",loaded["KN192-07_2008"]),
 pair("CD169_EM12_2005",loaded["CD169_EM12_2005"],"KNOX15RR_2008",loaded["KNOX15RR_2008"]),
 pair("KN192-07_2008",loaded["KN192-07_2008"],"KNOX15RR_2008",loaded["KNOX15RR_2008"])
]
three_year_any=any(p["gate_pass"] for p in pairs[:2]);triple=all(p["gate_pass"] for p in pairs)
out={
 "artifact_id":"JANUS-KUSTO-CAND003-THREE-YEAR-TRIPLE-LINEAGE-RUN-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-THREE-YEAR-TRIPLE-LINEAGE-PREREG-2026-09-22-v1.0.json",
 "target":C,
 "lineage_support":{k:support(v) for k,v in loaded.items()},
 "pairwise_results":pairs,
 "three_year_pass_at_least_one_2005_to_2008":three_year_any,
 "triple_lineage_all_pairs_pass":triple,
 "manifest":manifest,
 "claim_ceiling":"THREE_YEAR_PERSISTENT_LOCAL_SEAFLOOR_MORPHOLOGY_REPRODUCIBILITY_ONLY__NO_HYDROTHERMAL_IDENTITY__NO_ARTIFICIALITY"
}
p=OUT/"JANUS-KUSTO-CAND003-THREE-YEAR-TRIPLE-LINEAGE-RUN-2026-09-22-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps({"lineage_support":out["lineage_support"],"pairwise_results":pairs,"three_year_pass":three_year_any,"triple_lineage_all_pairs_pass":triple},indent=2))


# --- Support-matched posthoc preregistered sensitivity ---
from scipy.spatial import cKDTree

def thin_reference_2005(pts):
    cx,cy=cxy()
    rr=np.hypot(pts[:,0]-cx,pts[:,1]-cy)
    p=pts[rr<=1000].copy()
    cells={}
    for row in p:
        dx=row[0]-cx; dy=row[1]-cy
        ix=math.floor((dx+1000.0)/100.0)
        iy=math.floor((dy+1000.0)/100.0)
        xc=-1000.0+(ix+0.5)*100.0
        yc=-1000.0+(iy+0.5)*100.0
        key=(ix,iy)
        geom_d2=(dx-xc)**2+(dy-yc)**2
        cand=(geom_d2,float(row[0]),float(row[1]),float(row[2]))
        if key not in cells or cand[:3] < cells[key][:3]:
            cells[key]=cand
    rows=[]
    for key,cand in sorted(cells.items()):
        _,x,y,z=cand
        rows.append((x,y,z,key[0],key[1]))
    return np.asarray(rows,float)

def support_match(reference_thin, target_pts, estimator_radius_m):
    if len(reference_thin)==0 or len(target_pts)==0:
        return np.empty((0,5),float),np.empty((0,),float)
    tree=cKDTree(target_pts[:,:2])
    matched_ref=[]
    matched_target=[]
    for row in reference_thin:
        idx=tree.query_ball_point(row[:2], estimator_radius_m)
        if not idx:
            continue
        zt=float(np.median(target_pts[np.asarray(idx,dtype=int),2]))
        matched_ref.append(row)
        matched_target.append(zt)
    return np.asarray(matched_ref,float),np.asarray(matched_target,float)

def compare_support_matched(label, reference_thin, target_pts, estimator_radius_m):
    ref,zt=support_match(reference_thin,target_pts,estimator_radius_m)
    out={"target_lineage":label,"estimator_radius_m":estimator_radius_m,"matched_point_count":int(len(ref))}
    if len(ref)<50:
        out.update({"gate_pass":False,"verdict":"INSUFFICIENT_MATCHED_SUPPORT"})
        return out
    cx,cy=cxy()
    x=ref[:,0]-cx; y=ref[:,1]-cy
    zr=ref[:,2]
    _,rr=fit(x,y,zr)
    _,rt=fit(x,y,zt)
    pear=float(np.corrcoef(rr,rt)[0,1]) if np.std(rr)>0 and np.std(rt)>0 else None
    spear=float(spearmanr(rr,rt).statistic) if len(rr)>=3 else None
    sign=float(np.mean(np.sign(rr)==np.sign(rt)))
    p05r,p95r=np.percentile(zr,[5,95])
    p05t,p95t=np.percentile(zt,[5,95])
    relief_r=float(p95r-p05r); relief_t=float(p95t-p05t)
    gate=bool(pear is not None and spear is not None and pear>=0.60 and spear>=0.60 and sign>=0.65)
    out.update({
      "gate_pass":gate,
      "residual_pearson":pear,
      "residual_spearman":spear,
      "residual_sign_concordance":sign,
      "median_depth_offset_target_minus_2005_m":float(np.median(zt-zr)),
      "relief_2005_m":relief_r,
      "relief_target_m":relief_t,
      "relief_ratio_target_over_2005":None if relief_r==0 else relief_t/relief_r,
      "verdict":"SUPPORT_MATCHED_TEMPORAL_MORPHOLOGY_RECOVERY" if gate else "SUPPORT_MATCHED_TEMPORAL_MORPHOLOGY_STILL_FAILS"
    })
    return out

reference_thin=thin_reference_2005(loaded["CD169_EM12_2005"])
sens={}
for lab in ["KN192-07_2008","KNOX15RR_2008"]:
    sens[lab]={}
    for rad in [50,100,150]:
        sens[lab][str(rad)]=compare_support_matched(lab,reference_thin,loaded[lab],rad)

primary_pass_any=any(sens[lab]["100"].get("gate_pass",False) for lab in sens)
primary_pass_both=all(sens[lab]["100"].get("gate_pass",False) for lab in sens)
sout={
 "artifact_id":"JANUS-KUSTO-CAND003-SUPPORT-MATCHED-CROSSYEAR-SENSITIVITY-RUN-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-SUPPORT-MATCHED-CROSSYEAR-SENSITIVITY-PREREG-2026-09-22-v1.0.json",
 "target":C,
 "reference_2005_raw_points_within_1000m":int(np.sum(np.hypot(loaded["CD169_EM12_2005"][:,0]-cxy()[0],loaded["CD169_EM12_2005"][:,1]-cxy()[1])<=1000)),
 "reference_2005_thinned_support_points":int(len(reference_thin)),
 "selection_uses_depth":False,
 "results":sens,
 "primary_100m_recovery_any_2008_lineage":primary_pass_any,
 "primary_100m_recovery_both_2008_lineages":primary_pass_both,
 "strict_parent_fail_immutable":True,
 "claim_ceiling":"POSTHOC_PREREGISTERED_SUPPORT_MATCHED_SENSITIVITY_ONLY__STRICT_TEMPORAL_FAIL_UNCHANGED"
}
sp=OUT/"JANUS-KUSTO-CAND003-SUPPORT-MATCHED-CROSSYEAR-SENSITIVITY-RUN-2026-09-22-v1.0.json"
sp.write_text(json.dumps(sout,indent=2))
print("\nSUPPORT_MATCHED_SENSITIVITY")
print(json.dumps({
 "reference_2005_raw_points_within_1000m":sout["reference_2005_raw_points_within_1000m"],
 "reference_2005_thinned_support_points":sout["reference_2005_thinned_support_points"],
 "results":sens,
 "primary_100m_recovery_any_2008_lineage":primary_pass_any,
 "primary_100m_recovery_both_2008_lineages":primary_pass_both
},indent=2))


# --- preregistered processed-product demix diagnostic ---
PRODUCTS={
 "B1-81-1_Acceptl28-33.xyz.ascii":"0e8cac049b53dd6ba0e5013e8aaea2bfb48347db887edc3c4e3e9f697e96dcf0",
 "B1-81-1_Acceptl15-22.xyz.ascii":"ad2e4d60401f7ffaa5e14b3e43aa0e45da9d320c6b34fa7f3b873f4024be6007"
}

def load_cd_product(name,expected_sha):
    path=CD_ROOT+"/"+name
    f=ftp();h=hashlib.sha256();rows=[]
    cx,cy=cxy()
    try:
        sock=f.transfercmd("RETR "+path);stream=sock.makefile("rb")
        try:
            for raw in stream:
                h.update(raw);q=parse_xyz(raw)
                if q is None:continue
                lon,lat,z=q;x,y=gxy(lon,lat)
                if math.hypot(x-cx,y-cy)<=1000:
                    rows.append((x,y,z))
        finally:
            stream.close();sock.close()
    finally:
        f.close()
    sha=h.hexdigest()
    if sha!=expected_sha:
        raise RuntimeError(f"{name}: SHA mismatch {sha} != {expected_sha}")
    return np.asarray(rows,float),sha

def compare_product_support(label,reference_thin,target_pts):
    ref,zt=support_match(reference_thin,target_pts,100)
    out={"target_lineage":label,"estimator_radius_m":100,"matched_point_count":int(len(ref))}
    if len(ref)<30:
        out.update({"gate_pass":False,"verdict":"INSUFFICIENT_PRODUCT_SPECIFIC_SUPPORT"})
        return out
    cx,cy=cxy()
    x=ref[:,0]-cx;y=ref[:,1]-cy;zr=ref[:,2]
    _,rr=fit(x,y,zr);_,rt=fit(x,y,zt)
    pear=float(np.corrcoef(rr,rt)[0,1]) if np.std(rr)>0 and np.std(rt)>0 else None
    spear=float(spearmanr(rr,rt).statistic) if len(rr)>=3 else None
    sign=float(np.mean(np.sign(rr)==np.sign(rt)))
    p05r,p95r=np.percentile(zr,[5,95]);p05t,p95t=np.percentile(zt,[5,95])
    relr=float(p95r-p05r);relt=float(p95t-p05t)
    gate=bool(pear is not None and spear is not None and pear>=0.60 and spear>=0.60 and sign>=0.65)
    out.update({
      "gate_pass":gate,
      "residual_pearson":pear,
      "residual_spearman":spear,
      "residual_sign_concordance":sign,
      "median_depth_offset_target_minus_2005_m":float(np.median(zt-zr)),
      "relief_2005_m":relr,
      "relief_target_m":relt,
      "relief_ratio_target_over_2005":None if relr==0 else relt/relr,
      "verdict":"PASS_PRODUCT_SPECIFIC_TEMPORAL_MORPHOLOGY" if gate else "FAIL_PRODUCT_SPECIFIC_TEMPORAL_MORPHOLOGY"
    })
    return out

product_results={}
for name,shaexp in PRODUCTS.items():
    pts,sha=load_cd_product(name,shaexp)
    thin=thin_reference_2005(pts)
    product_results[name]={
      "sha256":sha,
      "raw_points_within_1000m":int(len(pts)),
      "geometry_thinned_support_points":int(len(thin)),
      "comparisons":{
        lab:compare_product_support(lab,thin,loaded[lab])
        for lab in ["KN192-07_2008","KNOX15RR_2008"]
      }
    }

passes_both={}
for name,r in product_results.items():
    comps=r["comparisons"]
    passes_both[name]=all(comps[lab].get("gate_pass",False) for lab in comps)
positive=(
    sum(bool(v) for v in passes_both.values())>=1 and
    any(not all(r["comparisons"][lab].get("gate_pass",False) for lab in r["comparisons"]) for r in product_results.values())
)
insufficient_any=any(
    any(c.get("verdict")=="INSUFFICIENT_PRODUCT_SPECIFIC_SUPPORT" for c in r["comparisons"].values())
    for r in product_results.values()
)
if positive:
    verdict="PRODUCT_MIXING_SIGNATURE_SUPPORTED"
elif insufficient_any:
    verdict="INSUFFICIENT_PRODUCT_SPECIFIC_SUPPORT"
else:
    verdict="NO_PRODUCT_MIXING_SIGNATURE"

pout={
 "artifact_id":"JANUS-KUSTO-CAND003-CD169-EM12-PRODUCT-DEMIX-RUN-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-CD169-EM12-PRODUCT-DEMIX-PREREG-2026-09-22-v1.0.json",
 "target":C,
 "products":product_results,
 "passes_both_2008_lineages":passes_both,
 "product_mixing_signature_pass":positive,
 "verdict":verdict,
 "strict_parent_fail_immutable":True,
 "claim_ceiling":"PROCESSED_PRODUCT_DEMIX_DIAGNOSTIC_ONLY"
}
pp=OUT/"JANUS-KUSTO-CAND003-CD169-EM12-PRODUCT-DEMIX-RUN-2026-09-22-v1.0.json"
pp.write_text(json.dumps(pout,indent=2))
print("\nPRODUCT_DEMIX_DIAGNOSTIC")
print(json.dumps(pout,indent=2))
