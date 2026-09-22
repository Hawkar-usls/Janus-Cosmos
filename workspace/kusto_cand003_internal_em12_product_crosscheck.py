#!/usr/bin/env python3
import ftplib, hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/EM12"
TARGET={"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432}
FILES={
 "ACCEPTL15_22":("B1-81-1_Acceptl15-22.xyz.ascii","ad2e4d60401f7ffaa5e14b3e43aa0e45da9d320c6b34fa7f3b873f4024be6007"),
 "ACCEPTL28_33":("B1-81-1_Acceptl28-33.xyz.ascii","0e8cac049b53dd6ba0e5013e8aaea2bfb48347db887edc3c4e3e9f697e96dcf0")
}
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
R=6371008.8;CELL=100.0;RAD=1000.0

def ftp():
    f=ftplib.FTP(timeout=120);f.connect(HOST,21);f.login("anonymous","janus-kusto@example.invalid");f.voidcmd("TYPE I");return f

def parse(raw):
    s=raw.decode("ascii","ignore").strip()
    if not s or s[0] in "#;!":return None
    p=s.replace(","," ").split()
    if len(p)<3:return None
    try:lon=float(p[0]);lat=float(p[1]);z=float(p[2])
    except:return None
    if not all(map(math.isfinite,(lon,lat,z))):return None
    return lon,lat,z

def en(lat,lon):
    north=math.radians(lat-TARGET["lat"])*R
    east=math.radians(lon-TARGET["lon"])*R*math.cos(math.radians((lat+TARGET["lat"])/2))
    return east,north

def load(name,expected):
    f=ftp();h=hashlib.sha256();pts=[]
    path=ROOT+"/"+name
    try:
        sock=f.transfercmd("RETR "+path);stream=sock.makefile("rb")
        try:
            for raw in stream:
                h.update(raw);q=parse(raw)
                if q is None:continue
                lon,lat,z=q;e,n=en(lat,lon)
                if math.hypot(e,n)<=RAD:pts.append((e,n,z))
        finally:stream.close();sock.close()
    finally:f.close()
    sha=h.hexdigest()
    if sha!=expected:raise RuntimeError(f"{name} SHA mismatch {sha}")
    return np.asarray(pts,float),sha

def grid(pts):
    d={}
    for e,n,z in pts:
        ix=math.floor((e+RAD)/CELL);iy=math.floor((n+RAD)/CELL)
        cx=-RAD+(ix+.5)*CELL;cy=-RAD+(iy+.5)*CELL
        if math.hypot(cx,cy)>RAD:continue
        d.setdefault((ix,iy),[]).append(float(z))
    return {k:{"x":-RAD+(k[0]+.5)*CELL,"y":-RAD+(k[1]+.5)*CELL,"z":float(np.median(v)),"n":len(v)} for k,v in d.items()}

def fit(x,y,z):
    A=np.column_stack([np.ones(len(x)),x,y]);coef=np.linalg.lstsq(A,z,rcond=None)[0]
    res=z-A@coef
    az=(math.degrees(math.atan2(coef[1],coef[2]))+360)%360
    return coef,res,az

pts={};sha={}
for k,(name,expected) in FILES.items():
    pts[k],sha[k]=load(name,expected)

g={k:grid(v) for k,v in pts.items()}
common=sorted(set(g["ACCEPTL15_22"]) & set(g["ACCEPTL28_33"]))
out={
 "artifact_id":"JANUS-KUSTO-CAND003-CD169-EM12-INTERNAL-PRODUCT-CROSSCHECK-RUN-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-CD169-EM12-INTERNAL-PRODUCT-CROSSCHECK-PREREG-2026-09-22-v1.0.json",
 "target":TARGET,
 "support":{
   k:{"raw_points_within_1000m":int(len(pts[k])),"grid_cells":len(g[k]),"sha256":sha[k]} for k in pts
 },
 "common_cell_count":len(common)
}
if len(common)<30:
    out.update({"gate_pass":False,"verdict":"INSUFFICIENT_COMMON_PRODUCT_SUPPORT"})
else:
    x=np.array([g["ACCEPTL15_22"][k]["x"] for k in common])
    y=np.array([g["ACCEPTL15_22"][k]["y"] for k in common])
    z1=np.array([g["ACCEPTL15_22"][k]["z"] for k in common])
    z2=np.array([g["ACCEPTL28_33"][k]["z"] for k in common])
    c1,r1,a1=fit(x,y,z1);c2,r2,a2=fit(x,y,z2)
    pear=float(np.corrcoef(r1,r2)[0,1]) if np.std(r1)>0 and np.std(r2)>0 else None
    spear=float(spearmanr(r1,r2).statistic)
    sign=float(np.mean(np.sign(r1)==np.sign(r2)))
    p051,p951=np.percentile(z1,[5,95]);p052,p952=np.percentile(z2,[5,95])
    ad=abs(a1-a2)%360;ad=min(ad,360-ad)
    gate=bool(pear is not None and pear>=.60 and spear>=.60 and sign>=.65)
    out.update({
      "residual_pearson":pear,
      "residual_spearman":spear,
      "residual_sign_concordance":sign,
      "median_depth_offset_28_33_minus_15_22_m":float(np.median(z2-z1)),
      "relief_15_22_m":float(p951-p051),
      "relief_28_33_m":float(p952-p052),
      "plane_gradient_azimuth_15_22_deg":a1,
      "plane_gradient_azimuth_28_33_deg":a2,
      "plane_gradient_azimuth_difference_deg":ad,
      "gate_pass":gate,
      "verdict":"INTERNAL_2005_PRODUCT_MORPHOLOGY_CONCORDANT" if gate else "INTERNAL_2005_PRODUCT_DISCORDANCE_CONFIRMED"
    })
out["strict_temporal_parent_fail_immutable"]=True
out["claim_ceiling"]="INTERNAL_PROCESSED_CD169_EM12_PRODUCT_CONCORDANCE_ONLY"
p=OUT/"JANUS-KUSTO-CAND003-CD169-EM12-INTERNAL-PRODUCT-CROSSCHECK-RUN-2026-09-22-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
