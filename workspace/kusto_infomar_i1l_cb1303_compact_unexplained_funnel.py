#!/usr/bin/env python3
from __future__ import annotations
import hashlib, io, json, math, os, shutil, subprocess, tempfile, zipfile
from pathlib import Path
import numpy as np
import requests
import rasterio
from rasterio.windows import Window
from scipy import ndimage

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1L-CB13_03-COMPACT-UNEXPLAINED-FUNNEL-PREREG-2026-09-25-v1.0.json").read_text())
IMP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1L-CB13_03-COMPACT-UNEXPLAINED-FUNNEL-IMPLEMENTATION-FREEZE-2026-09-25-v1.0.json").read_text())
I1K=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1K-C002-SOURCE-VALUE-IMPLEMENTATION-FREEZE-2026-09-25-v1.0.json").read_text())
TOKEN=os.environ.get("GITHUB_TOKEN");REPO=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
S=requests.Session();S.headers.update({"User-Agent":"JANUS-KUSTO-I1L/1.0","Accept-Encoding":"identity"})

def gh_artifact_json(spec):
    if not TOKEN: raise RuntimeError("GITHUB_TOKEN missing")
    u=f"https://api.github.com/repos/{REPO}/actions/artifacts/{spec['id']}/zip"
    r=S.get(u,headers={"Authorization":f"Bearer {TOKEN}","Accept":"application/vnd.github+json","X-GitHub-Api-Version":"2022-11-28"},timeout=180)
    r.raise_for_status();raw=r.content
    sha=hashlib.sha256(raw).hexdigest()
    if sha.lower()!=spec["sha256"].lower(): raise RuntimeError(f"artifact SHA drift {spec['id']} {sha}")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        js=[n for n in z.namelist() if n.lower().endswith(".json") and not n.startswith("__MACOSX/")]
        if len(js)!=1: raise RuntimeError(f"artifact {spec['id']} expected one JSON {js}")
        j=json.loads(z.read(js[0]))
    return j,sha

def download_source(spec,tmp):
    z=tmp/"bath.zip";h=hashlib.sha256();n=0
    with S.get(spec["url"],timeout=600,stream=True,allow_redirects=True) as r:
        r.raise_for_status()
        with open(z,"wb") as f:
            for b in r.iter_content(1024*1024):
                if not b: continue
                f.write(b);h.update(b);n+=len(b)
    if h.hexdigest()!=spec["archive_sha256"] or n!=spec["archive_bytes"]:raise RuntimeError("bath archive drift")
    od=tmp/"bath";od.mkdir();seven=shutil.which("7z")
    if not seven:raise RuntimeError("7z missing")
    subprocess.run([seven,"e","-y",f"-o{od}",str(z),spec["member"]],check=True,stdout=subprocess.DEVNULL)
    p=od/Path(spec["member"]).name
    mh=hashlib.sha256();mn=0
    with open(p,"rb") as f:
        while True:
            b=f.read(8*1024*1024)
            if not b:break
            mh.update(b);mn+=len(b)
    if mh.hexdigest()!=spec["member_sha256"] or mn!=spec["member_bytes"]:raise RuntimeError("bath member drift")
    return p,{"archive_sha256":h.hexdigest(),"member_sha256":mh.hexdigest()}

def read_win(ds,x,y,half,dx,dy):
    row,col=rasterio.transform.rowcol(ds.transform,x,y)
    hx=int(math.ceil(half/dx));hy=int(math.ceil(half/dy))
    ma=ds.read(1,window=Window(col-hx,row-hy,2*hx+1,2*hy+1),masked=True,boundless=True)
    a=np.asarray(ma.data,dtype=float);v=(~np.ma.getmaskarray(ma))&np.isfinite(a)
    rr=(np.arange(a.shape[0])-hy)*dy;cc=(np.arange(a.shape[1])-hx)*dx
    yy,xx=np.meshgrid(rr,cc,indexing="ij");d=np.hypot(xx,yy)
    return a,v,d,xx,yy,int(row),int(col)

def line_check(a,v):
    cr=a.shape[0]//2;cc=a.shape[1]//2
    hs=min(100,cc-2,a.shape[1]-cc-2);vs=min(100,cr-2,a.shape[0]-cr-2)
    rs=[];rgs=[]
    for off in range(-10,11):
        z=a[cr+off,cc-hs:cc+hs+1];ok=v[cr+off,cc-hs:cc+hs+1]
        g=np.abs(np.diff(z));gv=np.where(ok[:-1]&ok[1:],g,np.nan);rs.append(float(np.nanmax(gv)));rgs.append(gv)
    cs=[];cgs=[]
    for off in range(-10,11):
        z=a[cr-vs:cr+vs+1,cc+off];ok=v[cr-vs:cr+vs+1,cc+off]
        g=np.abs(np.diff(z));gv=np.where(ok[:-1]&ok[1:],g,np.nan);cs.append(float(np.nanmax(gv)));cgs.append(gv)
    def pct(z):return float(100*np.mean(np.asarray(z)<=z[10]))
    def mc(gs):
        t=gs[10];co=[]
        for i,g in enumerate(gs):
            if i==10:continue
            ok=np.isfinite(t)&np.isfinite(g)
            if ok.sum()>=5 and np.std(t[ok])>0 and np.std(g[ok])>0:co.append(float(np.corrcoef(t[ok],g[ok])[0,1]))
        return None if not co else float(np.median(co))
    rp,cp=pct(rs),pct(cs);rc,ccor=mc(rgs),mc(cgs)
    possible=bool((rp>=99 and (rc or -1)>=.90) or (cp>=99 and (ccor or -1)>=.90))
    return {"row_percentile":rp,"row_parallel_gradient_corr":rc,"col_percentile":cp,"col_parallel_gradient_corr":ccor,"possible_line_artefact":possible}

def characterize(ds,c):
    x=float(c["x_m"]);y=float(c["y_m"]);r=float(c["radius_m"])
    dx=abs(float(ds.transform.a));dy=abs(float(ds.transform.e))
    a50,v50,d50,x50,y50,row,col=read_win(ds,x,y,50,dx,dy)
    a100,v100,d100,x100,y100,_,_=read_win(ds,x,y,100,dx,dy)
    a250,v250,d250,x250,y250,_,_=read_win(ds,x,y,250,dx,dy)
    disk=(d50<=r)&v50;ann=(d50>1.25*r)&(d50<=2*r)&v50
    if disk.sum()<5 or ann.sum()<5:return {"status":"INSUFFICIENT_LOCAL_SUPPORT","valid_fraction_100m":float(v100.mean())}
    dm=float(np.median(a50[disk]));am=float(np.median(a50[ann]));sg=1 if dm-am>0 else -1 if dm-am<0 else 0
    extmask=(d50<=50)&v50
    idx=np.nanargmax(np.where(extmask,a50,np.nan)) if sg>=0 else np.nanargmin(np.where(extmask,a50,np.nan))
    er,ec=np.unravel_index(idx,a50.shape);ext=float(a50[er,ec]);ed=float(d50[er,ec]);thr=am+.5*(ext-am)
    mask=((a100>=thr) if sg>=0 else (a100<=thr))&v100
    lab,n=ndimage.label(mask,structure=np.ones((3,3),dtype=np.uint8))
    cr=a100.shape[0]//2;cc=a100.shape[1]//2
    de=float(x50[er,ec]);dn=float(y50[er,ec]);rr=int(round(cr+dn/dy));cl=int(round(cc+de/dx))
    L=int(lab[rr,cl]) if 0<=rr<lab.shape[0] and 0<=cl<lab.shape[1] else 0
    comp=(lab==L) if L>0 else np.zeros_like(lab,dtype=bool)
    area=float(comp.sum()*dx*dy);eq=float(math.sqrt(4*area/math.pi)) if area else 0.0
    touch=bool(np.any(comp[0]) or np.any(comp[-1]) or np.any(comp[:,0]) or np.any(comp[:,-1]))
    pts=np.column_stack([x100[comp],y100[comp]]) if comp.any() else np.empty((0,2))
    asp=None;spans=[0.0,0.0]
    if len(pts)>=2:
        q=pts-pts.mean(0);vals,vec=np.linalg.eigh(np.cov(q,rowvar=False));vec=vec[:,np.argsort(vals)[::-1]]
        pr=q@vec;ss=np.ptp(pr,axis=0);spans=[float(max(ss)),float(min(ss))];asp=None if min(ss)<=0 else float(max(ss)/min(ss))
    lc=line_check(a250,v250)
    broad=touch
    compact=bool(float(v100.mean())>=.95 and not touch and 4<=eq<=100 and ed<=max(20,2*r) and not lc["possible_line_artefact"])
    return {
      "status":"SCORABLE","target_pixel":[row,col],"native_relief_sign":"POSITIVE" if sg>0 else "NEGATIVE" if sg<0 else "ZERO",
      "disk_median":dm,"annulus_median":am,"signed_disk_minus_annulus_median":dm-am,
      "local_extremum_value":ext,"local_extremum_distance_m":ed,"half_relief_threshold":thr,
      "half_relief_component_area_m2":area,"half_relief_equivalent_diameter_m":eq,
      "principal_axis_span_m":spans,"aspect_ratio":asp,"component_touches_100m_window":touch,
      "valid_fraction_100m":float(v100.mean()),"valid_fraction_250m":float(v250.mean()),
      "line_artefact_check":lc,"broad_geologic_relief":broad,"compact_survivor":compact
    }

formal,fsha=gh_artifact_json(IMP["formal_artifact"])
triage,tsha=gh_artifact_json(IMP["triage_artifact"])
fm={c["candidate_id"]:c for c in formal["candidate_results"]}
order=triage["score_order_characterization"]
eligible=[]
for t in order:
    if t.get("gsi_wreck_within_500m") or t.get("nms_wreck_within_500m"):continue
    folk=[str(x).strip().lower() for x in t.get("folk_classes",[])]
    if "rock" in folk:continue
    c=fm[t["candidate_id"]].copy()
    c["folk_classes"]=t.get("folk_classes",[])
    c["nearest_gsi_wreck"]=t.get("nearest_gsi_wreck")
    c["nearest_nms_wreck"]=t.get("nearest_nms_wreck")
    eligible.append(c)

with tempfile.TemporaryDirectory(prefix="i1l_") as td:
    bp,bmeta=download_source(IMP["bathymetry_source"],Path(td))
    results=[]
    with rasterio.open(bp) as ds:
        for c in eligible:
            m=characterize(ds,c)
            results.append({
              "candidate_id":c["candidate_id"],"blind_rank":int(c["blind_rank"]),"aggregate_score":float(c["aggregate_score"]),
              "lon":float(c["lon"]),"lat":float(c["lat"]),"x_m":float(c["x_m"]),"y_m":float(c["y_m"]),"radius_m":float(c["radius_m"]),
              "folk_classes":c["folk_classes"],"morphology":m
            })

selected=next((r for r in results if r["candidate_id"]!="INFOMAR_CB13_03_C002" and r["morphology"].get("compact_survivor") is True),None)
status="COMPACT_UNEXPLAINED_SURVIVOR_SELECTED" if selected else "NO_COMPACT_UNEXPLAINED_SURVIVOR"
out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I1L-CB13_03-COMPACT-UNEXPLAINED-FUNNEL-RUN-2026-09-25-v1.0",
 "prereg":PRE["artifact_id"],"implementation_freeze":IMP["artifact_id"],
 "formal_artifact_sha_verified":fsha,"triage_artifact_sha_verified":tsha,"bathymetry_source_identity":bmeta,
 "supported_candidate_count":len(order),"eligible_nonrock_no_near_wreck_count":len(eligible),
 "results_in_frozen_score_order":results,"selected_target":selected,"status":status,
 "selection_rule_changed":False,"manual_visual_selection_used":False,
 "claim_ceiling":"POSTRESULT_COMPACT_UNEXPLAINED_TARGET_FUNNEL_ONLY"
}
raw=json.dumps(out,indent=2,ensure_ascii=False)
p=OUT/"JANUS-KUSTO-INFOMAR-I1L-CB13_03-COMPACT-UNEXPLAINED-FUNNEL-RUN-2026-09-25-v1.0.json";p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"status":status,"eligible_count":len(eligible),
 "selected":None if selected is None else {
   "candidate_id":selected["candidate_id"],"blind_rank":selected["blind_rank"],"score":selected["aggregate_score"],
   "lon":selected["lon"],"lat":selected["lat"],"radius_m":selected["radius_m"],"folk":selected["folk_classes"],
   "morphology":selected["morphology"]
 },
 "compact_survivor_ids":[r["candidate_id"] for r in results if r["morphology"].get("compact_survivor") is True],
 "output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2,ensure_ascii=False))
