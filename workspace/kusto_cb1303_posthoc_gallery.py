#!/usr/bin/env python3
import io,json,os,subprocess,tempfile,zipfile,requests,hashlib,math
from pathlib import Path
import numpy as np
import rasterio
from rasterio.windows import Window
from pyproj import Transformer
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"/"cb1303_gallery"; OUT.mkdir(parents=True,exist_ok=True)
TOKEN=os.environ["GITHUB_TOKEN"]; REPO=os.environ["GITHUB_REPOSITORY"]
ART=10835674072
FREEZE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1H-CB13_03-UNCHANGED-V2-IMPLEMENTATION-FREEZE-2026-09-25-v1.0.json").read_text())
H={"Authorization":f"Bearer {TOKEN}","Accept":"application/vnd.github+json","X-GitHub-Api-Version":"2022-11-28","User-Agent":"JANUS-KUSTO-posthoc-gallery"}

def dl(url,path=None):
 r=requests.get(url,headers=H if "api.github.com" in url else {"User-Agent":"JANUS-KUSTO-posthoc-gallery"},timeout=900,stream=bool(path))
 r.raise_for_status()
 if path:
  with open(path,"wb") as f:
   for b in r.iter_content(1024*1024):
    if b:f.write(b)
  return
 return r.content

raw=dl(f"https://api.github.com/repos/{REPO}/actions/artifacts/{ART}/zip")
with zipfile.ZipFile(io.BytesIO(raw)) as z:
 n=[x for x in z.namelist() if x.endswith(".json")][0]; run=json.loads(z.read(n))
cand={x["candidate_id"]:x for x in run["candidate_results"] if x.get("cross_channel_support")}
TARGETS=["INFOMAR_CB13_03_C064","INFOMAR_CB13_03_C062","INFOMAR_CB13_03_C024","INFOMAR_CB13_03_C020","INFOMAR_CB13_03_C013","INFOMAR_CB13_03_C016","INFOMAR_CB13_03_C087","INFOMAR_CB13_03_C002"]
to_utm=Transformer.from_crs("EPSG:4326","EPSG:32629",always_xy=True)
with tempfile.TemporaryDirectory() as td:
 td=Path(td)
 paths={}
 for label,key in [("bathymetry","bathymetry"),("backscatter","heldout_backscatter")]:
  spec=FREEZE["inputs"][key]; zp=td/f"{label}.zip"; outd=td/label; outd.mkdir()
  dl(spec["archive_url"],zp)
  subprocess.run(["7z","e","-y",f"-o{outd}",str(zp),spec["member"]],check=True,stdout=subprocess.DEVNULL)
  paths[label]=outd/Path(spec["member"]).name
 with rasterio.open(paths["bathymetry"]) as bd, rasterio.open(paths["backscatter"]) as ks:
  manifest=[]
  for cid in TARGETS:
   c=cand[cid]; lon,lat=float(c["lon"]),float(c["lat"]); x,y=to_utm.transform(lon,lat)
   # 180 m square bathymetry
   br,bc=bd.index(x,y); bhalf=max(40,int(round(90/abs(bd.transform.a))))
   bw=Window(max(0,bc-bhalf),max(0,br-bhalf),2*bhalf+1,2*bhalf+1)
   ba=bd.read(1,window=bw,masked=True).astype("float64")
   # 180 m square backscatter using reference approx metres/pixel
   kr,kc=ks.index(lon,lat); kh=60
   kw=Window(max(0,kc-kh),max(0,kr-kh),2*kh+1,2*kh+1)
   ka=ks.read(1,window=kw,masked=True).astype("float64")
   # robust display limits
   def lim(a):
    v=a.compressed()
    return (float(np.quantile(v,.02)),float(np.quantile(v,.98))) if len(v) else (0,1)
   bv=lim(ba); kv=lim(ka)
   fig,ax=plt.subplots(1,2,figsize=(10,5))
   ax[0].imshow(ba,cmap="gray",vmin=bv[0],vmax=bv[1]); ax[0].set_title(f"{cid} bathymetry"); ax[0].plot(ba.shape[1]/2,ba.shape[0]/2,"r+",ms=12,mew=2)
   ax[1].imshow(ka,cmap="gray",vmin=kv[0],vmax=kv[1]); ax[1].set_title("held-out backscatter"); ax[1].plot(ka.shape[1]/2,ka.shape[0]/2,"r+",ms=12,mew=2)
   for a in ax:a.set_axis_off()
   fig.suptitle(f"{lat:.7f}, {lon:.7f} | r={c['radius_m']:.1f} m | backscore={c['aggregate_score']:.2f}")
   fig.tight_layout()
   png=OUT/f"{cid}.png"; fig.savefig(png,dpi=160,bbox_inches="tight"); plt.close(fig)
   manifest.append({"candidate_id":cid,"lat":lat,"lon":lon,"radius_m":c["radius_m"],"aggregate_score":c["aggregate_score"],"q95":c["same_stratum_control_q95"],"png":png.name})
 (OUT/"manifest.json").write_text(json.dumps({"status":"POSTHOC_VISUAL_CHARACTERIZATION_ONLY","targets":manifest},indent=2))
 print(json.dumps({"status":"PASS","targets":len(manifest),"out":str(OUT)},indent=2))
