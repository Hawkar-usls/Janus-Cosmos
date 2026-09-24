#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, io, json, math, os, zipfile
from pathlib import Path
import requests
from pyproj import Geod

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1J-CB13_03-EXTREME-UNEXPLAINED-TRIAGE-PREREG-2026-09-25-v1.0.json").read_text())
FORM=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1I-CB13_03-FORMAL-UNSEEN-V2-PASS-RECEIPT-2026-09-25-v1.0.json").read_text())

REPO=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
TOKEN=os.environ.get("GITHUB_TOKEN")
ART_ID=int(FORM["run"]["artifact_id"])
ART_SHA=FORM["run"]["artifact_zip_sha256"].lower()
UA={"User-Agent":"JANUS-KUSTO-I1J-CB13_03/1.0"}
S=requests.Session(); S.headers.update(UA)
GEOD=Geod(ellps="WGS84")

GSI_WRECK="https://gsi.geodata.gov.ie/server/rest/services/Marine/IE_GSI_MI_Shipwrecks_IE_Waters_WGS84_LAT/FeatureServer/0"
SEDIMENT="https://gsi.geodata.gov.ie/server/rest/services/Marine/IE_GSI_MI_Seabed_Sediment_Classification_IE_Waters_WGS84/FeatureServer/1"
SAMPLES="https://gsi.geodata.gov.ie/server/rest/services/Marine/IE_GSI_MI_Seabed_Sediment_Samples_IE_Waters_WGS84_LAT/FeatureServer/0"
NMS_URL="https://www.arcgis.com/sharing/rest/content/items/d4b084c880b546fabe38345461b563d2/data"

def get_json(url,params=None,timeout=180,headers=None):
    h=dict(headers or {})
    r=S.get(url,params=params,headers=h,timeout=timeout,allow_redirects=True)
    r.raise_for_status()
    j=r.json()
    if isinstance(j,dict) and "error" in j: raise RuntimeError(f"{url}: {j['error']}")
    return j,r.url

def load_formal():
    if not TOKEN: raise RuntimeError("GITHUB_TOKEN required")
    url=f"https://api.github.com/repos/{REPO}/actions/artifacts/{ART_ID}/zip"
    r=S.get(url,headers={"Authorization":f"Bearer {TOKEN}","Accept":"application/vnd.github+json","X-GitHub-Api-Version":"2022-11-28"},timeout=180)
    r.raise_for_status()
    raw=r.content
    sha=hashlib.sha256(raw).hexdigest()
    if sha.lower()!=ART_SHA: raise RuntimeError(f"formal artifact SHA mismatch {sha} != {ART_SHA}")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        js=[n for n in z.namelist() if n.lower().endswith(".json") and not n.startswith("__MACOSX/")]
        if len(js)!=1: raise RuntimeError(f"expected one JSON in formal artifact, got {js}")
        j=json.loads(z.read(js[0]))
    if not j["population"]["primary_pass"]: raise RuntimeError("formal parent is not PASS")
    return j,sha

def env_params(lon,lat,meters,outfields="*"):
    dlat=meters/110540.0
    dlon=meters/(111320.0*max(0.05,math.cos(math.radians(lat))))
    env={"xmin":lon-dlon,"ymin":lat-dlat,"xmax":lon+dlon,"ymax":lat+dlat,"spatialReference":{"wkid":4326}}
    return {
      "where":"1=1","geometry":json.dumps(env,separators=(",",":")),"geometryType":"esriGeometryEnvelope",
      "inSR":"4326","spatialRel":"esriSpatialRelIntersects","outFields":outfields,
      "returnGeometry":"true","outSR":"4326","f":"json"
    }

def gsi_wrecks(lon,lat,radius,search_m=5000.0):
    j,url=get_json(GSI_WRECK+"/query",env_params(lon,lat,search_m))
    out=[]
    for f in j.get("features",[]) or []:
        g=f.get("geometry") or {}; x=g.get("x"); y=g.get("y")
        if x is None or y is None: continue
        _,_,d=GEOD.inv(lon,lat,float(x),float(y)); d=abs(float(d))
        if d>search_m: continue
        cls="DIRECT" if d<=radius else "NEAR" if d<=500 else "REGIONAL" if d<=2000 else "FAR"
        out.append({"distance_m":d,"distance_class":cls,"lon":float(x),"lat":float(y),"attributes":f.get("attributes") or {}})
    out.sort(key=lambda x:x["distance_m"])
    return out,url

def sediment_at(lon,lat):
    pt={"x":lon,"y":lat,"spatialReference":{"wkid":4326}}
    p={"where":"1=1","geometry":json.dumps(pt,separators=(",",":")),"geometryType":"esriGeometryPoint",
       "inSR":"4326","spatialRel":"esriSpatialRelIntersects","outFields":"OBJECTID,BIOZONE,SUBSTRATE,FOLK_5,EUNIS,MSFD_BBHT,DATASOURCE,RESOLUTION",
       "returnGeometry":"false","f":"json"}
    j,url=get_json(SEDIMENT+"/query",p)
    return [f.get("attributes") or {} for f in j.get("features",[]) or []],url

def nearest_samples(lon,lat,search_m=2000.0):
    j,url=get_json(SAMPLES+"/query",env_params(lon,lat,search_m))
    rows=[]
    for f in j.get("features",[]) or []:
        g=f.get("geometry") or {}; x=g.get("x"); y=g.get("y")
        if x is None or y is None: continue
        _,_,d=GEOD.inv(lon,lat,float(x),float(y)); d=abs(float(d))
        if d<=search_m:
            rows.append({"distance_m":d,"lon":float(x),"lat":float(y),"attributes":f.get("attributes") or {}})
    rows.sort(key=lambda x:x["distance_m"])
    return rows,url

def normalize(s):
    return "".join(ch.lower() for ch in str(s) if ch.isalnum())

def to_float(v):
    if v is None: return None
    s=str(v).strip().replace("−","-")
    if not s: return None
    try: return float(s)
    except: return None

def load_nms():
    r=S.get(NMS_URL,timeout=300,allow_redirects=True)
    r.raise_for_status()
    raw=r.content
    if raw[:2]==b"PK":
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            cs=[n for n in z.namelist() if n.lower().endswith(".csv")]
            if not cs: raise RuntimeError("NMS item zip has no CSV")
            raw=z.read(cs[0])
    text=raw.decode("utf-8-sig",errors="replace")
    rd=csv.DictReader(io.StringIO(text))
    fields=rd.fieldnames or []
    norms={f:normalize(f) for f in fields}
    lat_candidates=[f for f,n in norms.items() if n in ("latitude","lat","latdd","ddlatitude","ddlat") or "latitude" in n]
    lon_candidates=[f for f,n in norms.items() if n in ("longitude","long","lon","londd","ddlongitude","ddlong") or "longitude" in n]
    if not lat_candidates or not lon_candidates:
        raise RuntimeError(f"NMS coordinate columns not found; fields={fields}")
    rows=[]; total=0
    for row in rd:
        total+=1
        lat=next((to_float(row.get(f)) for f in lat_candidates if to_float(row.get(f)) is not None),None)
        lon=next((to_float(row.get(f)) for f in lon_candidates if to_float(row.get(f)) is not None),None)
        if lat is None or lon is None or not(-90<=lat<=90 and -180<=lon<=180): continue
        rows.append({"lon":lon,"lat":lat,"attributes":row})
    return rows,{"url":r.url,"sha256":hashlib.sha256(raw).hexdigest(),"total_records":total,"located_records":len(rows),"fields":fields,"lat_fields":lat_candidates,"lon_fields":lon_candidates}

def nearest_nms(lon,lat,nms,search_m=5000.0):
    rows=[]
    for q in nms:
        _,_,d=GEOD.inv(lon,lat,q["lon"],q["lat"]); d=abs(float(d))
        if d<=search_m:
            rows.append({"distance_m":d,"lon":q["lon"],"lat":q["lat"],"attributes":q["attributes"]})
    rows.sort(key=lambda x:x["distance_m"])
    return rows

formal,formal_zip_sha=load_formal()
supported=[c for c in formal["candidate_results"] if c.get("cross_channel_support") is True]
supported.sort(key=lambda x:(-float(x["aggregate_score"]),int(x["blind_rank"])))
if len(supported)!=35: raise RuntimeError(f"expected 35 supported candidates, got {len(supported)}")
nms,nms_meta=load_nms()

rows=[]
for c in supported:
    lon=float(c["lon"]); lat=float(c["lat"]); radius=float(c["radius_m"])
    gwr,gurl=gsi_wrecks(lon,lat,radius)
    sed,surl=sediment_at(lon,lat)
    nwr=nearest_nms(lon,lat,nms)
    folk=sorted({str(x.get("FOLK_5") or "").strip() for x in sed if str(x.get("FOLK_5") or "").strip()})
    nearest_gsi=gwr[0] if gwr else None
    nearest_nms_row=nwr[0] if nwr else None
    wreck500=bool((nearest_gsi and nearest_gsi["distance_m"]<=500) or (nearest_nms_row and nearest_nms_row["distance_m"]<=500))
    folk_ok=bool(folk) and all(x.lower() not in ("rock","unclassified") for x in folk)
    rows.append({
      "candidate_id":c["candidate_id"],"blind_rank":int(c["blind_rank"]),"lon":lon,"lat":lat,"radius_m":radius,
      "aggregate_score":float(c["aggregate_score"]),"control_q95":float(c["same_stratum_control_q95"]),
      "backscatter_metrics":c["metrics"],"metric_positive_robust_z":c["metric_positive_robust_z"],
      "folk_classes":folk,"sediment_records":sed,
      "nearest_gsi_wreck":nearest_gsi,"nearest_nms_wreck":nearest_nms_row,
      "gsi_wreck_within_500m":bool(nearest_gsi and nearest_gsi["distance_m"]<=500),
      "nms_wreck_within_500m":bool(nearest_nms_row and nearest_nms_row["distance_m"]<=500),
      "primary_eligible":bool((not wreck500) and folk_ok),
      "_gsi_query_url":gurl,"_sediment_query_url":surl
    })

selected=next((r for r in rows if r["primary_eligible"]),None)
selection_mode="PRIMARY_NONROCK_NO_WRECK_500M"
if selected is None:
    selected=next((r for r in rows if not (r["gsi_wreck_within_500m"] or r["nms_wreck_within_500m"])),None)
    selection_mode="FALLBACK_NO_WRECK_500M"
if selected is None: raise RuntimeError("no target survives frozen selection rule")

samples,sampurl=nearest_samples(selected["lon"],selected["lat"],2000.0)
selected["nearest_official_sediment_samples"]=samples[:10]
selected["_sample_query_url"]=sampurl

out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I1J-CB13_03-EXTREME-UNEXPLAINED-TRIAGE-RUN-2026-09-25-v1.0",
 "prereg":PRE["artifact_id"],"formal_parent":FORM["artifact_id"],
 "formal_artifact_zip_sha256_verified":formal_zip_sha,
 "supported_candidate_count":len(supported),
 "official_sources":{
   "gsi_mi_shipwreck_layer":GSI_WRECK,
   "nms_wreck_inventory":{"download":NMS_URL,**nms_meta},
   "gsi_mi_seabed_sediment_layer":SEDIMENT,
   "gsi_mi_sediment_samples_layer":SAMPLES
 },
 "selection_mode":selection_mode,
 "selected_target":selected,
 "score_order_characterization":rows,
 "formal_parent_changed":False,"candidate_coordinates_changed":False,"candidate_radii_changed":False,
 "promotion":False,
 "claim_ceiling":"POSTRESULT_EXTREME_UNEXPLAINED_TARGET_TRIAGE_ONLY"
}
p=OUT/"JANUS-KUSTO-INFOMAR-I1J-CB13_03-EXTREME-UNEXPLAINED-TRIAGE-RUN-2026-09-25-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False); p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"selection_mode":selection_mode,
 "selected_candidate":selected["candidate_id"],"blind_rank":selected["blind_rank"],
 "lon":selected["lon"],"lat":selected["lat"],"radius_m":selected["radius_m"],
 "aggregate_score":selected["aggregate_score"],"control_q95":selected["control_q95"],
 "folk_classes":selected["folk_classes"],
 "nearest_gsi_wreck_m":None if selected["nearest_gsi_wreck"] is None else selected["nearest_gsi_wreck"]["distance_m"],
 "nearest_nms_wreck_m":None if selected["nearest_nms_wreck"] is None else selected["nearest_nms_wreck"]["distance_m"],
 "nearest_sample_m":None if not samples else samples[0]["distance_m"],
 "nms_located_records":nms_meta["located_records"],
 "output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2,ensure_ascii=False))
