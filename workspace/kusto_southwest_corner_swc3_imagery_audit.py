#!/usr/bin/env python3
from __future__ import annotations

import csv, hashlib, io, json, math, os, re, zipfile
from collections import Counter
from pathlib import Path

import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC3-INDEPENDENT-IMAGERY-PREREG-2026-09-24-v1.0.json").read_text())
AUTH=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOUTHWEST-CORNER-SWC1-BLIND-BATHYMETRY-RECEIPT-2026-09-24-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-SouthwestCorner-SWC3-imagery-audit/1.0"}

def get(url,headers=None,timeout=180):
    h=dict(UA)
    if headers:h.update(headers)
    r=requests.get(url,headers=h,timeout=timeout,allow_redirects=True);r.raise_for_status();return r

def get_candidates():
    token=os.environ.get("GITHUB_TOKEN")
    repo=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
    aid=int(AUTH["run"]["artifact_id"])
    u=f"https://api.github.com/repos/{repo}/actions/artifacts/{aid}/zip"
    blob=get(u,headers={
      "Authorization":f"Bearer {token}",
      "Accept":"application/vnd.github+json",
      "X-GitHub-Api-Version":"2022-11-28"
    }).content
    zsha=hashlib.sha256(blob).hexdigest()
    if zsha!=AUTH["run"]["artifact_zip_sha256"]:raise RuntimeError("candidate artifact SHA mismatch")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        js=[n for n in zf.namelist() if n.lower().endswith(".json")]
        if len(js)!=1:raise RuntimeError(f"expected one candidate JSON, found {js}")
        raw=zf.read(js[0])
    jsha=hashlib.sha256(raw).hexdigest()
    if jsha!=AUTH["run"]["candidate_json_sha256"]:raise RuntimeError("candidate JSON SHA mismatch")
    j=json.loads(raw)
    return j,zsha,jsha

def norm(s):
    return re.sub(r"[^a-z0-9]+","_",str(s).strip().lower()).strip("_")

def parse_float(x):
    try:
        s=str(x).strip()
        if not s:return None
        return float(s)
    except Exception:return None

def discover_coord_pair(headers,rows):
    nh={h:norm(h) for h in headers}
    lat_heads=[h for h,n in nh.items() if n in ("lat","latitude","decimal_latitude","site_latitude","deployment_latitude") or n.endswith("_lat") or "latitude" in n]
    lon_heads=[h for h,n in nh.items() if n in ("lon","long","longitude","decimal_longitude","site_longitude","deployment_longitude") or n.endswith("_lon") or n.endswith("_long") or "longitude" in n]
    candidates=[]
    for la in lat_heads:
        for lo in lon_heads:
            vals=[]
            for row in rows[:min(5000,len(rows))]:
                a=parse_float(row.get(la));b=parse_float(row.get(lo))
                if a is None or b is None:continue
                if -90<=a<=90 and -180<=b<=180:vals.append((a,b))
            if not vals:continue
            frac=len(vals)/max(1,min(5000,len(rows)))
            in_region=sum(1 for a,b in vals if -35.5<=a<=-32.5 and 113.5<=b<=116.0)/len(vals)
            candidates.append({"lat_col":la,"lon_col":lo,"valid_rows_sample":len(vals),"valid_fraction_sample":frac,"regional_fraction_sample":in_region})
    good=[x for x in candidates if x["valid_rows_sample"]>=10 and x["regional_fraction_sample"]>=0.8]
    if len(good)!=1:
        raise RuntimeError(f"coordinate schema not uniquely resolved: all={candidates}, good={good}")
    return good[0],candidates

def hav(lat1,lon1,lat2,lon2):
    R=6371008.8
    p1=math.radians(lat1);p2=math.radians(lat2)
    dp=math.radians(lat2-lat1);dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(min(1.0,math.sqrt(a)))

candj,zsha,jsha=get_candidates()
cands=list(candj["frozen_candidates"])
url=PRE["independent_imagery_source"]["direct_csv"]
resp=get(url)
raw=resp.content
text=raw.decode("utf-8-sig","replace")
reader=csv.DictReader(io.StringIO(text))
headers=reader.fieldnames or []
rows=list(reader)
if not rows:raise RuntimeError("imagery CSV has no rows")
coord,coord_trials=discover_coord_pair(headers,rows)
la_col=coord["lat_col"];lo_col=coord["lon_col"]

pts=[]
for idx,row in enumerate(rows,1):
    lat=parse_float(row.get(la_col));lon=parse_float(row.get(lo_col))
    if lat is None or lon is None or not(-90<=lat<=90 and -180<=lon<=180):continue
    pts.append({"row_index":idx,"lat":float(lat),"lon":float(lon)})

if not pts:raise RuntimeError("no georeferenced imagery rows after schema binding")

# Deduplicate identical coordinate records for coverage counting while retaining raw row count.
uniq={}
for p in pts:
    key=(round(p["lat"],8),round(p["lon"],8))
    if key not in uniq:uniq[key]=p
unique_pts=list(uniq.values())

results=[]
for c in cands:
    lat=float(c["lat"]);lon=float(c["lon"]);radius=float(c["radius_m"])
    best=None
    for p in unique_pts:
        d=hav(lat,lon,p["lat"],p["lon"])
        if best is None or d<best["distance_m"]:
            best={"row_index":p["row_index"],"lat":p["lat"],"lon":p["lon"],"distance_m":d}
    d=best["distance_m"] if best else None
    if d is not None and d<=radius:state="DIRECT_IMAGE_SUPPORT"
    elif d is not None and d<=1000:state="NEAR_IMAGE_COVERAGE"
    elif d is not None and d<=5000:state="REGIONAL_IMAGE_COVERAGE"
    else:state="NO_IMAGE_COVERAGE"
    results.append({
      "candidate_id":c["candidate_id"],"blind_rank":c["blind_rank"],
      "candidate_lat":lat,"candidate_lon":lon,"radius_m":radius,
      "nearest_imagery_record":best,"coverage_state":state
    })

counts=Counter(x["coverage_state"] for x in results)
direct=[x for x in results if x["coverage_state"]=="DIRECT_IMAGE_SUPPORT"]
near=[x for x in results if x["coverage_state"]=="NEAR_IMAGE_COVERAGE"]
regional=[x for x in results if x["coverage_state"]=="REGIONAL_IMAGE_COVERAGE"]

# For direct matches only, preserve nearest raw CSV row fields for later semantic inspection.
direct_rows=[]
for q in direct:
    ri=int(q["nearest_imagery_record"]["row_index"])
    row=rows[ri-1]
    direct_rows.append({
      "candidate_id":q["candidate_id"],"blind_rank":q["blind_rank"],
      "distance_m":q["nearest_imagery_record"]["distance_m"],
      "radius_m":q["radius_m"],
      "source_row_index":ri,
      "source_row":row
    })

out={
 "artifact_id":"JANUS-KUSTO-SOUTHWEST-CORNER-SWC3-INDEPENDENT-IMAGERY-COVERAGE-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "candidate_authority":{"artifact_zip_sha256_verified":zsha,"candidate_json_sha256_verified":jsha,"frozen_candidate_count":len(cands)},
 "source":{"url":url,"http_status":resp.status_code,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),
   "headers":headers,"raw_rows":len(rows),"georeferenced_rows":len(pts),"unique_georeferenced_coordinates":len(unique_pts)},
 "coordinate_binding":{"selected":coord,"all_trials":coord_trials},
 "coverage_counts":{
   "DIRECT_IMAGE_SUPPORT":counts.get("DIRECT_IMAGE_SUPPORT",0),
   "NEAR_IMAGE_COVERAGE":counts.get("NEAR_IMAGE_COVERAGE",0),
   "REGIONAL_IMAGE_COVERAGE":counts.get("REGIONAL_IMAGE_COVERAGE",0),
   "NO_IMAGE_COVERAGE":counts.get("NO_IMAGE_COVERAGE",0)
 },
 "candidate_results":results,
 "direct_match_source_rows":direct_rows,
 "truth_firewall":{"candidate_coordinates_changed":False,"candidate_radii_changed":False,
   "prediction_map_used":False,"imagery_labels_used_to_select_matches":False},
 "interpretation":"INDEPENDENT_IMAGERY_COVERAGE_AUDIT_COMPLETE",
 "claim_ceiling":"DIRECT_GEOGRAPHIC_IMAGERY_RECORD_SUPPORT_ONLY__NOT_OBJECT_IDENTITY__NEAR_AND_REGIONAL_ARE_NOT_VISUAL_CONFIRMATION"
}
p=OUT/"JANUS-KUSTO-SOUTHWEST-CORNER-SWC3-INDEPENDENT-IMAGERY-COVERAGE-RUN-2026-09-24-v1.0.json"
rawout=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(rawout)
print(json.dumps({
  "artifact_id":out["artifact_id"],"csv_sha256":out["source"]["sha256"],
  "raw_rows":len(rows),"georeferenced_rows":len(pts),"unique_coords":len(unique_pts),
  "coord_binding":coord,"coverage_counts":out["coverage_counts"],
  "direct_matches":[{"candidate_id":x["candidate_id"],"blind_rank":x["blind_rank"],"distance_m":x["distance_m"],"radius_m":x["radius_m"],"source_row_index":x["source_row_index"]} for x in direct_rows],
  "nearest_overall":sorted([{"candidate_id":x["candidate_id"],"blind_rank":x["blind_rank"],"distance_m":x["nearest_imagery_record"]["distance_m"],"radius_m":x["radius_m"],"state":x["coverage_state"]} for x in results if x["nearest_imagery_record"]],key=lambda q:q["distance_m"])[:10],
  "output_json_sha256":hashlib.sha256(rawout.encode()).hexdigest()
},indent=2))
