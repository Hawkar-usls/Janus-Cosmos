#!/usr/bin/env python3
from __future__ import annotations
import csv, io, json, math, time
from pathlib import Path
from urllib.parse import quote

import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace/kusto_global_anomaly_out"; OUT.mkdir(parents=True,exist_ok=True)
RECEIPT=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-GLOBAL-STAGE1-BLIND-MORPHOLOGY-RECEIPT-2026-09-23-v1.0.json").read_text())
CANDS=RECEIPT["region2"]["candidates"]
BASE="https://www.ncei.noaa.gov/erddap/tabledap/deep_sea_corals.csv"
UA={"User-Agent":"JANUS-KUSTO-independent-photo-coverage/1.0"}
R=6371008.8

VARS=[
 "DatasetID","SurveyID","EventID","latitude","longitude","DepthInMeters",
 "ObservationDate","VehicleName","SamplingEquipment","RecordType","Locality",
 "ImageURL","HighlightImageURL"
]

def hav(lat1,lon1,lat2,lon2):
    p1,p2=math.radians(lat1),math.radians(lat2)
    dp=p2-p1;dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(min(1.0,math.sqrt(a)))

def get(url):
    last=None
    for i in range(7):
        try:
            r=requests.get(url,headers=UA,timeout=120)
            if r.status_code==429:
                time.sleep(min(30,2**i));continue
            r.raise_for_status();return r
        except Exception as e:
            last=e
            if i==6: raise
            time.sleep(min(30,2**i))
    raise last

query=",".join(VARS)
constraints=[
 'latitude>=20.9','latitude<=23.6',
 'longitude>=-162.6','longitude<=-157.4'
]
url=BASE+"?"+query+"&"+"&".join(quote(x,safe="=><-._") for x in constraints)
r=get(url)
text=r.text
rows=list(csv.DictReader(io.StringIO(text)))

qualified=[]
for row in rows:
    try:
        lat=float(row["latitude"]);lon=float(row["longitude"])
    except: continue
    img=(row.get("ImageURL") or "").strip()
    himg=(row.get("HighlightImageURL") or "").strip()
    if not img and not himg: continue
    sid=(row.get("SurveyID") or "").strip()
    if sid in {"FK140307","ZHNG09RR"}: continue
    rec=dict(row)
    rec["latitude"]=lat;rec["longitude"]=lon
    qualified.append(rec)

def classify(d,radius):
    if d<=radius:return "DIRECT_IMAGE_SUPPORT"
    if d<=5000:return "NEAR_IMAGE_COVERAGE"
    if d<=25000:return "REGIONAL_IMAGE_COVERAGE"
    return "NO_IMAGE_COVERAGE"

results=[]
for c in CANDS:
    hits=[]
    for rec in qualified:
        d=hav(c["lat"],c["lon"],rec["latitude"],rec["longitude"])
        if d<=25000:
            x=dict(rec);x["distance_m"]=d
            hits.append(x)
    hits.sort(key=lambda x:x["distance_m"])
    nearest=None
    nearest_class="NO_IMAGE_COVERAGE"
    if qualified:
        nearest=min(qualified,key=lambda rec:hav(c["lat"],c["lon"],rec["latitude"],rec["longitude"]))
        nearest=dict(nearest)
        nearest["distance_m"]=hav(c["lat"],c["lon"],nearest["latitude"],nearest["longitude"])
        nearest_class=classify(nearest["distance_m"],c["radius_m"])
    results.append({
      "candidate_id":c["id"],"lon":c["lon"],"lat":c["lat"],"radius_m":c["radius_m"],
      "nearest":nearest,
      "distance_class":nearest_class,
      "qualifying_images_within_25km_count":len(hits),
      "qualifying_images_within_25km":hits[:500]
    })

counts={}
for x in results: counts[x["distance_class"]]=counts.get(x["distance_class"],0)+1

out={
 "artifact_id":"JANUS-KUSTO-GLOBAL-20-INDEPENDENT-PHOTO-COVERAGE-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-GLOBAL-20-INDEPENDENT-PHOTO-COVERAGE-PREREG-2026-09-23-v1.0.json",
 "source":{
   "dataset":"NOAA NCEI ERDDAP deep_sea_corals",
   "endpoint":BASE,
   "query_extent":{"lat_min":20.9,"lat_max":23.6,"lon_min":-162.6,"lon_max":-157.4},
   "regional_rows_returned":len(rows),
   "qualifying_georeferenced_image_records":len(qualified)
 },
 "counts_by_distance_class":counts,
 "results":results,
 "claim_ceiling":"INDEPENDENT_VISUAL_COVERAGE_AUDIT_ONLY"
}
p=OUT/"JANUS-KUSTO-GLOBAL-20-INDEPENDENT-PHOTO-COVERAGE-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "regional_rows_returned":len(rows),
 "qualifying_georeferenced_image_records":len(qualified),
 "counts_by_distance_class":counts,
 "summary":[{
   "candidate_id":x["candidate_id"],"distance_class":x["distance_class"],
   "nearest_distance_m":None if x["nearest"] is None else x["nearest"]["distance_m"],
   "nearest_survey":None if x["nearest"] is None else x["nearest"].get("SurveyID"),
   "nearest_event":None if x["nearest"] is None else x["nearest"].get("EventID"),
   "nearest_vehicle":None if x["nearest"] is None else x["nearest"].get("VehicleName"),
   "nearest_image":None if x["nearest"] is None else (x["nearest"].get("HighlightImageURL") or x["nearest"].get("ImageURL")),
   "within25km":x["qualifying_images_within_25km_count"]
 } for x in results]
},indent=2))
