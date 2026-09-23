#!/usr/bin/env python3
from __future__ import annotations
import csv,hashlib,io,json,math,os,zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from pyproj import Transformer, Geod

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM4-AUV-COORDINATE-COVERAGE-PREREG-2026-09-24-v1.0.json").read_text())
EM3=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM3-FORMAL-DUAL-V2-V3-RESULT-RECEIPT-2026-09-24-v1.0.json").read_text())
CAND=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-ELIZABETH-MIDDLETON-EM1-BLIND-BATHYMETRY-RECEIPT-2026-09-24-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-ElizabethMiddleton-EM4-AUV-coordinate/1.0"}
SQUIDLE="https://squidle.org/api/annotation/export"
GEOD=Geod(ellps="WGS84")
TO_WGS=Transformer.from_crs("EPSG:32757","EPSG:4326",always_xy=True)

def get(url,headers=None,params=None,timeout=300):
    h=dict(UA)
    if headers:h.update(headers)
    r=requests.get(url,headers=h,params=params,timeout=timeout,allow_redirects=True)
    r.raise_for_status();return r

def github_artifact_json(aid,expected_zip_sha,expected_json_sha=None):
    token=os.environ.get("GITHUB_TOKEN")
    repo=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
    if not token:raise RuntimeError("GITHUB_TOKEN required")
    u=f"https://api.github.com/repos/{repo}/actions/artifacts/{aid}/zip"
    blob=get(u,headers={"Authorization":f"Bearer {token}","Accept":"application/vnd.github+json","X-GitHub-Api-Version":"2022-11-28"}).content
    zsha=hashlib.sha256(blob).hexdigest()
    if zsha!=expected_zip_sha:raise RuntimeError(f"artifact sha mismatch {aid}")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        js=[n for n in zf.namelist() if n.lower().endswith(".json")]
        if len(js)!=1:raise RuntimeError(f"artifact {aid}: expected one JSON {js}")
        raw=zf.read(js[0])
    jsha=hashlib.sha256(raw).hexdigest()
    if expected_json_sha and jsha!=expected_json_sha:raise RuntimeError(f"artifact json sha mismatch {aid}")
    return json.loads(raw),zsha,jsha

candj,cz,cj=github_artifact_json(
 int(CAND["run"]["artifact_id"]),
 CAND["run"]["artifact_zip_sha256"],
 CAND["run"]["candidate_json_sha256"]
)
em2j,ez,ej=github_artifact_json(
 int(EM3["run"]["artifact_id"]),
 EM3["run"]["artifact_zip_sha256"],
 EM3["run"]["output_json_sha256"]
)

v2_supported={
 x["candidate_id"] for x in em2j["methods"]["v2"]["candidate_results"]
 if x.get("cross_channel_support") is True
}

filters={"filters":[{"name":"annotation_set","op":"has","val":{"name":"usergroups","op":"any","val":{"name":"id","op":"eq","val":346}}}]}
columns=PRE["source"]["coordinate_only_columns"]
fops={"operations":[{"module":"pandas","method":"json_normalize"},{"method":"sort_index","kwargs":{"axis":1}}]}
params={
 "q":json.dumps(filters,separators=(",",":")),
 "template":"dataframe.csv",
 "disposition":"attachment",
 "nologin":"true",
 "include_columns":json.dumps(columns,separators=(",",":")),
 "f":json.dumps(fops,separators=(",",":"))
}
r=get(SQUIDLE,params=params,timeout=600)
raw=r.content
csv_sha=hashlib.sha256(raw).hexdigest()
df=pd.read_csv(io.BytesIO(raw),low_memory=False)
for col in ["point.pose.lat","point.pose.lon"]:
    if col not in df.columns:raise RuntimeError(f"missing coordinate column {col}; got {list(df.columns)}")
if any("label" in str(c).lower() for c in df.columns):
    raise RuntimeError(f"label column unexpectedly returned: {list(df.columns)}")

df["point.pose.lat"]=pd.to_numeric(df["point.pose.lat"],errors="coerce")
df["point.pose.lon"]=pd.to_numeric(df["point.pose.lon"],errors="coerce")
df=df[np.isfinite(df["point.pose.lat"]) & np.isfinite(df["point.pose.lon"])].copy()
df=df[(df["point.pose.lat"].between(-90,90)) & (df["point.pose.lon"].between(-180,180))].copy()

def sval(row,col):
    v=row[col] if col in row.index else ""
    return "" if pd.isna(v) else str(v)
keys=[]
for _,row in df.iterrows():
    mk=sval(row,"point.media.key");pt=sval(row,"point.pose.timestamp")
    if mk or pt:
        key="PRIMARY|"+mk+"|"+pt
    else:
        dep=sval(row,"point.media.deployment.name");mt=sval(row,"point.media.timestamp_start")
        key=f"FALLBACK|{dep}|{float(row['point.pose.lat']):.7f}|{float(row['point.pose.lon']):.7f}|{mt}"
    keys.append(key)
df["_dedup_key"]=keys
df=df.drop_duplicates("_dedup_key",keep="first").reset_index(drop=True)

lats=df["point.pose.lat"].to_numpy(dtype=float)
lons=df["point.pose.lon"].to_numpy(dtype=float)

def nearest_for_candidate(c):
    lon,lat=TO_WGS.transform(float(c["x_m"]),float(c["y_m"]))
    # vectorized geodesic distances
    _,_,dist=GEOD.inv(np.full_like(lons,lon),np.full_like(lats,lat),lons,lats)
    dist=np.asarray(dist,dtype=float)
    idx=int(np.argmin(dist));nearest=float(dist[idx])
    radius=float(c["radius_m"])
    direct_idx=np.where(dist<=radius)[0]
    near_idx=np.where(dist<=2*radius)[0]
    deps=[]
    for i in direct_idx:
        col="point.media.deployment.name"
        if col in df.columns and not pd.isna(df.iloc[i][col]):
            deps.append(str(df.iloc[i][col]))
    nearest_row=df.iloc[idx]
    return {
      "candidate_id":c["candidate_id"],"blind_rank":c["blind_rank"],"region_id":c["region_id"],
      "radius_m":radius,"candidate_lon":lon,"candidate_lat":lat,
      "v2_crosschannel_support":c["candidate_id"] in v2_supported,
      "nearest_annotation_pose_distance_m":nearest,
      "direct_auv_annotation_coverage":bool(len(direct_idx)),
      "near_auv_annotation_coverage_descriptive":bool(len(near_idx)),
      "unique_annotation_poses_within_radius":int(len(direct_idx)),
      "unique_annotation_poses_within_2r":int(len(near_idx)),
      "unique_deployments_within_radius":sorted(set(deps)),
      "nearest_pose":{
        "lat":float(nearest_row["point.pose.lat"]),
        "lon":float(nearest_row["point.pose.lon"]),
        "media_key":sval(nearest_row,"point.media.key"),
        "pose_timestamp":sval(nearest_row,"point.pose.timestamp"),
        "media_timestamp_start":sval(nearest_row,"point.media.timestamp_start"),
        "deployment":sval(nearest_row,"point.media.deployment.name"),
        "campaign":sval(nearest_row,"point.media.deployment.campaign.name")
      }
    }

results=[nearest_for_candidate(c) for c in candj["frozen_candidates"]]
direct=[x for x in results if x["direct_auv_annotation_coverage"]]
supported=[x for x in results if x["v2_crosschannel_support"]]
supported_direct=[x for x in supported if x["direct_auv_annotation_coverage"]]
near=[x for x in results if x["near_auv_annotation_coverage_descriptive"]]

out={
 "artifact_id":"JANUS-KUSTO-ELIZABETH-MIDDLETON-EM4-AUV-COORDINATE-COVERAGE-RUN-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "candidate_authority":{"artifact_id":CAND["run"]["artifact_id"],"artifact_zip_sha256_verified":cz,"candidate_json_sha256_verified":cj},
 "acoustic_authority":{"artifact_id":EM3["run"]["artifact_id"],"artifact_zip_sha256_verified":ez,"output_json_sha256_verified":ej},
 "auv_coordinate_source":{
   "doi":PRE["source"]["doi"],"endpoint":r.url,"http_status":r.status_code,
   "csv_bytes":len(raw),"csv_sha256":csv_sha,"returned_columns":list(df.columns.drop("_dedup_key",errors="ignore")),
   "raw_rows_with_valid_coordinates_before_dedup":None,
   "unique_georeferenced_media_poses":int(len(df)),
   "labels_requested":False,"label_columns_present":False
 },
 "coverage_summary":{
   "frozen_candidate_count":len(results),
   "direct_candidate_count":len(direct),
   "near_2r_candidate_count":len(near),
   "v2_crosschannel_supported_candidate_count":len(supported),
   "v2_supported_with_direct_auv_annotation_coverage":len(supported_direct),
   "direct_coverage_fraction_all":len(direct)/len(results) if results else None,
   "direct_coverage_fraction_v2_supported":len(supported_direct)/len(supported) if supported else None
 },
 "candidate_results":results,
 "interpretation":"DIRECT_GEOREFERENCED_AUV_ANNOTATION_COVERAGE_FOUND" if direct else "NO_DIRECT_AUV_ANNOTATION_POINT_INTERSECTIONS",
 "hard_limitations":[
   "ANNOTATION_POINTS_ARE_NOT_COMPLETE_AUV_TRACK_GEOMETRY",
   "NO_LABEL_OR_IMAGE_CONTENT_WAS_READ",
   "NO_INTERSECTION_DOES_NOT_IMPLY_NO_AUV_IMAGERY_EXISTS"
 ],
 "next_gate":"EM5_LABEL_AND_MEDIA_INTERPRETATION_ONLY_FOR_FROZEN_DIRECT_INTERSECTIONS" if direct else "RECOVER_FULL_AUV_TRACK_OR_MEDIA_FOOTPRINT_BEFORE_VISUAL_ABSENCE_CLAIM",
 "claim_ceiling":"DIRECT_GEOREFERENCED_AUV_ANNOTATION_COVERAGE_ONLY__NOT_VISUAL_OBJECT_CONFIRMATION"
}
# Record count before dedup separately after creating final df.
out["auv_coordinate_source"]["raw_rows_with_valid_coordinates_before_dedup"]=int(len(keys))
p=OUT/"JANUS-KUSTO-ELIZABETH-MIDDLETON-EM4-AUV-COORDINATE-COVERAGE-RUN-2026-09-24-v1.0.json"
rawout=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(rawout)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "unique_georeferenced_media_poses":len(df),
 "direct_candidate_count":len(direct),
 "near_2r_candidate_count":len(near),
 "v2_supported_candidate_count":len(supported),
 "v2_supported_direct_count":len(supported_direct),
 "direct_candidates":[{"candidate_id":x["candidate_id"],"rank":x["blind_rank"],"region":x["region_id"],"radius_m":x["radius_m"],"nearest_m":x["nearest_annotation_pose_distance_m"],"poses_inside":x["unique_annotation_poses_within_radius"],"v2_support":x["v2_crosschannel_support"],"deployments":x["unique_deployments_within_radius"]} for x in direct],
 "output_json_sha256":hashlib.sha256(rawout.encode()).hexdigest()
},indent=2,ensure_ascii=False))
