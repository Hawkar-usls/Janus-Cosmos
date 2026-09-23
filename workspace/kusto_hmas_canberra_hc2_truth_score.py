#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,math,os,zipfile
from pathlib import Path
import requests
from pyproj import CRS, Transformer

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-HMAS-CANBERRA-HC0-KNOWN-WRECK-RECOVERY-PREREG-2026-09-23-v1.0.json").read_text())
AUTH=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-HMAS-CANBERRA-HC1-BLIND-BATHYMETRY-RECEIPT-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-HMAS-Canberra-HC2-truth-score/1.0"}

def get_artifact():
    token=os.environ.get("GITHUB_TOKEN")
    repo=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
    aid=int(AUTH["run"]["artifact_id"])
    u=f"https://api.github.com/repos/{repo}/actions/artifacts/{aid}/zip"
    r=requests.get(u,headers={"User-Agent":UA["User-Agent"],"Authorization":f"Bearer {token}","Accept":"application/vnd.github+json"},timeout=120)
    r.raise_for_status();blob=r.content
    zsha=hashlib.sha256(blob).hexdigest()
    if zsha!=AUTH["run"]["artifact_zip_sha256"]:raise RuntimeError("artifact zip SHA mismatch")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        js=[n for n in zf.namelist() if n.lower().endswith(".json")]
        if len(js)!=1:raise RuntimeError(f"expected one JSON {js}")
        raw=zf.read(js[0])
    jsha=hashlib.sha256(raw).hexdigest()
    if jsha!=AUTH["run"]["candidate_json_sha256"]:raise RuntimeError("candidate JSON SHA mismatch")
    return json.loads(raw),zsha,jsha

cand,zsha,jsha=get_artifact()
crs_text=cand["input"]["crs"]
if not crs_text:raise RuntimeError("candidate raster CRS missing")
crs=CRS.from_user_input(crs_text)
tx=Transformer.from_crs("EPSG:4326",crs,always_xy=True)
lon=float(PRE["truth_unlock"]["anchor_lon"]);lat=float(PRE["truth_unlock"]["anchor_lat"])
truth_x,truth_y=tx.transform(lon,lat)
rows=[]
for c in cand["frozen_candidates"]:
    d=math.hypot(float(c["x_m"])-truth_x,float(c["y_m"])-truth_y)
    rows.append({"candidate_id":c["candidate_id"],"blind_rank":int(c["blind_rank"]),"radius_m":c["radius_m"],"distance_to_official_dive_site_anchor_m":d})
rows.sort(key=lambda q:q["distance_to_official_dive_site_anchor_m"])
nearest=rows[0] if rows else None
thr=float(PRE["truth_unlock"]["case_level_recovery_radius_m"])
inside=[q for q in rows if q["distance_to_official_dive_site_anchor_m"]<=thr]
out={
 "artifact_id":"JANUS-KUSTO-HMAS-CANBERRA-HC2-TRUTH-UNLOCK-RECOVERY-SCORE-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "candidate_authority":{"receipt":AUTH["artifact_id"],"artifact_zip_sha256_verified":zsha,"candidate_json_sha256_verified":jsha},
 "truth":{"authority":PRE["truth_unlock"]["authority"],"anchor_lat":lat,"anchor_lon":lon,"anchor_is_not_wreck_centroid":True,"projected_crs":crs.to_string(),"projected_x":truth_x,"projected_y":truth_y,"case_level_radius_m":thr},
 "score":{"nearest_candidate":nearest,"candidate_count_within_radius":len(inside),"candidates_within_radius":inside,"case_level_recovery":bool(inside)},
 "all_candidate_distances":rows,
 "interpretation":"CASE_LEVEL_KNOWN_OBJECT_RECOVERY_POSITIVE_CONTROL_PASS" if inside else "CASE_LEVEL_KNOWN_OBJECT_RECOVERY_POSITIVE_CONTROL_MISS",
 "limitations":["TARGET_CENTRIC_SURVEY","OFFICIAL_DIVE_SITE_ANCHOR_IS_NOT_WRECK_CENTROID","NO_POPULATION_LOCALIZATION_CLAIM"],
 "claim_ceiling":"CASE_LEVEL_KNOWN_OBJECT_RECOVERY_POSITIVE_CONTROL_ONLY"
}
p=OUT/"JANUS-KUSTO-HMAS-CANBERRA-HC2-TRUTH-UNLOCK-RECOVERY-SCORE-2026-09-23-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({"artifact_id":out["artifact_id"],"nearest":nearest,"within_200m":len(inside),"case_level_recovery":bool(inside),"interpretation":out["interpretation"],"output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()},indent=2))
