#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import zipfile
from pathlib import Path

import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHARK-ROCK-G3-TRUTH-UNLOCK-SCORING-PREREG-2026-09-23-v1.0.json").read_text())
G2=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHARK-ROCK-G2B-BLIND-MULTISURVEY-MORPHOLOGY-RECEIPT-2026-09-23-v1.0.json").read_text())

ARTIFACT_ID=int(G2["run"]["artifact_id"])
EXPECTED_ZIP_SHA=G2["run"]["artifact_zip_sha256"]
EXPECTED_JSON_SHA=G2["run"]["candidate_json_sha256"]
TRUTH=PRE["truth"]["coordinate"]
UNC=float(PRE["scoring_geometry"]["grid_center_uncertainty_m"])
UA={"User-Agent":"JANUS-KUSTO-SharkRock-G3/1.0"}

def get_g2b():
    token=os.environ.get("GITHUB_TOKEN")
    repo=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
    if not token:
        raise RuntimeError("GITHUB_TOKEN required")
    url=f"https://api.github.com/repos/{repo}/actions/artifacts/{ARTIFACT_ID}/zip"
    r=requests.get(url,headers={
      **UA,
      "Authorization":f"Bearer {token}",
      "Accept":"application/vnd.github+json",
      "X-GitHub-Api-Version":"2022-11-28"
    },timeout=180,allow_redirects=True)
    r.raise_for_status()
    blob=r.content
    zsha=hashlib.sha256(blob).hexdigest()
    if zsha!=EXPECTED_ZIP_SHA:
        raise RuntimeError(f"G2B artifact ZIP SHA mismatch {zsha} != {EXPECTED_ZIP_SHA}")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        names=[n for n in zf.namelist() if n.lower().endswith(".json") and not n.startswith("__MACOSX/")]
        if len(names)!=1:
            raise RuntimeError(f"Expected exactly one G2B JSON, found {names}")
        raw=zf.read(names[0])
    jsha=hashlib.sha256(raw).hexdigest()
    if jsha!=EXPECTED_JSON_SHA:
        raise RuntimeError(f"G2B JSON SHA mismatch {jsha} != {EXPECTED_JSON_SHA}")
    return json.loads(raw),zsha,jsha

def hav_m(lon1,lat1,lon2,lat2):
    R=6371008.8
    p1=math.radians(lat1);p2=math.radians(lat2)
    dp=p2-p1;dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(min(1.0,math.sqrt(a)))

g2,zsha,jsha=get_g2b()

per_survey={}
support_surveys=[]
for sid in sorted(g2["surveys"]):
    cand=g2["surveys"][sid]["frozen_candidates"]
    scored=[]
    for q in cand:
        d=hav_m(TRUTH["lon"],TRUTH["lat"],float(q["lon"]),float(q["lat"]))
        r=float(q["radius_m"])
        supported=bool(d<=r+UNC)
        scored.append((d,int(q["blind_rank"]),q,supported))
    scored.sort(key=lambda x:(x[0],x[1]))
    nearest=scored[0] if scored else None
    supporting=[x for x in scored if x[3]]
    supporting.sort(key=lambda x:(x[1],x[0]))
    if supporting:
        support_surveys.append(sid)
    n=len(cand)
    def pct(rank):
        if n<=1:return 1.0
        return 1.0-(rank-1)/(n-1)
    per_survey[sid]={
      "candidate_count":n,
      "nearest_candidate_distance_m":nearest[0] if nearest else None,
      "nearest_candidate_id":nearest[2]["candidate_id"] if nearest else None,
      "nearest_candidate_radius_m":nearest[2]["radius_m"] if nearest else None,
      "nearest_candidate_blind_rank":nearest[1] if nearest else None,
      "nearest_candidate_rank_percentile":pct(nearest[1]) if nearest else None,
      "truth_supported_by_any_candidate_patch":bool(supporting),
      "supporting_candidate_count":len(supporting),
      "best_supporting_candidate_id":supporting[0][2]["candidate_id"] if supporting else None,
      "best_supporting_candidate_distance_m":supporting[0][0] if supporting else None,
      "best_supporting_candidate_radius_m":supporting[0][2]["radius_m"] if supporting else None,
      "best_supporting_candidate_blind_rank":supporting[0][1] if supporting else None,
      "best_supporting_candidate_rank_percentile":pct(supporting[0][1]) if supporting else None
    }

clusters=g2.get("persistent_clusters") or []
cluster_scored=[]
for q in clusters:
    d=hav_m(TRUTH["lon"],TRUTH["lat"],float(q["lon"]),float(q["lat"]))
    supported=bool(
      int(q["distinct_survey_count"])>=2 and
      d<=float(q["radius_m"])+UNC
    )
    cluster_scored.append((d,q,supported))
cluster_scored.sort(key=lambda x:x[0])
nearest_cluster=cluster_scored[0] if cluster_scored else None
support_clusters=[x for x in cluster_scored if x[2]]
support_clusters.sort(key=lambda x:(-int(x[1]["distinct_survey_count"]),x[0],x[1]["cluster_id"]))

ns=len(support_surveys)
persistent=bool(support_clusters)
if persistent:
    interpretation="BLIND_CROSS_SURVEY_PERSISTENT_MORPHOLOGY_RECOVERY_OF_ROV_CONFIRMED_NATURAL_FEATURE"
elif ns>=2:
    interpretation="MULTI_EPOCH_LOCAL_RECOVERY_WITHOUT_FROZEN_CLUSTER_MATCH"
elif ns==1:
    interpretation="SINGLE_SURVEY_CASE_LEVEL_RECOVERY_ONLY"
else:
    interpretation="BLIND_DETECTOR_DID_NOT_RECOVER_THIS_KNOWN_NATURAL_ANOMALY_AT_FROZEN_SCALES"

out={
 "artifact_id":"JANUS-KUSTO-SHARK-ROCK-G3-TRUTH-UNLOCK-SCORING-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "g2b_authority":{
   "receipt":G2["artifact_id"],
   "artifact_id":ARTIFACT_ID,
   "artifact_zip_sha256_verified":zsha,
   "candidate_json_sha256_verified":jsha,
   "candidate_total":int(g2["candidate_total"]),
   "persistent_cluster_count":int(g2["persistent_cluster_count"])
 },
 "truth":{
   "benchmark_id":PRE["truth"]["benchmark_id"],
   "coordinate":TRUTH,
   "semantic_truth":PRE["truth"]["semantic_truth"],
   "truth_source_role":PRE["truth"]["truth_source_role"]
 },
 "scoring_geometry":{
   "grid_center_uncertainty_m":UNC,
   "candidate_support_rule":PRE["scoring_geometry"]["candidate_support_rule"],
   "persistent_cluster_support_rule":PRE["scoring_geometry"]["persistent_cluster_support_rule"]
 },
 "per_survey":per_survey,
 "cross_survey":{
   "number_of_surveys_with_truth_supported":ns,
   "survey_ids_with_truth_supported":sorted(support_surveys),
   "nearest_persistent_cluster_distance_m":nearest_cluster[0] if nearest_cluster else None,
   "nearest_persistent_cluster_id":nearest_cluster[1]["cluster_id"] if nearest_cluster else None,
   "nearest_persistent_cluster_radius_m":nearest_cluster[1]["radius_m"] if nearest_cluster else None,
   "nearest_persistent_cluster_distinct_survey_count":nearest_cluster[1]["distinct_survey_count"] if nearest_cluster else None,
   "truth_supported_by_any_persistent_cluster":persistent,
   "supporting_persistent_cluster_count":len(support_clusters),
   "best_supporting_persistent_cluster":({
      "cluster_id":support_clusters[0][1]["cluster_id"],
      "distance_m":support_clusters[0][0],
      "radius_m":support_clusters[0][1]["radius_m"],
      "distinct_survey_count":support_clusters[0][1]["distinct_survey_count"],
      "surveys":support_clusters[0][1]["surveys"]
   } if support_clusters else None)
 },
 "interpretation":interpretation,
 "semantic_firewall":{
   "ground_truth_is_natural":True,
   "artificiality_promotion_allowed":False,
   "morphology_recovery_if_any_means":"anomaly sensitivity and repeatability only",
   "pre_dive_shipwreck_hypothesis_status":"REJECTED_BY_ROV_GROUNDTRUTH"
 },
 "hard_rules_preserved":[
   "G2B_OUTPUT_HASH_VERIFIED_AND_IMMUTABLE",
   "NO_CANDIDATE_ADDITION_DELETION_OR_RECENTERING",
   "NO_RADIUS_CHANGE",
   "NO_THRESHOLD_RETUNING",
   "PREDECLARED_SUPPORT_GEOMETRY_USED"
 ],
 "claim_ceiling":PRE["claim_ceiling"]
}
p=OUT/"JANUS-KUSTO-SHARK-ROCK-G3-TRUTH-UNLOCK-SCORING-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "per_survey":{k:{
   "nearest_distance_m":v["nearest_candidate_distance_m"],
   "nearest_rank":v["nearest_candidate_blind_rank"],
   "supported":v["truth_supported_by_any_candidate_patch"],
   "best_support_rank":v["best_supporting_candidate_blind_rank"]
 } for k,v in per_survey.items()},
 "surveys_supporting_truth":sorted(support_surveys),
 "supporting_survey_count":ns,
 "nearest_persistent_cluster_distance_m":out["cross_survey"]["nearest_persistent_cluster_distance_m"],
 "truth_supported_by_persistent_cluster":persistent,
 "best_supporting_persistent_cluster":out["cross_survey"]["best_supporting_persistent_cluster"],
 "interpretation":interpretation
},indent=2))
