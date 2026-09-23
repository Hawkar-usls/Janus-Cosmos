#!/usr/bin/env python3
from __future__ import annotations
import hashlib, io, json, os, zipfile
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)

SWC2_ART=10779120528
SWC3_ART=10778503754
SWC2_SHA="a3ae4db4e153bb7a455998d4f21b7766b42ccb9e365da6490d5d71cd8a5f747b"
SWC3_SHA="1106e9adf6ef621507e83f1b184d88e99f3d2eca17bdd13079e6cfe1773ab36b"
UA={"User-Agent":"JANUS-KUSTO-SWC-posthoc-acoustic-imagery-intersection/1.0"}

def artifact_json(aid,expected_sha):
    token=os.environ["GITHUB_TOKEN"]
    repo=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
    r=requests.get(f"https://api.github.com/repos/{repo}/actions/artifacts/{aid}/zip",headers={
      **UA,"Authorization":f"Bearer {token}","Accept":"application/vnd.github+json","X-GitHub-Api-Version":"2022-11-28"
    },timeout=120);r.raise_for_status()
    blob=r.content;sha=hashlib.sha256(blob).hexdigest()
    if sha!=expected_sha:raise RuntimeError(f"artifact {aid} SHA mismatch")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        js=[n for n in zf.namelist() if n.lower().endswith(".json")]
        if len(js)!=1:raise RuntimeError(f"artifact {aid}: expected one JSON, got {js}")
        raw=zf.read(js[0])
    return json.loads(raw),hashlib.sha256(raw).hexdigest()

swc2,swc2_json_sha=artifact_json(SWC2_ART,SWC2_SHA)
swc3,swc3_json_sha=artifact_json(SWC3_ART,SWC3_SHA)
if swc2_json_sha!="68e5cdce48270443cb94ed299454eb6eadd672f4a8283316feb2ec67bcd85c4d":
    raise RuntimeError("SWC2 JSON SHA mismatch")
if swc3_json_sha!="6aa7ea1dccb0bd922cace0cb2137995fba5ead33ac2f3c668add1c9c943f7d44":
    raise RuntimeError("SWC3 JSON SHA mismatch")

coverage={x["candidate_id"]:x for x in swc3["candidate_results"]}
supported=[x for x in swc2["candidate_results"] if x.get("cross_channel_support") is True]
rows=[]
for x in supported:
    c=coverage[x["candidate_id"]]
    rows.append({
      "candidate_id":x["candidate_id"],
      "blind_rank":x["blind_rank"],
      "radius_m":x["radius_m"],
      "aggregate_score":x["aggregate_score"],
      "backscatter_support":True,
      "imagery_coverage_state":c["coverage_state"],
      "nearest_imagery_distance_m":c["nearest_imagery_record"]["distance_m"] if c.get("nearest_imagery_record") else None
    })
rows.sort(key=lambda q:q["blind_rank"])
counts={}
for x in rows:counts[x["imagery_coverage_state"]]=counts.get(x["imagery_coverage_state"],0)+1

out={
 "artifact_id":"JANUS-KUSTO-SOUTHWEST-CORNER-POSTHOC-SWC2-SUPPORTED-VS-IMAGERY-COVERAGE-2026-09-24-v1.0",
 "status":"READ_ONLY_DESCRIPTIVE_INTERSECTION",
 "authorities":{
   "swc2_artifact_id":SWC2_ART,"swc2_json_sha256":swc2_json_sha,
   "swc3_artifact_id":SWC3_ART,"swc3_json_sha256":swc3_json_sha
 },
 "supported_candidate_count":len(rows),
 "coverage_counts_among_supported":counts,
 "supported_candidates":rows,
 "scientific_promotion":False,
 "interpretation":"Descriptive only. This intersection was computed after both acoustic and imagery outcomes were frozen and cannot promote the primary result.",
 "claim_ceiling":"POSTHOC_DESCRIPTIVE_COVERAGE_OF_ALREADY_SUPPORTED_CANDIDATES_ONLY"
}
p=OUT/"JANUS-KUSTO-SOUTHWEST-CORNER-POSTHOC-SWC2-SUPPORTED-VS-IMAGERY-COVERAGE-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2);p.write_text(raw)
print(json.dumps({**{k:out[k] for k in ["artifact_id","supported_candidate_count","coverage_counts_among_supported"]},"supported_candidates":rows,"output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()},indent=2))
