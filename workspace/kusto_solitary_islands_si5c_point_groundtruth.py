#!/usr/bin/env python3
from __future__ import annotations
import collections, hashlib, io, json, math, os, tempfile, zipfile
from pathlib import Path

import requests, shapefile
from pyproj import CRS, Transformer

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI5-POSTRESULT-GROUNDTRUTH-CHARACTERIZATION-PREREG-2026-09-24-v1.0.json").read_text())
IMP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI5C-POINT-GROUNDTRUTH-IMPLEMENTATION-ADDENDUM-2026-09-24-v1.0.json").read_text())
AUTH=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI1-BLIND-BATHYMETRY-RECEIPT-2026-09-24-v1.0.json").read_text())
INV=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI5A-TRUTH-PACKAGE-INVENTORY-RECEIPT-2026-09-24-v1.0.json").read_text())

UA={"User-Agent":"JANUS-KUSTO-SolitaryIslands-SI5C-point-truth/1.0"}

def fetch_candidate_json():
    token=os.environ.get("GITHUB_TOKEN")
    repo=os.environ.get("GITHUB_REPOSITORY","hawkar-usls/Janus-Cosmos")
    aid=int(AUTH["run"]["artifact_id"])
    u=f"https://api.github.com/repos/{repo}/actions/artifacts/{aid}/zip"
    r=requests.get(u,headers={"User-Agent":UA["User-Agent"],"Authorization":f"Bearer {token}","Accept":"application/vnd.github+json"},timeout=120)
    r.raise_for_status(); blob=r.content
    zsha=hashlib.sha256(blob).hexdigest()
    if zsha!=AUTH["run"]["artifact_zip_sha256"]:raise RuntimeError("candidate artifact SHA mismatch")
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        js=[n for n in zf.namelist() if n.lower().endswith(".json")]
        if len(js)!=1:raise RuntimeError(f"expected one candidate json, got {js}")
        raw=zf.read(js[0])
    jsha=hashlib.sha256(raw).hexdigest()
    if jsha!=AUTH["run"]["candidate_json_sha256"]:raise RuntimeError("candidate JSON SHA mismatch")
    return json.loads(raw),zsha,jsha

def read_truth_layers():
    url=PRE["truth_source"]["data_package_url"]
    r=requests.get(url,headers=UA,timeout=240,allow_redirects=True);r.raise_for_status()
    blob=r.content
    sha=hashlib.sha256(blob).hexdigest()
    if sha!=INV["truth_package"]["zip_sha256"]:raise RuntimeError("truth package SHA mismatch")
    layers={
      "video":"SIMP_TowedVideoSubClass/SIMP_TowedVideo_SubstrateClass",
      "sediment":"SIMP_SedimentsMetadata/SIMP_Sediments_SamplesRetrieved"
    }
    out={}
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
      with tempfile.TemporaryDirectory() as td:
        for role,stem in layers.items():
          base=Path(td)/role
          for ext in (".shp",".shx",".dbf",".prj",".cpg"):
            n=stem+ext
            if n in zf.namelist():(base.with_suffix(ext)).write_bytes(zf.read(n))
          sf=shapefile.Reader(str(base.with_suffix(".shp")))
          prj=(base.with_suffix(".prj")).read_text(errors="replace")
          crs=CRS.from_wkt(prj)
          fields=[f[0] for f in sf.fields[1:]]
          rows=[]
          for sr in sf.iterShapeRecords():
            rec={fields[i]:v for i,v in enumerate(sr.record)}
            pt=sr.shape.points[0]
            rows.append({"x":float(pt[0]),"y":float(pt[1]),"record":rec})
          out[role]={"crs":crs,"rows":rows}
    return out,sha

candj,artifact_sha,candidate_sha=fetch_candidate_json()
truth,truth_sha=read_truth_layers()
candidates=candj["frozen_candidates"]
cand_crs=CRS.from_user_input(candj["input"]["crs"])
tv=Transformer.from_crs(truth["video"]["crs"],cand_crs,always_xy=True)
ts=Transformer.from_crs(truth["sediment"]["crs"],cand_crs,always_xy=True)

video=[]
for q in truth["video"]["rows"]:
    x,y=tv.transform(q["x"],q["y"])
    video.append((float(x),float(y),q["record"]))
sed=[]
for q in truth["sediment"]["rows"]:
    x,y=ts.transform(q["x"],q["y"])
    sed.append((float(x),float(y),q["record"]))

results=[]
video_covered=0;sed_covered=0
for c in candidates:
    cx,cy=float(c["x_m"]),float(c["y_m"]);rad=float(c["radius_m"])
    vd=[]
    for x,y,rec in video:
        d=math.hypot(x-cx,y-cy)
        vd.append((d,rec))
    vd.sort(key=lambda z:z[0])
    vinside=[(d,r) for d,r in vd if d<=rad]
    sd=[]
    for x,y,rec in sed:
        d=math.hypot(x-cx,y-cy)
        sd.append((d,rec))
    sd.sort(key=lambda z:z[0])
    sinside=[(d,r) for d,r in sd if d<=rad]

    if vinside:video_covered+=1
    if sinside:sed_covered+=1
    labels=collections.Counter(str(r.get("PrimarySub")) for _,r in vinside if r.get("PrimarySub") not in (None,""))
    transects=sorted({str(r.get("TransectID")) for _,r in vinside if r.get("TransectID") not in (None,"")})
    results.append({
      "candidate_id":c["candidate_id"],"blind_rank":c["blind_rank"],"radius_m":c["radius_m"],
      "video_observation_covered":bool(vinside),
      "video_observation_count_inside_radius":len(vinside),
      "unique_video_transects_inside_radius":len(transects),
      "transect_ids_inside_radius":transects,
      "video_primary_substrate_counts":dict(labels),
      "nearest_video_observation_distance_m":vd[0][0] if vd else None,
      "sediment_sample_covered":bool(sinside),
      "sediment_sample_count_inside_radius":len(sinside),
      "nearest_sediment_sample_distance_m":sd[0][0] if sd else None,
      "sediment_samples_inside_radius":[{"distance_m":d,"record":r} for d,r in sinside]
    })

video_dist=sorted([x["nearest_video_observation_distance_m"] for x in results if x["nearest_video_observation_distance_m"] is not None])
sed_dist=sorted([x["nearest_sediment_sample_distance_m"] for x in results if x["nearest_sediment_sample_distance_m"] is not None])
def quant(vals,p):
    if not vals:return None
    i=(len(vals)-1)*p;lo=int(math.floor(i));hi=int(math.ceil(i))
    if lo==hi:return vals[lo]
    return vals[lo]*(hi-i)+vals[hi]*(i-lo)

agg_labels=collections.Counter()
for x in results:
    agg_labels.update(x["video_primary_substrate_counts"])

out={
 "artifact_id":"JANUS-KUSTO-SOLITARY-ISLANDS-SI5C-FROZEN-CANDIDATE-POINT-GROUNDTRUTH-CHARACTERIZATION-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],"implementation_addendum":IMP["artifact_id"],
 "formal_parent":"JANUS-KUSTO-SOLITARY-ISLANDS-SI3-FORMAL-UNSEEN-V2-FAIL-RECEIPT-2026-09-24-v1.0",
 "candidate_authority":{"artifact_zip_sha256_verified":artifact_sha,"candidate_json_sha256_verified":candidate_sha,"candidate_count":len(candidates),"candidate_crs":cand_crs.to_string()},
 "truth_authority":{"truth_package_sha256_verified":truth_sha,"video_observation_count":len(video),"sediment_sample_count":len(sed)},
 "coverage":{
   "video_observation_covered_candidates":video_covered,
   "video_observation_coverage_fraction":video_covered/len(candidates),
   "sediment_sample_covered_candidates":sed_covered,
   "sediment_sample_coverage_fraction":sed_covered/len(candidates),
   "nearest_video_distance_m":{"min":min(video_dist),"median":quant(video_dist,.5),"q90":quant(video_dist,.9),"max":max(video_dist)} if video_dist else None,
   "nearest_sediment_distance_m":{"min":min(sed_dist),"median":quant(sed_dist,.5),"q90":quant(sed_dist,.9),"max":max(sed_dist)} if sed_dist else None,
   "video_primary_substrate_counts_across_covered_candidates":dict(agg_labels)
 },
 "candidate_results":results,
 "formal_si3_result_changed":False,
 "promotion":False,
 "interpretation":"POSTRESULT_SPATIAL_GROUNDTRUTH_CHARACTERIZATION_ONLY",
 "claim_ceiling":"POSTRESULT_POINT_GROUNDTRUTH_CHARACTERIZATION_ONLY__CANNOT_PROMOTE_SI3_FAIL"
}
p=OUT/"JANUS-KUSTO-SOLITARY-ISLANDS-SI5C-FROZEN-CANDIDATE-POINT-GROUNDTRUTH-CHARACTERIZATION-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False,default=str);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],"candidate_count":len(candidates),
 "video_observation_covered_candidates":video_covered,"video_coverage_fraction":video_covered/len(candidates),
 "sediment_sample_covered_candidates":sed_covered,"sediment_coverage_fraction":sed_covered/len(candidates),
 "nearest_video_distance_m":out["coverage"]["nearest_video_distance_m"],
 "nearest_sediment_distance_m":out["coverage"]["nearest_sediment_distance_m"],
 "video_primary_substrate_counts":dict(agg_labels),
 "formal_si3_result_changed":False,"promotion":False,
 "output_json_sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2,ensure_ascii=False,default=str))
