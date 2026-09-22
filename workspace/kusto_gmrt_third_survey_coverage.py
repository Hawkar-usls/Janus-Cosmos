#!/usr/bin/env python3
import hashlib, io, json, zipfile, re
from pathlib import Path
import requests, shapefile
from shapely.geometry import Point, shape

URL="https://www.gmrt.org/shapefiles/gmrt_swath_polygons.zip"
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
TARGETS=[
 {"id":"KN19207_CAND_002","lat":-3.9727527956056825,"lon":-12.272824298723462},
 {"id":"KN19207_CAND_003","lat":-4.015075679897318,"lon":-12.29915403590432}
]
KNOWN={"KN192-07","KNOX15RR"}
r=requests.get(URL,timeout=180,headers={"User-Agent":"JANUS-KUSTO-GMRT-third-survey-audit/1.0"})
r.raise_for_status()
b=r.content
sha=hashlib.sha256(b).hexdigest()
z=zipfile.ZipFile(io.BytesIO(b))
names=z.namelist()
work=OUT/"gmrt_swath_polygons_extract";work.mkdir(exist_ok=True)
z.extractall(work)
shps=list(work.rglob("*.shp"))
if not shps: raise RuntimeError("No .shp in GMRT swath zip")
if len(shps)>1: print("WARNING multiple shp:",[str(x) for x in shps])
shp=shps[0]
reader=shapefile.Reader(str(shp))
fields=[f[0] for f in reader.fields[1:]]
field_meta=[{"name":f[0],"type":f[1],"size":f[2],"decimal":f[3]} for f in reader.fields[1:]]

def norm(v):
    if isinstance(v,bytes):return v.decode("utf-8","replace")
    return v

def recdict(rec):
    return {fields[i]:norm(rec[i]) for i in range(len(fields))}

def possible_ids(d):
    vals=[]
    for k,v in d.items():
        if v is None:continue
        s=str(v).strip()
        ku=k.lower()
        if any(x in ku for x in ["survey","cruise","expedition","leg","id","name"]):
            if s and s.lower() not in {"none","null","nan"}:
                vals.append({"field":k,"value":s})
    return vals

def canonical_tokens(d):
    # Conservative extraction only; keep raw records authoritative.
    text=" ".join(str(v) for v in d.values() if v is not None)
    toks=sorted(set(re.findall(r'\b(?:KN192-07|KNOX15RR|[A-Z]{1,8}[0-9]{2,4}[-_A-Z0-9]*)\b',text)))
    return toks

hits={t["id"]:[] for t in TARGETS}
bbox_candidates={t["id"]:0 for t in TARGETS}
for sr in reader.iterShapeRecords():
    d=recdict(sr.record)
    bbox=sr.shape.bbox
    for t in TARGETS:
        x,y=t["lon"],t["lat"]
        if bbox[0] <= x <= bbox[2] and bbox[1] <= y <= bbox[3]:
            bbox_candidates[t["id"]]+=1
            geom=shape(sr.shape.__geo_interface__)
            p=Point(x,y)
            if geom.covers(p):
                hits[t["id"]].append({
                  "attributes":d,
                  "possible_identity_fields":possible_ids(d),
                  "canonical_tokens":canonical_tokens(d),
                  "bbox":bbox
                })

summary={}
for t in TARGETS:
    hh=hits[t["id"]]
    tokens=sorted(set(tok for h in hh for tok in h["canonical_tokens"]))
    known=sorted(set(x for x in tokens if x in KNOWN))
    new=sorted(set(x for x in tokens if x not in KNOWN))
    summary[t["id"]]={
      "target":t,
      "bbox_candidate_polygon_count":bbox_candidates[t["id"]],
      "exact_intersection_count":len(hh),
      "exact_intersections":hh,
      "canonical_tokens":tokens,
      "known_survey_tokens":known,
      "candidate_new_tokens":new
    }

out={
 "artifact_id":"JANUS-KUSTO-GMRT-CAND002003-THIRD-SURVEY-COVERAGE-RUN-2026-09-22-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-GMRT-CAND002003-THIRD-SURVEY-COVERAGE-PREREG-2026-09-22-v1.0.json",
 "source":{"url":URL,"bytes":len(b),"sha256":sha,"zip_members":names},
 "shapefile":{"path":str(shp),"shapeType":reader.shapeTypeName,"record_count":len(reader),"fields":field_meta},
 "targets":summary,
 "warning":"candidate_new_tokens are conservative text tokens only; raw exact_intersections and explicit survey attribution fields are authoritative. New survey independence requires manual/source review before any depth access.",
 "claim_ceiling":"EXACT_GMRT_SOURCE_COVERAGE_INVENTORY_ONLY"
}
p=OUT/"JANUS-KUSTO-GMRT-CAND002003-THIRD-SURVEY-COVERAGE-RUN-2026-09-22-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "archive_sha256":sha,
 "shapefile":out["shapefile"],
 "targets":{k:{
   "bbox_candidates":v["bbox_candidate_polygon_count"],
   "exact_count":v["exact_intersection_count"],
   "tokens":v["canonical_tokens"],
   "intersections":v["exact_intersections"]
 } for k,v in summary.items()}
},indent=2,ensure_ascii=False))
