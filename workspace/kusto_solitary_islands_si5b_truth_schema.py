#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json,zipfile,tempfile,os
from pathlib import Path
import requests, shapefile
from pyproj import CRS

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI5-POSTRESULT-GROUNDTRUTH-CHARACTERIZATION-PREREG-2026-09-24-v1.0.json").read_text())
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SOLITARY-ISLANDS-SI5A-TRUTH-PACKAGE-INVENTORY-RECEIPT-2026-09-24-v1.0.json").read_text())
URL=PRE["truth_source"]["data_package_url"]
r=requests.get(URL,headers={"User-Agent":"JANUS-KUSTO-SolitaryIslands-SI5B-schema/1.0"},timeout=240,allow_redirects=True);r.raise_for_status()
blob=r.content
sha=hashlib.sha256(blob).hexdigest()
if sha!=REC["truth_package"]["zip_sha256"]: raise RuntimeError("truth package SHA mismatch")

layers={
 "towed_video":"SIMP_TowedVideoSubClass/SIMP_TowedVideo_SubstrateClass",
 "sediment_samples":"SIMP_SedimentsMetadata/SIMP_Sediments_SamplesRetrieved"
}
audit={}
with zipfile.ZipFile(io.BytesIO(blob)) as zf:
    with tempfile.TemporaryDirectory() as td:
        for role,stem in layers.items():
            files={}
            for ext in (".shp",".shx",".dbf",".prj",".cpg"):
                name=stem+ext
                if name in zf.namelist():
                    outp=Path(td)/(role+ext)
                    outp.write_bytes(zf.read(name)); files[ext]=str(outp)
            if ".shp" not in files or ".dbf" not in files:
                raise RuntimeError(f"{role}: missing required shp/dbf")
            sf=shapefile.Reader(files[".shp"])
            fields=[{"name":f[0],"type":f[1],"size":f[2],"decimal":f[3]} for f in sf.fields[1:]]
            shape_types={}
            for s in sf.iterShapes():
                shape_types[str(s.shapeType)]=shape_types.get(str(s.shapeType),0)+1
            prj=None;crs=None
            if ".prj" in files:
                prj=Path(files[".prj"]).read_text(errors="replace")
                try: crs=CRS.from_wkt(prj).to_string()
                except Exception: crs=None
            samples=[]
            for i,sr in enumerate(sf.iterShapeRecords()):
                if i>=3:break
                recd={fields[k]["name"]:v for k,v in enumerate(sr.record)}
                samples.append({"shape_type":sr.shape.shapeType,"bbox":getattr(sr.shape,"bbox",None),"record":recd})
            audit[role]={
              "feature_count":len(sf),
              "shape_type_header":sf.shapeType,
              "shape_type_counts":shape_types,
              "fields":fields,
              "prj_wkt":prj,
              "crs":crs,
              "first_three_records":samples
            }

out={
 "artifact_id":"JANUS-KUSTO-SOLITARY-ISLANDS-SI5B-TRUTH-LAYER-SCHEMA-AUDIT-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],
 "truth_inventory_receipt":REC["artifact_id"],
 "truth_package_sha256_verified":sha,
 "layers":audit,
 "candidate_geometry_scoring_run":False,
 "formal_si3_result_changed":False,
 "claim_ceiling":"POSTRESULT_TRUTH_LAYER_SCHEMA_CRS_GEOMETRY_AUDIT_ONLY"
}
p=OUT/"JANUS-KUSTO-SOLITARY-ISLANDS-SI5B-TRUTH-LAYER-SCHEMA-AUDIT-2026-09-24-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False,default=str);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "layers":{k:{"feature_count":v["feature_count"],"shape_type_header":v["shape_type_header"],"shape_type_counts":v["shape_type_counts"],"crs":v["crs"],"fields":[f["name"] for f in v["fields"]],"first_three_records":v["first_three_records"]} for k,v in audit.items()},
 "formal_si3_result_changed":False
},indent=2,ensure_ascii=False,default=str))
