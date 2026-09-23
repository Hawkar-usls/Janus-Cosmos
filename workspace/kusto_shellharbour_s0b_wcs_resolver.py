#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,xml.etree.ElementTree as ET
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHELLHARBOUR-S0B-WCS-COVERAGE-RESOLUTION-PREREG-2026-09-23-v1.0.json").read_text())
WCS=PRE["wcs_endpoint"]
UA={"User-Agent":"JANUS-KUSTO-Shellharbour-S0B-WCS/1.0"}
S=requests.Session();S.headers.update(UA)

def get(params,timeout=180):
    r=S.get(WCS,params=params,timeout=timeout,allow_redirects=True);r.raise_for_status();return r

def lname(tag):return tag.split("}",1)[-1]

cap=get({"SERVICE":"WCS","REQUEST":"GetCapabilities","VERSION":"2.0.1"})
root=ET.fromstring(cap.content)
rows=[]
for el in root.iter():
    if lname(el.tag) in ("CoverageSummary","CoverageOfferingBrief"):
        ident=title=abstract=""
        for c in el.iter():
            n=lname(c.tag);txt=(c.text or "").strip()
            if n in ("CoverageId","Identifier","name") and not ident and txt:ident=txt
            elif n in ("Title","label") and not title and txt:title=txt
            elif n=="Abstract" and not abstract and txt:abstract=txt
        if ident:rows.append({"id":ident,"title":title,"abstract":abstract})

bindings={}
for role,tokens in PRE["required_roles"].items():
    matches=[x for x in rows if all(t in (" ".join(x.values())).lower() for t in tokens)]
    if len(matches)!=1:
        bindings[role]={"status":"NOT_UNIQUELY_RESOLVED","tokens":tokens,"matches":matches}
        continue
    m=matches[0]
    # WCS 2.0 DescribeCoverage
    d=get({"SERVICE":"WCS","REQUEST":"DescribeCoverage","VERSION":"2.0.1","COVERAGEID":m["id"]})
    bindings[role]={
      "status":"RESOLVED",
      "coverage":m,
      "describe_url":d.url,
      "describe_sha256":hashlib.sha256(d.content).hexdigest(),
      "describe_bytes":len(d.content),
      "describe_xml":d.content.decode("utf-8","replace")[:50000]
    }

out={
 "artifact_id":"JANUS-KUSTO-SHELLHARBOUR-S0B-WCS-COVERAGE-RESOLUTION-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "capabilities":{"url":cap.url,"sha256":hashlib.sha256(cap.content).hexdigest(),"bytes":len(cap.content),"coverage_count":len(rows)},
 "shellharbour_coverages":[x for x in rows if "shellharbour" in (" ".join(x.values())).lower()],
 "bindings":bindings,
 "all_required_coverages_resolved":all(v.get("status")=="RESOLVED" for v in bindings.values()),
 "getcoverage_called":False,
 "raster_values_read":False,
 "landform_truth_read":False,
 "next_gate":"S1_BATHYMETRY_ONLY_WCS_GETCOVERAGE" if all(v.get("status")=="RESOLVED" for v in bindings.values()) else "WCS_COVERAGE_NAME_REPAIR_ONLY",
 "claim_ceiling":"WCS_COVERAGE_IDENTITY_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-SHELLHARBOUR-S0B-WCS-COVERAGE-RESOLUTION-RUN-2026-09-23-v1.0.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "shellharbour_coverages":out["shellharbour_coverages"],
 "resolved":{k:v.get("status") for k,v in bindings.items()},
 "getcoverage_called":False
},indent=2,ensure_ascii=False))
