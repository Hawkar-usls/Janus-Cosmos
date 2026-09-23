#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,time
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHELLHARBOUR-S0B-WCS-COVERAGE-RESOLUTION-PREREG-2026-09-23-v1.0.json").read_text())
REP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHELLHARBOUR-S0B-WCS-TRANSPORT-REPAIR-2026-09-23-v1.0.json").read_text())
WCS=PRE["wcs_endpoint"]
UA={"User-Agent":"JANUS-KUSTO-Shellharbour-S0B-direct-describe/1.0","Accept-Encoding":"identity"}
S=requests.Session();S.headers.update(UA)

def describe(cid):
    last=None
    for a in range(5):
        try:
            r=S.get(WCS,params={"SERVICE":"WCS","REQUEST":"DescribeCoverage","VERSION":"2.0.1","COVERAGEID":cid},
                    timeout=120,allow_redirects=True)
            r.raise_for_status()
            content=r.content
            text=content.decode("utf-8","replace")
            informative=("CoverageDescription" in text or "boundedBy" in text or "domainSet" in text)
            return {"status":"RESOLVED" if informative else "NONINFORMATIVE",
                    "coverage_id":cid,"url":r.url,"http_status":r.status_code,
                    "sha256":hashlib.sha256(content).hexdigest(),"bytes":len(content),
                    "content_type":r.headers.get("content-type"),"describe_xml":text[:50000]}
        except Exception as e:
            last=e;time.sleep(min(12,2**a))
    return {"status":"TRANSPORT_BLOCKED","coverage_id":cid,"error":type(last).__name__+": "+str(last)}

bindings={role:describe(cid) for role,cid in REP["repair"]["coverage_ids_to_probe"].items()}
ok=all(v["status"]=="RESOLVED" for v in bindings.values())
out={
 "artifact_id":"JANUS-KUSTO-SHELLHARBOUR-S0B-WCS-COVERAGE-RESOLUTION-RUN-2026-09-23-v1.1",
 "prereg":PRE["artifact_id"],
 "transport_repair":REP["artifact_id"],
 "bindings":bindings,
 "all_required_coverages_resolved":ok,
 "getcoverage_called":False,
 "raster_values_read":False,
 "landform_truth_read":False,
 "next_gate":"S1_BATHYMETRY_ONLY_INPUT_BINDING" if ok else "WCS_COVERAGE_ROUTE_BLOCKED",
 "claim_ceiling":"WCS_COVERAGE_IDENTITY_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-SHELLHARBOUR-S0B-WCS-COVERAGE-RESOLUTION-RUN-2026-09-23-v1.1.json"
raw=json.dumps(out,indent=2,ensure_ascii=False);p.write_text(raw)
print(json.dumps({"artifact_id":out["artifact_id"],"resolved":{k:v["status"] for k,v in bindings.items()},
 "getcoverage_called":False,"raster_values_read":False},indent=2))
