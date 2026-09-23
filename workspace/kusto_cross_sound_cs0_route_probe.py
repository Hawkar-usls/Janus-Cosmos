#!/usr/bin/env python3
from __future__ import annotations
import json, time
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-CROSS-SOUND-ALASKA-UNSEEN-DUALCHANNEL-PREREG-2026-09-23-v1.0.json").read_text())
REP=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-CROSS-SOUND-CS0-SCIENCEBASE-ITEM-DOWNLOAD-ROUTE-REPAIR-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-CrossSound-CS0-route-probe/1.0"}

def probe(item_id):
    u=f"https://www.sciencebase.gov/catalog/file/get/{item_id}"
    # HEAD first. If unsupported/opaque, use streamed GET and close without reading body.
    attempts=[]
    try:
        r=requests.head(u,headers=UA,timeout=90,allow_redirects=True)
        attempts.append({"method":"HEAD","status":r.status_code,"url":r.url,"headers":dict(r.headers)})
        if 200 <= r.status_code < 400:
            return attempts[-1]
    except Exception as e:
        attempts.append({"method":"HEAD","error":type(e).__name__+": "+str(e)})
    try:
        r=requests.get(u,headers=UA,timeout=90,allow_redirects=True,stream=True)
        rec={"method":"GET_STREAM_HEADERS_ONLY","status":r.status_code,"url":r.url,"headers":dict(r.headers)}
        attempts.append(rec)
        r.close()
        if 200 <= rec["status"] < 400:
            return rec
    except Exception as e:
        attempts.append({"method":"GET_STREAM_HEADERS_ONLY","error":type(e).__name__+": "+str(e)})
    return {"attempts":attempts,"status":None}

results={}
for kind,spec in REP["authoritative_child_items"].items():
    rec=probe(spec["item_id"])
    h=rec.get("headers") or {}
    cd=h.get("content-disposition") or h.get("Content-Disposition")
    ct=h.get("content-type") or h.get("Content-Type")
    cl=h.get("content-length") or h.get("Content-Length")
    results[kind]={
      "item_id":spec["item_id"],"title":spec["title"],
      "download_endpoint":f"https://www.sciencebase.gov/catalog/file/get/{spec['item_id']}",
      "probe_method":rec.get("method"),"status_code":rec.get("status"),
      "resolved_url":rec.get("url"),"content_disposition":cd,
      "content_type":ct,"content_length":cl,
      "payload_body_read":False
    }

def eligible(x):
    if x.get("status_code") is None or not (200 <= int(x["status_code"]) < 400):return False
    cd=(x.get("content_disposition") or "").lower()
    ct=(x.get("content_type") or "").lower()
    # Item-level route may serve octet-stream/zip or named attachment.
    return ("attachment" in cd or "zip" in ct or "tiff" in ct or "octet-stream" in ct)

status="PASS_PRIMARY_DOWNLOAD_ROUTES_BOUND__NO_PAYLOAD_READ" if all(eligible(x) for x in results.values()) else "BLOCKED_ITEM_DOWNLOAD_ROUTE_NOT_CONFIRMED"
out={
 "artifact_id":"JANUS-KUSTO-CROSS-SOUND-CS0-PRIMARY-DOWNLOAD-ROUTE-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],"repair":REP["artifact_id"],"status":status,
 "products":results,
 "numeric_payload_read":False,"raster_values_read":False,
 "next_gate":"CS1_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if status.startswith("PASS") else "SOURCE_ROUTE_ALTERNATIVE_REQUIRED",
 "hard_rules_preserved":["CHILD_ITEM_IDENTITY_EXACT","NO_PAYLOAD_BODY_READ","NO_RASTER_VALUE_READ","NO_METHOD_CHANGE"],
 "claim_ceiling":"PRIMARY_DOWNLOAD_ROUTE_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-CROSS-SOUND-CS0-PRIMARY-DOWNLOAD-ROUTE-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
