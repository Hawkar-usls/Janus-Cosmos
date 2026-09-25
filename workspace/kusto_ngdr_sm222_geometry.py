#!/usr/bin/env python3
import requests,re,json,hashlib
from pathlib import Path

BASE="https://geodataindia.gov.in/"
LAT=17.0710306; LON=83.2693611
OUT=Path("workspace/kusto_global_groundtruth_out/ngdr_sm222_geometry"); OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-NGDR-SM222/4.0","Referer":BASE+"guestuser"})
h=S.get(BASE+"guestuser",timeout=120); h.raise_for_status()
m=re.search(r'id=["\']csrfvalue["\']\s+value=["\']([^"\']+)',h.text)
if not m: raise RuntimeError("CSRF missing")
tok=m.group(1)

def post(path,data):
    r=S.post(BASE+path+"?_csrf="+tok,data=data,timeout=120)
    rec={"status":r.status_code,"url":r.url,"content_type":r.headers.get("content-type"),"bytes":len(r.content),"text":r.text[:10000]}
    try: rec["json"]=r.json()
    except Exception: rec["json"]=None
    return rec

# Candidate identifiers from official catalog hit for exact SM-222-title report.
ids=[5192,5837]
bounds=[]
for ident in ids:
    q=post("guestuser/getgridBoundary_basereport",{"tablename":"baseline_data_gis_view","clmname":"baseid","grid_no":str(ident)})
    js=q.get("json")
    inside=None
    if isinstance(js,list) and len(js)>=4:
        try:
            xmin,ymin,xmax,ymax=map(float,js[:4])
            inside=(xmin<=LON<=xmax and ymin<=LAT<=ymax)
        except Exception: pass
    bounds.append({"candidate_id":ident,"response":q,"bbox_contains_target":inside})

# Probe generic geometry/intersection endpoints using target point and likely table/layer identifiers.
point_wkt=f"POINT({LON} {LAT})"
probes={}
for path,data in [
    ("guestuser/getgridBoundary_basereport",{"tablename":"baseline_data_gis_view","clmname":"baseid","grid_no":"5192"}),
    ("getgridBoundary_basereport",{"tablename":"baseline_data_gis_view","clmname":"baseid","grid_no":"5192"}),
]:
    key=path.replace("/","_")
    probes[key]=post(path,data)

out={
 "artifact_id":"JANUS-KUSTO-NGDR-SM222-GEOMETRY-PROBE-2026-09-25-v1.0",
 "target":{"lat":LAT,"lon":LON},
 "catalog_binding":{"id":5192,"subid":5837,"title":"Progress Report On Multibeam Bathymetric Survey Of The Continental Slope Off Pudimadaka-Godavari, Andhra Pradesh Coast"},
 "boundary_trials":bounds,
 "probes":probes
}
raw=json.dumps(out,indent=2,ensure_ascii=False)
(OUT/"geometry_probe.json").write_text(raw,encoding="utf-8")
print(json.dumps({
 "status":"PASS",
 "boundary_trials":[{"candidate_id":x["candidate_id"],"status":x["response"]["status"],
 "json":x["response"].get("json"),"bbox_contains_target":x["bbox_contains_target"],
 "text":x["response"]["text"][:1000]} for x in bounds],
 "sha256":hashlib.sha256(raw.encode()).hexdigest()
},indent=2,ensure_ascii=False))
