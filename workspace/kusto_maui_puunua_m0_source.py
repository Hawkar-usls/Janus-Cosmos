#!/usr/bin/env python3
from __future__ import annotations
import json,re,time
from pathlib import Path
from urllib.parse import urljoin
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-MAUI-PUUNUA-UNSEEN-DUALCHANNEL-V2-PREREG-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-Maui-Puunua-M0-source/1.0"}
PAGES=[
 "https://cmgds.marine.usgs.gov/data-releases/datarelease/10.5066-P9XJLL6C/",
 "https://cmgds.marine.usgs.gov/catalog/pcmsc/DataReleases/CMGDS_DR_tool/DR_P9XJLL6C/A-01-13-HW_backscatter_1m_UTM04_WGS84_metadata.html"
]

def get(u):
    last=None
    for a in range(6):
        try:
            r=requests.get(u,headers=UA,timeout=90,allow_redirects=True)
            if r.status_code==429:
                time.sleep(min(20,2**a));continue
            r.raise_for_status();return r
        except Exception as e:
            last=e
            if a==5:raise
            time.sleep(min(20,2**a))
    raise last

hrefs=set();page_meta=[]
for p in PAGES:
    r=get(p);page_meta.append({"url":p,"resolved":r.url,"bytes":len(r.content),"status":r.status_code})
    for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I):
        hrefs.add(urljoin(r.url,h))
    for u in re.findall(r'https?://[^\s<>"\']+',r.text,re.I):
        hrefs.add(u.rstrip(").,;"))

# Exact selected site was frozen before this discovery.
site=[u for u in hrefs if ("puunua" in u.lower() or "puunoa" in u.lower())]
# Canonical data release routes may omit .zip in the terminal URL but carry name=<file>.zip.
bathy=sorted({u for u in site if "bathy" in u.lower() and "backscatter" not in u.lower() and ".zip" in u.lower()})
back_all=sorted({u for u in site if "backscatter" in u.lower() and ".zip" in u.lower()})\n# Prefer the unique direct media ZIP as canonical file identity; retain landing/download aliases as provenance.\nback_media=[u for u in back_all if "/data-releases/media/" in u.lower()]\nback=back_media if len(back_media)==1 else back_all

# If landing did not expose bathy, probe the deterministic paired naming pattern derived only from
# the published backscatter filename. HEAD/stream headers only: no ZIP body is read here.
probes=[]
if len(back)==1 and not bathy:
    cand=back[0].replace("backscatter","bathymetry")
    for alt in [cand,cand.replace("_bathymetry_","_bathy_")]:
        try:
            h=requests.head(alt,headers=UA,timeout=60,allow_redirects=True)
            probes.append({"url":alt,"status":h.status_code,"resolved":h.url,"content_type":h.headers.get("content-type"),"content_length":h.headers.get("content-length")})
            if 200<=h.status_code<400:
                bathy=[alt];break
        except Exception as e:
            probes.append({"url":alt,"error":type(e).__name__+": "+str(e)})

status="PASS_EXACT_PUUNUA_PAIR_BOUND__NO_RASTER_VALUES_READ" if len(bathy)==1 and len(back)==1 else "BLOCKED_PUUNUA_PAIR_BINDING"
out={
 "artifact_id":"JANUS-KUSTO-MAUI-PUUNUA-M0-SOURCE-BINDING-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],"status":status,
 "pages_read":page_meta,
 "selected_site":"PUUNUA_POINT",
 "candidate_urls":{"bathymetry":bathy,"backscatter":back},\n "backscatter_all_routes":back_all,
 "site_related_urls":sorted(site),
 "header_only_pair_probes":probes,
 "raster_values_read":False,"numeric_zip_body_read":False,
 "next_gate":"M1_BATHYMETRY_ONLY_BLIND_CANDIDATE_FREEZE" if status.startswith("PASS") else "SOURCE_BINDING_REPAIR_ONLY",
 "claim_ceiling":"SOURCE_FILE_BINDING_ONLY"
}
p=OUT/"JANUS-KUSTO-MAUI-PUUNUA-M0-SOURCE-BINDING-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
