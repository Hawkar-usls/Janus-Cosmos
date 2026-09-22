#!/usr/bin/env python3
import json, requests, re
from urllib.parse import urlencode

pts=[
 ("KN19207_CAND_002",-3.9727527956056825,-12.272824298723462),
 ("KN19207_CAND_003",-4.015075679897318,-12.29915403590432),
]
FP="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-KNOX15RR-probe/1.0"})
for cid,lat,lon in pts:
    p={"geometry":f"{lon},{lat}","geometryType":"esriGeometryPoint","inSR":"4326","outSR":"4326",
       "spatialRel":"esriSpatialRelIntersects","outFields":"*","returnGeometry":"false","f":"json"}
    d=S.get(FP,params=p,timeout=60).json()
    print("\nPOINT",cid)
    for f in d.get("features",[]):
        a=f.get("attributes",{})
        if a.get("SURVEY_ID")=="KNOX15RR":
            print(json.dumps(a,indent=2,default=str))
            url=a.get("DOWNLOAD_URL")
            if url:
                try:
                    r=S.get(url,timeout=60,allow_redirects=True)
                    print("DOWNLOAD_URL_STATUS",r.status_code,r.url,r.headers.get("content-type"),len(r.content))
                    print(r.text[:5000] if "text" in r.headers.get("content-type","") or "html" in r.headers.get("content-type","") else "")
                except Exception as e: print("DOWNLOAD_URL_ERROR",repr(e))
bases=[
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/data/version1/MB/",
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/roger-revelle/KNOX15RR/multibeam/data/version1/MB/",
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/revelle/KNOX15RR/multibeam/data/version1/MB/",
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/",
]
for b in bases:
    try:
        r=S.get(b,timeout=60,allow_redirects=True)
        print("\nBASE",b,"STATUS",r.status_code,"FINAL",r.url,"CT",r.headers.get("content-type"),"LEN",len(r.content))
        if r.status_code==200: print(r.text[:6000])
    except Exception as e: print("BASE_ERR",b,repr(e))
