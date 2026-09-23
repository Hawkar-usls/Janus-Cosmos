#!/usr/bin/env python3
from __future__ import annotations
import hashlib,io,json
from pathlib import Path
import requests
from rasterio.io import MemoryFile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHELLHARBOUR-S1-WMS-FORMAT-PROBE-PREREG-2026-09-23-v1.0.json").read_text())
P=PRE["probe"]
UA={"User-Agent":"JANUS-KUSTO-Shellharbour-S1-WMS-format-probe/1.0","Accept-Encoding":"identity"}

rows=[]
for fmt in P["formats"]:
    params={"SERVICE":"WMS","VERSION":"1.1.1","REQUEST":"GetMap","LAYERS":P["layer"],"STYLES":"",
            "SRS":P["srs"],"BBOX":",".join(str(x) for x in P["bbox"]),
            "WIDTH":str(P["width"]),"HEIGHT":str(P["height"]),"FORMAT":fmt,"TRANSPARENT":"TRUE"}
    try:
        r=requests.get(P["endpoint"],params=params,headers=UA,timeout=120,allow_redirects=True)
        raw=r.content
        rec={"format":fmt,"http_status":r.status_code,"content_type":r.headers.get("content-type"),
             "bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),"url":r.url}
        try:
            with MemoryFile(raw) as mf:
                with mf.open() as src:
                    rec.update({"gdal_open":True,"driver":src.driver,"shape":[src.height,src.width],
                                "count":src.count,"crs":str(src.crs) if src.crs else None,
                                "transform":list(src.transform)[:6],"nodata":src.nodata,
                                "mask_flags":[str(x) for x in src.mask_flag_enums[0]]})
        except Exception as e:
            rec.update({"gdal_open":False,"open_error":type(e).__name__+": "+str(e),
                        "text_prefix":raw[:1200].decode("utf-8","replace")})
        rows.append(rec)
    except Exception as e:
        rows.append({"format":fmt,"transport_error":type(e).__name__+": "+str(e)})

out={"artifact_id":"JANUS-KUSTO-SHELLHARBOUR-S1-WMS-FORMAT-PROBE-RUN-2026-09-23-v1.0",
     "prereg":PRE["artifact_id"],"results":rows,"pixel_values_read":False,"backscatter_read":False,
     "claim_ceiling":"TRANSPORT_CAPABILITY_DIAGNOSTIC_ONLY"}
p=OUT/"JANUS-KUSTO-SHELLHARBOUR-S1-WMS-FORMAT-PROBE-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
