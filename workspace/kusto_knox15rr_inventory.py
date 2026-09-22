#!/usr/bin/env python3
import requests,re,json
base="https://data.ngdc.noaa.gov/platforms/ocean/ships/roger_revelle/KNOX15RR/multibeam/data/version1/MB/"
s=requests.Session(); s.headers.update({"User-Agent":"JANUS-KUSTO-KNOX15RR-inventory/1.0"})
for url in [base,base+"generated/",base+"processed/",base+"products/"]:
    r=s.get(url,timeout=60)
    print("\nURL",url,"STATUS",r.status_code,"CT",r.headers.get("content-type"),"LEN",len(r.content))
    if r.status_code==200:
        hrefs=re.findall(r'href="([^"]+)"',r.text,re.I)
        print("HREF_COUNT",len(hrefs))
        print(json.dumps(hrefs[:500],indent=2))
