#!/usr/bin/env python3
import json,re,requests,datetime as dt
from urllib.parse import urljoin,urlparse
from pathlib import Path
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-unseen-route-repair-v2/1.0"}

def get(u): return requests.get(u,headers=UA,timeout=60)
def hrefs(r,base): return sorted(set(urljoin(base,h) for h in re.findall(r'href=["\']([^"\']+)["\']',r.text,re.I))) if r.ok else []
def data_lines(txt): return [x for x in txt.splitlines() if x.strip() and not x.lstrip().startswith("#")]

pfbase="https://data.ngdc.noaa.gov/platforms/ocean/ships/pathfinder/PF0501/multibeam/data/version1/MB/generated/"
pr=get(pfbase); pfnv=[u for u in hrefs(pr,pfbase) if u.lower().endswith(".fnv")]
start=dt.datetime(2005,4,30,tzinfo=dt.timezone.utc).timestamp()
end=dt.datetime(2005,5,31,tzinfo=dt.timezone.utc).timestamp()
rows=[];files=[]
for u in pfnv:
    rr=get(u);kept=0;total=0
    if rr.ok:
        for line in data_lines(rr.text):
            p=line.split()
            if len(p)<19:continue
            try:
                epoch=float(p[6]);portlon=float(p[15]);portlat=float(p[16]);stbdlon=float(p[17]);stbdlat=float(p[18])
            except:continue
            total+=1
            if start<=epoch<end:
                rows.append((epoch,portlon,portlat,stbdlon,stbdlat));kept+=1
    files.append({"url":u,"parsed_rows":total,"rows_inside_official_cruise_window":kept})
lons=[x for r in rows for x in (r[1],r[3])];lats=[x for r in rows for x in (r[2],r[4])]
pf={
 "fnv_files":len(pfnv),
 "rows_in_official_cruise_window":len(rows),
 "files_with_rows_in_window":sum(x["rows_inside_official_cruise_window"]>0 for x in files),
 "envelope":None if not rows else {"west":min(lons),"east":max(lons),"south":min(lats),"north":max(lats)},
 "official_metadata_bounds":{"west":-76.973175,"east":-72.47391,"south":28.144133,"north":35.14623},
 "out_of_bounds_endpoint_fraction":None if not rows else sum(
   not(-77.5<=lon<=-72.0 and 27.5<=lat<=36.0)
   for r in rows for lon,lat in ((r[1],r[2]),(r[3],r[4]))
 )/(2*len(rows))
}

exroots=[
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/okeanos_explorer/EX2105/multibeam/data/version1/MB/",
 "https://data.ngdc.noaa.gov/platforms/ocean/ships/okeanos_explorer/EX2105/multibeam/data/version2/MB/"
]
ex=[]
for root in exroots:
    r=get(root);ls=hrefs(r,root)
    ex.append({
      "root":root,"status":r.status_code,
      "links":[{"url":u,"basename":urlparse(u).path.rstrip("/").split("/")[-1],"is_dir":u.endswith("/")} for u in ls]
    })
out={"artifact_id":"JANUS-KUSTO-UNSEEN-PF0501-EX2105-SOURCE-ROUTE-REPAIR-V2-2026-09-23-v1.0","depth_values_read":False,"PF0501_time_filtered":pf,"EX2105_root_links":ex,"claim_ceiling":"SOURCE_ROUTE_AND_METADATA_TIME_FILTER_PROBE_ONLY"}
p=OUT/"JANUS-KUSTO-UNSEEN-PF0501-EX2105-SOURCE-ROUTE-REPAIR-V2-2026-09-23-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
