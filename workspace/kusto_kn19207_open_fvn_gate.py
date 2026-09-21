#!/usr/bin/env python3
import concurrent.futures, hashlib, json, math, re
from pathlib import Path
from urllib.parse import urljoin
import requests

TARGET={"id":"CD169_TOBI_LOCUS","lat":-3.865418,"lon":-12.14244}
BASE="https://data.ngdc.noaa.gov/platforms/ocean/ships/knorr/KN192-07/multibeam/data/version1/MB/generated/"
OUT=Path("workspace/kusto_open_seafloor_out")
OUT.mkdir(parents=True,exist_ok=True)
S=requests.Session()
S.headers.update({"User-Agent":"JANUS-KUSTO-open-fnv-gate/1.0"})

def xy(lon,lat):
    lat0=TARGET["lat"]
    return ((lon-TARGET["lon"])*111320*math.cos(math.radians(lat0)), (lat-TARGET["lat"])*111320)

def point_in_poly(poly):
    # Target is origin in local coordinates
    inside=False
    n=len(poly)
    for i in range(n):
        x1,y1=poly[i]; x2,y2=poly[(i+1)%n]
        if ((y1>0)!=(y2>0)):
            xin=x1+(x2-x1)*(-y1)/(y2-y1)
            if xin>0: inside=not inside
    return inside

def segdist(a,b):
    ax,ay=a; bx,by=b
    vx,vy=bx-ax,by-ay
    d=vx*vx+vy*vy
    if d==0:return math.hypot(ax,ay)
    t=max(0,min(1,-(ax*vx+ay*vy)/d))
    return math.hypot(ax+t*vx,ay+t*vy)

def parse_line(line):
    p=line.split()
    if len(p)<19:return None
    try:
        return {
            "iso":f"{p[0]}-{p[1]}-{p[2]}T{p[3]}:{p[4]}:{float(p[5]):09.6f}Z",
            "epoch":float(p[6]),
            "navlon":float(p[7]),"navlat":float(p[8]),
            "heading":float(p[9]),"speed":float(p[10]),"draft":float(p[11]),
            "roll":float(p[12]),"pitch":float(p[13]),"heave":float(p[14]),
            "portlon":float(p[15]),"portlat":float(p[16]),
            "stbdlon":float(p[17]),"stbdlat":float(p[18]),
        }
    except Exception:return None

def fetch(url):
    r=requests.get(url,timeout=60,headers={"User-Agent":"JANUS-KUSTO-open-fnv-gate/1.0"})
    r.raise_for_status()
    return url,r.text,hashlib.sha256(r.content).hexdigest()

html=S.get(BASE,timeout=60).text
names=sorted(set(re.findall(r'href="([^"]+\.fnv)"',html,re.I)))
urls=[urljoin(BASE,n) for n in names]
print(f"FNV files discovered: {len(urls)}")

results=[]
nearest={"distance_m":float("inf")}
covered=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=20) as ex:
    futs=[ex.submit(fetch,u) for u in urls]
    for k,fut in enumerate(concurrent.futures.as_completed(futs),1):
        url,text,sha=fut.result()
        rows=[r for r in (parse_line(x) for x in text.splitlines()) if r]
        file_hits=[]
        for i,r in enumerate(rows):
            center=xy(r["navlon"],r["navlat"])
            cross=(xy(r["portlon"],r["portlat"]),xy(r["stbdlon"],r["stbdlat"]))
            d=segdist(*cross)
            if d<nearest["distance_m"]:
                nearest={"distance_m":d,"file":url.split("/")[-1],"row":i,"record":r,
                         "center_distance_m":math.hypot(*center),
                         "port_distance_m":math.hypot(*cross[0]),"stbd_distance_m":math.hypot(*cross[1])}
            if i+1<len(rows):
                q=[
                    xy(r["portlon"],r["portlat"]),
                    xy(r["stbdlon"],r["stbdlat"]),
                    xy(rows[i+1]["stbdlon"],rows[i+1]["stbdlat"]),
                    xy(rows[i+1]["portlon"],rows[i+1]["portlat"]),
                ]
                if point_in_poly(q):
                    # preserve both records for exact raw-file narrowing
                    file_hits.append({
                        "row_i":i,
                        "record_i":r,
                        "record_j":rows[i+1],
                        "quad_local_m":q,
                        "center_distance_i_m":math.hypot(*center),
                        "cross_track_segment_distance_i_m":d
                    })
        if file_hits:
            covered.append({"fnv_file":url.split("/")[-1],"sha256":sha,"hits":file_hits})
        results.append({"fnv_file":url.split("/")[-1],"sha256":sha,"rows":len(rows),"hit_count":len(file_hits)})
        if k%50==0: print(f"scanned {k}/{len(urls)}")

summary={
 "artifact_id":"JANUS-KUSTO-KN19207-FNV-EXACT-SWATH-GATE-2026-09-21-v1.0",
 "target":TARGET,
 "source":{
   "survey_id":"KN192-07",
   "provider":"NOAA NCEI",
   "instrument":"SeaBeam 3012",
   "generated_fvn_index":BASE,
   "fnv_semantics":"MB-System fast navigation: tMXYHScRPr=X=Y+X+Y; = port-most and + starboard-most values"
 },
 "files_discovered":len(urls),
 "files_scanned":len(results),
 "files_with_target_inside_adjacent_ping_swath_quad":len(covered),
 "target_swath_covered":bool(covered),
 "coverage_hits":covered,
 "nearest_ping_cross_swath_segment":nearest,
 "file_manifest":sorted(results,key=lambda x:x["fnv_file"]),
 "claim_ceiling":"RAW_SURVEY_SWATH_GEOMETRY_COVERAGE_ONLY__NOT_PROCESSED_DEPTH_TRUTH__NOT_MORPHOLOGY_IDENTITY",
 "next_if_pass":"DOWNLOAD_ONLY_MATCHING_RAW_XSE_AND_GENERATED_FBT_INF__VERIFY_SOUNDINGS_AROUND_TARGET",
 "next_if_fail":"PRESERVE_NEGATIVE_AND_SEARCH_OTHER_GEBCO_TID11_SOURCE_LINEAGES"
}
p=OUT/"JANUS-KUSTO-KN19207-FNV-EXACT-SWATH-GATE-2026-09-21-v1.0.json"
p.write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps({k:summary[k] for k in ["files_discovered","files_scanned","files_with_target_inside_adjacent_ping_swath_quad","target_swath_covered"]},indent=2))
print("nearest:",json.dumps(nearest,indent=2)[:5000])
print("hits:",json.dumps(covered,indent=2)[:12000])
