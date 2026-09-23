#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures, json, math, re, statistics, struct, time
from pathlib import Path
from urllib.parse import urljoin, urlparse
import requests
from shapely.geometry import LineString, box

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-SHARK-ROCK-BLIND-TILE-PREREG-2026-09-23-v1.0.json").read_text())
UA={"User-Agent":"JANUS-KUSTO-shark-rock-G2A/1.0"}
LAYER="https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_footprints/MapServer/0/query"
T=PRE["blind_tile"]; TILE=box(T["lon_min"],T["lat_min"],T["lon_max"],T["lat_max"])
M_PER_DEG=111320.0
MAX_FBT_SAMPLE_PER_SURVEY=12

def get(url,params=None,stream=False):
    last=None
    for a in range(7):
        try:
            r=requests.get(url,params=params,headers=UA,timeout=120,stream=stream)
            if r.status_code==429:
                time.sleep(min(30,2**a)); continue
            r.raise_for_status(); return r
        except Exception as e:
            last=e
            if a==6: raise
            time.sleep(min(30,2**a))
    raise last

def hrefs(u):
    return sorted(set(urljoin(u,h) for h in re.findall(r'href=["\']([^"\']+)["\']',get(u).text,re.I)))

def survey_features(sid):
    d=get(LAYER,params={"where":f"SURVEY_ID='{sid}'","outFields":"SURVEY_ID,PLATFORM,INSTRUMENT,SURVEY_YEAR,DOWNLOAD_URL","returnGeometry":"false","f":"json"}).json()
    return [f.get("attributes") or {} for f in d.get("features",[])]

def derive_base(download_url):
    p=urlparse(download_url); parts=[x for x in p.path.split("/") if x]
    ship=parts[-2]; survey=parts[-1]
    if survey.endswith("_mb.html"): survey=survey[:-8]
    else: survey=survey.split(".")[0]
    return f"https://data.ngdc.noaa.gov/platforms/ocean/ships/{ship}/{survey}/multibeam/data/"

def list_fnv(base,maxdepth=4):
    seen=set(); out=[]
    def walk(u,d):
        if u in seen or d>maxdepth:return
        seen.add(u)
        for x in hrefs(u):
            if not x.startswith(base) or x.rstrip("/")==u.rstrip("/"): continue
            if x.endswith("/"): walk(x,d+1)
            elif x.lower().endswith(".fnv"): out.append(x)
    walk(base,0)
    return sorted(set(out))

def parse_fnv(text):
    rows=[]
    for line in text.splitlines():
        p=line.split()
        if len(p)<19:continue
        try: rows.append((float(p[6]),float(p[15]),float(p[16]),float(p[17]),float(p[18])))
        except: pass
    return rows

def fnv_tile_hit(url):
    rows=parse_fnv(get(url).text); n=0
    for _,plon,plat,slon,slat in rows:
        if LineString([(plon,plat),(slon,slat)]).intersects(TILE): n+=1
    return url if n else None

def fbt_url(fnv): return fnv[:-4]+".fbt"

def xy(lon,lat):
    lat0=(T["lat_min"]+T["lat_max"])/2
    return lon*M_PER_DEG*math.cos(math.radians(lat0)), lat*M_PER_DEG

def ll(x,y):
    lat0=(T["lat_min"]+T["lat_max"])/2
    return x/(M_PER_DEG*math.cos(math.radians(lat0))), y/M_PER_DEG

def parse_spacing_fbt(data):
    off=0; cross=[]; ping_xy=[]; records=0; good_beams=0; in_tile_records=0
    ids={b"V4":(4,90),b"V5":(5,98),b"cc":("comment",2),b"##":("oldcomment",30)}
    while off+2<=len(data):
        tag=data[off:off+2]
        if tag not in ids: break
        ver,hs=ids[tag]
        if isinstance(ver,str):
            off+=hs+128; continue
        h=data[off:off+hs]
        if len(h)<hs: break
        t,lon,lat,sd,alt=struct.unpack_from(">5d",h,2)
        heading,speed,roll,pitch,heave,bx,bl=struct.unpack_from(">7f",h,42)
        if ver==4:
            nb,na,ns,shd=struct.unpack_from(">4h",h,70); dscale,xscale=struct.unpack_from(">2f",h,78)
        else:
            nb,na,ns,shd=struct.unpack_from(">4i",h,70); dscale,xscale=struct.unpack_from(">2f",h,86)
        p=off+hs
        if min(nb,na,ns)<0 or nb>100000: break
        flags=data[p:p+nb]; p+=nb
        bath_start=p; p+=2*nb
        across_start=p; p+=2*nb
        along_start=p; p+=2*nb
        p+=2*na+6*ns
        if p>len(data): break
        if T["lon_min"]<=lon<=T["lon_max"] and T["lat_min"]<=lat<=T["lat_max"]:
            in_tile_records+=1
            ping_xy.append(xy(lon,lat))
            vals=[]
            for i in range(nb):
                if flags[i]!=0: continue
                a=struct.unpack_from(">h",data,across_start+2*i)[0]
                vals.append(xscale*a); good_beams+=1
            vals.sort()
            for a,b in zip(vals[:-1],vals[1:]):
                d=abs(b-a)
                if 0<d<1000: cross.append(float(d))
        records+=1; off=p
    along=[]
    for (x1,y1),(x2,y2) in zip(ping_xy[:-1],ping_xy[1:]):
        d=math.hypot(x2-x1,y2-y1)
        if 0<d<5000: along.append(float(d))
    return {
      "records_total":records,
      "records_nav_inside_tile":in_tile_records,
      "good_beams_nav_inside_tile":good_beams,
      "cross_track_adjacent_spacing_m":cross,
      "along_track_ping_spacing_m":along
    }

def q(vals,p):
    if not vals:return None
    v=sorted(vals); i=(len(v)-1)*p; a=int(math.floor(i)); b=int(math.ceil(i))
    if a==b:return v[a]
    return v[a]*(b-i)+v[b]*(i-a)

surveys={}
for sid in PRE["survey_ids"]:
    feats=survey_features(sid)
    bases=sorted(set(derive_base(f["DOWNLOAD_URL"]) for f in feats if f.get("DOWNLOAD_URL")))
    all_fnv=[]
    for b in bases:
        try: all_fnv.extend(list_fnv(b))
        except: pass
    all_fnv=sorted(set(all_fnv))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        hits=[x for x in ex.map(fnv_tile_hit,all_fnv) if x]
    # Deterministic spread over lexicographically ordered tile-hit files.
    hits=sorted(hits)
    if len(hits)<=MAX_FBT_SAMPLE_PER_SURVEY:
        sample=hits
    else:
        idx=sorted(set(round(i*(len(hits)-1)/(MAX_FBT_SAMPLE_PER_SURVEY-1)) for i in range(MAX_FBT_SAMPLE_PER_SURVEY)))
        sample=[hits[i] for i in idx]
    file_stats=[]; cross=[]; along=[]
    for fnv in sample:
        fu=fbt_url(fnv)
        try:
            r=get(fu); s=parse_spacing_fbt(r.content)
            cross.extend(s["cross_track_adjacent_spacing_m"]); along.extend(s["along_track_ping_spacing_m"])
            file_stats.append({"fbt":fu,"bytes":len(r.content),"records_total":s["records_total"],"records_nav_inside_tile":s["records_nav_inside_tile"],"good_beams_nav_inside_tile":s["good_beams_nav_inside_tile"]})
        except Exception as e:
            file_stats.append({"fbt":fu,"error":type(e).__name__+": "+str(e)})
    cross_med=q(cross,.5); along_med=q(along,.5)
    vals=[v for v in [cross_med,along_med] if v is not None and v>0]
    characteristic=max(vals) if vals else None
    surveys[sid]={
      "metadata_features":feats,
      "fnv_tile_hit_count":len(hits),
      "sampled_fbt_count":len(sample),
      "sample_rule":"deterministic_even_spread_over_sorted_tile_hit_fnv",
      "file_stats":file_stats,
      "cross_track_spacing_m":{"n":len(cross),"p10":q(cross,.1),"median":cross_med,"p90":q(cross,.9)},
      "along_track_spacing_m":{"n":len(along),"p10":q(along,.1),"median":along_med,"p90":q(along,.9)},
      "characteristic_spacing_m":characteristic
    }
    print(sid, "tile_fnv",len(hits),"sample_fbt",len(sample),"spacing",characteristic,flush=True)

eligible={k:v["characteristic_spacing_m"] for k,v in surveys.items() if v["characteristic_spacing_m"] is not None}
if not eligible: raise RuntimeError("No survey yielded blind sampling-density estimate")
common_base=max(eligible.values())
# Freeze physically interpretable common radii before truth unlock.
raw=[max(25.0,common_base*m) for m in PRE["stage_g2_scale_freeze"]["radius_multiplier_candidates"]]
radii=[]
for x in raw:
    z=min(float(PRE["stage_g2_scale_freeze"]["max_radius_m"]), max(float(PRE["stage_g2_scale_freeze"]["min_radius_m"]), round(x/25.0)*25.0))
    if z not in radii:radii.append(z)
grid=max(25.0, round((common_base*2.0)/25.0)*25.0)
grid=min(grid,250.0)

out={
 "artifact_id":"JANUS-KUSTO-SHARK-ROCK-G2A-SAMPLING-SCALE-FREEZE-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "truth_coordinate_read":False,
 "rov_truth_read":False,
 "depth_values_used_for_scale_selection":False,
 "surveys":surveys,
 "common_characteristic_spacing_m":common_base,
 "frozen_common_grid_m":grid,
 "frozen_common_radii_m":radii,
 "scale_logic":"max survey characteristic spacing where each characteristic=max(median adjacent good-beam cross-track spacing, median ping-nav along-track spacing); radii=4/8/16/32x rounded to 25m and capped by prereg",
 "next_gate":"G2B_BLIND_MULTISURVEY_MORPHOLOGY_DISCOVERY",
 "claim_ceiling":"SAMPLING_DENSITY_AND_SCALE_FREEZE_ONLY"
}
p=OUT/"JANUS-KUSTO-SHARK-ROCK-G2A-SAMPLING-SCALE-FREEZE-RUN-2026-09-23-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "characteristic_spacing_m":eligible,
 "common_characteristic_spacing_m":common_base,
 "frozen_common_grid_m":grid,
 "frozen_common_radii_m":radii,
 "truth_coordinate_read":False,
 "depth_values_used_for_scale_selection":False
},indent=2))
