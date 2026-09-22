#!/usr/bin/env python3
from __future__ import annotations
import ftplib, hashlib, json, math, struct, zipfile
from datetime import datetime, timezone
from pathlib import Path

HOST="livftp.noc.ac.uk"
ZIP_PATH="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11280/EM12/EM12raw.zip"
ZIP_SHA="060b65865312330e866f842959a89674b744ff44569dcb1fa13daa7b5b0cefb0"
LOCAL=Path("/tmp/CD169_EM12raw.zip")
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)

EM_SVP=0x029A
EM_12S_BATH=0x0296
SVP_SIZE=419
BATH_SIZE=926

FROZEN={
 "LINE18_ERA":{
   "file":"0018_240205_000335_raw.all",
   "pings":{
      60720:[80],
      60721:[79,80],
      60722:[80],
      60723:[79,80],
      60724:[80],
      60725:[79,80]
   }
 },
 "LINE32_ERA":{
   "file":"0032_270205_144041_raw.all",
   "pings":{21395:[0,1]}
 }
}

def get_zip():
    h=hashlib.sha256()
    f=ftplib.FTP(timeout=180); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I")
    with LOCAL.open("wb") as out:
        def cb(b): out.write(b); h.update(b)
        f.retrbinary("RETR "+ZIP_PATH,cb,blocksize=1024*1024)
    try:f.quit()
    except:f.close()
    if h.hexdigest()!=ZIP_SHA: raise RuntimeError("ZIP SHA mismatch")

def digits(b,a,z):
    s=b[a:z]
    if len(s)!=(z-a) or not all(48<=x<=57 for x in s): return None
    return int(s.decode("ascii"))

def parse_time(payload,pos_mode=False):
    try:
        dd=digits(payload,0,2); mm=digits(payload,2,4); yy=digits(payload,4,6)
        hh=digits(payload,6,8); mi=digits(payload,8,10); ss=digits(payload,10,12); cs=digits(payload,12,14)
        if None in (dd,mm,yy,hh,mi,ss,cs): return None
        year=2000+yy if yy<80 else 1900+yy
        dt=datetime(year,mm,dd,hh,mi,ss,cs*10000,tzinfo=timezone.utc)
        return dt.timestamp(),dt.isoformat().replace("+00:00","Z")
    except Exception:
        return None

def signed8(x): return x-256 if x>=128 else x

def parse_svp(payload,file_name,offset):
    tt=parse_time(payload)
    if tt is None or len(payload)<SVP_SIZE: return None
    n=struct.unpack_from("<h",payload,14)[0]
    if n<0 or n>100: return None
    prof=[]
    for i in range(n):
        d=struct.unpack_from("<h",payload,16+4*i)[0]
        v=struct.unpack_from("<h",payload,18+4*i)[0]
        prof.append({"depth_m":float(d),"velocity_m_s":0.1*float(v),"stored_velocity":int(v)})
    return {
      "file":file_name,"offset":offset,"epoch":tt[0],"utc":tt[1],
      "point_count":n,"profile":prof,
      "first_velocity_m_s":prof[0]["velocity_m_s"] if prof else None,
      "deepest_profile_depth_m":prof[-1]["depth_m"] if prof else None
    }

def parse_selected_bath(payload,file_name,offset,wanted_beams):
    tt=parse_time(payload)
    if tt is None or len(payload)<BATH_SIZE: return None
    ping=struct.unpack_from("<H",payload,14)[0]
    bath_res=payload[16]
    bath_quality_signed=signed8(payload[17])
    keel_raw=struct.unpack_from("<H",payload,18)[0]
    heading_raw=struct.unpack_from("<H",payload,20)[0]
    roll_raw=struct.unpack_from("<h",payload,22)[0]
    pitch_raw=struct.unpack_from("<h",payload,24)[0]
    ping_heave_raw=struct.unpack_from("<h",payload,26)[0]
    sound_vel_raw=struct.unpack_from("<H",payload,28)[0]
    bath_mode=payload[30]
    if bath_res==1:
        depth_scale=0.1; off_scale=0.2; tt_scale_ms=0.2; keel_scale=0.1
    elif bath_res==2:
        depth_scale=0.2; off_scale=0.5; tt_scale_ms=0.8; keel_scale=0.2
    else:
        depth_scale=off_scale=tt_scale_ms=keel_scale=None
    beams=[]
    for i in wanted_beams:
        base=32+11*i
        stored_bath=struct.unpack_from("<H",payload,base)[0]
        across=struct.unpack_from("<h",payload,base+2)[0]
        along=struct.unpack_from("<h",payload,base+4)[0]
        tt_raw=struct.unpack_from("<h",payload,base+6)[0]
        amp_raw=signed8(payload[base+8])
        qual=payload[base+9]
        heave_raw=signed8(payload[base+10])
        beams.append({
          "beam_index":i,
          "stored_bath":stored_bath,
          "raw_depth_m":None if stored_bath==0 or depth_scale is None else depth_scale*stored_bath,
          "across_track_m":None if off_scale is None else off_scale*across,
          "along_track_m":None if off_scale is None else off_scale*along,
          "travel_time_ms":None if tt_scale_ms is None else tt_scale_ms*tt_raw,
          "amplitude_db":0.5*amp_raw,
          "quality_byte":qual,
          "beam_heave_m":0.1*heave_raw,
          "stored_travel_time":tt_raw,
          "stored_amplitude":amp_raw,
          "stored_heave":heave_raw
        })
    return {
      "file":file_name,"offset":offset,"epoch":tt[0],"utc":tt[1],"ping_number":ping,
      "bath_res":bath_res,"bath_quality_signed":bath_quality_signed,
      "bath_quality_raw_byte":payload[17],
      "keel_depth_m":None if keel_scale is None else keel_scale*keel_raw,
      "heading_deg":0.1*heading_raw,
      "roll_deg":0.01*roll_raw,
      "pitch_deg":0.01*pitch_raw,
      "ping_heave_m":0.01*ping_heave_raw,
      "sound_velocity_m_s":0.1*sound_vel_raw,
      "bath_mode":bath_mode,
      "beams":beams
    }

def scan_records(data,file_name):
    svps=[]; selected=[]
    wanted_by_ping={}
    family=None
    for fam,spec in FROZEN.items():
        if spec["file"]==file_name:
            family=fam
            wanted_by_ping={int(k):v for k,v in spec["pings"].items()}
            break
    n=len(data)
    seen=set()
    for i in range(n-2):
        if data[i]!=0x02: continue
        typ=(data[i]<<8)|data[i+1]
        if typ==EM_SVP and i+2+SVP_SIZE<=n:
            key=(i,"SVP")
            if key in seen: continue
            seen.add(key)
            x=parse_svp(data[i+2:i+2+SVP_SIZE],file_name,i)
            if x: svps.append(x)
        elif family and typ==EM_12S_BATH and i+2+BATH_SIZE<=n:
            payload=data[i+2:i+2+BATH_SIZE]
            if len(payload)<16: continue
            try: ping=struct.unpack_from("<H",payload,14)[0]
            except: continue
            if ping in wanted_by_ping:
                key=(i,"BATH")
                if key in seen: continue
                seen.add(key)
                x=parse_selected_bath(payload,file_name,i,wanted_by_ping[ping])
                if x:
                    x["family"]=family
                    selected.append(x)
    return svps,selected

get_zip()
all_svps=[]; selected=[]
entry_inventory=[]
with zipfile.ZipFile(LOCAL,"r") as z:
    for info in z.infolist():
        data=z.read(info.filename)
        svps,pings=scan_records(data,info.filename)
        all_svps.extend(svps); selected.extend(pings)
        entry_inventory.append({"file":info.filename,"svp_count":len(svps),"selected_ping_count":len(pings)})

all_svps.sort(key=lambda x:(x["epoch"],x["file"],x["offset"]))
selected.sort(key=lambda x:(x["epoch"],x["file"],x["ping_number"]))

def bind_svp(p):
    prev=[s for s in all_svps if s["epoch"]<=p["epoch"]]
    if not prev: return None
    s=prev[-1]
    return {
      "file":s["file"],"utc":s["utc"],"epoch":s["epoch"],
      "age_seconds_at_ping":p["epoch"]-s["epoch"],
      "point_count":s["point_count"],
      "first_velocity_m_s":s["first_velocity_m_s"],
      "deepest_profile_depth_m":s["deepest_profile_depth_m"],
      "profile":s["profile"]
    }

for p in selected:
    p["latest_preserved_svp_at_or_before_ping"]=bind_svp(p)

by_family={}
for fam in FROZEN:
    pp=[p for p in selected if p["family"]==fam]
    by_family[fam]={
      "expected_file":FROZEN[fam]["file"],
      "selected_pings_found":len(pp),
      "pings":pp,
      "unique_recorded_surface_sound_velocity_m_s":sorted(set(p["sound_velocity_m_s"] for p in pp)),
      "unique_bath_res":sorted(set(p["bath_res"] for p in pp)),
      "unique_bath_mode":sorted(set(p["bath_mode"] for p in pp)),
      "unique_svp_bindings":[
        {"file":x[0],"utc":x[1],"age_seconds_min":min(x[2]),"age_seconds_max":max(x[2])}
        for x in sorted({
          (p["latest_preserved_svp_at_or_before_ping"]["file"],
           p["latest_preserved_svp_at_or_before_ping"]["utc"],
           tuple([p["latest_preserved_svp_at_or_before_ping"]["age_seconds_at_ping"]]))
          for p in pp if p["latest_preserved_svp_at_or_before_ping"]
        }, key=lambda q:(q[1],q[0]))
      ] if False else []
    }
    # Build proper compact unique binding summary
    tmp={}
    for p in pp:
        s=p["latest_preserved_svp_at_or_before_ping"]
        if not s: continue
        key=(s["file"],s["utc"])
        tmp.setdefault(key,[]).append(s["age_seconds_at_ping"])
    by_family[fam]["unique_svp_bindings"]=[
      {"file":k[0],"utc":k[1],"age_seconds_min":min(v),"age_seconds_max":max(v)}
      for k,v in sorted(tmp.items(),key=lambda kv:kv[0][1])
    ]

# Diagnostic comparisons only; no causal promotion.
l18=by_family["LINE18_ERA"]; l32=by_family["LINE32_ERA"]
svp18=l18["pings"][0]["latest_preserved_svp_at_or_before_ping"] if l18["pings"] else None
svp32=l32["pings"][0]["latest_preserved_svp_at_or_before_ping"] if l32["pings"] else None
same_svp = bool(svp18 and svp32 and svp18["file"]==svp32["file"] and svp18["utc"]==svp32["utc"])
same_profile = bool(svp18 and svp32 and svp18["profile"]==svp32["profile"])
surface18=l18["unique_recorded_surface_sound_velocity_m_s"]
surface32=l32["unique_recorded_surface_sound_velocity_m_s"]
diagnostic={
  "same_latest_preserved_svp_record":same_svp,
  "same_latest_preserved_svp_profile_values":same_profile,
  "line18_recorded_surface_sound_velocity_m_s":surface18,
  "line32_recorded_surface_sound_velocity_m_s":surface32,
  "recorded_surface_sound_velocity_sets_equal":surface18==surface32
}

out={
 "artifact_id":"JANUS-KUSTO-CAND003-RAW-LINE18-VS-LINE32-SVP-ANCILLARY-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CAND003-RAW-LINE18-VS-LINE32-SVP-ANCILLARY-PREREG-2026-09-23-v1.0.json",
 "source":{"path":ZIP_PATH,"sha256":ZIP_SHA,"entries":66},
 "all_preserved_svp_count":len(all_svps),
 "all_preserved_svp_inventory":all_svps,
 "families":by_family,
 "diagnostic":diagnostic,
 "entry_inventory":entry_inventory,
 "depth_used_for_selection":False,
 "quality_used_for_selection":False,
 "interpretation_ceiling":"RAW_EM12_SVP_AND_ANCILLARY_MECHANISM_DIAGNOSTIC_ONLY__NO_CAUSAL_REFRACTION_CLAIM"
}
p=OUT/"JANUS-KUSTO-CAND003-RAW-LINE18-VS-LINE32-SVP-ANCILLARY-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
  "all_preserved_svp_count":out["all_preserved_svp_count"],
  "svp_inventory":[{"file":s["file"],"utc":s["utc"],"n":s["point_count"],"first_velocity_m_s":s["first_velocity_m_s"],"deepest_profile_depth_m":s["deepest_profile_depth_m"]} for s in all_svps],
  "families":by_family,
  "diagnostic":diagnostic
},indent=2))
