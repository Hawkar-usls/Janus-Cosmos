#!/usr/bin/env python3
"""Source-exact reader for MB-System MBF_MBLDEOIH FBT records, versions 1-5.

Implementation authority:
  MB-System src/mbio/mbr_mbldeoih.c
"""
from __future__ import annotations
import datetime as _dt
import math
import struct
from collections import Counter

_TAGS={
    b"dd":(1,38),
    b"nn":(2,44),
    b"DD":(3,48),
    b"V4":(4,90),
    b"V5":(5,98),
}
_COMMENTS={b"##":38,b"cc":2}

def _old_time(year, day, minute, sec, msec):
    try:
        dt=_dt.datetime(int(year),1,1,tzinfo=_dt.timezone.utc)
        dt += _dt.timedelta(days=int(day)-1, minutes=int(minute), seconds=int(sec), milliseconds=int(msec))
        return dt.timestamp()
    except Exception:
        return float("nan")

def parse_records(data: bytes, *, emit_arrays=True):
    off=0
    n=len(data)
    while off<n:
        if off+2>n:
            raise RuntimeError(f"trailing byte at {off}/{n}")
        tag=data[off:off+2]
        start=off
        if tag in _COMMENTS:
            hs=_COMMENTS[tag]
            end=off+hs+128
            if end>n:
                raise RuntimeError(f"truncated comment {tag!r} at {off}")
            yield {"kind":"comment","tag":tag.decode("ascii"),"version":0,"offset":start,"next_offset":end}
            off=end
            continue
        if tag not in _TAGS:
            raise RuntimeError(f"unknown MBLDEOIH tag {tag!r} at byte {off}")
        version,hs=_TAGS[tag]
        if off+hs>n:
            raise RuntimeError(f"truncated header {tag!r} at {off}")
        h=data[off:off+hs]
        if version>=4:
            time_d,lon,lat,sensordepth,altitude=struct.unpack_from(">5d",h,2)
            heading,speed,roll,pitch,heave,bx,bl=struct.unpack_from(">7f",h,42)
            if version==4:
                nb,na,ns,sensorhead=struct.unpack_from(">4h",h,70)
                depth_scale,distance_scale=struct.unpack_from(">2f",h,78)
            else:
                nb,na,ns,sensorhead=struct.unpack_from(">4i",h,70)
                depth_scale,distance_scale=struct.unpack_from(">2f",h,86)
            heading=float(heading)%360.0
            bath_is_absolute=False
        else:
            vals=struct.unpack_from(">14h",h,2)
            year,day,minute,sec,msec,lon2u,lon2b,lat2u,lat2b,heading_raw,speed_raw,nb,na,ns=vals
            p=30
            if version==1:
                depth_scale_i,distance_scale_i,td_short,alt_short=struct.unpack_from(">4h",h,p)
                sensor_depth_mm=int(depth_scale_i)*int(td_short)
                altitude_mm=int(depth_scale_i)*int(alt_short)
                bx_i=bl_i=ss_type=0
            elif version==2:
                depth_scale_i,distance_scale_i,td_short,alt_short,bx_i,bl_i,ss_type=struct.unpack_from(">7h",h,p)
                sensor_depth_mm=int(depth_scale_i)*int(td_short)
                altitude_mm=int(depth_scale_i)*int(alt_short)
            else:
                depth_scale_i,distance_scale_i=struct.unpack_from(">2h",h,p); p+=4
                sensor_depth_mm,altitude_mm=struct.unpack_from(">2i",h,p); p+=8
                bx_i,bl_i,ss_type=struct.unpack_from(">3h",h,p)
            time_d=_old_time(year,day,minute,sec,msec)
            lon=float(lon2u)/60.0 + float(lon2b)/600000.0
            lat=float(lat2u)/60.0 + float(lat2b)/600000.0 - 90.0
            sensordepth=0.001*float(sensor_depth_mm)
            altitude=0.001*float(altitude_mm)
            heading=(0.0054932*float(heading_raw))%360.0
            speed=0.01*float(speed_raw)
            roll=pitch=heave=0.0
            bx=0.01*bx_i if bx_i>0 else 2.0
            bl=0.01*bl_i if bl_i>0 else 2.0
            sensorhead=0
            depth_scale=0.001*float(depth_scale_i)
            distance_scale=0.001*float(distance_scale_i)
            bath_is_absolute=(version<3)
        if nb<0 or na<0 or ns<0 or nb>200000 or na>200000 or ns>2000000:
            raise RuntimeError(f"implausible dimensions v{version} @{off}: {nb},{na},{ns}")
        p=off+hs
        need=int(nb) + 6*int(nb) + 2*int(na) + 6*int(ns)
        end=p+need
        if end>n:
            raise RuntimeError(f"truncated arrays v{version} @{off}: need {need}, remain {n-p}")
        flags=memoryview(data)[p:p+nb]; p+=nb
        def shorts(count):
            nonlocal p
            if count<=0:
                return ()
            vals=struct.unpack_from(f">{count}h",data,p)
            p+=2*count
            return vals
        bath=shorts(nb)
        across=shorts(nb)
        along=shorts(nb)
        amp=shorts(na)
        ss=shorts(ns); ssx=shorts(ns); ssy=shorts(ns)
        rec={
            "kind":"data","tag":tag.decode("ascii"),"version":version,
            "offset":start,"next_offset":end,
            "time_d":float(time_d),"longitude":float(lon),"latitude":float(lat),
            "sensordepth_m":float(sensordepth),"altitude_m":float(altitude),
            "heading_deg":float(heading),"speed":float(speed),
            "roll":float(roll),"pitch":float(pitch),"heave":float(heave),
            "beam_xwidth_deg":float(bx),"beam_lwidth_deg":float(bl),
            "beams_bath":int(nb),"beams_amp":int(na),"pixels_ss":int(ns),
            "sensorhead":int(sensorhead),"depth_scale":float(depth_scale),
            "distance_scale":float(distance_scale),
            "bath_is_absolute":bool(bath_is_absolute),
        }
        if emit_arrays:
            rec["beamflag"]=bytes(flags)
            rec["bath"]=bath
            rec["bath_acrosstrack"]=across
            rec["bath_alongtrack"]=along
        off=end
        yield rec

def absolute_depth_m(rec, stored_bath):
    if rec["bath_is_absolute"]:
        return rec["depth_scale"]*stored_bath
    return rec["depth_scale"]*stored_bath + rec["sensordepth_m"]

def summarize(data: bytes):
    versions=Counter()
    tags=Counter()
    data_count=comment_count=0
    first=None; last=None
    nav_bounds=[math.inf,math.inf,-math.inf,-math.inf]
    bath_min=math.inf; bath_max=-math.inf; good_count=0
    last_offset=0
    for rec in parse_records(data,emit_arrays=True):
        tags[rec["tag"]]+=1
        last_offset=rec["next_offset"]
        if rec["kind"]=="comment":
            comment_count+=1
            continue
        data_count+=1;versions[str(rec["version"])]+=1
        if first is None:first={k:rec[k] for k in ("tag","version","time_d","longitude","latitude","heading_deg","beams_bath")}
        last={k:rec[k] for k in ("tag","version","time_d","longitude","latitude","heading_deg","beams_bath")}
        lon=rec["longitude"];lat=rec["latitude"]
        nav_bounds[0]=min(nav_bounds[0],lon);nav_bounds[1]=min(nav_bounds[1],lat)
        nav_bounds[2]=max(nav_bounds[2],lon);nav_bounds[3]=max(nav_bounds[3],lat)
        for fl,b in zip(rec["beamflag"],rec["bath"]):
            if fl==0:
                z=absolute_depth_m(rec,b)
                if math.isfinite(z):
                    bath_min=min(bath_min,z);bath_max=max(bath_max,z);good_count+=1
    return {
        "file_bytes":len(data),"bytes_consumed":last_offset,"eof_exact":last_offset==len(data),
        "tags":dict(tags),"versions":dict(versions),"data_records":data_count,"comment_records":comment_count,
        "first_data":first,"last_data":last,
        "nav_bounds":None if data_count==0 else nav_bounds,
        "good_beam_count":good_count,
        "good_depth_range_m":None if good_count==0 else [bath_min,bath_max],
    }
