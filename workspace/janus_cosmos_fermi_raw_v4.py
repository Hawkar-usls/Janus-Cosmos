#!/usr/bin/env python3
"""JANUS Cosmos V4: byte-pinned Fermi-LAT Extended event/spacecraft FITS gate.

No external FITS library is needed for the small scalar columns consumed here.
Only a strictly bounded subset of FITS BINTABLE is implemented; unsupported
types fail closed. This is NOT an IRF-forward-folded photon likelihood.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from janus_cosmos_fermi_observation_gate import verify_history
from janus_cosmos_simulation_audit import AuditError, ROOT, digest

SCHEMA = "janus.cosmos.fermi-lat-extended-raw-ingest.v4"
MANIFEST_SCHEMA = "janus.cosmos.fermi-lat-fits-manifest.v4"
MANIFEST = ROOT / "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-FITS-MANIFEST-2026-10-18.json"
BLOCK = 2880
SCALES = {"B":(1,">B"),"I":(2,">h"),"J":(4,">i"),"K":(8,">q"),
          "E":(4,">f"),"D":(8,">d"),"L":(1,None),"A":(1,None)}


@dataclass(frozen=True)
class HDU:
    header:dict[str,Any]
    content:bytes


def sha256(path:Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def _card_value(card:bytes) -> Any:
    s=card.decode("ascii")
    raw=s[10:].lstrip()
    if raw.startswith("'"):
        match=re.match(r"'((?:[^']|'')*)'",raw)
        if not match: raise AuditError("malformed FITS quoted card")
        return match.group(1).replace("''", "'").strip()
    raw=raw.split("/",1)[0].strip()
    if raw=="T": return True
    if raw=="F": return False
    try: return float(raw.replace("D","E")) if any(c in raw for c in ".eED") else int(raw)
    except ValueError: return raw


def read_fits(path:Path) -> list[HDU]:
    """Conservative FITS reader, refusing truncated/overlarge/malformed HDUs."""
    if path.stat().st_size > 50_000_000 or path.stat().st_size % BLOCK:
        raise AuditError("FITS size/alignment invalid")
    blob=path.read_bytes()
    if not blob.startswith(b"SIMPLE  ="):
        raise AuditError("FITS SIMPLE header missing")
    pos=0
    hdus=[]
    while pos<len(blob):
        start=pos
        h={}
        found_end=False
        for i in range(10000):
            off=start+80*i
            if off+80>len(blob): raise AuditError("truncated FITS header")
            card=blob[off:off+80]
            key=card[:8].decode("ascii",errors="strict").strip()
            if key=="END":
                found_end=True
                pos=start+((i+1)*80+BLOCK-1)//BLOCK*BLOCK
                break
            if card[8:10]==b"= ": h[key]=_card_value(card)
        if not found_end: raise AuditError("FITS END card missing")
        if pos>len(blob): raise AuditError("truncated padded FITS header")
        naxis=int(h.get("NAXIS",-1))
        if not 0<=naxis<=5: raise AuditError("invalid FITS axis count")
        if h.get("XTENSION")=="BINTABLE":
            nrow=int(h.get("NAXIS2",-1)); width=int(h.get("NAXIS1",-1))
            if nrow<0 or width<0 or nrow>1_000_000 or width>65536:
                raise AuditError("unsafe FITS table dimensions")
            data_len=nrow*width+int(h.get("PCOUNT",0))
        elif naxis==0: data_len=0
        else:
            dims=[int(h.get(f"NAXIS{i}",-1)) for i in range(1,naxis+1)]
            if any(n<0 for n in dims): raise AuditError("invalid FITS image dimensions")
            data_len=(abs(int(h.get("BITPIX",0)))//8)*math.prod(dims)+int(h.get("PCOUNT",0))
        if data_len<0 or pos+data_len>len(blob):
            raise AuditError("truncated FITS data")
        hdus.append(HDU(h,blob[pos:pos+data_len]))
        pos+=((data_len+BLOCK-1)//BLOCK)*BLOCK
    if len(hdus)<2 or pos!=len(blob):
        raise AuditError("FITS extensions missing or misaligned")
    return hdus


def _column_info(h:dict[str,Any]) -> dict[str,tuple[int,int,str,int]]:
    n=int(h.get("TFIELDS",-1))
    if not 1<=n<=200: raise AuditError("invalid FITS column count")
    col={}; off=0
    for i in range(1,n+1):
        name=h.get(f"TTYPE{i}")
        fmt=h.get(f"TFORM{i}")
        m=re.fullmatch(r"(\d*)([A-Z])",str(fmt or ""))
        if not isinstance(name,str) or not m:
            raise AuditError("unsupported FITS TFORM")
        repeat=int(m.group(1) or 1); kind=m.group(2)
        if kind not in SCALES and kind!="X":
            raise AuditError("unsupported variable/complex FITS data type")
        if repeat<=0 or repeat>8192: raise AuditError("unsafe FITS repeat count")
        size=(repeat+7)//8 if kind=="X" else repeat*SCALES[kind][0]
        if name in col: raise AuditError("duplicate FITS column")
        col[name]=(off,size,kind,repeat);off+=size
    if off!=int(h.get("NAXIS1",-1)):
        raise AuditError("FITS field-width mismatch")
    return col


def rows(hdu:HDU,fields:set[str]) -> list[dict[str,Any]]:
    h=hdu.header
    if h.get("XTENSION")!="BINTABLE": raise AuditError("expected FITS BINTABLE")
    meta=_column_info(h)
    if not fields.issubset(meta): raise AuditError(f"FITS required fields missing: {fields-set(meta)}")
    stride=int(h["NAXIS1"])
    out=[]
    for row in range(int(h["NAXIS2"])):
        record={}
        start=row*stride
        for name in fields:
            off,size,kind,rep=meta[name]
            b=hdu.content[start+off:start+off+size]
            if kind=="X": val=int.from_bytes(b,"big") >> (size*8-rep)
            elif kind=="L": val=tuple(x==ord("T") for x in b)
            elif kind=="A": val=b.decode("ascii").strip()
            else:
                form=SCALES[kind][1]
                val=tuple(v[0] for v in struct.iter_unpack(form,b))
                if rep==1: val=val[0]
            if kind=="L" and rep==1: val=val[0]
            record[name]=val
        out.append(record)
    return out


def angular_sep_deg(ra:float,dec:float,ra0:float,dec0:float) -> float:
    a,b,c,d=map(math.radians,(ra,dec,ra0,dec0))
    v=math.sin(b)*math.sin(d)+math.cos(b)*math.cos(d)*math.cos(a-c)
    return math.degrees(math.acos(max(-1.,min(1.,v))))


def _check_times(items:list[dict[str,Any]],start:str,stop:str) -> None:
    previous=-math.inf
    for item in items:
        a=float(item[start]);b=float(item[stop])
        if not math.isfinite(a) or not math.isfinite(b) or a>=b or a<previous:
            raise AuditError("invalid/unsorted interval or overlapping GTI/SC")
        previous=b


def ingest_one(config:dict[str,Any],root:Path) -> dict[str,Any]:
    grb=config["grb"]
    hdus={}
    provenance={}
    for type_ in ("event","spacecraft"):
        entry=config[type_]
        p=root/entry["filename"]
        if not p.is_file() or p.stat().st_size!=entry["bytes"] or sha256(p)!=entry["sha256"]:
            raise AuditError(f"RAW_FITS_SHA256_OR_SIZE_MISMATCH:{grb}/{type_}")
        parsed=read_fits(p)
        if parsed[0].header.get("INSTRUME")!="LAT" or parsed[0].header.get("ORIGIN")!="FSSC":
            raise AuditError("NASA FSSC LAT primary FITS identity missing")
        hdus[type_]=parsed
        provenance[type_] = {"source_url":entry["url"],"sha256":entry["sha256"],
                             "bytes":entry["bytes"],"rows_by_hdu":{str(h.header.get("EXTNAME",i)):h.header.get("NAXIS2",0) for i,h in enumerate(parsed)}}
    evh={h.header.get("EXTNAME"):h for h in hdus["event"]}
    sch={h.header.get("EXTNAME"):h for h in hdus["spacecraft"]}
    if "EVENTS" not in evh or "GTI" not in evh or "SC_DATA" not in sch:
        raise AuditError("LAT EVENTS/GTI/SC_DATA HDUs required")
    eh=evh["EVENTS"].header
    eh0=hdus["event"][0].header
    sh0=hdus["spacecraft"][0].header
    for key in ("TSTART","TSTOP"):
        if abs(float(eh0[key])-float(sh0[key]))>1e-6:
            raise AuditError("LAT event/spacecraft time windows differ")
    if eh.get("DSTYP3")!="POS(RA,DEC)" or eh.get("DSTYP4")!="TIME" or eh.get("DSTYP5")!="ENERGY":
        raise AuditError("FSSC Extended query provenance filters not in FITS")
    sky=re.fullmatch(r"circle\(([-\d.]+),([-\d.]+),([-\d.]+)\)",str(eh.get("DSVAL3","")))
    if not sky or any(abs(float(v)-float(config[field]))>1e-5 for v,field in zip(sky.groups(),("ra_deg","dec_deg","query_radius_deg"))):
        raise AuditError("FITS spatial query differs from preregistration")
    query_time=[float(x) for x in str(eh.get("DSVAL4","")).split(":")]
    query_energy=[float(x) for x in str(eh.get("DSVAL5","")).split(":")]
    if len(query_time)!=2 or any(abs(x-float(eh0[k]))>1e-4 for x,k in zip(query_time,("TSTART","TSTOP"))):
        raise AuditError("FITS temporal query provenance invalid")
    if query_energy!=config["query_energy_MeV"]:
        raise AuditError("FITS energy query differs from preregistration")
    if _column_info(eh)["EVENT_CLASS"][2:] != ("X",32) or _column_info(eh)["EVENT_TYPE"][2:] != ("X",32):
        raise AuditError("Pass8 event masks must be 32-bit FITS bit strings")
    e=rows(evh["EVENTS"],{"TIME","ENERGY","RA","DEC","THETA","ZENITH_ANGLE","EVENT_CLASS","EVENT_TYPE","RUN_ID","EVENT_ID"})
    gti=rows(evh["GTI"],{"START","STOP"})
    sc=rows(sch["SC_DATA"],{"START","STOP","DATA_QUAL","LAT_CONFIG","IN_SAA","LIVETIME","RA_SCZ","DEC_SCZ"})
    _check_times(gti,"START","STOP");_check_times(sc,"START","STOP")
    if not gti or not sc: raise AuditError("empty GTI or spacecraft attitude table")
    selection=config["selection"]
    if selection["source_bit_mask"]!=128 or selection["spacecraft_data_qual"]!=1 or selection["spacecraft_lat_config"]!=1 or selection["spacecraft_in_saa"] is not False or selection["require_gti"] is not True:
        raise AuditError("preregistered Pass8 SOURCE spacecraft/GTI cuts changed")
    trigger_met=float(config["trigger_met_estimated"])
    ra0=float(config["ra_deg"]);dec0=float(config["dec_deg"])
    accepted=[]; seen=set(); rejected={"not_gti":0,"bad_spacecraft_or_gap":0,"class":0,
                                   "energy":0,"zenith":0,"theta":0,"roi":0}
    for row in e:
        t=float(row["TIME"]);en=float(row["ENERGY"])
        if not math.isfinite(t) or not math.isfinite(en) or en<=0:raise AuditError("nonfinite raw event")
        key=(row["RUN_ID"],row["EVENT_ID"])
        if key in seen:raise AuditError("duplicate LAT event identity")
        seen.add(key)
        if not any(x["START"]<=t<=x["STOP"] for x in gti):rejected["not_gti"]+=1;continue
        matching=[x for x in sc if x["START"]<=t<=x["STOP"]]
        if len(matching)!=1 or matching[0]["DATA_QUAL"]!=1 or matching[0]["LAT_CONFIG"]!=1 or matching[0]["IN_SAA"]:
            rejected["bad_spacecraft_or_gap"]+=1;continue
        if not row["EVENT_CLASS"] & selection["source_bit_mask"]:rejected["class"]+=1;continue
        if not selection["min_energy_MeV"]<=en<=selection["max_energy_MeV"]:rejected["energy"]+=1;continue
        if float(row["ZENITH_ANGLE"])>selection["max_zenith_deg"]:rejected["zenith"]+=1;continue
        if float(row["THETA"])>selection["max_theta_deg"]:rejected["theta"]+=1;continue
        angle=angular_sep_deg(float(row["RA"]),float(row["DEC"]),ra0,dec0)
        if angle>selection["max_roi_deg"]:rejected["roi"]+=1;continue
        accepted.append({"met_s":t,"after_trigger_s":round(t-trigger_met,6),
            "energy_MeV":round(en,6),"angular_sep_deg":round(angle,6),
            "theta_deg":round(float(row["THETA"]),6),"zenith_deg":round(float(row["ZENITH_ANGLE"]),6),
            "event_class_mask":row["EVENT_CLASS"],"event_type_mask":row["EVENT_TYPE"],
            "event_id":row["EVENT_ID"],"run_id":row["RUN_ID"]})
    accepted.sort(key=lambda x:(x["met_s"],x["run_id"],x["event_id"]))
    live=sum(float(x["LIVETIME"]) for x in sc if x["DATA_QUAL"]==1 and x["LAT_CONFIG"]==1 and not x["IN_SAA"])
    return {"grb":grb,"redshift":config["redshift"],"class":config["class"],
        "selection":selection,"trigger_met_estimated":config["trigger_met_estimated"],
        "trigger_time_not_emission_time":True,"query_id":config["query_id"],
        "fits_provenance":provenance,"instrument_metadata":{"pass8_event_data_proc_ver":eh0.get("PROC_VER"),
             "spacecraft_data_proc_ver":sh0.get("PROC_VER"),"lat_event_timesys":eh0.get("TIMESYS"),
             "lat_event_timeref":eh0.get("TIMEREF"),"event_query_filters":{k:eh.get(k) for k in ("DSTYP2","DSVAL2","DSVAL3","DSVAL4","DSVAL5")}},
        "counts":{"events_in_archive_roi":len(e),"after_diagnostic_selection":len(accepted),
        "rejected_by_first_failed_cut":rejected,"gti_intervals":len(gti),"spacecraft_intervals":len(sc)},
        "spacecraft_nominal_livetime_s":round(live,6),"events":accepted}


def run(manifest:dict[str,Any],raw_root:Path) -> dict[str,Any]:
    if manifest.get("schema")!=MANIFEST_SCHEMA or manifest.get("experiment_id")!="JANUS-COSMOS-ALI-V4-FERMI-LAT-RAW-JOINT-2026-10-18":
        raise AuditError("V4 FITS manifest identity invalid")
    if manifest.get("claim_ceiling")!="NO_OBSERVATIONAL_LIV_INFERENCE_WITHOUT_IRF_GBM_SOURCE_LAG_MODEL":
        raise AuditError("V4 science ceiling must be explicit")
    names=[x.get("grb") for x in manifest.get("bursts",[])]
    if names!=["090510","090902B","130427A"]: raise AuditError("three preregistered bursts required")
    if any(x.get("selection")!=manifest.get("common_selection") for x in manifest["bursts"]):
        raise AuditError("common diagnostic selection must be fixed")
    history=verify_history()
    bursts=[ingest_one(x,raw_root) for x in manifest["bursts"]]
    result={"schema":SCHEMA,"experiment_id":manifest["experiment_id"],
        "observation_type":"INDIVIDUAL_REAL_FERMI_LAT_PASS8_EXTENDED_FITS_EVENTS_AND_SPACECRAFT",
        "frozen_history":history,"manifest_sha256":hashlib.sha256((json.dumps(manifest,sort_keys=True,separators=(',',':'))).encode()).hexdigest(),
        "sources":bursts,"counts":{"bursts":len(bursts),
         "raw_events_in_query_rois":sum(b["counts"]["events_in_archive_roi"] for b in bursts),
         "diagnostic_source_region_candidates":sum(b["counts"]["after_diagnostic_selection"] for b in bursts)},
        "instrument_caveat":"Event masks, spacecraft quality, GTI, incidence/zenith cuts are verified; no energy-dependent effective area, PSF, EDISP, background, event-by-event source probability or GBM calibration is forward folded.",
        "source_caveat":"The observed trigger is not an emission timestamp; arbitrary burst-dependent emission lags can absorb propagation effects.",
        "gates":{"RAW_EXTENDED_FITS":"PASS_SHA256_VERIFIED","SPACECRAFT_AND_GTI":"PASS_PARSED_DIAGNOSTIC_QUALITY_CUTS",
                 "IRF_FORWARD_FOLDING":"NOT_IMPLEMENTED","LAT_EVENT_ASSOCIATION":"NOT_ESTIMATED",
                 "GBM_PROMPT_SPECTROTEMPORAL_TEMPLATE":"NOT_INGESTED",
                 "INTRINSIC_SOURCE_LAG":"FREE_UNBOUNDED",
                 "MULTIBURST_PHYSICAL_LIKELIHOOD":"BLOCKED_SOURCE_AND_INSTRUMENT_DEGENERACY",
                 "LIV_OBSERVED":"NOT_ESTABLISHED","LATTICE_SPACING_CONSTRAINT":"NOT_DERIVED",
                 "SIMULATION_HYPOTHESIS":"NOT_ESTABLISHED","AMPS_FIREWALL":"OPEN_UNCONFIRMED","BHB_BOOM":"NOT_ESTABLISHED","OSIRIS_PROOF_AUTHORITY":"NO_CHANGE"},
        "claim_ceiling":"Real, SHA-pinned FITS events only support event-level and instrument-quality diagnostics; no physical lattice constraint or LIV preference is computed."}
    result["result_sha256"]=digest(result)
    return result


def main(argv:list[str]|None=None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest",type=Path,default=MANIFEST)
    ap.add_argument("--raw-dir",type=Path,required=True)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args(argv)
    try:
        raw=json.loads(args.manifest.read_text(encoding="utf-8"))
        out=json.dumps(run(raw,args.raw_dir),sort_keys=True,indent=2,ensure_ascii=False)+"\n"
        if args.output:args.output.write_text(out,encoding="utf-8")
        else:sys.stdout.write(out)
        return 0
    except (AuditError,OSError,ValueError,KeyError,struct.error) as exc:
        print(f"JANUS-FERMI-V4-RAW-ERROR: {exc}",file=sys.stderr)
        return 2


if __name__=="__main__":
    raise SystemExit(main())
