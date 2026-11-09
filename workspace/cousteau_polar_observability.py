#!/usr/bin/env python3
"""Scientific firewall: audit Argo observation visibility, geolocation and data gaps.

This tool NEVER controls a float, changes satellite measurements, or manufactures
observations for unsampled cells. It is a QC/coverage companion to Cousteau.
"""
import argparse
import collections
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

BANDS = [(-90, -80), (-80, -70), (-70, -60), (-60, -50)]
SECTOR_EDGES = [-180, -120, -60, 0, 60, 120, 180]
GOOD_POSITION_QC = {"1"}
ESTIMATED_POSITION_QC = {"8"}
QC_DEFINITION = {
    "1": "good position",
    "8": "estimated/interpolated position; under-ice is ONE possible cause",
    "9": "missing position",
}

def latitude_band(lat):
    if lat == -50:
        return "[-60,-50]"
    for low, high in BANDS:
        if low <= lat < high:
            return f"[{low},{high})"
    return None

def longitude_sector(lon):
    if lon == 180:
        lon = -180
    if lon < -180 or lon >= 180:
        return None
    for a,b in zip(SECTOR_EDGES, SECTOR_EDGES[1:]):
        if a <= lon < b:
            return f"[{a},{b})"
    return None

def safe_float(value):
    try:
        v = float(value)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None

def source_integrity(profile):
    sha=profile.get("source_sha256")
    url=profile.get("source_url")
    if not isinstance(sha,str) or len(sha)!=64 or any(c not in "0123456789abcdef" for c in sha.lower()):
        return False
    return isinstance(url,str) and url.startswith(("https://data-argo.ifremer.fr/","https://nrlgodae1.nrlmry.navy.mil/","https://usgodae.org/"))

def profile_gate(profile):
    mode=str(profile.get("data_mode","")).strip()
    pos=str(profile.get("position_qc","")).strip()
    lat=safe_float(profile.get("latitude_deg"))
    lon=safe_float(profile.get("longitude_deg"))
    geo_valid = lat is not None and lon is not None and -90<=lat<=-50 and -180<=lon<=360
    if geo_valid and lon>180:lon-=360
    has_lineage=source_integrity(profile)
    samples=profile.get("sample_levels")
    physical_ok=bool(isinstance(samples,list) and samples and
                    all(isinstance(v,dict) and
                        safe_float(v.get("pressure_dbar")) is not None and
                        safe_float(v.get("temp_c")) is not None and
                        safe_float(v.get("psal_psu")) is not None
                        for v in samples))
    if not has_lineage or mode!="D" or not physical_ok:
        tier="REJECT_PHYSICAL_QC_OR_PROVENANCE"
    elif pos=="1" and geo_valid:
        tier="ADMIT_GEOGRAPHIC_QC1"
    elif pos=="8" and geo_valid:
        tier="RETAIN_PHYSICS__ESTIMATED_LOCATION__NO_FINE_MAP"
    else:
        tier="RETAIN_PHYSICS__LOCATION_UNKNOWN__NO_MAP"
    # This is a diagnostic. It does not infer ice from estimated coordinates.
    return {"tier":tier,"position_qc":pos or "MISSING","latitude_deg":lat if geo_valid else None,
            "longitude_deg":lon if geo_valid else None,
            "has_source_sha256":has_lineage,
            "profile_data_mode":mode,
            "physical_samples_present":physical_ok,
            "ice_presence":"UNKNOWN_WITHOUT_INDEPENDENT_SEA_ICE_PRODUCT",
            "geography_usable":tier=="ADMIT_GEOGRAPHIC_QC1"}

def lag_diagnostic(profile):
    # GDAC profile date and JULD_LOCATION are not receipt timestamps.
    sample=profile.get("observed_at_utc")
    receipt=profile.get("received_at_utc")
    if not sample or not receipt:
        return {"state":"TELEMETRY_DELAY_UNKNOWN__NOT_ZERO","lag_days":None}
    try:
        a=datetime.fromisoformat(sample.replace("Z","+00:00"))
        b=datetime.fromisoformat(receipt.replace("Z","+00:00"))
        if a.tzinfo is None or b.tzinfo is None:
            raise ValueError("timestamps require timezone")
        lag=(b-a).total_seconds()/86400
        if lag<0:
            return {"state":"CLOCK_OR_PROVENANCE_CONFLICT","lag_days":None}
        return {"state":"EXPLICIT_UTC_TIMESTAMPS","lag_days":round(lag,4)}
    except (ValueError, TypeError):
        return {"state":"INVALID_RECEIPT_TIME__NOT_ZERO","lag_days":None}

def polar_audit(observations):
    profiles=observations.get("profiles",[])
    if not isinstance(profiles,list):
        raise ValueError("profiles must be a list")
    tiers=collections.Counter()
    qc=collections.Counter()
    seasons=collections.Counter()
    cells=collections.Counter()
    float_ids=set()
    delayed_explicit=0
    lag_unknown=0
    for prof in profiles:
        if not isinstance(prof,dict):
            tiers["REJECT_INVALID_RECORD"]+=1
            continue
        check=profile_gate(prof)
        tiers[check["tier"]]+=1
        qc[check["position_qc"]]+=1
        if prof.get("wmo"):
            float_ids.add(str(prof["wmo"]))
        d=str(prof.get("index_date_utc",""))
        if len(d)>=6 and d[4:6].isdigit():
            month=int(d[4:6])
            if 1<=month<=12:
                seasons["AUSTRAL_WINTER_APR_SEP" if 4<=month<=9 else "AUSTRAL_SUMMER_OCT_MAR"]+=1
            else:
                seasons["UNKNOWN_MONTH"]+=1
        else:
            seasons["UNKNOWN_MONTH"]+=1
        lag=lag_diagnostic(prof)
        if lag["lag_days"] is None:lag_unknown+=1
        else:delayed_explicit+=1
        if check["geography_usable"]:
            band=latitude_band(check["latitude_deg"])
            sector=longitude_sector(check["longitude_deg"])
            if band and sector:
                cells[band+"_"+sector]+=1
    # Coverage of the OBSERVED sample only. Do not equate zero counts with no ocean data.
    return {"schema":"janus.cousteau.polar_observability_audit.v1",
            "source_run_id":observations.get("run_id","not_in_source"),
            "profile_count":len(profiles),"distinct_wmo":len(float_ids),
            "sample_biased":True,
            "position_qc_counts":dict(sorted(qc.items())),
            "admission_tiers":dict(sorted(tiers.items())),
            "observation_seasons":dict(sorted(seasons.items())),
            "telemetry_delay_provenance":{"explicit_lag_count":delayed_explicit,
                "unknown_lag_count":lag_unknown,
                "never_interpret_JULD_LOCATION_as_delivery_time":True},
            "observed_geographic_cells":dict(sorted(cells.items())),
            "geographic_cells_with_no_qc1_sample_in_subset":24-len(cells),
            "absence_interpretation":"UNSAMPLED_IN_THIS_SUBSET__NOT_PROOF_OF_NO_OCEAN_OBSERVATIONS",
            "under_ice_inference":"POSITION_QC_8_ALONE_DOES_NOT_PROVE_UNDER_ICE",
            "model_promotion_authorized":False}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path)
    parser.add_argument("--out",type=Path)
    args=parser.parse_args()
    if args.input:
        src=args.input.read_bytes()
        obj=json.loads(src)
        audit=polar_audit(obj)
        audit["source_json_sha256"]=hashlib.sha256(src).hexdigest()
        audit["source_path"]=str(args.input)
    else:
        audit=polar_audit({"profiles":[]})
        audit["status"]="NO_REAL_ARGO_INPUT_AVAILABLE__NO_REAL_PROFILE_CLAIM"
    txt=json.dumps(audit,indent=2,sort_keys=True)+"\n"
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(txt,encoding="utf-8")
    print(txt)

if __name__=="__main__":
    main()
