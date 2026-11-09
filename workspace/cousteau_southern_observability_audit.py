#!/usr/bin/env python3
"""JANUS Cousteau observational coverage firewall.

Audits an already collected source-pinned Argo sample; does not infer conditions
in empty cells and never labels gridded bathymetry as direct soundings.
Stdlib only; works offline on archived source evidence.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

DEPTHS = (10, 100, 300, 600, 1000)
LAT_BANDS = ("50-60S", "60-70S", "70-80S", "80-90S")

def digest(data):
    return hashlib.sha256(data).hexdigest()

def band(lat):
    if not -90 <= lat <= -50:
        raise ValueError("outside preregistered Southern Ocean domain")
    return LAT_BANDS[min(3, int((-lat - 50)//10))]

def polar_season(month):
    if month in (12, 1, 2): return "austral_summer_DJF"
    if month in (3, 4, 5): return "austral_autumn_MAM"
    if month in (6, 7, 8): return "austral_winter_JJA"
    return "austral_spring_SON"

def build(observations, receipt, original_bytes, *, source_run):
    if digest(original_bytes) != receipt["observations_sha256"]:
        raise ValueError("observation SHA256 mismatch: REFUSE")
    if observations["index_sha256"] != receipt["source_index_sha256"]:
        raise ValueError("source index hash differs: REFUSE")
    if receipt["status"] != "REAL_ARGO_DATA_TRAINED_AND_EVALUATED":
        raise ValueError("no verified Argo training evidence")
    profiles = observations["profiles"]
    if len(profiles) != observations["admissible_profiles"] or not profiles:
        raise ValueError("missing profiles")
    if len(set(p["wmo"] for p in profiles)) != len(profiles):
        raise ValueError("not WMO-disjoint sampled profiles")
    month=Counter(); year=Counter(); lat=Counter(); lng=Counter()
    depth=Counter(); season=Counter(); lat_month=Counter()
    for p in profiles:
        if p["data_mode"] != "D" or p["position_qc"] != "1":
            raise ValueError("untrusted profile QC entered geographically precise audit")
        b=band(float(p["latitude_deg"]))
        m=int(p["index_date_utc"][4:6])
        if not 1<=m<=12: raise ValueError("bad profile month")
        l=float(p["longitude_deg"])
        if not -180<=l<=360: raise ValueError("longitude invalid")
        year[p["index_date_utc"][:4]]+=1
        month[m]+=1; lat[b]+=1; season[polar_season(m)]+=1
        lon_180=((l+180)%360)-180
        sector=int((lon_180+180)//60)
        lng["%d_to_%dE" % (sector*60-180,sector*60-120)]+=1
        lat_month[(b,m)]+=1
        for sample in p["sample_levels"]:
            dep=int(sample["pressure_dbar"])
            if dep not in DEPTHS: raise ValueError("unexpected preregistered pressure")
            depth[dep]+=1
    n=len(profiles)
    sparse=[{"latitude_band":b,"month":m,"count":lat_month[(b,m)],
             "status":"OBSERVED_IN_SELECTED_SAMPLE" if lat_month[(b,m)] else "NO_OBSERVATION_IN_SELECTED_SAMPLE"}
            for b in LAT_BANDS for m in range(1,13)]
    zero=sum(x["count"]==0 for x in sparse)
    checks={
        "source_receipt_provenance":"VERIFIED_SHA256",
        "fraction_lat_month_cells_empty_in_sample":round(zero/48,4),
        "sample_representative_of_antarctic_ocean":False,
        "winter_JJA_count":season["austral_winter_JJA"],
        "austral_summer_DJF_count":season["austral_summer_DJF"],
        "under_ice_profile_geolocation_completeness":"UNKNOWN_FROM_QC1_ONLY_SAMPLE",
        "satellite_to_float_communication_latency":"UNKNOWN_NO_JULD_LOCATION_OR_TRANSMISSION_TIMES",
        "bathymetry_direct_sounding_share":"UNKNOWN_REQUIRES_IBCSO_TID_AND_RID",
        "ice_shelf_melt_estimation_allowed":False,
        "under_ice_navigation_accuracy_estimation_allowed":False,
        "automated_model_promotion_allowed":False
    }
    warnings=[
        "Profiles are a WMO-disjoint, 90-float convenience sample, NOT a random sample of the ocean.",
        "Empty lat-month cells mean no observations in THIS SAMPLE, not no ocean process.",
        "Geolocated data exclude POSITION_QC=8; zero such records in accepted sample cannot establish zero under-ice observations.",
        "Reported pressures 10/100/300/600/1000 dbar are linearly interpolated from quality-controlled measured levels, not necessarily direct instrument samples at those exact pressures.",
        "Index date is measurement date, not necessarily transmission/receipt time.",
        "Ridge holdout used distinct floats but is retrospective same-era interpolation, not prospective forecasting.",
        "Argo cannot directly measure submerged ice-shelf undersides, bathymetry, melt flux or robot navigation.",
    ]
    if season["austral_summer_DJF"] > n*0.4:
        warnings.append("SUMMER_DOMINATED_SELECTION: Southern Ocean year-round conditions not established.")
    if lat["80-90S"] == 0:
        warnings.append("NO_PROFILES_80_90S: high-latitude gap in selected sample.")
    model=receipt["result"].get("model_benchmark",{})
    return {
        "schema":"janus.cousteau.southern_ocean.observability_firewall.v1",
        "source_run_id":source_run, "status":"SOURCE_VERIFIED_COVERAGE_AUDIT__NOT_GLOBAL_OCEAN_MODEL",
        "source_observations_sha256":digest(original_bytes),
        "source_gdac_index_sha256":observations["index_sha256"],
        "source_doi":"10.17882/42182",
        "sample":{"unique_floats":n, "download_attempts":observations["profiles_download_attempted"],
                  "rejections":observations["rejections"],
                  "by_year":dict(sorted(year.items())),
                  "by_month":{str(m).zfill(2):month[m] for m in range(1,13)},
                  "by_season":dict(sorted(season.items())),
                  "by_latitude_band":{b:lat[b] for b in LAT_BANDS},
                  "by_60deg_longitude":dict(sorted(lng.items())),
                  "pressure_coverage_dbar":{str(d):depth[d] for d in DEPTHS},
                  "latitude_month_cells":sparse},
        "review":checks,"limits_and_warnings":warnings,
        "existing_float_holdout_metrics":model.get("results",{}),
        "research_origin":{"video_url":"https://youtu.be/oMPHmMnS0Mc",
                           "video_date":"2026-11-06",
                           "relevant_timecodes":["15:09-17:50","19:29-21:54","23:17-25:25","32:08-34:48","45:15-47:43"],
                           "video_is_inspiration_not_measurement_source":True},
        "external_controls":[
            {"dataset":"IBCSO_2025","required_grids":["bed","tid","rid"],"state":"NOT_INGESTED"},
            {"dataset":"MEOP_CTD","source":"https://www.ncei.noaa.gov/archive/archive-management-system/OAS/bin/prd/jquery/accession/details/297525","state":"NOT_INGESTED_REQUIRES_SEPARATE_QC"},
            {"dataset":"SOCCOM_BGC_ARGO","state":"NOT_INGESTED_REQUIRES_BGC_ADJUSTED_QC"}
        ]
    }

def selftest():
    p=[]
    for j in range(8):
        p.append({"wmo":str(1000+j),"data_mode":"D","position_qc":"1",
                  "latitude_deg":-53-j*4, "longitude_deg":-100+j*20,
                  "index_date_utc":"2023%02d01000000"%(j+1),
                  "sample_levels":[{"pressure_dbar":d,"temp_c":1,"psal_psu":34} for d in DEPTHS]})
    obs={"profiles":p,"admissible_profiles":8,"profiles_download_attempted":10,
         "rejections":{"NO_ADMISSIBLE_D_MODE_QC1_PROFILE":2},
         "index_sha256":"hash"}
    raw=json.dumps(obs).encode()
    receipt={"observations_sha256":digest(raw),"source_index_sha256":"hash",
             "status":"REAL_ARGO_DATA_TRAINED_AND_EVALUATED","result":{}}
    report=build(obs,receipt,raw,source_run="selftest_synthetic")
    assert sum(c["count"] for c in report["sample"]["latitude_month_cells"])==8
    assert report["review"]["automated_model_promotion_allowed"] is False
    assert report["sample"]["pressure_coverage_dbar"]["1000"]==8
    try: build(obs,receipt,raw+b"\n",source_run="tamper")
    except ValueError: pass
    else: raise AssertionError("tampering was not rejected")
    print("OBSERVABILITY_SELFTEST_PASS_SYNTHETIC_ONLY")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--self-test",action="store_true")
    parser.add_argument("--observations")
    parser.add_argument("--receipt")
    parser.add_argument("--out")
    args=parser.parse_args()
    if args.self_test: selftest(); return
    if not (args.observations and args.receipt and args.out):
        parser.error("--observations, --receipt and --out required")
    raw=Path(args.observations).read_bytes()
    obs=json.loads(raw)
    receipt=json.loads(Path(args.receipt).read_text())
    report=build(obs,receipt,raw,source_run=receipt["run_id"])
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print(json.dumps({"status":report["status"],"floats":report["sample"]["unique_floats"],
                      "winter":report["review"]["winter_JJA_count"],
                      "summer":report["review"]["austral_summer_DJF_count"],
                      "unobserved_bins":sum(x["count"]==0 for x in report["sample"]["latitude_month_cells"])},
                     sort_keys=True))

if __name__=="__main__": main()
