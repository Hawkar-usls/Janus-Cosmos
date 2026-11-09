#!/usr/bin/env python3
"""Bounded, evidence-first Southern Ocean Argo-GDAC D-mode CTD profile benchmark."""
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import math
import os
from pathlib import Path
import time

import numpy as np
import requests
from netCDF4 import Dataset

OUT = Path("data/cousteau/southern_ocean")
RUN_ID = os.environ.get("GITHUB_RUN_ID", "local")
RAW_DIR = Path(os.environ.get("COUSTEAU_ARGO_RAW_DIR", "/tmp/cousteau-argo-core"))
URLS = (
    "https://data-argo.ifremer.fr/",
    "https://nrlgodae1.nrlmry.navy.mil/pub/outgoing/argo/",
    "https://usgodae.org/ftp/outgoing/argo/",
)
TARGET_P = (10, 100, 300, 600, 1000)
HEADERS = {"User-Agent": "JANUS-Cousteau-Argo-SouthernOcean/1.0 (research)"}

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def sha(b):
    return hashlib.sha256(b).hexdigest()

def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

def download(url, maximum=250000000):
    with requests.get(url, headers=HEADERS, timeout=(20, 120), stream=True) as resp:
        resp.raise_for_status()
        chunks = []
        n = 0
        for chunk in resp.iter_content(262144):
            if not chunk: continue
            n += len(chunk)
            if n > maximum: raise ValueError("download over byte cap")
            chunks.append(chunk)
        return b"".join(chunks)

def get_index():
    errors = []
    for root in URLS:
        url = root + "argo_profile_index.txt"
        try:
            payload = download(url, maximum=250000000)
            if b"file,date,latitude,longitude" not in payload[:6000].lower():
                raise ValueError("missing Argo GDAC index header")
            if len(payload) < 100000: raise ValueError("incomplete Argo index")
            return payload, url, errors
        except (requests.RequestException, ValueError) as e:
            errors.append({"url":url, "error":type(e).__name__, "detail":str(e)[:150]})
    raise RuntimeError("All three GDAC index endpoints unavailable: "+repr(errors))

def find_candidates(payload):
    txt = io.StringIO(payload.decode("utf-8-sig", errors="replace"))
    non_comment = (line for line in txt if not line.lstrip().startswith("#") and line.strip())
    reader = csv.DictReader(non_comment)
    expected = {"file", "date", "latitude", "longitude"}
    if not expected.issubset(set(reader.fieldnames or [])):
        raise ValueError("GDAC index schema mismatch: %s" % reader.fieldnames)
    seen = {}
    count_region = 0
    count_d = 0
    for row in reader:
        try:
            lat=float(row["latitude"])
            lon=float(row["longitude"])
            date=row["date"]
            path=row["file"]
            if not (-90 <= lat <= -50 and -180 <= lon <= 180 and
                    date[:8].isdigit() and "20220101" <= date[:8] < "20250101"):
                continue
            count_region += 1
            parts=path.split("/")
            if len(parts)<4 or parts[-2]!="profiles" or ".." in parts: continue
            name=parts[-1]
            if not (name.startswith("D") and name.endswith(".nc")): continue
            count_d += 1
            wmo=parts[-3]
            if not wmo.isdigit(): continue
            item={"file":path, "wmo":wmo, "date_utc":date, "lat":lat, "lon":lon}
            # Keep the most recent candidate from each WMO. One profile per WMO.
            if wmo not in seen or (date, path) > (seen[wmo]["date_utc"], seen[wmo]["file"]):
                seen[wmo]=item
        except (ValueError, TypeError, KeyError): continue
    arr=sorted(seen.values(),key=lambda z:sha((z["wmo"]+"|20261109").encode()))
    return arr, {"regional_2022_2024_index_entries":count_region,
                 "d_mode_file_entries":count_d,"unique_d_mode_wmo_candidates":len(arr)}

def normalize_path(path):
    if path.startswith("dac/"): return path
    return "dac/"+path

def extract_qc(ds, key, i):
    raw=np.ma.filled(np.ma.array(ds.variables[key][i]), b"9")
    return np.asarray(raw).astype("U1").ravel()

def extract_data(ds, key, i):
    return np.asarray(np.ma.filled(np.ma.array(ds.variables[key][i]), np.nan),dtype=float).ravel()

def extract_float(ds, key, i):
    return float(np.asarray(ds.variables[key][i]).ravel()[0])

def profile_from_nc(blob, item):
    if not (blob.startswith(b"CDF") or blob.startswith(bytes.fromhex("89484446"))):
        return None, "INVALID_NETCDF_SIGNATURE"
    with Dataset("argo-profile-inmemory.nc", memory=blob) as ds:
        required=("PRES_ADJUSTED","TEMP_ADJUSTED","PSAL_ADJUSTED",
                  "PRES_ADJUSTED_QC","TEMP_ADJUSTED_QC","PSAL_ADJUSTED_QC",
                  "DATA_MODE","POSITION_QC","LATITUDE","LONGITUDE")
        if any(k not in ds.variables for k in required):
            return None, "MISSING_CORE_ADJUSTED_VARIABLES"
        n=len(ds.dimensions["N_PROF"])
        excluded_pos=0
        for i in range(n):
            mode=str(np.asarray(ds.variables["DATA_MODE"][i]).astype("U1").ravel()[0])
            pos=str(np.asarray(ds.variables["POSITION_QC"][i]).astype("U1").ravel()[0])
            if pos=="8": excluded_pos+=1
            if mode!="D" or pos!="1": continue
            lat=extract_float(ds,"LATITUDE",i)
            lon=extract_float(ds,"LONGITUDE",i)
            if not (-90<=lat<=-50 and -180<=lon<=360): continue
            qc=(extract_qc(ds,"PRES_ADJUSTED_QC",i)=="1") & (
                extract_qc(ds,"TEMP_ADJUSTED_QC",i)=="1") & (
                extract_qc(ds,"PSAL_ADJUSTED_QC",i)=="1")
            p=extract_data(ds,"PRES_ADJUSTED",i)
            t=extract_data(ds,"TEMP_ADJUSTED",i)
            s=extract_data(ds,"PSAL_ADJUSTED",i)
            valid=(qc & np.isfinite(p) & np.isfinite(t) & np.isfinite(s) &
                   (p>=0) & (p<=2400) & (t>=-2.5) & (t<=35) &
                   (s>=20) & (s<=40))
            if np.count_nonzero(valid)<12: continue
            pres=p[valid]; temp=t[valid]; psal=s[valid]
            order=np.argsort(pres)
            pres=pres[order];temp=temp[order];psal=psal[order]
            unique,idx=np.unique(pres,return_index=True)
            pres=unique;temp=temp[idx];psal=psal[idx]
            samples=[]
            for depth in TARGET_P:
                k=int(np.searchsorted(pres,depth))
                if k==0 or k==len(pres): continue
                if pres[k]-pres[k-1]>120: continue
                if depth-pres[k-1]>65 or pres[k]-depth>65: continue
                samples.append({"pressure_dbar":depth,
                    "temp_c":round(float(np.interp(depth,pres,temp)),6),
                    "psal_psu":round(float(np.interp(depth,pres,psal)),6)})
            if len(samples)>=3:
                return {"wmo":item["wmo"],"source_index_path":item["file"],
                        "index_date_utc":item["date_utc"],"latitude_deg":lat,
                        "longitude_deg":lon,"data_mode":"D","position_qc":"1",
                        "n_good_levels":int(np.count_nonzero(valid)),
                        "sample_levels":samples,"under_ice_qc8_excluded_in_file":excluded_pos}, None
        return None, "NO_ADMISSIBLE_D_MODE_QC1_PROFILE"

def features(row, pressure):
    lat=row["latitude_deg"]
    lon=math.radians(row["longitude_deg"])
    try: month=int(row["index_date_utc"][4:6])
    except Exception: month=1
    p=math.log1p(pressure)/7
    x=(lat+60)/20
    return [p,x,math.sin(lon),math.cos(lon),math.sin(month*math.tau/12),
            math.cos(month*math.tau/12),p*x,p*math.sin(lon),p*math.cos(lon)]

def benchmark(accepted):
    wmos=sorted(set(r["wmo"] for r in accepted),
                key=lambda w:sha(("testsplit|"+w+"|20261109").encode()))
    if len(wmos)<24: raise ValueError("Not enough distinct QC1 D-mode floats")
    test=set(wmos[:max(6,len(wmos)//4)])
    split={"train_wmos":len(wmos)-len(test),"test_wmos":len(test)}
    all_scores={}
    for variable in ("temp_c","psal_psu"):
        train=[];holdout=[]
        for row in accepted:
            for item in row["sample_levels"]:
                rec=(row,item)
                (holdout if row["wmo"] in test else train).append(rec)
        if len(train)<50 or len(holdout)<12: raise ValueError("Insufficient valid pressure samples")
        depth_average={}
        for d in TARGET_P:
            vals=[p[variable] for r,p in train if p["pressure_dbar"]==d]
            if vals:depth_average[d]=float(np.mean(vals))
        holdout=[(r,p) for r,p in holdout if p["pressure_dbar"] in depth_average]
        X=np.array([features(r,p["pressure_dbar"]) for r,p in train])
        Y=np.array([p[variable] for r,p in train])
        T=np.array([features(r,p["pressure_dbar"]) for r,p in holdout])
        truth=np.array([p[variable] for r,p in holdout])
        mean=X.mean(axis=0)
        std=X.std(axis=0)
        std[std<1e-9]=1
        Z=(X-mean)/std
        # Fixed alpha, no holdout-informed tuning.
        alpha=10
        coeff=np.linalg.solve(Z.T@Z+alpha*np.eye(Z.shape[1]),Z.T@(Y-Y.mean()))
        predicted=float(Y.mean()) + ((T-mean)/std)@coeff
        naive=np.array([depth_average[p["pressure_dbar"]] for r,p in holdout])
        def metric(estimate):
            e=estimate-truth
            return {"rmse":round(float(np.sqrt(np.mean(e*e))),6),
                    "mae":round(float(np.mean(np.abs(e))),6)}
        all_scores[variable]={
            "n_train_pressure_samples":len(train),
            "n_test_pressure_samples":len(holdout),
            "ridge":metric(predicted),"depth_mean_baseline":metric(naive),
            "ridge_beats_depth_mean_baseline":metric(predicted)["rmse"]<metric(naive)["rmse"],
            "ridge_parameters":{"alpha":alpha,"feature_mean":mean.tolist(),
                                "feature_scale":std.tolist(),"coefficients":coeff.tolist(),
                                "intercept":float(Y.mean())},
            "holdout_errors_by_pressure":[
                {"pressure_dbar":d, "n":sum(p["pressure_dbar"]==d for r,p in holdout),
                 "ridge_rmse":metric(predicted[np.array([p["pressure_dbar"]==d for r,p in holdout])])["rmse"]}
                for d in TARGET_P if any(p["pressure_dbar"]==d for r,p in holdout)
            ]
        }
    return {"split":split,"mode":"float-disjoint retrospective same-era interpolation",
            "results":all_scores,"automatic_promotion":False,
            "claim_ceiling":"LIMITED_CONVENIENCE_SAMPLE__NOT_PROSPECTIVE_FORECAST__NO_ACOUSTIC_CLAIMS"}

def execute():
    run_start=utc()
    receipt=OUT/"argo_runs"/("run-"+RUN_ID+".json")
    try:
        index,index_url,fetch_errors=get_index()
        candidates,stats=find_candidates(index)
        RAW_DIR.mkdir(parents=True,exist_ok=True)
        idx_path=RAW_DIR/"argo_profile_index.txt"
        idx_path.write_bytes(index)
        accepted=[]
        manifest=[]
        rejection={}
        attempted=0
        for item in candidates[:160]:
            if len(accepted)>=90:break
            attempted+=1
            rel=normalize_path(item["file"])
            file_data=None
            last_error=""
            url_used=None
            for root in URLS:
                try:
                    candidate_url=root+rel
                    file_data=download(candidate_url,maximum=10000000)
                    url_used=candidate_url
                    break
                except (requests.RequestException,ValueError) as ex:
                    last_error=str(ex)[:120]
            if file_data is None:
                rejection["DOWNLOAD_FAILURE"]=rejection.get("DOWNLOAD_FAILURE",0)+1
                continue
            try:
                profile,reason=profile_from_nc(file_data,item)
            except Exception as ex:
                profile,reason=None, "PARSE_"+type(ex).__name__
            manifest.append({"index_path":item["file"],"wmo":item["wmo"],
                             "url":url_used,"sha256":sha(file_data),
                             "bytes":len(file_data),"outcome":reason or "ACCEPT"})
            if reason:
                rejection[reason]=rejection.get(reason,0)+1
                continue
            rawpath=RAW_DIR/(item["wmo"]+"-"+Path(rel).name)
            rawpath.write_bytes(file_data)
            profile["source_url"]=url_used
            profile["source_sha256"]=sha(file_data)
            accepted.append(profile)
        profiles_path=OUT/"argo_observations"/("argo_core_D_QC1_run-"+RUN_ID+".json")
        obs={"schema":"janus.cousteau.argo_core_qc1_observations.v1",
             "generated_utc":utc(),"index_url":index_url,
             "index_sha256":sha(index),"index_bytes":len(index),
             "index_selection_stats":stats,"index_fetch_errors":fetch_errors,
             "requested_profiles_max":90,"profiles_download_attempted":attempted,
             "admissible_profiles":len(accepted),"rejections":rejection,
             "download_manifest":manifest,"profiles":accepted}
        save(profiles_path,obs)
        if len(accepted)<24:
            result={"state":"DATA_GATE_UNMET_NO_TRAINING","accepted_profiles":len(accepted)}
        else:
            result={"state":"REAL_ARGO_DATA_TRAINED_AND_EVALUATED",
                    "model_benchmark":benchmark(accepted)}
        status=result["state"]
        save(receipt,{"schema":"janus.cousteau.argo_train_receipt.v1",
             "run_id":RUN_ID,"started_utc":run_start,"finished_utc":utc(),
             "status":status,"source_index_sha256":sha(index),
             "observations_path":str(profiles_path),
             "observations_sha256":sha(profiles_path.read_bytes()),
             "quality_gate":{"accepted":len(accepted),"attempted":attempted,"rejects":rejection},
             "result":result,"source_doi":"https://doi.org/10.17882/42182",
             "provenance_ceiling":"D-mode QC1 adjusted profile convenience sample, not the full Argo archive"})
        print(json.dumps({"status":status,"accepted":len(accepted),
                          "attempted":attempted,"rejections":rejection,
                          "benchmark":result.get("model_benchmark",{}).get("results",{})},
                          indent=2))
    except BaseException as ex:
        save(receipt,{"schema":"janus.cousteau.argo_train_receipt.v1",
                      "run_id":RUN_ID,"started_utc":run_start,
                      "finished_utc":utc(),"status":"FAILURE_NO_TRAINING_CLAIM",
                      "error":type(ex).__name__,"message":str(ex)[:1800]})
        raise

def selftest():
    x=[{"wmo":str(1000000+i),"latitude_deg":-58-(i%15),
        "longitude_deg":i*21-180,"index_date_utc":"20230615000000",
        "sample_levels":[{"pressure_dbar":d,"temp_c":2-.002*d+(i%10)*.03,
                          "psal_psu":34+.0003*d+(i%7)*.01} for d in TARGET_P]}
       for i in range(32)]
    result=benchmark(x)
    assert result["split"]["test_wmos"]>0
    assert all(result["results"][var]["ridge"]["rmse"]>=0 for var in ("temp_c","psal_psu"))
    assert len(list(TARGET_P))==5
    print("ARGO_SELF_TEST_PASS__SYNTHETIC_ONLY")

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--self-test",action="store_true")
    args=a.parse_args()
    selftest() if args.self_test else execute()
