#!/usr/bin/env python3
"""Freeze *already downloaded* official NASA GBM and Fermitools CALDB assets.

This manifest is a single auditable source of byte checksums. It does not
download anything or interpret detector counts as a physical LIV detection.
"""

import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PACKAGE_SHA="f4ac356ee8e0ff621eb7d5acd39d069cafd8d9f0b2cfa8e533780e1004ae0e10"
PACKAGE_URL="https://api.anaconda.org/download/fermi/fermitools/2.5.8/linux-64/fermitools-2.5.8-py312h410524c_0.conda"
IRF_BASE="share/fermitools/data/caldb/data/glast/lat/bcf"
V4_OBSERVED_SHA="09d087ff9e19ae545b26870fda403541144b6a6e62db45ed35a55c8d2010b1be"
V4_MANIFEST_SHA="4c3ba242e8b0b718938c5bf5d264e21d57ef7e3f0293ebf9b72fbf55b3ca112b"
BURSTS=[("090510","2009","bn090510016",("n6","n7","b1"),"v00"),
        ("090902B","2009","bn090902462",("n0","n1","b0"),"v00"),
        ("130427A","2013","bn130427324",("n6","n7","b1"),"v01")]


def file_info(path:Path,name:str,url:str) -> dict:
    if not path.is_file():raise ValueError(f"required source absent: {path}")
    h=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1024*1024),b""):h.update(part)
    return {"name":name,"bytes":path.stat().st_size,"sha256":h.hexdigest(),"origin_url":url}


def make_manifest(caldb_dir:Path,gbm_dir:Path) -> dict:
    cal=[]
    for kind in ("ea","psf","edisp"):
        for partition in ("FB","PSF","EDISP"):
            name=f"{('aeff' if kind=='ea' else kind)}_P8R3_SOURCE_V3_{partition}.fits"
            p=f"{IRF_BASE}/{kind}/{name}"
            item=file_info(caldb_dir/p,p,PACKAGE_URL)
            item["kind"]=kind;item["partition"]=partition
            item["extraction_path_inside_package"]=p
            cal.append(item)
    bursts=[]
    for grb,year,bn,dets,cv in BURSTS:
        detectors=[]
        for det in dets:
            base=f"https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/{year}/{bn}/current/"
            names={"tte":f"glg_tte_{det}_{bn}_v00.fit",
                   "cspec":f"glg_cspec_{det}_{bn}_{cv}.pha",
                   "rsp2":f"glg_cspec_{det}_{bn}_v00.rsp2"}
            files={role:file_info(gbm_dir/name,name,base+name) for role,name in names.items()}
            detectors.append({"detector":det,"event_and_response_files":files})
        bursts.append({"grb":grb,"gbm_trigger_id":bn,"detectors":detectors,
              "detector_selection":"PILOT_FIXED_TWO_NAI_AND_ONE_BGO_NOT_OPTIMAL_RESPONSE_SELECTION",
              "data_scope":"TTE_SOURCE_AND_BACKGROUND_EVENTS_AND_CSPEC_RSP2_NOT_A_JOINT_SPECTRAL_FIT"})
    return {"schema":"janus.cosmos.fermi-irf-gbm-manifest.v5",
        "experiment_id":"JANUS-COSMOS-ALI-V5-GBM-LAT-IRF-2026-10-19",
        "retrieval_utc":"2026-10-18", "lat_caldb":{
           "irf_name":"P8R3_SOURCE_V3", "package_version":"fermitools-2.5.8",
           "package_url":PACKAGE_URL,"package_sha256":PACKAGE_SHA,
           "installation_source":"OFFICIAL_FERMI_CHANNEL_ANACONDA_PACKAGE",
           "caldb_files":cal,
           "document_url":"https://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_LAT_IRFs/IRF_overview.html"},
        "gbm_source_archive":"NASA_HEASARC_TRIGGER_DATA",
        "gbm_bursts":bursts,
        "v4_frozen_observed_json_sha256":V4_OBSERVED_SHA,
        "v4_frozen_manifest_sha256":V4_MANIFEST_SHA,
        "predeclared_diagnostic_bands_keV":{
            "nai_soft":[25,100],"nai_hard":[100,300],"bgo":[300,3000]},
        "prompt_windows_s_relative_trigger":{"090510":[0,2],"090902B":[0,40],"130427A":[0,140]},
        "background_window_relative_trigger_s":[-20,-5],
        "pulse_bin_width_s":{"090510":0.032,"090902B":0.256,"130427A":0.256},
        "claim_ceiling":"IRF_GRIDS_AND_REAL_GBM_TTE_ACCEPTED_BUT_NO_TRUE_PSF_EDISP_GBM_RSP2_FORWARD_FOLDED_LAT_GBM_LIV_LIKELIHOOD"}


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caldb-dir",type=Path,required=True)
    parser.add_argument("--gbm-dir",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    manifest=make_manifest(args.caldb_dir,args.gbm_dir)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(manifest,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":main()
