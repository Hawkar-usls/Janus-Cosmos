#!/usr/bin/env python3
"""Second independently implemented JSON cross-check of a V7 view against V4.

This is implementation-level independence, not independent telescope
observation, external lab replication, source class verification, or LIV test.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
V4=ROOT/"data/cosmos/COSMOS-ALI-V4-FERMI-RAW-EVENT-DIAGNOSTICS-2026-10-18.json"
MANIFEST=ROOT/"data/cosmos/COSMOS-ALI-V4-FERMI-RAW-FITS-MANIFEST-2026-10-18.json"


def reject(problem: str) -> None:
    raise ValueError("INDEPENDENT_IMPLEMENTATION_CONTROL_FAIL: "+problem)


def bytesha(path:Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(archive:dict,view:dict,manifest:dict,v4_file_sha:str,manifest_file_sha:str) -> dict:
    source=view.get("source",{})
    if source.get("v4_event_diagnostics_file_sha256") != v4_file_sha:
        reject("archived V4 file identity")
    if source.get("v4_manifest_file_sha256") != manifest_file_sha:
        reject("FITS manifest file identity")
    original=archive.get("sources",[])
    displayed=view.get("sources",[])
    if len(original)!=3 or len(displayed)!=3 or len(view.get("events",[]))!=175:
        reject("cohort size")
    original_bursts={x["grb"]:x for x in original}
    displayed_bursts={x["grb"]:x for x in displayed}
    if set(original_bursts)!=set(displayed_bursts):
        reject("burst labels")
    if {x["grb"] for x in manifest["bursts"]}!=set(original_bursts):
        reject("FITS source cohort")
    # Compare event multiset independently, using immutable numerical tuples.
    def identity(burst,e):
        return (burst,int(e["run_id"]),int(e["event_id"]),
                float(e["met_s"]),float(e["energy_MeV"]),
                float(e["after_trigger_s"]),float(e["angular_sep_deg"]),
                float(e["theta_deg"]),float(e["zenith_deg"]),
                int(e["event_class_mask"]),int(e["event_type_mask"]))
    expected=Counter(identity(s["grb"],e) for s in original for e in s["events"])
    actual=Counter(identity(e["grb"],e) for e in view["events"])
    if expected != actual:
        reject("event-level content changed or duplicated")
    if len(expected)!=175 or any(m!=1 for m in expected.values()):
        reject("duplicate archived event rows")
    if archive["counts"]["diagnostic_source_region_candidates"]!=len(expected):
        reject("V4 summary does not match events")
    roi=0
    for grb in sorted(original_bursts):
        s=original_bursts[grb]; v=displayed_bursts[grb]
        c=s["counts"]
        if v["counts"]["diagnostic_selected"]!=len(s["events"]) or v["counts"]["events_in_query_roi"]!=c["events_in_archive_roi"]:
            reject("count per GRB")
        if float(v["trigger_met_estimated"])!=float(s["trigger_met_estimated"]):
            reject("trigger estimate")
        for e in s["events"]:
            if not math.isclose(float(e["met_s"])-float(s["trigger_met_estimated"]),
                                float(e["after_trigger_s"]),abs_tol=3e-6):
                reject("MET/trigger mismatch")
        roi+=c["events_in_archive_roi"]
    if roi!=view["counts"]["archive_query_roi_events"] or roi!=6520:
        reject("query ROI denominator")
    if view.get("gates",{}).get("observational_LIV_p_value") is not None or view["gates"].get("lattice_spacing_fit") is not None:
        reject("forbidden scientific inference")
    if view["gates"].get("independent_external_replication")!="NOT_DONE":
        reject("false independent validation claim")
    return {"internal_second_implementation":"PASS","events_compared":len(expected),
            "query_roi_events_accounted":roi,"external_replication":"NOT_DONE"}


def main()->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--view",type=Path,required=True)
    args=parser.parse_args()
    try:
        archive=json.loads(V4.read_text(encoding="utf-8"))
        manifest=json.loads(MANIFEST.read_text(encoding="utf-8"))
        view=json.loads(args.view.read_text(encoding="utf-8"))
        print(json.dumps(compare(archive,view,manifest,bytesha(V4),bytesha(MANIFEST)),
                         ensure_ascii=False,sort_keys=True))
        return 0
    except (ValueError,KeyError,OSError,TypeError) as exc:
        print(str(exc));return 2

if __name__=="__main__":
    raise SystemExit(main())
