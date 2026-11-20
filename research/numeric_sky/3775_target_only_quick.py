#!/usr/bin/env python3
"""Quick, independently attestable R3775 NEOWISE + Pan-STARRS target lane.
No control claims: matched controls are checked by the full R3775 workflow.
"""
import importlib
import json
import pathlib
import sys
sys.path.insert(0,"research/numeric_sky")
m=importlib.import_module("3775_independent_epochs")
def main():
    q=m.fixed_source()
    anchor=m.query_allwise()
    p=m.P["preselected_allwise_neighbor"]
    target={"id":p["id"],"ra_deg":p["ra_deg"],"dec_deg":p["dec_deg"]}
    if anchor["anchor_status"]!="PASS_FROZEN_NEIGHBOR_MATCH":
        status="HOLD_FROZEN_ARCHIVE_SOURCE_MISSING";epoch=None;optical=None
    else:
        epoch=m.neowise_one(target)
        optical=m.panstarrs_one(target)
        status="TARGET_ONLY_ARCHIVED_WITH_NO_CONTROL_CLAIM"
    out={"status":status,"coordinate":q["coordinate"],"anchor_status":anchor["anchor_status"],
        "source_offset_arcsec":q["angular_separation_arcsec"],
        "neowise":{k:v for k,v in epoch.items() if k!="observations"} if epoch else None,
        "panstarrs":optical,"not_anomaly_test_without_matched_controls":True,
        "no_new_telescope_observation":True,"production_unchanged":True}
    m.save("target_only_quick.json",out)
    print("R3775_QUICK",json.dumps(out,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
