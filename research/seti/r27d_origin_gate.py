#!/usr/bin/env python3
"""Real-control provenance and leave-source-out feasibility gate, not invented origin likelihoods."""
import json,pathlib
import r26c_reference_alignment as r
import donor_eye_r4_4 as eye
OUT=pathlib.Path("runtime/r27d");OUT.mkdir(parents=True,exist_ok=True)
def main():
    r.BLC_OUT=OUT
    sources={"blc":{"id":14366830512,"digest":"bf50c2db96dad4dae126921272e06b7c3a41eb3657a9d011325d524bc630e8a1"},"voyager":{"id":14367546526,"digest":"106f34686ad948a1eef91ea14ad60fb3939cd1d0e9baf34bdb88c921e1b7648e"}}
    records=[]
    for k,v in sources.items():
        r.BLC_P["artifacts"][k+"_r27d"]=v;z,n,p=r.archive(k+"_r27d")
        d=r.get(z,n,"blc_result.json" if k=="blc" else "result.json")
        if k=="blc":
            assert d["anchor_raw_recovered"]
            records.append({"group":"BLC1_2019_Parkes","label":"TERRESTRIAL_RFI","reference":"https://www.nature.com/articles/s41550-021-01508-8","measurements":len(d["rows"]),"known_positive_detections":1,"receiver":"Parkes","label_basis":"published intermodulation study; related scans not independent causes"})
        else:
            assert d["status"]=="KNOWN_CARRIER_RECOVERED"
            records.append({"group":"Voyager1_GBT_sample","label":"HUMAN_SPACECRAFT","reference":"https://github.com/UCBerkeleySETI/turbo_seti","measurements":1,"known_positive_detections":1,"receiver":"GBT","label_basis":"published known transmitter; receiver and source confounded"})
        z.close();p.unlink()
    splits=[]
    for q in records:
        train=[x for x in records if x["group"]!=q["group"]];labels={x["label"] for x in train};splits.append({"held_out":q["group"],"training_labels":sorted(labels),"status":"HOLD_MISSING_CLASS" if q["label"] not in labels else "POSSIBLE"})
    model={"id":"JANUS-R27D-REAL-PROVENANCE-GATE","cause_classes":[{"id":"terrestrial_rfi"},{"id":"human_spacecraft"},{"id":"other_unknown"}],"tests":[{"id":"morphology","available":True,"decision_usable":False,"calibrated":False},{"id":"receiver_independent_repeat","available":False,"decision_usable":False,"calibrated":False}],"model_set":[{"id":"equal","provenance":"declared sensitivity, not empirical prior"},{"id":"rfi_dominant","provenance":"declared sensitivity, not empirical prior","prior_by_cause":{"terrestrial_rfi":.8,"human_spacecraft":.1,"other_unknown":.1}}]}
    receipt=eye.solve(model);assert receipt["status"]=="ROBUST_NON_IDENTIFIABLE_UNDER_DECLARED_MODEL_SET"
    d={"status":"REAL_CALIBRATION_HOLD_INSUFFICIENT_SOURCE_GROUPS","records":records,"source_group_count":len(records),"leave_source_out":splits,"eye_receipt":receipt,"scope":"Two known source groups from different receiver domains do not support empirical origin likelihoods. Repeated windows do not constitute independent labels. No true extraterrestrial positive. No source classification or production adoption.","production_adoption":False}
    (OUT/"result.json").write_text(json.dumps(d,indent=2)+"\n");print("R27D_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
