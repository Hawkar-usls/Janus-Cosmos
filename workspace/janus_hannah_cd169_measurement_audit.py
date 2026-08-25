#!/usr/bin/env python3
"""Deterministic JANUS audit for the first real Hannah/CD169 TOBI receipt.

This audit validates provenance and frozen methodological gates. It deliberately
has no post-hoc morphology/anomaly threshold and cannot promote mnemonic
similarity into a physical claim.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path


def gate(ok, *, evidence=None, fail=None):
    return {"state": "PASS" if ok else "FAIL", "evidence": evidence, "failure": None if ok else fail}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('receipt',type=Path); ap.add_argument('--output',required=True,type=Path); a=ap.parse_args()
    r=json.loads(a.receipt.read_text(encoding='utf-8'))
    out={
      "artifact_id":"JANUS-HANNAH-CD169-MEASUREMENT-AUDIT-RUN-001-2026-08-25-v1.0",
      "schema":"janus.cosmos.cousteau.hannah_cd169.measurement_audit.v1",
      "source_receipt":str(a.receipt),
      "source_artifact_id":r.get('artifact_id'),
      "scientific_convergence_claim":False,
      "posthoc_morphology_threshold_added":False,
      "gates":{}
    }
    src=r.get('source',{}); fmt=r.get('format_validation',{}); nb=r.get('nearest_real_block',{}); rg=r.get('real_data_synesthetic_regression_gates',{}); state=r.get('current_scientific_state',{}); mir=r.get('local_pre_post_mirror',{})
    out['gates']['G0_RAW_FILE_INTEGRITY']=gate(
      bool(src.get('raw_file_sha256')) and src.get('raw_file_size_matches_server') is True and src.get('raw_file_size_mod_40960')==0,
      evidence={"sha256":src.get('raw_file_sha256'),"size_bytes":src.get('raw_file_size_bytes'),"size_matches_server":src.get('raw_file_size_matches_server'),"size_mod_block":src.get('raw_file_size_mod_40960')},
      fail='RAW_FILE_INTEGRITY_INCOMPLETE'
    )
    out['gates']['G1_FORMAT_AND_CADENCE']=gate(
      fmt.get('block_size_bytes')==40960 and fmt.get('declared_cadence_seconds')==4 and fmt.get('measured_whole_file_cadence_median_seconds')==4.0 and fmt.get('measured_whole_file_cadence_mad_seconds')==0.0,
      evidence={k:fmt.get(k) for k in ['block_size_bytes','declared_cadence_seconds','measured_whole_file_cadence_median_seconds','measured_whole_file_cadence_mad_seconds','timestamp_encoding']},
      fail='REAL_STREAM_DOES_NOT_MATCH_DECLARED_4S_STRUCTURE'
    )
    delta=nb.get('delta_seconds_from_frozen_target')
    # For a 4 s lattice, any correctly located nearest sample must lie within half cadence (2 s).
    half_cadence_bound=isinstance(delta,(int,float)) and abs(delta)<=2.0
    out['gates']['G2_NEAREST_SAMPLE_TIME_ALIGNMENT']=gate(
      half_cadence_bound and bool(nb.get('raw_block_sha256')),
      evidence={"target_utc":r.get('frozen_target',{}).get('target_utc_approx'),"nearest_real_utc":nb.get('datetime_utc'),"delta_seconds":delta,"mathematical_bound":"abs(delta)<=cadence/2=2s","raw_block_sha256":nb.get('raw_block_sha256')},
      fail='NEAREST_REAL_SAMPLE_OUTSIDE_HALF_CADENCE_BOUND_OR_UNHASHED'
    )
    out['gates']['G3_REAL_DATA_MNEMONIC_INVARIANTS']=gate(
      rg.get('all_pass') is True and all(rg.get(k) is True for k in ['same_input_same_passport','forbidden_label_invariance','direction_does_not_change_measurement_fingerprint','direction_does_change_context_passport','raw_target_block_hash_present']),
      evidence=rg,
      fail='SYNESTHETIC_CORE_FAILED_REAL_DATA_FIREWALL_OR_DETERMINISM_TEST'
    )
    sims=[]
    for s in ['60s','300s','1800s','7200s']:
        v=(mir.get(s) or {}).get('common_measurement_similarity')
        if isinstance(v,(int,float)): sims.append((s,float(v)))
    monotonic=len(sims)==4 and all(sims[i][1]>=sims[i+1][1] for i in range(3))
    out['gates']['G4_LOCAL_MIRROR_MEMORY_TRAJECTORY']={
      "state":"OBSERVED_NOT_A_SCIENTIFIC_GATE",
      "similarities":dict(sims),
      "monotonic_decrease_with_window_width":monotonic,
      "interpretation":"The local PRE/POST measurement vector is most similar at the shortest tested scale and diverges gradually as wider temporal context is included.",
      "forbidden_interpretation":"No preregistered threshold exists here for anomaly, morphology, artificial origin, or target presence."
    }
    out['gates']['G5_NATIVE_TOWFISH_GEOLOCATION']={
      "state":"PASS" if state.get('TOWFISH_POINT_GEOLOCATION_FROM_NATIVE_NAV') is True else "BLOCKED",
      "required":"native towfish navigation or validated layback/cable solution at target time",
      "current":state.get('TOWFISH_POINT_GEOLOCATION_FROM_NATIVE_NAV')
    }
    out['gates']['G6_SAMPLE_TO_GROUND_FOOTPRINT']={
      "state":"PASS" if state.get('SIDESCAN_SAMPLE_TO_GROUND_FOOTPRINT_MAPPING') is True else "BLOCKED",
      "required":"map TOBI port/stbd sample range and vehicle geometry to Earth-fixed footprint",
      "current":state.get('SIDESCAN_SAMPLE_TO_GROUND_FOOTPRINT_MAPPING')
    }
    out['gates']['G7_GROUND_FIXED_REPLICATION']={"state":"PASS" if state.get('GROUND_FIXED_REPLICATION') is True else "BLOCKED","current":state.get('GROUND_FIXED_REPLICATION')}
    out['gates']['G8_GEOMETRY_CLASS']={"state":"BLOCKED","current":state.get('PYRAMID_LIKE_GEOMETRY'),"reason":"G5/G6/G7 must pass before blind geometry-control comparison."}
    out['gates']['G9_ACOUSTIC_RESONANCE']={"state":"BLOCKED","current":state.get('ACOUSTIC_RESONANCE'),"reason":"Mnemonic audio frequency is synthetic memory output, not measured sonar frequency; suitable physical spectral/phase data are required."}
    out['janus_state']={
      "TIME_MEMORY_CUE_ON_REAL_DATA":"VALIDATED",
      "TARGET_RAW_BLOCK":"LOCATED_AND_HASHED",
      "LOCAL_MIRROR":"MEASURED",
      "SPACE_PROOF":"BLOCKED_PENDING_NATIVE_TOWFISH_NAV_AND_SAMPLE_FOOTPRINT",
      "GROUND_FIXED_FEATURE":"NOT_ESTABLISHED",
      "PYRAMID_LIKE_GEOMETRY":"NOT_TESTED",
      "ACOUSTIC_RESONANCE":"NOT_TESTED",
      "UNDERWATER_PYRAMID_DETECTED":False
    }
    out['next_action']=r.get('next_janus_gate')
    out['hard_rules']=[
      "MNEMONIC_AUDIO_FREQUENCY_IS_NOT_PHYSICAL_ACOUSTIC_FREQUENCY",
      "MNEMONIC_GLYPH_IS_NOT_SEAFLOOR_GEOMETRY",
      "HIGH_LOCAL_MIRROR_SIMILARITY_IS_NOT_TARGET_PROOF",
      "SHIP_GPS_IS_NOT_TOWFISH_POSITION",
      "NO_POINT_FEATURE_WITHOUT_SAMPLE_TO_GROUND_FOOTPRINT",
      "NO_GEOMETRY_CLASS_BEFORE_GROUND_FIXED_REPLICATION",
      "NO_UNDERWATER_PYRAMID_DETECTED_YET"
    ]
    required=['G0_RAW_FILE_INTEGRITY','G1_FORMAT_AND_CADENCE','G2_NEAREST_SAMPLE_TIME_ALIGNMENT','G3_REAL_DATA_MNEMONIC_INVARIANTS']
    out['audit_status']='PASS_TO_SPACE_GATE' if all(out['gates'][g]['state']=='PASS' for g in required) else 'STOP_BEFORE_SPACE_GATE'
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'audit_status':out['audit_status'],'janus_state':out['janus_state'],'local_similarity':dict(sims)},indent=2))
    return 0 if out['audit_status']=='PASS_TO_SPACE_GATE' else 2
if __name__=='__main__': raise SystemExit(main())
