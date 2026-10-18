#!/usr/bin/env python3
"""JANUS Cosmos V3: real Fermi-LAT published summaries, NO LIV inference.

This only overlays specified V2 lattice forward predictions on a documented
transcription of published 2FLGC Table 7 high-energy photons. It explicitly
blocks physical speed or spacing estimation from trigger-relative times:
emission lags, full LAT event data, instrument responses and a joint source
likelihood have not been supplied. No synthetic event is called observed.
"""

from __future__ import annotations

import argparse
import json
import sys
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any

from janus_cosmos_lattice_tof import AuditError, PRECISION, ROOT, file_sha256, frozen_v1_gate, lattice_delay_seconds
from janus_cosmos_simulation_audit import digest, num, show

SCHEMA = "janus.cosmos.fermi-published-tof-method-gate.v3"
CATALOG_SCHEMA = "janus.cosmos.published-fermi-photon-summary.v3"
CATALOG_PATH = ROOT / "data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-HE-PHOTONS-2026-10-18.json"
CATALOG_SHA256 = "77893a8cb8367d8b78f4d0bde60cfe790c3d67c2b4291c8bbdfaf34ba7b7982a"
FROZEN_V2_PATH = ROOT / "data/cosmos/COSMOS-ALI-V2-LATTICE-TOF-REVEAL-2026-10-18.json"
FROZEN_V2_SHA256 = "06558fae6e8fbc63c86ad19f3ed1f3cd102d3ecbe3f6da43275468b9419a5ff6"
PAPER_HTML_URL = "https://arxiv.org/html/1906.11403v1"
OFFICIAL_COLUMNS_URL = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/lat_data_columns.html"

EXPECTED_GRBS = {
    "080916C":2,"090510":1,"090902B":7,"090926A":2,"100414A":2,
    "130427A":17,"131231A":2,"160509A":2,"160623A":1,
    "160625B":1,"171010A":1,
}


def verify_history(v2_file:Path=FROZEN_V2_PATH) -> dict[str,str]:
    prior=frozen_v1_gate()
    if file_sha256(v2_file)!=FROZEN_V2_SHA256:
        raise AuditError("FROZEN_COSMOS_V2_REVEAL_SHA256_MISMATCH")
    v2=json.loads(v2_file.read_text(encoding="utf-8"))
    if v2.get("gates",{}).get("ACTUAL_LIV_OBSERVED")!="NOT_ESTABLISHED" or v2.get("gates",{}).get("QUANTUM_FIREWALL")!="OPEN_UNCONFIRMED":
        raise AuditError("V2_OBSERVATIONAL_OR_AMPS_GATE_PROMOTED")
    return {"frozen_v1_file_sha256":prior["frozen_v1_file_sha256"],
            "frozen_v2_reveal_file_sha256":FROZEN_V2_SHA256,
            "liv_status":"NOT_ESTABLISHED","amps_status":"OPEN_UNCONFIRMED"}


def validate_catalog(raw:dict[str,Any]) -> dict[str,Any]:
    if not isinstance(raw,dict) or raw.get("schema")!=CATALOG_SCHEMA or raw.get("experiment_id")!="JANUS-COSMOS-ALI-V3-FERMI-PHOTON-TOF-2026-10-18":
        raise AuditError("V3 catalog schema and experiment identity required")
    if raw.get("provenance_status")!="TRANSCRIBED_FROM_PRIMARY_PUBLISHED_FERMI_LAT_SECOND_GRB_CATALOG_TABLE_7":
        raise AuditError("no fabricated real-data provenance")
    src=raw.get("source",{})
    if not isinstance(src,dict) or src.get("versioned_html_table_url")!=PAPER_HTML_URL or src.get("official_event_columns")!=OFFICIAL_COLUMNS_URL:
        raise AuditError("primary paper and NASA event-column source must match")
    if src.get("table")!="Table 7: Fermi-LAT GRBs with photon energies E>10 GeV" or src.get("event_catalog_fits")!="NOT_DOWNLOADED_OR_VERIFIED":
        raise AuditError("raw LAT FITS files have not been retrieved or verified")
    if src.get("time_zero")!="SECONDS_AFTER_FERMI_GBM_TRIGGER_NOT_SOURCE_EMISSION_TIME" or src.get("precision")!="PUBLISHED_ROUNDED_TABLE_VALUES_NOT_FULL_PRECISION_TIMESTAMPS":
        raise AuditError("trigger-relative and rounded nature of times must be declared")
    if src.get("association_selection")!="P_ASSOCIATED_GREATER_THAN_0_90_PER_PUBLISHED_TABLE_CRITERION_NOT_PER_EVENT_PROBABILITY":
        raise AuditError("catalog association selection cannot be forged as individual probabilities")
    p=raw.get("analysis_protocol",{})
    if not isinstance(p,dict) or p.get("reference_energy_status")!="HYPOTHETICAL_ONLY_NO_MATCHED_LOW_ENERGY_PHOTON_MEASURED" or p.get("inference_ceiling")!="PREDICTION_OVERLAYS_AND_SOURCE_ORDER_CONTROLS_ONLY_NO_ASTROPHYSICAL_LIV_BOUND":
        raise AuditError("hypothetical reference energy and no-LIV-limit ceiling required")
    if p.get("true_event_pair_lag_available") is not False or p.get("source_intrinsic_lag_calibrated") is not False or p.get("instrument_response_reconstructed") is not False or p.get("raw_event_files_ingested") is not False or p.get("statistical_trials_calibrated") is not False:
        raise AuditError("missing astrophysical data/control must not be marked available")
    if p.get("prohibited_operations")!=["FIT_LATTICE_SPACING_TO_TRIGGER_RELATIVE_SINGLE_PHOTON_TIME",
                                      "ASSUME_GBM_TRIGGER_IS_EMISSION_OF_A_LOW_ENERGY_PHOTON",
                                      "REPORT_LIV_OR_SIMULATION_DETECTION"]:
        raise AuditError("non-identifiable fit may not be enabled")
    ref=num(p.get("reference_energy_eV"),"reference_energy_eV")
    if ref!=Decimal("1e9") or p.get("lattice_spacings_m")!=["4e-27","1.616255e-35"]:
        raise AuditError("V3 preregistered diagnostic parameters must match")
    cos=p.get("cosmology",{})
    if not isinstance(cos,dict) or {k:cos.get(k) for k in ("H0_km_s_Mpc","Omega_m","Omega_lambda")}!={"H0_km_s_Mpc":"70","Omega_m":"0.3","Omega_lambda":"0.7"}:
        raise AuditError("cosmology must retain the V2 illustrative model")
    bursts=raw.get("bursts")
    if not isinstance(bursts,list) or len(bursts)!=11:
        raise AuditError("expected 11 redshift-known primary-table GRBs")
    found={}
    for burst in bursts:
        if not isinstance(burst,dict) or not isinstance(burst.get("grb"),str):
            raise AuditError("invalid GRB identifier")
        name=burst["grb"]
        if name in found or name not in EXPECTED_GRBS:
            raise AuditError("duplicate or non-preregistered GRB")
        if burst.get("class")!=("SHORT" if name=="090510" else "LONG"):
            raise AuditError("unexpected GRB source classification")
        z=num(burst.get("redshift"),"redshift")
        if not Decimal(0)<z<=Decimal(5):
            raise AuditError("published GRB redshift must be positive and finite")
        events=burst.get("events")
        if not isinstance(events,list) or len(events)!=EXPECTED_GRBS[name]:
            raise AuditError("truncated or fabricated high-energy event counts")
        previous_energy=None
        for item in events:
            if not isinstance(item,dict) or set(item)!={"energy_GeV","time_s"}:
                raise AuditError("only primary table energy and trigger-relative time may be claimed")
            energy=num(item["energy_GeV"],"energy_GeV")
            t=num(item["time_s"],"time_s")
            if energy<=10 or t<0 or (previous_energy is not None and energy>=previous_energy):
                raise AuditError("GRB photons must be above 10 GeV, sorted by descending energy, and time>=0")
            previous_energy=energy
        found[name]=True
    if set(found)!=set(EXPECTED_GRBS):
        raise AuditError("nonmatching set of published bursts")
    return raw


def evaluate_catalog(raw:dict[str,Any],catalog_hash:str=CATALOG_SHA256,v2_file:Path=FROZEN_V2_PATH) -> dict[str,Any]:
    validate_catalog(raw)
    history=verify_history(v2_file)
    cos=raw["analysis_protocol"]["cosmology"]
    ref=num(raw["analysis_protocol"]["reference_energy_eV"],"reference_energy_eV")
    result_bursts=[]
    total_photons=total_pairs=earlier=later=0
    for burst in raw["bursts"]:
        z=num(burst["redshift"],"redshift")
        events=[]
        for index,e in enumerate(burst["events"]):
            energy=num(e["energy_GeV"],"energy_GeV")
            time=num(e["time_s"],"time_s")
            predictions={}
            with localcontext() as ctx:
                ctx.prec=PRECISION
                for spacing in raw["analysis_protocol"]["lattice_spacings_m"]:
                    c={"redshift":burst["redshift"],"lattice_spacing_m":spacing,
                       "energy_low_eV":show(ref),"energy_high_eV":show(energy*Decimal("1e9"))}
                    predictions[spacing]=show(lattice_delay_seconds(c,cos))
            events.append({"index_in_published_energy_descending_order":index,
                           "measured_energy_GeV":e["energy_GeV"],
                           "measured_arrival_after_GBM_trigger_s":e["time_s"],
                           "published_event_class":"HIGH_ENERGY_EVENT_P_GT_0_90_SELECTION",
                           "lattice_predicted_delay_vs_unobserved_1GeV_reference_s":predictions,
                           "fitted_lattice_spacing_m":None,
                           "physical_speed_inference":"BLOCKED_UNKNOWN_SOURCE_EMISSION_TIME"})
        inversions=0
        ordered=0
        for i in range(len(events)):
            for j in range(i+1,len(events)):
                # High-energy photon at i, lower-energy photon at j.
                delta=num(events[i]["measured_arrival_after_GBM_trigger_s"],"t_hi") - num(events[j]["measured_arrival_after_GBM_trigger_s"],"t_lo")
                if delta<0: inversions+=1
                elif delta>0: ordered+=1
                else: raise AuditError("identical published times preclude a strict order")
        pairs=inversions+ordered
        total_pairs+=pairs
        earlier+=inversions
        later+=ordered
        total_photons+=len(events)
        result_bursts.append({"grb":burst["grb"],"published_type":burst["class"],
           "published_redshift":burst["redshift"],"photon_count":len(events),
           "source_table":"AJELLO_ET_AL_2019_TABLE_7_ROUNDED",
           "all_photon_times_reference":"GBM_TRIGGER_NOT_EMISSION_TIME",
           "high_vs_low_energy_time_order_pairs":{
               "higher_energy_arrived_earlier":inversions,
               "higher_energy_arrived_later":ordered,
               "non_simultaneity_confounder":"NO_COMMON_EMISSION_TIME_INFERRED"},
           "photons":events,
           "source_delay_estimated":False,
           "physical_model_preference":"UNIDENTIFIABLE_WITH_UNBOUNDED_SOURCE_LAGS"})
    result={"schema":SCHEMA,"experiment_id":raw["experiment_id"],
       "observation_type":"REAL_PUBLISHED_FERMI_LAT_2FLGC_ROUNDED_TABLE_7_VALUES_NOT_RAW_EVENT_FITS",
       "source_paper_table":raw["source"],
       "source_input_file_sha256":catalog_hash,
       "frozen_history":history,
       "model":"V2_1D_LATTICE_RADIAL_FLRW_ANSATZ_CONDITIONAL_PROJECTION",
       "cosmology":cos,
       "reference_energy_eV":show(ref),"reference_photon":"HYPOTHETICAL_NOT_AN_OBSERVED_EVENT",
       "selection":"CURATED_11_GRBS_NOT_COMPLETENESS_CORRECTED_POPULATION",
       "grbs":result_bursts,
       "counts":{"bursts":len(result_bursts),"published_high_energy_photons":total_photons,
                 "within_burst_energy_order_pairs":total_pairs,
                 "higher_energy_arrived_earlier":earlier,
                 "higher_energy_arrived_later":later,
                 "events_with_verified_emission_time":0,"matched_low_energy_photons":0},
       "negative_control_example":{"grb":"090902B","higher_energy_GeV":"40.0",
          "higher_energy_time_s":"81.7","lower_energy_GeV":"21.7",
          "lower_energy_time_s":"331.9","observed_pair_time_difference_high_minus_low_s":"-250.2",
          "interpretation":"A HIGHER-ENERGY photon arrived earlier; this neither proves nor falsifies dispersion without source emission lags"},
       "gates":{"PUBLISHED_PRIMARY_TABLE_NUMBERS":"PASS_CURATED_TRANSCRIPTION_WITH_SOURCE_URL",
                "REAL_ASTROPHYSICAL_MEASUREMENTS":"PRESENT_AS_PUBLISHED_ROUNDED_SUMMARIES",
                "FULL_RAW_LAT_EVENT_FITS":"NOT_ACQUIRED",
                "LOW_ENERGY_CORRELATED_EVENT":"NOT_ACQUIRED_HYPOTHETICAL_1GEV_REFERENCE",
                "INTRINSIC_SOURCE_LAG":"UNBOUNDED_UNKNOWN",
                "INSTRUMENT_RESPONSE":"NOT_RECONSTRUCTED",
                "ASSOCIATION_PROBABILITY":"TABLE_THRESHOLD_ONLY_NOT_PER_PHOTON",
                "LAT_DETECTION_SELECTION":"NOT_MODELED",
                "REAL_LIV_LIKELIHOOD_FIT":"BLOCKED_NOT_IDENTIFIABLE",
                "LATTICE_SPACING_CONSTRAINT":"NOT_DERIVED",
                "LIV_OBSERVED":"NOT_ESTABLISHED",
                "SIMULATION_HYPOTHESIS":"NOT_ESTABLISHED",
                "AMPS_FIREWALL":"OPEN_UNCONFIRMED",
                "BHB_BOOM":"NOT_ESTABLISHED",
                "OSIRIS_PROOF_AUTHORITY":"NO_CHANGE"},
       "next_gate":"RETRIEVE_SHA_PINNED_LAT_EXTENDED_EVENT_AND_SPACECRAFT_FITS_PLUS_GBM_PROMPT_DATA_FOR_MULTIBURST_FORWARD_FOLDED_SOURCE_AND_INSTRUMENT_LIKELIHOOD",
       "claim_ceiling":"These are real cited catalog values but only forward-prediction overlays. No LIV test statistic, spacetime lattice bound, simulation proof, quantum firewall or black-hole explosion is inferred."}
    result["result_sha256"]=digest(result)
    return result


def main(argv:list[str]|None=None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--batch",type=Path,default=CATALOG_PATH)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args(argv)
    try:
        if file_sha256(args.batch)!=CATALOG_SHA256:
            raise AuditError("V3_PRIMARY_SOURCE_TRANSCRIPTION_SHA256_MISMATCH")
        raw=json.loads(args.batch.read_text(encoding="utf-8"))
        out=evaluate_catalog(raw,file_sha256(args.batch))
        body=json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n"
        if args.output: args.output.write_text(body,encoding="utf-8")
        else: sys.stdout.write(body)
        return 0
    except (OSError,ValueError,KeyError) as exc:
        print(f"JANUS-FERMI-V3-ERROR: {exc}",file=sys.stderr)
        return 2


if __name__=="__main__":
    raise SystemExit(main())
