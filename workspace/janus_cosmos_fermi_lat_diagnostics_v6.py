#!/usr/bin/env python3
"""V6 LAT source/control ROI and spacecraft/GTI exposure-proxy pilot.

The proxy combines nearest-bin FRONT/BACK effective areas with real SC_DATA
livetime and boresight pointings. It omits azimuth/efficiency response terms,
PSF, energy dispersion, time-varying diffuse background and source assignment.
No event is counted as a proven GRB photon; no observed LIV p-value is made.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

from janus_cosmos_fermi_calibration_v5 import MANIFEST_PATH as CALDB_MANIFEST, area_lookup, audit_caldb
from janus_cosmos_fermi_raw_v4 import (
    MANIFEST as V4_MANIFEST, angular_sep_deg, ingest_one, read_fits, rows, sha256,
)
from janus_cosmos_simulation_audit import AuditError, ROOT, digest

SCHEMA = "janus.cosmos.fermi-lat-sc-gtis-approx-exposure-controls.v6"
OUT = ROOT / "data/cosmos/COSMOS-ALI-V6-LAT-SOURCE-CONTROL-EXPOSURE-2026-10-19.json"
PROMPT = {"090510": 2., "090902B": 40., "130427A": 140.}
ENERGY_SAMPLES_MEV = (1000., 10000.)


def intersection(*intervals: tuple[float, float]) -> float:
    return max(0., min(x[1] for x in intervals) - max(x[0] for x in intervals))


def exposure_proxy(sc: list[dict], gti: list[dict], interval: tuple[float, float],
                   ra: float, dec: float, aeff: dict) -> dict:
    """Real livetime but simplified boresight/IRF (NOT calibrated exposure)."""
    if interval[1] <= interval[0]:
        raise AuditError("nonpositive spacecraft integration interval")
    exposure = {str(int(energy)): 0. for energy in ENERGY_SAMPLES_MEV}
    gtis = [(float(x["START"]), float(x["STOP"])) for x in gti]
    exposure_sc_intervals = 0
    nominal_live = 0.
    for s in sc:
        t0, t1 = float(s["START"]), float(s["STOP"])
        if t0 >= t1 or s["DATA_QUAL"] != 1 or s["LAT_CONFIG"] != 1 or s["IN_SAA"]:
            continue
        inc = angular_sep_deg(float(s["RA_SCZ"]), float(s["DEC_SCZ"]), ra, dec)
        if inc > 65:
            continue
        span = sum(intersection(interval, (t0, t1), good) for good in gtis)
        if span <= 0:
            continue
        fraction = float(s["LIVETIME"]) / (t1 - t0)
        if not 0 <= fraction <= 1 + 1e-4:
            raise AuditError("spacecraft nominal livetime exceeds interval")
        nominal_live += span * min(1., fraction)
        exposure_sc_intervals += 1
        for energy in ENERGY_SAMPLES_MEV:
            area = area_lookup(aeff["FRONT"], energy, inc) + area_lookup(aeff["BACK"], energy, inc)
            exposure[str(int(energy))] += area * 10000 * span * min(1., fraction)
    return {"nearest_FB_aeff_times_livetime_cm2_s_by_energy_MeV": {k: round(v, 4) for k, v in exposure.items()},
            "quality_and_theta_selected_sc_intervals": exposure_sc_intervals,
            "quality_selected_on_axis_approx_livetime_s": round(nominal_live, 6),
            "is_real_calibrated_lat_exposure": False}


def run(caldb_manifest: dict, v4_manifest: dict, caldb_dir: Path, raw_dir: Path) -> dict:
    if v4_manifest["schema"] != "janus.cosmos.fermi-lat-fits-manifest.v4":
        raise AuditError("incorrect frozen V4 raw selection")
    if caldb_manifest.get("v4_frozen_manifest_sha256") != sha256(V4_MANIFEST):
        raise AuditError("V5 CALDB pin not tied to frozen V4")
    cal, aeff = audit_caldb(caldb_manifest, caldb_dir)
    base = json.loads((ROOT / "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-EVENT-DIAGNOSTICS-2026-10-18.json").read_text())
    out = []
    for config, prior in zip(v4_manifest["bursts"], base["sources"]):
        grb = config["grb"]
        if grb != prior["grb"]:
            raise AuditError("V4 source ordering changed")
        source = ingest_one(config, raw_dir)
        if source["events"] != prior["events"]:
            raise AuditError("V4 source-region photon events not exactly reproduced")
        extended = copy.deepcopy(config)
        extended["selection"]["max_roi_deg"] = 10
        whole = ingest_one(extended, raw_dir)
        inner_ids = {(v["run_id"], v["event_id"]) for v in source["events"]}
        wider_ids = {(v["run_id"], v["event_id"]) for v in whole["events"]}
        if not inner_ids.issubset(wider_ids):
            raise AuditError("5deg source region not a subset of 10deg archive")
        annulus = [ev for ev in whole["events"] if (ev["run_id"], ev["event_id"]) not in inner_ids]
        if not all(5 < float(ev["angular_sep_deg"]) <= 10 for ev in annulus):
            raise AuditError("LAT control annulus geometrically malformed")
        sc_path = raw_dir / config["spacecraft"]["filename"]
        ev_path = raw_dir / config["event"]["filename"]
        sc = rows(next(h for h in read_fits(sc_path) if h.header.get("EXTNAME") == "SC_DATA"),
                  {"START", "STOP", "LIVETIME", "DATA_QUAL", "LAT_CONFIG", "IN_SAA", "RA_SCZ", "DEC_SCZ"})
        gti = rows(next(h for h in read_fits(ev_path) if h.header.get("EXTNAME") == "GTI"), {"START", "STOP"})
        trigger = float(next(b["exact_gbm_trigger_met_s"] for b in json.loads((ROOT / "data/cosmos/COSMOS-ALI-V5-CALDB-GBM-LAT-INSTRUMENT-PULSE-2026-10-19.json").read_text())["gbm_bursts"] if b["grb"] == grb))
        tstop = PROMPT[grb]
        windows = {"pre_source_control": (-100., -20.), "prompt": (0., tstop), "post_source_control": (tstop + 20., tstop + 120.)}
        summaries = []
        for name, (lo, hi) in windows.items():
            proxy = exposure_proxy(sc, gti, (trigger + lo, trigger + hi), config["ra_deg"], config["dec_deg"], aeff)
            def count(items):
                return sum(trigger + lo <= float(e["met_s"]) < trigger + hi for e in items)
            summaries.append({"time_window": name, "relative_to_exact_gbm_trigger_s": [lo, hi],
                              "source_circle_5deg_count": count(source["events"]),
                              "control_annulus_5to10deg_count": count(annulus),
                              **proxy})
        types = []
        for event in source["events"]:
            mask = int(event["event_type_mask"])
            if sum(bool(mask & bit) for bit in (1, 2)) != 1 or sum(bool(mask & bit) for bit in (4, 8, 16, 32)) != 1 or sum(bool(mask & bit) for bit in (64, 128, 256, 512)) != 1:
                raise AuditError("LAT Pass8 EVENT_TYPE partition invalid")
            types.append(mask)
        out.append({"grb": grb, "redshift": config["redshift"], "exact_gbm_trigger_met_s": trigger,
                    "five_degree_source_region_candidate_total": len(inner_ids),
                    "five_to_ten_degree_control_annulus_candidate_total": len(annulus),
                    "control_annulus_over_source_solid_angle_ratio": round((math.cos(math.radians(5)) - math.cos(math.radians(10))) / (1 - math.cos(math.radians(5))), 6),
                    "three_time_windows": summaries,
                    "event_type_3_partitions_decoded_count": len(types),
                    "LAT_PSF_EDISP_SOURCE_POSTERIOR_NOT_COMPUTED": True,
                    "post_source_control_can_contain_afterglow": True})
    result = {"schema": SCHEMA, "experiment_id": "JANUS-COSMOS-ALI-V6-LAT-GBM-JOINT-2026-10-19",
              "counts": {"v4_lat_source_region_events_recovered_exactly": sum(x["five_degree_source_region_candidate_total"] for x in out),
                         "num_bursts": 3},
              "lat_irf_set": cal["irf_name"], "SC_AND_GTI_REAL": True,
              "approx_exposure_source": "REAL_LIVETIME_AND_POINTING_NEAREST_FB_EFFECTIVE_AREA_ONLY_NO_PHI_OR_DEADTIME_EFFICIENCY_CORRECTIONS",
              "source_control_background": "5DEG_CIRCLE_VS_5TO10DEG_ANNULUS_3_WINDOWS_NO_DIFFUSE_TEMPLATE",
              "photon_direction_energy_irf": "NOT_EVALUATED_PSF_AND_ENERGY_DISPERSION_PENDING",
              "GRB_sources": out,
              "scientific_gates": {"FROZEN_V4_REAL_LAT_EVENTS": "PASS_EXACT_REPRODUCTION",
                                   "LAT_SC_GTI_LIVETIME_BORESIGHT": "PASS_APPROX_PROXY_ONLY",
                                   "LAT_TRUE_FORWARD_FOLDED_SOURCE_LIKELIHOOD": "BLOCKED_PSF_EDISP_AND_BACKGROUND",
                                   "OBSERVED_JOINT_LIV_TEST": "BLOCKED_INSTRUMENT_AND_SOURCE_LAG_MODEL",
                                   "AMPS_FIREWALL": "OPEN_UNCONFIRMED",
                                   "BHB_BOOM": "NOT_ESTABLISHED"},
              "real_joint_observed_liv_p_value": None,
              "claim_ceiling": "REAL_EXPOSURE_CONTROL_DIAGNOSTICS_NOT_LAT_PSF_EDISP_CALIBRATED_LIKELIHOOD"}
    result["result_sha256"] = digest(result)
    return result


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--caldb-dir", type=Path, required=True)
    p.add_argument("--raw-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, default=OUT)
    a = p.parse_args()
    m = json.loads(CALDB_MANIFEST.read_text())
    v4 = json.loads(V4_MANIFEST.read_text())
    result = run(m, v4, a.caldb_dir, a.raw_dir)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n")
    print(json.dumps({"output": str(a.output), "counts": result["counts"], "gates": result["scientific_gates"]}))


if __name__ == "__main__":
    main()
