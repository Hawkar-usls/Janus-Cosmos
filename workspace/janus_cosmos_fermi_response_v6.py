#!/usr/bin/env python3
"""JANUS Cosmos V6: real GBM RSP2/Poisson on-off forward-folding pilot.

This is a response-aware *spectral diagnostic*, not a calibrated LAT/GBM LIV
likelihood. The GBM photon law (single cutoff power law) is falsifiable and
usually imperfect for GRBs with multiple spectral components. No p-value is
reported from uncalibrated model selection or source-intrinsic lag priors.
"""

from __future__ import annotations

import argparse
import json
import math
import struct
from pathlib import Path

from janus_cosmos_fermi_calibration_v5 import (
    MANIFEST_PATH, _tte_events, _valid_ebounds, checked_file, fits_large_tte,
    get_hdu, read_fits,
)
from janus_cosmos_simulation_audit import AuditError, ROOT, digest

SCHEMA = "janus.cosmos.fermi-gbm-rsp2-onoff-folded-spectral-pilot.v6"
OUT = ROOT / "data/cosmos/COSMOS-ALI-V6-GBM-RSP2-ONOFF-SPECTRAL-2026-10-19.json"

# Coarse windows are fixed before examining the spectral-fit statistics.
WINDOWS = {
    "090510": [(0., .25), (.25, .50), (.50, .75), (.75, 1.), (1., 1.5), (1.5, 2.)],
    "090902B": [(0., 5.), (5., 10.), (10., 20.), (20., 30.), (30., 40.)],
    "130427A": [(0., 10.), (10., 20.), (20., 40.), (40., 70.), (70., 100.), (100., 140.)],
}
ALPHAS = [-1.6, -1.2, -0.8, -0.4]
EPEAK_KEV = [150., 300., 600., 1200., 2400., 4800.]
BG_WINDOWS = [(-20., -5.)]
MIN_BACK_EXPOSURE_S = 10.


def decode_drm(hdu) -> dict:
    """Decode OGIP variable-length PI/PE arrays and the RSP2 FITS heap.

    The table has exactly 140 true-energy rows, with start-channel and
    contiguous n-channel groups. This bounded decoder refuses all other
    variable array layouts, malformed heap pointers, or duplicate channels.
    """
    h = hdu.header
    if h.get("EXTNAME") != "SPECRESP MATRIX" or int(h.get("NAXIS2", -1)) != 140:
        raise AuditError("unexpected GBM RSP2 DRM identity/energy row count")
    forms = [h.get(f"TFORM{i}") for i in range(1, 7)]
    if forms[:3] != ["1E", "1E", "1I"] or forms[3] not in {"PI(1)", "PI(128)"} or forms[4] not in {"PI(1)", "PI(128)"} or forms[5] != "PE(128)":
        raise AuditError("unexpected GBM RSP2 variable-length array format")
    width = int(h.get("NAXIS1", -1))
    nrows = int(h["NAXIS2"])
    if width != 34 or int(h.get("DETCHANS", -1)) != 128:
        raise AuditError("wrong GBM RSP2 channel/row extent")
    data = hdu.content
    heap_start = int(h.get("THEAP", width * nrows))
    if heap_start < width * nrows or heap_start > len(data):
        raise AuditError("invalid RSP2 heap base")
    matrix = []
    energies = []
    prev = 0.

    def vector(row: int, off: int, fmt: str, elt_size: int) -> tuple:
        length, displacement = struct.unpack_from(">ii", data, row * width + off)
        start = heap_start + displacement
        if length < 0 or length > 16384 or displacement < 0 or start + length * elt_size > len(data):
            raise AuditError("RSP2 variable array pointer escaped FITS heap")
        return tuple(v[0] for v in struct.iter_unpack(fmt, data[start:start + length * elt_size]))

    for row in range(nrows):
        offset = row * width
        elo, ehi, groups = struct.unpack_from(">ffh", data, offset)
        if not (math.isfinite(elo) and math.isfinite(ehi) and 0 < elo < ehi and elo >= prev) or not 0 <= groups <= 128:
            raise AuditError("RSP2 true-energy row invalid or out of order")
        prev = ehi
        starts = vector(row, 10, ">h", 2)
        sizes = vector(row, 18, ">h", 2)
        values = vector(row, 26, ">f", 4)
        if len(starts) != groups or len(sizes) != groups or len(values) != sum(sizes):
            raise AuditError("RSP2 DRM channel-group heap mismatch")
        counts = [0.] * 128
        touched = set()
        cursor = 0
        for first, num in zip(starts, sizes):
            # GBM OGIP RSP2 F_CHAN is 1-based (1..128) although GBM
            # TTE PHA/EBOUNDS CHANNEL is 0-based (0..127).
            if first < 1 or num < 1 or first + num - 1 > 128:
                raise AuditError("RSP2 DRM group outside detector channel range")
            for i in range(num):
                channel = first + i - 1
                value = float(values[cursor + i])
                if channel in touched or not math.isfinite(value) or value < 0:
                    raise AuditError("RSP2 overlapping/negative or nonfinite matrix entry")
                touched.add(channel)
                counts[channel] = value
            cursor += num
        energies.append((float(elo), float(ehi)))
        matrix.append(counts)
    return {"true_energy_edges_keV": energies, "response_cm2": matrix,
            "time_span_met": (float(h["TSTART"]), float(h["TSTOP"])),
            "true_energy_bins": len(energies), "channels": 128}


def weighted_drm(hdus, start_met: float, stop_met: float) -> tuple[list[tuple[float, float]], list[list[float]], dict]:
    if start_met >= stop_met:
        raise AuditError("invalid GBM RSP2 window")
    mats = [decode_drm(h) for h in hdus if h.header.get("EXTNAME") == "SPECRESP MATRIX"]
    if not mats:
        raise AuditError("no response matrices in GBM RSP2")
    energy = mats[0]["true_energy_edges_keV"]
    response = [[0.] * 128 for _ in energy]
    duration = stop_met - start_met
    used = []
    coverage = 0.
    for i, item in enumerate(mats):
        if item["true_energy_edges_keV"] != energy:
            raise AuditError("RSP2 true-energy grids change between attitude DRMs")
        left, right = item["time_span_met"]
        overlap = max(0., min(right, stop_met) - max(left, start_met))
        if overlap <= 0:
            continue
        coverage += overlap
        weight = overlap / duration
        used.append({"matrix_number": i + 1, "weight": round(weight, 8)})
        for ei, source in enumerate(item["response_cm2"]):
            target = response[ei]
            for c, v in enumerate(source):
                target[c] += weight * v
    if abs(coverage - duration) > max(1e-5, duration * 1e-5):
        raise AuditError("GBM prompt segment not entirely covered by physical RSP2 matrices")
    return energy, response, {"rsp2_matrices_used": used, "rsp2_total_weight": round(coverage / duration, 8)}


def channel_counts(times: list[float], channels: list[int], lo: float, hi: float) -> list[int]:
    if lo >= hi or len(times) != len(channels):
        raise AuditError("invalid TTE segment")
    counts = [0] * 128
    for time, channel in zip(times, channels):
        if lo <= time < hi:
            counts[channel] += 1
    return counts


def photon_bin_flux(energies: list[tuple[float, float]], alpha: float, epeak: float) -> list[float]:
    if alpha <= -2 or epeak <= 0:
        raise AuditError("invalid cutoff powerlaw shape")
    return [((math.sqrt(lo * hi) / 100.) ** alpha) *
            math.exp(-(alpha + 2.) * math.sqrt(lo * hi) / epeak) * (hi - lo)
            for lo, hi in energies]


def folded_rate(energy: list[tuple[float, float]], response: list[list[float]], alpha: float, epeak: float) -> list[float]:
    if len(energy) != len(response):
        raise AuditError("DRM true energy axis mismatch")
    photons = photon_bin_flux(energy, alpha, epeak)
    out = [0.] * 128
    for profile, row in zip(photons, response):
        for c, effective_area in enumerate(row):
            out[c] += profile * effective_area
    return out


def profile_background(on: int, off: int, predicted_source: float, tau: float) -> float:
    """Analytic profile of on-region background b in Poisson(on,s+b)*Poisson(off,tau*b)."""
    if on < 0 or off < 0 or predicted_source < 0 or tau <= 0 or not math.isfinite(tau):
        raise AuditError("invalid Poisson on-off nuisance input")
    k = 1. + tau
    linear = k * predicted_source - on - off
    return max(0., (-linear + math.sqrt(linear * linear + 4. * k * off * predicted_source)) / (2. * k))


def poisson_deviance(observed: int, expected: float) -> float:
    if expected < 0 or observed < 0 or not math.isfinite(expected):
        raise AuditError("invalid Poisson expected count")
    if observed == 0:
        return 2 * expected
    if expected <= 0:
        return math.inf
    return 2 * (expected - observed + observed * math.log(observed / expected))


def cstat_onoff(on: list[int], off: list[int], shape: list[float], amplitude: float, duration: float, bg_duration: float, allowed: list[bool]) -> float:
    if not (len(on) == len(off) == len(shape) == len(allowed) == 128):
        raise AuditError("on-off channel vector mismatch")
    tau = bg_duration / duration
    total = 0.
    for n, m, rate, use in zip(on, off, shape, allowed):
        if not use:
            continue
        s = max(0., amplitude * rate * duration)
        b = profile_background(n, m, s, tau)
        total += poisson_deviance(n, s + b) + poisson_deviance(m, tau * b)
    return total


def profiled_shape(detectors: list[dict], alpha: float, epeak: float) -> dict:
    rows = []
    for d in detectors:
        rate = folded_rate(d["energy"], d["response"], alpha, epeak)
        rows.append((d, rate))

    def objective(log_amplitude: float) -> float:
        amp = math.exp(log_amplitude)
        return sum(cstat_onoff(d["on"], d["off"], rate, amp, d["duration"], d["bg_duration"], d["allowed"]) for d, rate in rows)

    # One-dimensional bounded golden search; values outside the broad
    # normalization range deliberately trigger a boundary warning.
    left, right = math.log(1e-9), math.log(1e3)
    ratio = (math.sqrt(5.) - 1.) / 2.
    c, e = right - ratio * (right - left), left + ratio * (right - left)
    fc, fe = objective(c), objective(e)
    for _ in range(48):
        if fc <= fe:
            right, e, fe = e, c, fc
            c = right - ratio * (right - left)
            fc = objective(c)
        else:
            left, c, fc = c, e, fe
            e = left + ratio * (right - left)
            fe = objective(e)
    middle = (left + right) / 2.
    return {"alpha": alpha, "epeak_keV": epeak,
            "amplitude_ph_cm2_s_keV_at_100keV": math.exp(middle),
            "onoff_poisson_deviance": objective(middle),
            "amplitude_boundary_warning": abs(middle - math.log(1e-9)) < .03 or abs(middle - math.log(1e3)) < .03}


def fit_segment(detectors: list[dict]) -> dict:
    if len(detectors) != 3:
        raise AuditError("expected two GBM NaI and one BGO detector")
    models = [profiled_shape(detectors, alpha, ep) for alpha in ALPHAS for ep in EPEAK_KEV]
    best = min(models, key=lambda model: model["onoff_poisson_deviance"])
    null = sum(cstat_onoff(d["on"], d["off"], [0.] * 128, 0., d["duration"], d["bg_duration"], d["allowed"]) for d in detectors)
    fitted_channels = sum(sum(d["allowed"]) for d in detectors)
    stress_ratio = best["onoff_poisson_deviance"] / fitted_channels
    return {**best,
            "null_no_photon_source_deviance": round(null, 6),
            "fit_deviance": round(best["onoff_poisson_deviance"], 6),
            "source_vs_no_source_diagnostic_delta_deviance": round(max(0., null - best["onoff_poisson_deviance"]), 6),
            "unscaled_deviance_per_fitted_channel_diagnostic": round(stress_ratio, 6),
            "spectral_mismatch_stress_gate": "HIGH_RESIDUAL_MODEL_REVISE" if stress_ratio > 5 else "INCONCLUSIVE_WITHOUT_CALIBRATED_PREDICTIVE_CHECKS",
            "stress_gate_threshold_descriptive_not_p_value": 5.,
            "not_a_discovery_significance_or_calibrated_p_value": True,
            "profiled_background_per_channel": True,
            "fixed_template_search_alpha": ALPHAS,
            "fixed_template_search_epeak_keV": EPEAK_KEV,
            "photon_model": "SINGLE_CUTOFF_POWERLAW_MIDPOINT_BIN_QUADRATURE_NOT_BAND_PLUS_POWERLAW",
            "no_interpolation_between_template_grid_nodes": True}


def run(manifest: dict, gbm_dir: Path) -> dict:
    if manifest["schema"] != "janus.cosmos.fermi-irf-gbm-manifest.v5" or [b["grb"] for b in manifest["gbm_bursts"]] != list(WINDOWS):
        raise AuditError("V6 GBM data cohort differs from frozen V5")
    bursts = []
    for burst in manifest["gbm_bursts"]:
        name = burst["grb"]
        prepared = []
        for detector in burst["detectors"]:
            label = detector["detector"]
            files = detector["event_and_response_files"]
            tte = fits_large_tte(checked_file(gbm_dir, files["tte"]))
            rsp = read_fits(checked_file(gbm_dir, files["rsp2"]))
            checked_file(gbm_dir, files["cspec"])
            times, channels = _tte_events(get_hdu(tte, "EVENTS"))
            emin, emax = _valid_ebounds(get_hdu(tte, "EBOUNDS"))
            low, high = (8., 900.) if label.startswith("n") else (300., 30000.)
            allowed = [a < high and b > low for a, b in zip(emin, emax)]
            if not any(allowed):
                raise AuditError("no calibrated GBM channels")
            back = [0] * 128
            for lo, hi in BG_WINDOWS:
                segment = channel_counts(times, channels, lo, hi)
                back = [a + b for a, b in zip(back, segment)]
            bg_t = sum(hi - lo for lo, hi in BG_WINDOWS)
            if bg_t < MIN_BACK_EXPOSURE_S:
                raise AuditError("insufficient GBM off-source exposure")
            trig = float(tte[0].header["TRIGTIME"])
            if abs(trig - float(rsp[0].header["TRIGTIME"])) > 1e-5:
                raise AuditError("GBM TTE vs RSP2 trigger mismatch")
            prepared.append({"name": label, "trigger": trig, "times": times, "channels": channels,
                             "off": back, "bg_duration": bg_t, "allowed": allowed, "rsp": rsp})
        time_rows = []
        for lo, hi in WINDOWS[name]:
            detector_inputs = []
            metadata = []
            for d in prepared:
                energy, response, drm_meta = weighted_drm(d["rsp"], d["trigger"] + lo, d["trigger"] + hi)
                on = channel_counts(d["times"], d["channels"], lo, hi)
                detector_inputs.append({"name": d["name"], "energy": energy, "response": response,
                                        "on": on, "off": d["off"], "duration": hi - lo,
                                        "bg_duration": d["bg_duration"], "allowed": d["allowed"]})
                metadata.append({"detector": d["name"], "on_selected_counts": sum(v for v, allow in zip(on, d["allowed"]) if allow),
                                 "off_selected_counts": sum(v for v, allow in zip(d["off"], d["allowed"]) if allow),
                                 "fitted_channel_count": sum(d["allowed"]), **drm_meta})
            spectral = fit_segment(detector_inputs)
            time_rows.append({"start_relative_trigger_s": lo, "stop_relative_trigger_s": hi,
                              "detector_counts_and_time_weighted_rsp2": metadata,
                              "onoff_spectral_grid_fit": spectral})
        bursts.append({"grb": name, "gbm_trigger_id": burst["gbm_trigger_id"],
                       "exact_trigger_met": prepared[0]["trigger"],
                       "time_resolved_spectral_fit_segments": time_rows,
                       "spectral_goodness_of_fit_calibrated": False,
                       "background_model_time_variable": False,
                       "response_ebounds_reconciliation_calibrated": False})
    result = {"schema": SCHEMA, "experiment_id": "JANUS-COSMOS-ALI-V6-LAT-GBM-JOINT-2026-10-19",
              "data_kind": "ACTUAL_GBM_COUNTS_OGIP_RSP2_FORWARD_FOLDED_PILOT_NOT_VALIDATED_GLOBAL_MODEL",
              "source_manifest": str(MANIFEST_PATH.relative_to(ROOT)),
              "source_manifest_sha256": __import__("hashlib").sha256(MANIFEST_PATH.read_bytes()).hexdigest(),
              "response_units": "cm2", "photon_model_units": "ph_cm_minus2_s_minus1_keV_minus1",
              "background_treatment": "ONOFF_POISSON_PROFILE_CHANNEL_NUISANCE_WITH_FIXED_-20_TO_-5_s_OFFWINDOW",
              "time_variable_background_inference": "NOT_IMPLEMENTED",
              "pulse_spectral_evolution": "MULTIPLE_FIXED_TIME_WINDOWS_INDEPENDENT_CUTOFF_POWERLAW_GRID",
              "spectral_template_limitation": "SINGLE_COMPONENT_CUTOFF_POWERLAW_IGNORES_POSSIBLE_ADDITIONAL_GE_V_OR_HARD_COMPONENT",
              "time_weighted_RSP2": "ATTITUDE_DEPENDENT_DRMS_EXPOSURE_WEIGHTED_OVER_SEGMENT",
              "GBM_detector_coupling": "TWO_NAI_ONE_BGO_SHARED_ALPHA_EPEAK_AND_NORMALIZATION_PER_TIME_SEGMENT",
              "gbm_bursts": bursts,
              "scientific_gates": {
                  "REAL_GBM_RSP2_OGIP_HEAP_MATRICES": "PASS_BYTE_PINNED_AND_DECODED",
                  "REAL_GBM_ONOFF_SPECTRAL_COUNT_FORWARD_FOLD": "PILOT_GRID_WITH_SIMPLE_SPECTRAL_LAW",
                  "GBM_TIME_VARYING_BACKGROUND_AND_RESPONSE_CHANNEL_SYSTEMATICS": "UNVALIDATED",
                  "GBM_MODEL_COMPLETENESS": "NOT_ESTABLISHED_SINGLE_CPL_NOT_BAND_PLUS_PL",
                  "LAT_FULL_PSF_EDISP_EXPOSURE_BACKGROUND": "NOT_YET_VALIDATED",
                  "INDEPENDENT_GE_V_SOURCE_LAG_PRIOR": "NOT_ESTABLISHED",
                  "OBSERVED_LAT_GBM_LIV_HYPOTHESIS_TEST": "BLOCKED",
                  "LORENTZ_VIOLATION": "NOT_ESTABLISHED",
                  "AMPS_FIREWALL": "OPEN_UNCONFIRMED",
                  "BHB_BOOM": "NOT_ESTABLISHED"},
              "observational_liv_eta_estimate": None,
              "observational_joint_p_value": None,
              "black_hole_firewall_detection": False,
              "claim_ceiling": "REAL_GBM_RSP2_ONOFF_FORWARD_FOLDED_SPECTRAL_DIAGNOSTIC_ONLY_NOT_PHYSICAL_JOINT_LAT_GBM_LIV_TEST"}
    result["result_sha256"] = digest(result)
    return result


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gbm-dir", type=Path, required=True)
    p.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    p.add_argument("--output", type=Path, default=OUT)
    a = p.parse_args()
    value = run(json.loads(a.manifest.read_text()), a.gbm_dir)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n")
    print(json.dumps({"output": str(a.output), "burst_count": len(value["gbm_bursts"]),
                      "spectral_segments": sum(len(x["time_resolved_spectral_fit_segments"]) for x in value["gbm_bursts"]),
                      "observational_p_value": value["observational_joint_p_value"]}))


if __name__ == "__main__":
    main()
