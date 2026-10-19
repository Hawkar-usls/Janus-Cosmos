#!/usr/bin/env python3
"""V6-R1 GBM spectral *stress test*: Band + power law and off-time drift.

The actual GBM TTE and RSP2 bytes are inherited from the frozen V5 manifest.
This is a reproducible count-domain model-comparison *pilot*, not a validated
physical LAT/GBM likelihood, a LAT energy measurement, or a LIV detection.

SciPy/NumPy are optional runtime dependencies; the pure functions and unit
tests have no dependency beyond the Python standard library.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from janus_cosmos_fermi_calibration_v5 import (
    MANIFEST_PATH, _tte_events, _valid_ebounds, checked_file,
    fits_large_tte, get_hdu, read_fits,
)
from janus_cosmos_fermi_response_v6 import WINDOWS, channel_counts, weighted_drm
from janus_cosmos_simulation_audit import AuditError, ROOT, digest

SCHEMA = "janus.cosmos.fermi-gbm-band-plus-powerlaw-drift-independent-holdouts.v6r1"
OUT = ROOT / "data/cosmos/COSMOS-ALI-V6R1-GBM-BAND-POWERLAW-DRIFT-2026-10-19.json"

# Fixed before V6-R1 fits. The latest off interval is held out from drift
# fitting, and never reused to estimate the on-window background rate.
BACKGROUND = {
    "090510": ((-23., -17.), (-17., -11.), (-11., -5.)),
    "090902B": ((-26., -19.), (-19., -12.), (-12., -5.)),
    "130427A": ((-130., -90.), (-90., -50.), (-50., -15.)),
}
ALPHAS = (-1.4, -1.0, -0.6)
BETAS = (-3.0, -2.4)
EPEAK_KEV = (150., 400., 1000., 2500.)
EXTRA_PL_INDEX = -1.8
EXTRA_PL_RELATIVE_AT_1MEV = (0., 0.5, 2., 8.)


def band_shape(energy_keV: float, alpha: float, beta: float, epeak_keV: float) -> float:
    """Standard Band photon spectrum, amplitude 1 at nominal 100-keV pivot."""
    if not (energy_keV > 0 and epeak_keV > 0 and alpha > -2 and beta < alpha):
        raise AuditError("invalid Band shape")
    e0 = epeak_keV / (alpha + 2.)
    cross = (alpha - beta) * e0
    if energy_keV <= cross:
        return (energy_keV / 100.) ** alpha * math.exp(-energy_keV / e0)
    return ((cross / 100.) ** (alpha - beta) * math.exp(beta - alpha)
            * (energy_keV / 100.) ** beta)


def band_plus_pl_shape(energy_keV: float, alpha: float, beta: float,
                       epeak_keV: float, extra_fraction: float,
                       extra_index: float = EXTRA_PL_INDEX) -> float:
    if not (0 <= extra_fraction <= 100 and -4 < extra_index < -1):
        raise AuditError("invalid extra power-law component")
    band = band_shape(energy_keV, alpha, beta, epeak_keV)
    powerlaw = extra_fraction * band_shape(1000., alpha, beta, epeak_keV)
    return band + powerlaw * (energy_keV / 1000.) ** extra_index


def integrate_photon_bins(edges: list[tuple[float, float]], shape) -> list[float]:
    """Four-subinterval Simpson quadrature of dN/dE in each true-energy bin."""
    result = []
    for lo, hi in edges:
        if not (0 < lo < hi):
            raise AuditError("invalid true-energy bin")
        step = (hi - lo) / 4.
        values = [shape(lo + step * j) for j in range(5)]
        result.append(step / 3. * (values[0] + 4 * values[1] + 2 * values[2] + 4 * values[3] + values[4]))
    return result


def exponential_time_kernel(start: float, stop: float, origin: float, slope: float) -> float:
    """Stable integral of exp(slope * (t - origin)) dt, in seconds."""
    if not (math.isfinite(start) and math.isfinite(stop) and math.isfinite(slope) and start < stop and abs(slope) <= .04):
        raise AuditError("background time kernel outside bounds")
    width = stop - start
    midpoint = (start + stop) / 2.
    half = slope * width / 2.
    if abs(half) < 1e-5:
        sinhc = 1. + half ** 2 / 6. + half ** 4 / 120.
    else:
        sinhc = math.sinh(half) / half
    return width * math.exp(slope * (midpoint - origin)) * sinhc


def fitted_background_slope(train_counts: tuple[int, int],
                            train_windows: tuple[tuple[float, float], tuple[float, float]]) -> float:
    """Shared per-detector log-linear drift from disjoint off-source windows.

    Shape is profiled from total counts; per-channel Poisson nuisance counts
    are *not* replaced by this point estimate during spectral fitting.
    """
    if len(train_counts) != 2 or min(train_counts) < 0 or not sum(train_counts):
        return 0.
    import scipy.optimize  # optional dependency needed only for actual fitting
    a, b = train_counts
    origin = (train_windows[0][0] + train_windows[-1][1]) / 2.

    def negative_log_likelihood(slope: float) -> float:
        k0, k1 = (exponential_time_kernel(*win, origin, slope) for win in train_windows)
        norm = k0 + k1
        return -(a * math.log(k0 / norm) + b * math.log(k1 / norm))

    fitted = scipy.optimize.minimize_scalar(negative_log_likelihood, bounds=(-.039, .039), method="bounded",
                                           options={"xatol": 1e-6})
    if not fitted.success:
        raise AuditError("background drift optimization failed")
    return float(fitted.x)


def onoff_deviance_variable(on, off, expected_source, on_kernel: float, off_kernel: float) -> float:
    """Vectorized profiled ON/OFF Poisson deviance; same off count reused by windows.

    This *per-window* diagnostic is not a globally factorized likelihood:
    the shared off sample induces dependence between adjacent prompt windows.
    """
    import numpy as np
    on, off, source = (np.asarray(v, dtype=float) for v in (on, off, expected_source))
    if on.shape != off.shape or on.shape != source.shape or on_kernel <= 0 or off_kernel <= 0:
        raise AuditError("invalid variable-background on/off shapes or weights")
    if not (np.all(np.isfinite(on)) and np.all(np.isfinite(off)) and np.all(np.isfinite(source))):
        raise AuditError("nonfinite on/off input")
    if np.any(on < 0) or np.any(off < 0) or np.any(source < 0):
        raise AuditError("negative on/off value")
    tau = off_kernel / on_kernel
    k = 1. + tau
    linear = k * source - on - off
    back = np.maximum(0., (-linear + np.sqrt(linear * linear + 4. * k * off * source)) / (2. * k))
    mu_on = source + back
    mu_off = tau * back

    def dev(n, mu):
        with np.errstate(divide="ignore", invalid="ignore"):
            term = np.where(n > 0, n * np.log(np.where(mu > 0, n / mu, 1.)), 0.)
        if np.any((n > 0) & (mu <= 0)):
            return math.inf
        return float(2. * np.sum(mu - n + term))

    return dev(on, mu_on) + dev(off, mu_off)


def fit_templates(detectors: list[dict]) -> dict:
    """Fit a frozen, finite Band(+PL) grid and one common amplitude.

    Real GBM channels and actual time-overlap RSP2 matrices enter the source
    prediction. A nominal 1-MeV power-law ratio is shape-only; its significance
    is *not* inferred from choosing the best template post hoc.
    """
    import numpy as np
    from scipy.optimize import minimize_scalar

    for d in detectors:
        if len(d["energy"]) != len(d["response"]) or len(d["allowed"]) != 128:
            raise AuditError("GBM detector response/channel mismatch")
    arrays = [{**d, "matrix": np.asarray(d["response"], dtype=float),
               "on_arr": np.asarray(d["on"], dtype=float)[d["allowed"]],
               "off_arr": np.asarray(d["off"], dtype=float)[d["allowed"]]}
              for d in detectors]
    models = []
    for alpha in ALPHAS:
        for beta in BETAS:
            for epeak in EPEAK_KEV:
                for q in EXTRA_PL_RELATIVE_AT_1MEV:
                    folded = []
                    for d in arrays:
                        photons = np.asarray(integrate_photon_bins(
                            d["energy"], lambda E: band_plus_pl_shape(E, alpha, beta, epeak, q)), dtype=float)
                        folded.append((photons @ d["matrix"])[d["allowed"]] * d["duration"])

                    def objective(logamp):
                        amp = math.exp(logamp)
                        return sum(onoff_deviance_variable(d["on_arr"], d["off_arr"], amp * signal,
                                                            d["on_bg_kernel"], d["off_bg_kernel"])
                                   for d, signal in zip(arrays, folded))

                    candidate = minimize_scalar(objective, bounds=(-23., 6.9), method="bounded",
                                                options={"xatol": 2e-3, "maxiter": 100})
                    if not candidate.success or not math.isfinite(candidate.fun):
                        raise AuditError("Band+PL amplitude fit failed")
                    models.append({"alpha": alpha, "beta": beta, "epeak_keV": epeak,
                                   "extra_index": EXTRA_PL_INDEX, "extra_relative_flux_at_1MeV": q,
                                   "amplitude_at_100_keV": round(math.exp(candidate.x), 8),
                                   "onoff_deviance": round(float(candidate.fun), 6),
                                   "amplitude_fit_boundary_warning": bool(candidate.x < -22.8 or candidate.x > 6.7)})
    band = min((r for r in models if r["extra_relative_flux_at_1MeV"] == 0.), key=lambda x: x["onoff_deviance"])
    combined = min(models, key=lambda x: x["onoff_deviance"])
    return {"band_only": band, "band_plus_powerlaw": combined,
            "descriptive_deviance_improvement_not_lrt": round(max(0., band["onoff_deviance"] - combined["onoff_deviance"]), 6),
            "num_predeclared_shape_templates": len(models),
            "model_selection_significance": None,
            "grid_is_coarse_and_fit_parameter_coverage_unvalidated": True}


def _prepare_detector(item: dict, grb: str, gbm_dir: Path) -> dict:
    files = item["event_and_response_files"]
    tte = fits_large_tte(checked_file(gbm_dir, files["tte"]))
    rsp = read_fits(checked_file(gbm_dir, files["rsp2"]))
    checked_file(gbm_dir, files["cspec"])
    times, channels = _tte_events(get_hdu(tte, "EVENTS"))
    low, high = ((8., 900.) if item["detector"].startswith("n") else (300., 30000.))
    elo, ehi = _valid_ebounds(get_hdu(tte, "EBOUNDS"))
    allowed = [a < high and b > low for a, b in zip(elo, ehi)]
    if not allowed or not any(allowed):
        raise AuditError("GBM detector has no accepted channel")
    train1, train2, heldout = BACKGROUND[grb]
    train_counts = [channel_counts(times, channels, *w) for w in (train1, train2)]
    hold_counts = channel_counts(times, channels, *heldout)
    selected_totals = tuple(sum(v for v, ok in zip(counts, allowed) if ok) for counts in train_counts)
    slope = fitted_background_slope(selected_totals, (train1, train2))
    origin = (train1[0] + train2[1]) / 2.
    off_kernel = sum(exponential_time_kernel(*w, origin, slope) for w in (train1, train2))
    hold_kernel = exponential_time_kernel(*heldout, origin, slope)
    off_by_channel = [a + b for a, b in zip(*train_counts)]
    total_off = sum(v for v, ok in zip(off_by_channel, allowed) if ok)
    predicted_holdout = total_off * hold_kernel / off_kernel
    observed_holdout = sum(v for v, ok in zip(hold_counts, allowed) if ok)
    return {"name": item["detector"], "trigger": float(tte[0].header["TRIGTIME"]), "times": times,
            "channels": channels, "rsp": rsp, "allowed": allowed, "off": off_by_channel,
            "train_counts": selected_totals, "train_intervals": (train1, train2),
            "heldout_interval": heldout, "holdout_observed": observed_holdout,
            "holdout_expected": predicted_holdout, "background_slope_per_second": slope,
            "kernel_origin": origin, "off_bg_kernel": off_kernel}


def run(manifest: dict, gbm_dir: Path) -> dict:
    if manifest["schema"] != "janus.cosmos.fermi-irf-gbm-manifest.v5" or set(BACKGROUND) != {x["grb"] for x in manifest["gbm_bursts"]}:
        raise AuditError("not the frozen V5 GBM sample")
    bursts = []
    for burst in manifest["gbm_bursts"]:
        grb = burst["grb"]
        prepared = [_prepare_detector(det, grb, gbm_dir) for det in burst["detectors"]]
        if len(prepared) != 3 or not (prepared[0]["name"].startswith("n") and prepared[1]["name"].startswith("n") and prepared[2]["name"].startswith("b")):
            raise AuditError("independent BGO holdout requires two NaI and one BGO")
        segment_results = []
        for lo, hi in WINDOWS[grb]:
            per_detector = []
            for d in prepared:
                energy, response, _ = weighted_drm(d["rsp"], d["trigger"] + lo, d["trigger"] + hi)
                per_detector.append({"energy": energy, "response": response,
                                     "on": channel_counts(d["times"], d["channels"], lo, hi),
                                     "off": d["off"], "duration": hi - lo,
                                     "on_bg_kernel": exponential_time_kernel(lo, hi, d["kernel_origin"], d["background_slope_per_second"]),
                                     "off_bg_kernel": d["off_bg_kernel"], "allowed": d["allowed"]})
            all_fit = fit_templates(per_detector)
            nai_fit = fit_templates(per_detector[:2])
            # Detector holdout: *never* use BGO on-source counts to select the
            # template/amplitude. OFF background for BGO remains its own data.
            bgo = per_detector[2]
            import numpy as np
            def heldout_bgo_deviance(model):
                photon = np.asarray(integrate_photon_bins(
                    bgo["energy"], lambda E: band_plus_pl_shape(E, model["alpha"], model["beta"],
                                                                 model["epeak_keV"], model["extra_relative_flux_at_1MeV"])))
                rate = photon @ np.asarray(bgo["response"])
                a = np.asarray(bgo["allowed"], dtype=bool)
                return round(onoff_deviance_variable(np.asarray(bgo["on"])[a], np.asarray(bgo["off"])[a],
                                                      np.asarray(rate)[a] * model["amplitude_at_100_keV"] * bgo["duration"],
                                                      bgo["on_bg_kernel"], bgo["off_bg_kernel"]), 6)
            bgo_band = heldout_bgo_deviance(nai_fit["band_only"])
            bgo_combo = heldout_bgo_deviance(nai_fit["band_plus_powerlaw"])
            segment_results.append({"relative_trigger_interval_s": [lo, hi], "three_detector_fit": all_fit,
                                    "nai_train_bgo_holdout": {
                                        "nai_only_best_band": nai_fit["band_only"],
                                        "nai_only_best_band_plus_powerlaw": nai_fit["band_plus_powerlaw"],
                                        "bgo_holdout_band_deviance": bgo_band,
                                        "bgo_holdout_band_plus_powerlaw_deviance": bgo_combo,
                                        "bgo_predeclared_independent_detector_evaluation": True,
                                        "bgo_absolute_fit_p_value": None},
                                    "background_extrapolation_over_30s": hi - prepared[0]["train_intervals"][-1][1] > 30.})
        bursts.append({"grb": grb, "gbm_trigger_id": burst["gbm_trigger_id"],
                       "detector_pretrigger_train_and_holdout_background": [
                           {"detector": d["name"], "train_intervals": d["train_intervals"],
                            "heldout_interval": d["heldout_interval"], "train_counts": d["train_counts"],
                            "log_linear_drift_per_second": round(d["background_slope_per_second"], 8),
                            "heldout_observed_counts": d["holdout_observed"],
                            "heldout_predicted_counts_plugin_not_predictive_interval": round(d["holdout_expected"], 4)}
                           for d in prepared], "spectral_segments": segment_results})
    out = {"schema": SCHEMA, "input_manifest": str(MANIFEST_PATH.relative_to(ROOT)),
           "input_v5_manifest_sha256": __import__("hashlib").sha256(MANIFEST_PATH.read_bytes()).hexdigest(),
           "data_kind": "REAL_TTE_RSP2_ONOFF_SPECTRAL_DIAGNOSTIC_NOT_CALIBRATED_JOINT_LIV_LIKELIHOOD",
           "photons_per_energy_units": "shape_x_amplitude_ph_cm_minus2_s_minus1_keV_minus1",
           "spectral_photon_model": "BAND_PLUS_POWERLAW_AT_1MEV_RATIO_FIXED_GRID_SIMPSON_TRUE_ENERGY_QUADRATURE",
           "background_model": "PER_DETECTOR_SHARED_LOG_LINEAR_DRIFT_PRETRIGGER_OFF_TRAIN_TWO_WINDOWS_FINAL_PRETRIGGER_WINDOW_HELD_OUT",
           "background_and_fit_important_limitations": [
               "two off intervals give weak slope constraints; fit ignores slope uncertainty",
               "prompt windows share the same OFF counts and are statistically dependent",
               "late prompt intervals extrapolate beyond short pretrigger background control",
               "GBM TTE/CSPEC/RSP2 channel-edge systematics uncalibrated",
               "no LAT response or intrinsic GeV emission-lag prior is in this file",
               "predeclared coarse grid fit lacks calibrated search and predictive coverage"],
           "gbm_bursts": bursts,
           "observed_liv_p_value": None,
           "real_grb_instrument_aware_liv_lrt": None,
           "scientific_gates": {"REAL_GBM_TTE_RSP2_FOLD": "PASS_SHA256_PINNED_BAND_PLUS_PL_PILOT",
                                "GBM_TIME_VARIABLE_BACKGROUND": "LOG_LINEAR_FIT_AND_INDEPENDENT_PRETRIGGER_HOLDOUT_PILOT",
                                "GBM_INDEPENDENT_DETECTOR_VALIDATION": "BGO_HOLDOUT_DIAGNOSTIC_NOT_CALIBRATED_P",
                                "LAT_FULL_PSF_EDISP_DIFFUSE": "NOT_IMPLEMENTED_HERE",
                                "SOURCE_INTRINSIC_LAG_PRIOR": "NOT_ESTABLISHED",
                                "OBSERVED_MULTIBURST_LIV_P_VALUE": "BLOCKED",
                                "BHB_FIREWALL": "OPEN_UNCONFIRMED",
                                "BHB_BOOM": "NOT_ESTABLISHED"},
           "claim_ceiling": "REAL_GBM_RSP2_BAND_PLUS_PL_DRIFT_STRESS_TEST_ONLY_NO_PHYSICAL_LIV_DETECTION"}
    out["result_sha256"] = digest(out)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gbm-dir", type=Path, required=True)
    p.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    p.add_argument("--output", type=Path, default=OUT)
    args = p.parse_args()
    result = run(json.loads(args.manifest.read_text()), args.gbm_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output), "grbs": len(result["gbm_bursts"]),
                      "segments": sum(len(x["spectral_segments"]) for x in result["gbm_bursts"]),
                      "observed_p": result["observed_liv_p_value"]}))


if __name__ == "__main__":
    main()
