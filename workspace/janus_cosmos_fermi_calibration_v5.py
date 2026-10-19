#!/usr/bin/env python3
"""JANUS Cosmos V5: pinned NASA GBM timing and LAT Pass8 IRF table audit.

This is instrument *data preparation*, not a validated LAT/GBM spectral or
time-of-flight likelihood. In particular, nearest-bin Aeff is illustrative;
neither PSF/EDISP kernels nor GBM RSP2 matrices are forward folded.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import struct
from pathlib import Path

from janus_cosmos_fermi_raw_v4 import (
    BLOCK, HDU, _card_value, _column_info, read_fits, rows, sha256,
)
from janus_cosmos_simulation_audit import AuditError, ROOT, digest

SCHEMA = "janus.cosmos.fermi-gbm-irf-audit.v5"
MANIFEST_PATH = ROOT / "data/cosmos/COSMOS-ALI-V5-FERMI-LAT-GBM-CALDB-MANIFEST-2026-10-19.json"
V4_PATH = ROOT / "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-EVENT-DIAGNOSTICS-2026-10-18.json"
V4_MANIFEST_PATH = ROOT / "data/cosmos/COSMOS-ALI-V4-FERMI-RAW-FITS-MANIFEST-2026-10-18.json"
EXPECTED_BURSTS = ("090510", "090902B", "130427A")
EXPECTED_IRF = "P8R3_SOURCE_V3"
KNOWN_PARTITIONS = ("FB", "PSF", "EDISP")


def checked_file(root: Path, item: dict) -> Path:
    """Strictly reject absolute/traversing names and unpinned bytes."""
    name = item["name"]
    if not isinstance(name, str) or not name or Path(name).is_absolute() or ".." in Path(name).parts:
        raise AuditError("unsafe or absent manifest name")
    path = root / name
    if not path.is_file() or path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
        raise AuditError(f"V5_INPUT_BYTES_OR_SHA256_MISMATCH:{name}")
    return path


def fits_large_tte(path: Path) -> list[HDU]:
    """Bounded FITS subset supporting >1M GBM TTE events (not general FITS)."""
    if not 0 < path.stat().st_size <= 100_000_000 or path.stat().st_size % BLOCK:
        raise AuditError("invalid GBM TTE byte extent")
    blob = path.read_bytes()
    if not blob.startswith(b"SIMPLE  ="):
        raise AuditError("GBM FITS primary header absent")
    out = []
    position = 0
    while position < len(blob):
        beginning = position
        header = {}
        for i in range(10000):
            offset = beginning + 80 * i
            if offset + 80 > len(blob):
                raise AuditError("truncated GBM FITS header")
            card = blob[offset:offset + 80]
            key = card[:8].decode("ascii", "strict").strip()
            if key == "END":
                position = beginning + ((i + 1) * 80 + BLOCK - 1) // BLOCK * BLOCK
                break
            if card[8:10] == b"= ":
                header[key] = _card_value(card)
        else:
            raise AuditError("no GBM FITS END card")
        axis = int(header.get("NAXIS", -1))
        if header.get("XTENSION") == "BINTABLE":
            nrows = int(header.get("NAXIS2", -1))
            stride = int(header.get("NAXIS1", -1))
            if not (0 <= nrows <= 5_000_000 and 0 <= stride <= 65536):
                raise AuditError("unsafe GBM FITS table dimensions")
            extent = nrows * stride + int(header.get("PCOUNT", 0))
        elif axis == 0:
            extent = 0
        else:
            raise AuditError("unexpected GBM TTE image extension")
        if extent < 0 or position + extent > len(blob):
            raise AuditError("truncated GBM FITS table")
        out.append(HDU(header, blob[position:position + extent]))
        position += ((extent + BLOCK - 1) // BLOCK) * BLOCK
    if position != len(blob) or len(out) != 4:
        raise AuditError("unexpected GBM TTE FITS structure")
    return out


def get_hdu(hdus: list[HDU], extension: str) -> HDU:
    found = [h for h in hdus if h.header.get("EXTNAME") == extension]
    if len(found) != 1:
        raise AuditError(f"required unique extension missing: {extension}")
    return found[0]


def _grid_and_finiteness(hdu: HDU, grid_name: str, cols: set[str], nenergy: int, ntheta: int) -> dict:
    table = rows(hdu, cols | {"ENERG_LO", "ENERG_HI", "CTHETA_LO", "CTHETA_HI"})
    if len(table) != 1:
        raise AuditError(f"unexpected {grid_name} rows")
    row = table[0]
    elo, ehi = row["ENERG_LO"], row["ENERG_HI"]
    clo, chi = row["CTHETA_LO"], row["CTHETA_HI"]
    if len(elo) != nenergy or len(ehi) != nenergy or len(clo) != ntheta or len(chi) != ntheta:
        raise AuditError(f"invalid {grid_name} axis shape")
    for lo, hi in ((elo, ehi), (clo, chi)):
        if any(not (math.isfinite(a) and math.isfinite(b) and a < b) for a, b in zip(lo, hi)):
            raise AuditError(f"invalid {grid_name} bin range")
    for col in cols:
        if len(row[col]) != nenergy * ntheta or not all(math.isfinite(float(x)) for x in row[col]):
            raise AuditError(f"invalid {grid_name} grid {col}")
    return row


def audit_caldb(manifest: dict, root: Path) -> tuple[dict, dict]:
    cal = manifest["lat_caldb"]
    if cal["irf_name"] != EXPECTED_IRF or cal["package_version"] != "fermitools-2.5.8":
        raise AuditError("incorrect LAT CALDB version; V2 must not masquerade as V3")
    files = cal["caldb_files"]
    if len(files) != 9 or {(x["kind"], x["partition"]) for x in files} != {
        (a, b) for a in ("ea", "psf", "edisp") for b in KNOWN_PARTITIONS
    }:
        raise AuditError("nine P8R3_SOURCE_V3 IRF families required")
    aeff_lookup = {}
    grids = []
    for item in files:
        path = checked_file(root, item)
        if EXPECTED_IRF not in item["name"] or not item["name"].endswith(f"_{item['partition']}.fits"):
            raise AuditError("inconsistent CALDB filename and version")
        hdus = read_fits(path)
        kind, partition = item["kind"], item["partition"]
        sides = ("FRONT", "BACK") if partition == "FB" else tuple(f"{partition}{i}" for i in range(4))
        # The SOURCE_PSFn/EDISPn event types require their own response
        # families: validate their naming and grid structure separately.
        if partition == "PSF":
            sides = tuple(f"PSF{i}" for i in range(4))
        if partition == "EDISP":
            sides = tuple(f"EDISP{i}" for i in range(4))
        expected_cols = {
            "ea": {"EFFAREA"},
            "psf": {"NCORE", "NTAIL", "SCORE", "STAIL", "GCORE", "GTAIL"},
            "edisp": {"F", "S1", "K1", "BIAS", "BIAS2", "S2", "K2", "PINDEX1", "PINDEX2"},
        }[kind]
        samples = []
        for side in sides:
            ext = (f"EFFECTIVE AREA_{side}" if kind == "ea" else
                   f"RPSF_{side}" if kind == "psf" else f"ENERGY DISPERSION_{side}")
            nenergy, ntheta = (74, 32) if kind == "ea" else (23, 8)
            grid = _grid_and_finiteness(get_hdu(hdus, ext), ext, expected_cols, nenergy, ntheta)
            if kind == "ea":
                if any(x < 0 for x in grid["EFFAREA"]):
                    raise AuditError("negative LAT effective area")
                aeff_lookup[side] = grid
                samples.append({"event_type": side, "sample_effective_area_m2_1GeV_30deg": area_lookup(grid, 1000, 30)})
        grids.append({"kind": kind, "partition": partition, "file_sha256": item["sha256"],
                      "event_type_tables_validated": list(sides), "row_values_finite": True,
                      "sample_nearest_bin": samples})
    return {"irf_name": EXPECTED_IRF, "caldb_grid_files_checked": len(grids),
            "caldb_tabular_checks": grids, "fb_aeff_units": "m2",
            "phidependence_livetime_efficiency_corrections_applied": False,
            "psf_edisp_kernel_normalization_and_forward_folding_verified": False}, aeff_lookup


def area_lookup(grid: dict, energy_mev: float, theta_deg: float) -> float:
    """Nearest (energy, cos(theta)) bin, not a full LAT effective exposure."""
    e = grid["ENERG_LO"]
    ce = math.cos(math.radians(theta_deg))
    energybin = min(range(len(e)), key=lambda i: abs(math.log(max(energy_mev, 1e-30)) - math.log(math.sqrt(e[i] * grid["ENERG_HI"][i]))))
    cb = min(range(len(grid["CTHETA_LO"])), key=lambda j: abs(ce - (grid["CTHETA_LO"][j] + grid["CTHETA_HI"][j]) / 2))
    return round(float(grid["EFFAREA"][cb * len(e) + energybin]), 6)


def _tte_events(hdu: HDU) -> tuple[list[float], list[int]]:
    if int(hdu.header.get("NAXIS2", -1)) <= 0 or int(hdu.header.get("NAXIS1", -1)) != 10:
        raise AuditError("invalid TTE EVENTS row width/count")
    meta = _column_info(hdu.header)
    if meta != {"TIME": (0, 8, "D", 1), "PHA": (8, 2, "I", 1)}:
        raise AuditError("unexpected GBM TTE event columns")
    times, channels = [], []
    last = -math.inf
    for time, pha in struct.iter_unpack(">dh", hdu.content):
        if not math.isfinite(time) or time < last:
            raise AuditError("unsorted or nonfinite GBM TTE time")
        if pha < 0 or pha > 127:
            raise AuditError("GBM pulse-height channels outside 0..127")
        times.append(time)
        channels.append(pha)
        last = time
    if len(times) != hdu.header["NAXIS2"]:
        raise AuditError("GBM TTE row count inconsistent")
    return times, channels


def _valid_ebounds(hdu: HDU) -> tuple[list[float], list[float]]:
    bounds = rows(hdu, {"CHANNEL", "E_MIN", "E_MAX"})
    if len(bounds) != 128 or any(int(x["CHANNEL"]) != i for i, x in enumerate(bounds)):
        raise AuditError("GBM 128-channel EBOUNDS mismatch")
    emin = [float(r["E_MIN"]) for r in bounds]
    emax = [float(r["E_MAX"]) for r in bounds]
    if not all(math.isfinite(a) and math.isfinite(b) and 0 <= a < b for a, b in zip(emin, emax)):
        raise AuditError("GBM EBOUNDS energies invalid")
    return emin, emax


def audit_detector(det: dict, burst: dict, root: Path, manifest: dict) -> dict:
    label = det["detector"]
    if label not in {"n0", "n1", "n6", "n7", "b0", "b1"}:
        raise AuditError("unexpected selected GBM detector")
    f = det["event_and_response_files"]
    if set(f) != {"tte", "cspec", "rsp2"}:
        raise AuditError("TTE/CSPEC/RSP2 must all be present")
    paths = {role: checked_file(root, meta) for role, meta in f.items()}
    tte = fits_large_tte(paths["tte"])
    cspec = read_fits(paths["cspec"])
    rsp = read_fits(paths["rsp2"])
    trig_times = []
    for hdus in (tte, cspec, rsp):
        h0 = hdus[0].header
        if h0.get("INSTRUME") != "GBM" or h0.get("TIMESYS") != "TT":
            raise AuditError("GBM instrument/time identity mismatch")
        if h0.get("DETNAM") != (f"NAI_{int(label[1:]):02d}" if label[0] == "n" else f"BGO_{int(label[1:]):02d}"):
            raise AuditError("GBM detector identity mismatch")
        trig_times.append(float(h0["TRIGTIME"]))
    if max(trig_times) - min(trig_times) > 1e-5:
        raise AuditError("GBM TTE/CSPEC/RSP2 trigger-time mismatch")
    ebmin, ebmax = _valid_ebounds(get_hdu(tte, "EBOUNDS"))
    a, b = _valid_ebounds(get_hdu(cspec, "EBOUNDS"))
    # The v01 CSPEC for GRB 130427A was generated with channel edges that
    # differ from v00 TTE. Record, do not silently relabel channels.
    cspec_edge_relative_mismatch = max(abs(x - y) / max(y, 1) for one, two in ((a, ebmin), (b, ebmax)) for x, y in zip(one, two))
    ra, rb = _valid_ebounds(get_hdu(rsp, "EBOUNDS"))
    # The genuine archive RSP2 EBOUNDS are sometimes *not* identical to
    # TTE/CSPEC (e.g. NAI_06 on 090510). Record this measured difference,
    # and explicitly block any attempt to call it a forward-folded fit.
    rsp_edge_relative_mismatch = max(abs(x - y) / max(y, 1) for one, two in ((ra, ebmin), (rb, ebmax)) for x, y in zip(one, two))
    evh = get_hdu(tte, "EVENTS")
    times, pha = _tte_events(evh)
    gti = rows(get_hdu(tte, "GTI"), {"START", "STOP"})
    if not gti or any(x["START"] >= x["STOP"] for x in gti):
        raise AuditError("GBM GTI absent/invalid")
    start, stop = manifest["prompt_windows_s_relative_trigger"][burst["grb"]]
    bgstart, bgstop = manifest["background_window_relative_trigger_s"]
    bw = float(manifest["pulse_bin_width_s"][burst["grb"]])
    count = round((stop - start) / bw)
    if count < 1 or count > 2000 or abs(start + count * bw - stop) > bw * 0.51:
        raise AuditError("GBM time binning invalid")
    edges = [start + i * bw for i in range(count + 1)]
    # The GBM TTE EVENTS.TIME and GTI START/STOP *both* use seconds relative
    # to TRIGTIME, while FITS TSTART/TSTOP and RSP2 attitude windows use MET.
    if not (min(g["START"] for g in gti) < 0 < max(g["STOP"] for g in gti)):
        raise AuditError("GBM relative GTI not straddling trigger")
    prompt_in_gti = any(g["START"] <= start and g["STOP"] >= stop for g in gti)
    background_in_gti = any(g["START"] <= bgstart and g["STOP"] >= bgstop for g in gti)
    if not prompt_in_gti:
        raise AuditError("GBM prompt time window outside TTE GTI")
    if times[0] < min(g["START"] for g in gti) - 1e-3 or times[-1] > max(g["STOP"] for g in gti) + 1e-3:
        raise AuditError("GBM relative EVENTS.TIME inconsistent with relative GTI")
    reltime = times
    raw_energies = {"soft": (25, 100), "hard": (100, 300)} if label.startswith("n") else {"bgo": (300, 3000)}
    bands = {}
    for band, (lo, hi) in raw_energies.items():
        # Response-overlap selection, not exact photon energy reconstruction.
        channel_mask = [a < hi and b > lo for a, b in zip(ebmin, ebmax)]
        hist = [0] * count
        bgcount = 0 if background_in_gti else None
        for time, pulse in zip(reltime, pha):
            if not channel_mask[pulse]:
                continue
            if bgcount is not None and bgstart <= time < bgstop:
                bgcount += 1
            if start <= time <= stop:
                idx = min(count - 1, int((time - start) / bw))
                hist[idx] += 1
        bgrate = bgcount / (bgstop - bgstart) if bgcount is not None else None
        excess = [x - (bgrate * bw if bgrate is not None else 0.) for x in hist]
        peak_bin = max(range(len(excess)), key=excess.__getitem__)
        bands[band] = {"channel_overlap_selection": "GBM_PHA_CHANNEL_ENERGY_BOUNDS_INTERSECT_BAND",
                       "channels_selected": sum(channel_mask), "prompt_counts": sum(hist),
                       "offsource_counts": bgcount,
                       "offsource_rate_counts_per_s": round(bgrate, 6) if bgrate is not None else None,
                       "background_subtraction_available": background_in_gti,
                       "peak_time_relative_trigger_s": round(float((edges[peak_bin] + edges[peak_bin + 1]) / 2), 6),
                       "time_bin_counts": hist}
    matrices = [h for h in rsp if h.header.get("EXTNAME") == "SPECRESP MATRIX"]
    if not matrices or any(h.header.get("NAXIS2") != 140 for h in matrices):
        raise AuditError("GBM RSP2 energy matrices missing/invalid")
    matrix_spans = [(float(h.header["TSTART"]), float(h.header["TSTOP"])) for h in matrices]
    zero_duration_rsp = sum(a == b for a, b in matrix_spans)
    # GRB 090902B official RSP2 has a terminal zero-duration extension.
    # Retain its provenance, never count it as usable exposure.
    if zero_duration_rsp > 1 or any(a > b for a, b in matrix_spans) or any(a[1] > b[0] + 1e-3 for a, b in zip(matrix_spans, matrix_spans[1:])):
        raise AuditError("GBM RSP2 interval ordering invalid")
    if get_hdu(cspec, "SPECTRUM").header.get("NAXIS2", 0) < 1:
        raise AuditError("GBM CSPEC counts absent")
    return {"detector": label, "trigger_met_s": round(trig_times[0], 6),
            "tte_event_rows": int(len(times)), "tte_gti_intervals": len(gti),
            "tte_time_relative_trigger_start_stop_s": [round(float(times[0]), 3), round(float(times[-1]), 3)],
            "cspec_spectrum_rows": get_hdu(cspec, "SPECTRUM").header["NAXIS2"],
            "cspec_vs_tte_ebounds_max_relative_discrepancy": round(cspec_edge_relative_mismatch, 6),
            "cspec_tte_ebounds_exactly_aligned": cspec_edge_relative_mismatch <= 1e-5,
            "rsp2_attitude_response_matrices": len(matrices),
            "rsp2_zero_duration_terminal_matrices": zero_duration_rsp,
            "rsp2_vs_tte_ebounds_max_relative_discrepancy": round(rsp_edge_relative_mismatch, 6),
            "rsp2_tte_ebounds_exactly_aligned": rsp_edge_relative_mismatch <= 1e-5,
            "rsp2_covers_prompt": matrix_spans[0][0] <= trig_times[0] + start and matrix_spans[-1][1] >= trig_times[0] + stop,
            "rsp2_variable_length_matrix_read_and_forward_fold": False,
            "time_bin_width_s": bw, "pulse_bands": bands,
            "source_sha256": {r: f[r]["sha256"] for r in ("tte", "cspec", "rsp2")}}


def audit_lat_candidates(v4: dict, triggers: dict, aeff: dict) -> list[dict]:
    summary = []
    for src in v4["sources"]:
        grb = src["grb"]
        events = src["events"]
        areas = []
        for ev in events:
            # Pass8 EVENT_TYPE bit 0 FRONT, bit 1 BACK.
            typ = int(ev["event_type_mask"])
            side = "FRONT" if typ & 1 else "BACK" if typ & 2 else None
            if side is None:
                raise AuditError("LAT event cannot be assigned FB event type")
            areas.append(area_lookup(aeff[side], float(ev["energy_MeV"]), float(ev["theta_deg"])))
        times = [round(float(ev["met_s"]) - triggers[grb], 6) for ev in events]
        summary.append({"grb": grb, "real_lat_roi_source_candidate_count": len(events),
                        "gbm_exact_trigger_referenced_first_last_lat_candidate_s":
                        [min(times), max(times)] if times else [],
                        "nominal_nearest_fb_aeff_area_m2_min_max": [min(areas), max(areas)] if areas else [],
                        "area_units": "m2", "area_method": "NEAREST_BIN_FB_ONLY_UNCORRECTED_NOT_EXPOSURE",
                        "no_lat_background_or_source_membership_fit": True})
    return summary


def run(manifest: dict, caldb_dir: Path, gbm_dir: Path) -> dict:
    if manifest.get("schema") != "janus.cosmos.fermi-irf-gbm-manifest.v5" or manifest.get("experiment_id") != "JANUS-COSMOS-ALI-V5-GBM-LAT-IRF-2026-10-19":
        raise AuditError("incorrect V5 manifest identity")
    if manifest.get("claim_ceiling") != "IRF_GRIDS_AND_REAL_GBM_TTE_ACCEPTED_BUT_NO_TRUE_PSF_EDISP_GBM_RSP2_FORWARD_FOLDED_LAT_GBM_LIV_LIKELIHOOD":
        raise AuditError("V5 science claim ceiling weakened")
    if [x.get("grb") for x in manifest.get("gbm_bursts", [])] != list(EXPECTED_BURSTS):
        raise AuditError("not the three fixed V5 GRBs")
    if sha256(V4_PATH) != manifest["v4_frozen_observed_json_sha256"] or sha256(V4_MANIFEST_PATH) != manifest["v4_frozen_manifest_sha256"]:
        raise AuditError("V4 frozen lineage checksum mismatch")
    v4 = json.loads(V4_PATH.read_text(encoding="utf-8"))
    cal_summary, aeff = audit_caldb(manifest, caldb_dir)
    bursts = []
    triggers = {}
    for burst in manifest["gbm_bursts"]:
        grb = burst["grb"]
        if len(burst["detectors"]) != 3:
            raise AuditError("three fixed GBM detectors required per GRB")
        detectors = [audit_detector(det, burst, gbm_dir, manifest) for det in burst["detectors"]]
        if max(d["trigger_met_s"] for d in detectors) - min(d["trigger_met_s"] for d in detectors) > 1e-5:
            raise AuditError("GBM detector trigger timestamps differ")
        trig = detectors[0]["trigger_met_s"]
        triggers[grb] = trig
        bursts.append({"grb": grb, "gbm_trigger_id": burst["gbm_trigger_id"],
                       "exact_gbm_trigger_met_s": trig,
                       "naive_liv_geV_source_lag_prior_not_inferred": True,
                       "source_emission_lag_remains_unidentified": True,
                       "detectors": detectors})
    obs = {"schema": SCHEMA, "experiment_id": manifest["experiment_id"],
           "instrument_status": "REAL_GBM_TTE_CSPEC_RSP2_AND_LAT_CALDB_BYTES_AND_TABLE_AXES_VERIFIED",
           "lat_caldb": cal_summary, "gbm_bursts": bursts,
           "lat_real_candidate_aeff_diagnostics": audit_lat_candidates(v4, triggers, aeff),
           "counts": {"gbm_tte_events_all_selected_detectors": sum(d["tte_event_rows"] for b in bursts for d in b["detectors"]),
                      "gbm_tte_files": 9, "cspec_files": 9, "rsp2_files": 9, "lat_irf_files": 9,
                      "lat_v4_selected_roi_candidates": sum(len(s["events"]) for s in v4["sources"])},
           "observational_joint_liv_p_value": None,
           "observational_joint_liv_confidence_interval": None,
           "space_discreteness_detection": False,
           "black_hole_firewall_detection": False,
           "bhb_boom_detection": False,
           "scientific_gates": {
               "real_data_ingest": "PASS_SHA256_IDENTITIES_GTI_IRF_GRIDS_AND_GBM_RSP2_HEADERS",
               "lat_psf_edisp_forward_folding": "NOT_IMPLEMENTED",
               "lat_effective_exposure_and_background_forward_folding": "NOT_IMPLEMENTED",
               "gbm_rsp2_joint_spectral_forward_folding": "NOT_IMPLEMENTED",
               "relative_emission_lag_credible_physics": "NOT_IDENTIFIED",
               "joint_observational_liv_hypothesis_test": "BLOCKED_NO_VALIDATED_LIKELIHOOD",
               "bhb_firewall_quantum_model": "OPEN_UNCONFIRMED"},
           "claim_ceiling": "REAL_INSTRUMENT_TABLES_AND_PROMPT_LIGHTCURVE_DIAGNOSTICS_ONLY_NO_LIV_FIREWALL_OR_BHB_DETECTION"}
    obs["result_sha256"] = digest(obs)
    return obs


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    p.add_argument("--caldb-dir", type=Path, required=True)
    p.add_argument("--gbm-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    result = run(manifest, args.caldb_dir, args.gbm_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"result": str(args.output), "counts": result["counts"], "gates": result["scientific_gates"]}))


if __name__ == "__main__":
    main()
