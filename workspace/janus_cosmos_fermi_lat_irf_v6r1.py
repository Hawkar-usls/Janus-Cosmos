#!/usr/bin/env python3
"""V6-R1 real LAT CALDB parameter audit, NOT an event-level LAT likelihood.

Validate all nine pinned P8R3_SOURCE_V3 families, including PSF/EDISP
event-type scales, azimuth dependence, and efficiency parameter tables.
Do not evaluate a homemade PSF/EDISP probability density as calibrated when
the complete response convolution and diffuse sky model are absent.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from janus_cosmos_fermi_calibration_v5 import MANIFEST_PATH, checked_file, get_hdu
from janus_cosmos_fermi_raw_v4 import read_fits, rows
from janus_cosmos_simulation_audit import AuditError, ROOT, digest

SCHEMA = "janus.cosmos.fermi-lat-all-caldb-irf-scale-efficiency-audit.v6r1"
OUT = ROOT / "data/cosmos/COSMOS-ALI-V6R1-LAT-FULL-IRF-READINESS-2026-10-19.json"
TYPES = {"FB": ("FRONT", "BACK"),
         "PSF": ("PSF0", "PSF1", "PSF2", "PSF3"),
         "EDISP": ("EDISP0", "EDISP1", "EDISP2", "EDISP3")}
PSF_NAMES = ("NCORE", "NTAIL", "SCORE", "STAIL", "GCORE", "GTAIL")
EDISP_NAMES = ("F", "S1", "K1", "BIAS", "BIAS2", "S2", "K2", "PINDEX1", "PINDEX2")


def nearest_irf_parameter(grid: dict, col: str, energy_mev: float, theta_deg: float) -> float:
    """Nearest CALDB grid node (one diagnostic value, not IRF interpolation)."""
    e0, e1 = grid["ENERG_LO"], grid["ENERG_HI"]
    c0, c1 = grid["CTHETA_LO"], grid["CTHETA_HI"]
    if len(e0) != 23 or len(c0) != 8 or len(grid[col]) != 184 or energy_mev <= 0:
        raise AuditError("unexpected PSF/EDISP 23x8 CALDB grid")
    eidx = min(range(23), key=lambda i: abs(math.log(energy_mev / math.sqrt(e0[i] * e1[i]))))
    cos_theta = math.cos(math.radians(theta_deg))
    tidx = min(range(8), key=lambda i: abs(cos_theta - (c0[i] + c1[i]) / 2.))
    value = float(grid[col][tidx * 23 + eidx])
    if not math.isfinite(value):
        raise AuditError("nonfinite LAT CALDB parameter")
    return value


def audit_physical_parameter_bounds(kind: str, grid: dict) -> None:
    """Check physically required constraints; amplitude tables are not PDFs."""
    param_names = PSF_NAMES if kind == "psf" else EDISP_NAMES
    if any(len(grid[name]) != 184 for name in param_names):
        raise AuditError("incorrect LAT response parameter count")
    for name in param_names:
        if not all(math.isfinite(float(v)) for v in grid[name]):
            raise AuditError("nonfinite LAT IRF table parameter")
    if kind == "psf":
        positive = ("NCORE", "NTAIL", "SCORE", "STAIL")
        king = ("GCORE", "GTAIL")
        if any(min(grid[n]) < 0 for n in positive) or any(min(grid[n]) <= 1 for n in king):
            raise AuditError("invalid PSF Moffat amplitude/width/tail gamma")
        if any(min(grid[n]) <= 0 for n in ("SCORE", "STAIL")):
            raise AuditError("nonpositive PSF Moffat width")
    elif kind == "edisp":
        if min(grid["F"]) < 0 or max(grid["F"]) > 1:
            raise AuditError("EDISP mixture fraction outside 0..1")
        if any(min(grid[n]) <= 0 for n in ("S1", "S2", "K1", "K2", "PINDEX1", "PINDEX2")):
            raise AuditError("EDISP asymmetric exponential shape not positive")
    else:
        raise AuditError("unknown LAT IRF family")


def audit(manifest: dict, caldb_dir: Path) -> dict:
    if manifest["schema"] != "janus.cosmos.fermi-irf-gbm-manifest.v5" or manifest["lat_caldb"]["irf_name"] != "P8R3_SOURCE_V3":
        raise AuditError("wrong frozen LAT IRF manifest")
    files = manifest["lat_caldb"]["caldb_files"]
    if {(x["kind"], x["partition"]) for x in files} != {(k, p) for k in ("ea", "psf", "edisp") for p in TYPES} or len(files) != 9:
        raise AuditError("nine IRF families required")
    response_families = []
    for file in files:
        path = checked_file(caldb_dir, file)
        hdus = read_fits(path)
        type_rows = []
        for event_type in TYPES[file["partition"]]:
            kind = file["kind"]
            if kind == "ea":
                area = rows(get_hdu(hdus, "EFFECTIVE AREA_" + event_type), {"ENERG_LO", "ENERG_HI", "CTHETA_LO", "CTHETA_HI", "EFFAREA"})[0]
                phi = rows(get_hdu(hdus, "PHI_DEPENDENCE_" + event_type), {"PHIDEP0", "PHIDEP1"})[0]
                efficiency = rows(get_hdu(hdus, "EFFICIENCY_PARAMS_" + event_type), {"EFFICIENCY_PARS"})
                if len(area["EFFAREA"]) != 74 * 32 or len(phi["PHIDEP0"]) != 23 * 8 or len(phi["PHIDEP1"]) != 23 * 8 or len(efficiency) != 2:
                    raise AuditError("missing phi or efficiency correction parameter axes")
                if any(not math.isfinite(float(x)) for vector in (phi["PHIDEP0"], phi["PHIDEP1"]) for x in vector):
                    raise AuditError("nonfinite LAT phi dependence")
                if any(len(entry["EFFICIENCY_PARS"]) != 6 or any(not math.isfinite(float(x)) for x in entry["EFFICIENCY_PARS"]) for entry in efficiency):
                    raise AuditError("nonfinite LAT efficiency correction")
                type_rows.append({"event_type": event_type, "effective_area_cells": len(area["EFFAREA"]),
                                  "azimuth_phi_parameter_cells_per_term": len(phi["PHIDEP0"]),
                                  "efficiency_parameter_rows": len(efficiency),
                                  "phi_correction_applied_to_exposure": False,
                                  "efficiency_correction_applied_to_exposure": False})
            else:
                param_names = PSF_NAMES if kind == "psf" else EDISP_NAMES
                extension = ("RPSF_" if kind == "psf" else "ENERGY DISPERSION_") + event_type
                grid = rows(get_hdu(hdus, extension), set(param_names) | {"ENERG_LO", "ENERG_HI", "CTHETA_LO", "CTHETA_HI"})[0]
                audit_physical_parameter_bounds(kind, grid)
                scale_col = "PSFSCALE" if kind == "psf" else "EDISPSCALE"
                scale_extension = ("PSF_SCALING_PARAMS_" if kind == "psf" else "EDISP_SCALING_PARAMS_") + event_type
                scales = rows(get_hdu(hdus, scale_extension), {scale_col})[0][scale_col]
                if len(scales) != (3 if kind == "psf" else 6) or not all(math.isfinite(x) for x in scales):
                    raise AuditError("invalid IRF scaling parameters")
                anchor = {n: round(nearest_irf_parameter(grid, n, 1000., 30.), 7) for n in param_names}
                entry = {"event_type": event_type, "caldb_grid_energy_bins": 23,
                         "caldb_grid_cos_theta_bins": 8,
                         "caldb_scale_coefficients_in_native_units": [round(float(v), 9) for v in scales],
                         "nearest_grid_parameters_1GeV_theta30deg": anchor,
                         "parameter_bounds_validated": True,
                         "normalized_physical_kernel_evaluated": False}
                if kind == "psf":
                    c0, c1, power = scales
                    # This is the published LAT *angular scaling factor*
                    # sqrt((c0*(E/100 MeV)**power)**2+c1**2) in radians,
                    # NOT a 68%-containment angle or integrated point-source PSF.
                    if not (c0 > 0 and c1 > 0 and -.95 < power < -.2):
                        raise AuditError("invalid PSF angular scaling")
                    entry["psf_angular_scale_1GeV_radians_not_containment"] = round(math.hypot(c0 * (1000. / 100.) ** power, c1), 8)
                type_rows.append(entry)
        response_families.append({"kind": file["kind"], "partition": file["partition"],
                                  "byte_sha256": file["sha256"], "event_type_parameter_audits": type_rows})
    out = {"schema": SCHEMA, "input_manifest_sha256": hashlib.sha256(MANIFEST_PATH.read_bytes()).hexdigest(),
           "lat_irf": "P8R3_SOURCE_V3", "all_9_pinned_caldb_file_families_audited": True,
           "event_type_nonexclusive_families": ["FB", "PSF", "EDISP"],
           "source_event_type_use": "USE_EXACT_EVENT_TYPE_BIT_PARTITION_AND_CORRESPONDING_IRF_PER_EVENT_FOR_FUTURE_LAT_LIKELIHOOD",
           "response_families": response_families,
           "cannot_promote_to_calibrated_lat_forward_model_without": [
               "normalized event-by-event PSF angular kernels with calibrated core-tail weights",
               "normalized energy-dispersion Erec-given-Etrue kernels including scaling",
               "azimuth dependent effective area, livetime efficiency and theta/phi time convolution",
               "official Galactic and isotropic P8R3_SOURCE_V3 diffuse templates and source catalog ROI model",
               "unbinned or binned Poisson likelihood with source+background and benchmark against Fermitools",
               "physical source-emission time lag constraints and coverage tests"],
           "observed_liv_p_value": None,
           "scientific_gates": {"ACTUAL_P8R3_SOURCE_V3_PSF_EDISP_EA_TABLES": "PASS_BYTE_PINNED_PARAMETER_AUDIT",
                                "AZIMUTH_AND_LIVETIME_EFFICIENCY_TABLES": "AUDITED_NOT_APPLIED",
                                "LAT_NORMALIZED_PSF_AND_EDISP": "BLOCKED_NOT_EVALUATED",
                                "LAT_DIFFUSE_SKY_MODEL": "BLOCKED_NOT_FETCHED",
                                "LAT_PHOTON_SOURCE_PROBABILITIES": "NOT_COMPUTED",
                                "OBSERVED_LIV_STATISTIC": "BLOCKED",
                                "BHB_FIREWALL": "OPEN_UNCONFIRMED",
                                "BHB_BOOM": "NOT_ESTABLISHED"},
           "claim_ceiling": "REAL_IRF_PARAMETER_COVERAGE_AUDIT_ONLY_NO_LAT_CALIBRATED_PHOTON_LIKELIHOOD"}
    out["result_sha256"] = digest(out)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--caldb-dir", type=Path, required=True)
    p.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    p.add_argument("--output", type=Path, default=OUT)
    args = p.parse_args()
    result = audit(json.loads(args.manifest.read_text()), args.caldb_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output), "families": len(result["response_families"]),
                      "observed_p": result["observed_liv_p_value"]}))


if __name__ == "__main__":
    main()
