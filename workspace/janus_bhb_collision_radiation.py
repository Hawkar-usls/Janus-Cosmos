#!/usr/bin/env python3
"""BHB-R3A: conservation-exact e-/e+ -> two radial photons by Schwarzschild BH.

Selected kinematic channel, NOT a cross-section, angular distribution,
firewall microphysics, steady shell, radiation flux or astrophysical event.
One electron falls inward from infinity. One positron is held locally static
by an unspecified external support at finite r>r_s. All local energy and
radial momentum are conserved; the outward photon follows a null geodesic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any

from janus_bhb import BHBError, number, show
from janus_bhb_kerr import PRECISION, ROOT, geometry, photon


SCHEMA = "janus.cosmos.bhb.schwarzschild-qed-two-photon.v1"
DEFAULT_MASS = "0.00051099895"  # illustrative rounded electron rest energy, GeV
FW1_PATH = ROOT / "data/bhb/BHB-FIREWALL-FW1-RESULT-2026-10-18.json"
FW1_SHA256 = "b16a3f086f63a2df06423c2f399dd847d40e4857b0c2bc04c5f7a6b0d9923384"
FW1_INTERNAL_SHA256 = "66ac78102e190ba09f55733a294c4160b15e4aa18f1a18d2554346c6d61022b8"
EPS = Decimal("1e-57")


def digest(doc: Any) -> str:
    return hashlib.sha256(json.dumps(doc, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def verify_fw1(path: Path = FW1_PATH) -> dict[str, Any]:
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != FW1_SHA256:
        raise BHBError("FW1_FROZEN_RECEIPT_SHA256_MISMATCH")
    doc = json.loads(raw)
    if doc.get("result_sha256") != FW1_INTERNAL_SHA256:
        raise BHBError("FW1_INTERNAL_RECEIPT_SHA256_MISMATCH")
    if doc.get("gates", {}).get("FIREWALL_PHYSICAL_EXISTENCE") != "OPEN_UNCONFIRMED":
        raise BHBError("FW1_UNCONFIRMED_GATE_WAS_CHANGED")
    return {"fw1_sha256": FW1_SHA256, "fw1_internal_sha256": FW1_INTERNAL_SHA256}


def evaluate_case(raw: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise BHBError("case must be object")
    if raw.get("spacetime", "schwarzschild") != "schwarzschild":
        raise BHBError("R3A supports only Schwarzschild; not Kerr or AMPS dynamics")
    if raw.get("process") != "electron_positron_to_two_radial_photons":
        raise BHBError("only explicitly defined QED two-photon kinematic channel supported")
    if raw.get("target") != "externally_supported_static_positron":
        raise BHBError("static positron target needs explicit unsupported external support")
    with localcontext() as ctx:
        ctx.prec = PRECISION
        x = number(raw.get("r_over_rs"), "r_over_rs")
        e = number(raw.get("electron_specific_energy"), "electron_specific_energy")
        m = number(raw.get("electron_mass_energy_GeV", DEFAULT_MASS), "electron_mass_energy_GeV")
        if x <= 1 or e < 1 or m <= 0:
            raise BHBError("r/r_s>1, infall from infinity e>=1, mass>0 required")
        f = (x - 1) / x
        alpha = f.sqrt()
        q = (e * e - f).sqrt()
        gamma = e / alpha
        v_e = -q / e
        one_minus_abs_v = f / (e * (e + q))
        total_local_E = m * (gamma + 1)
        total_local_Pr = -m * q / alpha

        # Radial null photons have E=+p_out, E=-p_in in the static tetrad.
        # Rationalized positive outward solution avoids gamma-sqrt(gamma²-1)
        # catastrophic cancellation near the horizon.
        out_local = m / 2 * (1 + alpha / (e + q))
        in_local = m / 2 * (1 + (e + q) / alpha)
        ecm2 = 2 * m * m * (1 + gamma)
        ecm = ecm2.sqrt()
        balance_E = abs(total_local_E - out_local - in_local)
        balance_P = abs(total_local_Pr - (out_local - in_local))
        if max(balance_E,balance_P) > EPS * max(Decimal(1), total_local_E):
            raise BHBError("LOCAL_FOUR_MOMENTUM_CONSERVATION_FAILED")
        # Killing energies at infinity: the static support target carries
        # alpha*m, electron carries e*m. The inward photon never escapes.
        out_inf = alpha * out_local
        in_inf = alpha * in_local
        initial_inf = m * (e + alpha)
        budget_error = abs(initial_inf - out_inf - in_inf)
        if budget_error > EPS * max(Decimal(1),initial_inf):
            raise BHBError("INFINITY_KILLING_ENERGY_BUDGET_FAILED")

        radial = photon(geometry(2 * x, Decimal(0)), 2 * x, Decimal(0),
                        {"source":"hypothetical_zamo_emission","azimuth_sine":"0"})
        if radial["escape_gate"] != "EQUATORIAL_NULL_GEODESIC_REACHES_INFINITY":
            raise BHBError("radial outgoing Schwarzschild null geodesic did not escape")
        if abs(number(radial["E_inf_over_E_zamo"], "ray") - alpha) > Decimal("1e-23") * max(Decimal(1),alpha):
            raise BHBError("R2_NULL_GEODESIC_REDSHIFT_DISAGREEMENT")

        r_over_M = 2 * x
        static_accel_M = 1 / (r_over_M * r_over_M * alpha)
        row = {
            "schema": SCHEMA,
            "scenario_id": raw.get("scenario_id", "BHB-R3A-QED-CASE"),
            "process": "SELECTED_KINEMATIC_e_MINUS_e_PLUS_TO_TWO_GAMMA",
            "spacetime": "SCHWARZSCHILD_TEST_PARTICLE_EXTERIOR",
            "r_over_rs": show(x,55),
            "electron_specific_energy_from_infinity": show(e),
            "electron_rest_energy_GeV": show(m),
            "static_lapse": show(alpha),
            "electron_local_speed_over_c_signed": show(v_e,55),
            "one_minus_absolute_electron_speed_over_c": show(one_minus_abs_v,40),
            "electron_relative_gamma_to_static_positron": show(gamma),
            "E_cm_local_GeV": show(ecm),
            "E_cm_squared_GeV2": show(ecm2),
            "incoming_four_momentum_tetrad_GeV": {"energy":show(total_local_E),"radial_momentum":show(total_local_Pr)},
            "outgoing_photon": {"local_energy_GeV":show(out_local),"radial_direction":"out","E_infinity_GeV":show(out_inf),
                                "escape_gate":"PASS_RADIAL_NULL_GEODESIC_TO_INFINITY"},
            "ingoing_photon": {"local_energy_GeV":show(in_local),"radial_direction":"in","conserved_Killing_energy_GeV":show(in_inf),
                               "escape_gate":"CAPTURED_INTO_BLACK_HOLE"},
            "initial_Killing_energy_GeV": show(initial_inf),
            "fraction_initial_Killing_energy_in_escaping_photon": show(out_inf / initial_inf),
            "local_four_momentum_balance_max_residual_GeV": show(max(balance_E,balance_P), 8),
            "infinity_energy_balance_residual_GeV": show(budget_error,8),
            "local_qed_kinematic_gate": "PASS_EXACT_1_PLUS_1_FOUR_MOMENTUM",
            "null_transport_gate": "PASS_R2_SCHWARZSCHILD_LIMIT",
            "static_target_required_proper_acceleration_times_M": show(static_accel_M),
            "target_source_gate": "EXTERNAL_STATIC_SUPPORT_REQUIRED_NOT_ESTABLISHED",
            "qed_differential_cross_section": None,
            "annihilation_event_rate": None,
            "firewall_quantum_interaction_gate": "NOT_MODELED_PROXY_NOT_AMPS_FIREWALL",
            "radiation_observation_gate": "ONE_POSSIBLE_ESCAPING_PHOTON_ENERGY_NO_FLUX_OR_DETECTION",
            "claim_ceiling": "CONDITIONAL_E_PLUS_E_MINUS_COLLISION_PRODUCT_RAY_ONLY_NOT_BLACK_HOLE_BOOM",
        }
        row["result_sha256"] = digest(row)
        return row


def evaluate_batch(doc: dict[str, Any], fw1_path: Path = FW1_PATH) -> dict[str, Any]:
    if not isinstance(doc,dict) or not isinstance(doc.get("cases"),list) or not 1 <= len(doc["cases"]) <= 100:
        raise BHBError("batch requires 1..100 cases")
    source=verify_fw1(fw1_path)
    rows=[evaluate_case(c) for c in doc["cases"]]
    out={"schema":"janus.cosmos.bhb.r3a-qed-two-photon-batch.v1",
         "experiment_id":doc.get("experiment_id","BHB-R3A-ELECTRON-POSITRON"),
         "model_class":"CONDITIONAL_CLASSICAL_SCHWARZSCHILD_PLUS_LOCAL_QED_KINEMATICS",
         "fw1_source":source,"precision_digits":PRECISION,
         "cases":rows,
         "gates":{"LOCAL_4MOMENTUM_CONSERVATION":"PASS",
                  "ESCAPING_SELECTED_QED_PHOTON_KINEMATICS":"PASS",
                  "STATIONARY_SHELL_SUPPORT":"OPEN_UNESTABLISHED",
                  "PHOTON_YIELD_AND_ANGULAR_PROBABILITY":"NOT_MODELED",
                  "AMPS_FIREWALL_MICROPHYSICS":"NOT_MODELED",
                  "OBSERVED_BHB_RADIATION":"NOT_ESTABLISHED"}}
    out["result_sha256"]=digest(out)
    return out


def main(argv: list[str] | None = None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--batch",type=Path,required=True)
    ap.add_argument("--fw1",type=Path,default=FW1_PATH)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args(argv)
    try:
        raw=json.loads(args.batch.read_text(encoding="utf-8"))
        out=evaluate_batch(raw,args.fw1)
        txt=json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False)+"\n"
        if args.output:
            args.output.write_text(txt,encoding="utf-8")
        else:
            sys.stdout.write(txt)
        return 0
    except (BHBError,ValueError,OSError,KeyError) as exc:
        print(f"BHB-R3A-ERROR: {exc}",file=sys.stderr)
        return 2


if __name__=="__main__":
    raise SystemExit(main())
