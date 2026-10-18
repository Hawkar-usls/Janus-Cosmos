#!/usr/bin/env python3
"""JANUS Cosmos V2: preregistered synthetic photon time-of-flight lattice gate.

Model: one specified isotropic local nearest-neighbor photon dispersion in
a preferred comoving frame, against a Lorentz-invariant photon reference.
The FLRW travel-time expression is a small-dispersion, leading-order ansatz,
NOT a derived complete quantum-gravity or discrete-spacetime cosmology.

Three strict phases: --prepare creates synthetic observations from a separately
stored holdout truth, --score sees only observations and protocol, --reveal
checks the truth SHA and reports injection recovery. All inputs synthetic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any

from janus_cosmos_simulation_audit import AuditError, PI, PRECISION, ROOT, digest, num, show

SCHEMA = "janus.cosmos.lattice-photon-tof.v2"
PROTOCOL_SCHEMA = "janus.cosmos.lattice-tof-protocol.v2"
TRUTH_SCHEMA = "janus.cosmos.lattice-tof-synthetic-truth.v2"
OBSERVATIONS_SCHEMA = "janus.cosmos.lattice-tof-synthetic-observations.v2"
SCORE_SCHEMA = "janus.cosmos.lattice-tof-blind-score.v2"
REVEAL_SCHEMA = "janus.cosmos.lattice-tof-reveal.v2"
V1_FILE = ROOT / "data/cosmos/COSMOS-ALI-SIMULATION-AUDIT-RESULT-2026-10-18.json"
V1_SHA256 = "3b78e75eb93933e25755f501b45db1d78c4385356e0778e251055fd4353cba49"

C = Decimal("299792458")
H_SI = Decimal("6.62607015e-34")
EV_J = Decimal("1.602176634e-19")
MPC_M = Decimal("3.0856775814913673e22")  # frozen adopted SI conversion
N_SIMPSON = 384  # fixed even count, deterministic Decimal evaluation
LIV_THRESHOLD_DELTA_CHI2 = Decimal("25")  # engineering threshold, not calibrated p-value
Q_MAX = Decimal("0.01")  # enforce perturbative local dispersion across light path


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_v1_gate(path: Path = V1_FILE) -> dict[str, str]:
    if file_sha256(path) != V1_SHA256:
        raise AuditError("FROZEN_COSMOS_V1_SHA256_MISMATCH")
    item = json.loads(path.read_text(encoding="utf-8"))
    gates = item.get("gates", {})
    if gates.get("SIMULATION_HYPOTHESIS_OBSERVED") != "NOT_ESTABLISHED" or gates.get("AMPS_FIREWALL_PHYSICAL_EXISTENCE") != "OPEN_UNCONFIRMED":
        raise AuditError("FROZEN_COSMOS_V1_SCIENCE_GATE_PROMOTED")
    return {"frozen_v1_file_sha256":V1_SHA256, "simulation_gate":"NOT_ESTABLISHED", "amps_gate":"OPEN_UNCONFIRMED"}


def validate_protocol(raw: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(raw, dict) or raw.get("schema") != PROTOCOL_SCHEMA:
        raise AuditError("invalid versioned V2 TOF protocol")
    if raw.get("source_video_id") != "QBpuorUxuSM" or raw.get("all_observations") != "SYNTHETIC_ONLY":
        raise AuditError("video identity and synthetic-only input gate required")
    if raw.get("scientific_ceiling") != "NO_OBSERVED_LIV_NO_SIMULATION_DETECTION":
        raise AuditError("science ceiling cannot be changed")
    cos = raw.get("cosmology", {})
    if cos.get("assumption") != "ILLUSTRATIVE_FLAT_LAMBDA_CDM_SMALL_DISPERSION":
        raise AuditError("cosmological approximation must be declared")
    h0 = num(cos.get("H0_km_s_Mpc"), "H0")
    om = num(cos.get("Omega_m"), "Omega_m")
    ol = num(cos.get("Omega_lambda"), "Omega_lambda")
    if h0 <= 0 or not 0 <= om <= 1 or not 0 <= ol <= 1 or om + ol != 1:
        raise AuditError("positive H0 and flat 0<=Omega<=1 sum=1 required")
    if raw.get("quadrature", {}).get("method") != "SIMPSON_FIXED_DECIMAL_384":
        raise AuditError("frozen quadrature declaration required")
    sources = raw.get("sources")
    if not isinstance(sources, dict) or not {"nist","magic","beane","prior_v1"} <= sources.keys():
        raise AuditError("original research and continuity citations required")
    cases = raw.get("cases")
    if not isinstance(cases, list) or not 8 <= len(cases) <= 100:
        raise AuditError("preregister 8..100 synthetic cases")
    ids = []
    for case in cases:
        if not isinstance(case, dict) or not isinstance(case.get("case_id"), str) or not case["case_id"]:
            raise AuditError("case IDs required")
        ids.append(case["case_id"])
        z, a, lo, hi, sigma = (num(case.get(k),k) for k in
            ("redshift","lattice_spacing_m","energy_low_eV","energy_high_eV","timing_sigma_s"))
        if not 0 <= z <= 5 or a <= 0 or lo <= 0 or hi < lo or sigma <= 0:
            raise AuditError("invalid redshift, spacing, photon energy, or timing sigma")
        lag = case.get("source_lag_policy")
        if lag not in {"KNOWN_ZERO","KNOWN_VALUE","FREE_UNBOUNDED","BOUNDED_SYMMETRIC"}:
            raise AuditError("explicit source-lag treatment required")
        if lag == "KNOWN_VALUE":
            num(case.get("known_intrinsic_source_lag_s"), "known_intrinsic_source_lag_s")
        if lag == "BOUNDED_SYMMETRIC" and num(case.get("max_abs_intrinsic_source_lag_s"), "max_abs_intrinsic_source_lag_s") < 0:
            raise AuditError("source-lag bound must be nonnegative")
        with localcontext() as ctx:
            ctx.prec = PRECISION
            q_at_source = a * hi * EV_J / (2 * (H_SI/(2*PI)) * C) * (1+z)
            if q_at_source > Q_MAX:
                raise AuditError("perturbative q<=0.01 violated along source-to-observer path")
    if len(set(ids)) != len(ids):
        raise AuditError("duplicate V2 scenario IDs")
    return raw


def lattice_delay_seconds(case: dict[str, Any], cosmology: dict[str, Any]) -> Decimal:
    """Leading-order FLRW delay for the *postulated* small-q dispersion.

    omega(k)=2c/a sin(ka/2), E=hbar*omega; v_g/c=sqrt(1-q^2).
    q(E,z)=a*E_obs*(1+z)/(2*hbar*c) is an explicitly assumed adiabatic
    leading-order photon-energy redshift rule. The exact nonlinear lattice
    E(k,z) and its microscopic cosmological embedding are NOT derived.

    dt_high-low=integral_0^zs [(1-v_high/c)-(1-v_low/c)] dz/H(z),
    valid only for small absolute velocity defects and chosen flat LCDM.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        zmax = num(case["redshift"], "redshift")
        if zmax == 0 or num(case["energy_high_eV"], "energy_high_eV") == num(case["energy_low_eV"], "energy_low_eV"):
            return Decimal(0)
        a = num(case["lattice_spacing_m"], "lattice_spacing_m")
        hi = num(case["energy_high_eV"], "energy_high_eV")
        lo = num(case["energy_low_eV"], "energy_low_eV")
        hbar = H_SI / (2*PI)
        qhi = a*hi*EV_J / (2*hbar*C)
        qlo = a*lo*EV_J / (2*hbar*C)
        h0 = num(cosmology["H0_km_s_Mpc"], "H0_km_s_Mpc")*1000/MPC_M
        om = num(cosmology["Omega_m"], "Omega_m")
        ol = num(cosmology["Omega_lambda"], "Omega_lambda")

        def velocity_defect(q: Decimal) -> Decimal:
            q2 = q*q
            if q2 > Q_MAX*Q_MAX:
                raise AuditError("dispersion q out of perturbative domain")
            # Stable 1-sqrt(1-q^2) avoids catastrophic cancellation.
            return q2 / (1+(1-q2).sqrt())

        dz = zmax/N_SIMPSON
        total = Decimal(0)
        for i in range(N_SIMPSON+1):
            z = dz*i
            factor = 1+z
            h = h0*(om*factor**3+ol).sqrt()
            term = (velocity_defect(qhi*factor) - velocity_defect(qlo*factor))/h
            total += (1 if i in {0,N_SIMPSON} else (4 if i%2 else 2))*term
        return +(total*dz/3)


def _truth_lag(raw: Any, case: dict[str, Any], cosmology: dict[str, Any]) -> Decimal:
    if raw == "MIMIC_LATTICE":
        return lattice_delay_seconds(case,cosmology)/(1+num(case["redshift"],"redshift"))
    return num(raw,"intrinsic_source_lag_s")


def prepare(protocol: dict[str, Any], truth: dict[str, Any], protocol_hash: str, truth_hash: str) -> dict[str, Any]:
    validate_protocol(protocol)
    if not isinstance(truth, dict) or truth.get("schema") != TRUTH_SCHEMA or truth.get("origin") != "DETERMINISTIC_SYNTHETIC_INJECTIONS":
        raise AuditError("sealed truth must be explicitly synthetic")
    labels = truth.get("injections",[])
    if not isinstance(labels,list) or len(labels)!=len(protocol["cases"]):
        raise AuditError("sealed case matrix mismatch")
    label_by_id = {v.get("case_id"):v for v in labels if isinstance(v,dict)}
    if len(label_by_id)!=len(labels) or set(label_by_id)!={c["case_id"] for c in protocol["cases"]}:
        raise AuditError("truth ID mismatch or duplicate")
    observations = []
    for c in protocol["cases"]:
        t=label_by_id[c["case_id"]]
        if t.get("injected_model") not in {"LORENTZ_INVARIANT","LATTICE_DISPERSION"}:
            raise AuditError("unrecognized synthetic injection")
        with localcontext() as ctx:
            ctx.prec = PRECISION
            z=num(c["redshift"],"redshift")
            intrinsic=_truth_lag(t.get("intrinsic_source_lag_s"),c,protocol["cosmology"])
            noise=num(t.get("noise_s"),"noise_s")
            propagation=lattice_delay_seconds(c,protocol["cosmology"]) if t["injected_model"]=="LATTICE_DISPERSION" else Decimal(0)
            observed=(1+z)*intrinsic+propagation+noise
            observations.append({"case_id":c["case_id"],"observed_high_minus_low_arrival_s":show(observed),
                                 "timing_sigma_s":c["timing_sigma_s"], "origin":"SYNTHETIC_ONLY"})
    result={"schema":OBSERVATIONS_SCHEMA,"protocol_file_sha256":protocol_hash,"sealed_truth_file_sha256":truth_hash,
            "origin":"SYNTHETIC_ONLY","case_count":len(observations),"observations":observations,
            "preparation_gate":"SEALED_LABELS_NOT_COPIED_TO_OBSERVATIONS"}
    result["result_sha256"]=digest(result)
    return result


def _profile_chisquare(residual: Decimal, sigma: Decimal, case: dict[str, Any]) -> Decimal | None:
    with localcontext() as ctx:
        ctx.prec=PRECISION
        p=case["source_lag_policy"]
        z=num(case["redshift"],"redshift")
        if p=="FREE_UNBOUNDED":
            return None  # both hypotheses exactly profiled; no identifiability
        if p=="KNOWN_ZERO":
            fitted=Decimal(0)
        elif p=="KNOWN_VALUE":
            fitted=(1+z)*num(case["known_intrinsic_source_lag_s"],"known_intrinsic_source_lag_s")
        else:
            bound=(1+z)*num(case["max_abs_intrinsic_source_lag_s"],"max_abs_intrinsic_source_lag_s")
            fitted=min(max(residual,-bound),bound)
        return +((residual-fitted)/sigma)**2


def score(protocol: dict[str, Any], observations: dict[str, Any], protocol_hash: str,
          v1_path: Path=V1_FILE) -> dict[str, Any]:
    validate_protocol(protocol)
    frozen=frozen_v1_gate(v1_path)
    if not isinstance(observations,dict) or observations.get("schema")!=OBSERVATIONS_SCHEMA or observations.get("origin")!="SYNTHETIC_ONLY":
        raise AuditError("only known synthetic observation schema may be scored")
    if observations.get("protocol_file_sha256")!=protocol_hash or observations.get("preparation_gate")!="SEALED_LABELS_NOT_COPIED_TO_OBSERVATIONS":
        raise AuditError("observation protocol provenance mismatch")
    chk=observations.get("result_sha256")
    if not isinstance(chk,str) or chk!=digest({k:v for k,v in observations.items() if k!="result_sha256"}):
        raise AuditError("observation payload integrity mismatch")
    obs=observations.get("observations",[])
    if not isinstance(obs,list) or [o.get("case_id") for o in obs if isinstance(o,dict)] != [c["case_id"] for c in protocol["cases"]]:
        raise AuditError("observation case list mismatches preregistered protocol")
    rows=[]
    for c,o in zip(protocol["cases"],obs):
        with localcontext() as ctx:
            ctx.prec=PRECISION
            if o.get("origin")!="SYNTHETIC_ONLY" or num(o.get("timing_sigma_s"),"sigma")!=num(c["timing_sigma_s"],"sigma"):
                raise AuditError("injected observation origin or sigma changed")
            measured=num(o.get("observed_high_minus_low_arrival_s"),"arrival_lag")
            sigma=num(c["timing_sigma_s"],"sigma")
            predicted=lattice_delay_seconds(c,protocol["cosmology"])
            chi_li=_profile_chisquare(measured,sigma,c)
            chi_lattice=_profile_chisquare(measured-predicted,sigma,c)
            separation=predicted/sigma
            if chi_li is None or chi_lattice is None:
                decision="UNIDENTIFIABLE_FREE_INTRINSIC_LAG"
                delta=None
            else:
                delta=chi_li-chi_lattice
                if delta >= LIV_THRESHOLD_DELTA_CHI2:
                    decision="LATTICE_PREFERRED_SYNTHETIC"
                elif delta <= -LIV_THRESHOLD_DELTA_CHI2:
                    decision="LORENTZ_PREFERRED_SYNTHETIC"
                else:
                    decision="INCONCLUSIVE_SYNTHETIC"
            rows.append({"case_id":c["case_id"], "model":"LOCAL_NEAREST_NEIGHBOR_1D_ISOTROPIC_FLRW_ANSATZ",
                         "observed_lag_s":show(measured),"predicted_lattice_minus_lorentz_lag_s":show(predicted),
                         "timing_sigma_s":show(sigma),"ideal_separation_sigma_if_SOURCE_KNOWN":show(separation),
                         "source_lag_policy":c["source_lag_policy"],
                         "chi2_lorentz":show(chi_li) if chi_li is not None else None,
                         "chi2_lattice":show(chi_lattice) if chi_lattice is not None else None,
                         "delta_chi2_lorentz_minus_lattice":show(delta) if delta is not None else None,
                         "model_preference":decision,
                         "evidence_gate":"SYNTHETIC_CALIBRATION_ONLY_NO_OBSERVED_LIV"})
    result={"schema":SCORE_SCHEMA,"experiment_id":protocol.get("experiment_id"),
            "model_class":"ONE_SPECIFIED_QUADRATIC_LIV_LATTICE_TOY_NOT_ALL_DISCRETE_SPACETIME",
            "protocol_file_sha256":protocol_hash,"observations_internal_sha256":observations["result_sha256"],
            "blind_to_injection_labels":True,"truth_file_read_during_scoring":False,"frozen_cosmos_v1":frozen,
            "constant_definition":{"c_m_per_s":"299792458_exact_SI","h_J_s":"6.62607015e-34_exact_SI",
                                   "eV_in_J":"1.602176634e-19_exact_SI","Mpc_m":show(MPC_M),
                                   "speed":"v_g/c=sqrt(1-(a*E/(2*hbar*c))^2)",
                                   "delta_t":"integral_0^z [delta_high(zeta)-delta_low(zeta)]/H(zeta) dzeta",
                                   "delta":"q^2/(1+sqrt(1-q^2))","q":"a*E_obs*(1+z)/(2*hbar*c)",
                                   "expansion":"small q<=0.01; leading-order FLRW and adiabatic energy redshift"},
            "threshold":{"absolute_delta_chi2":show(LIV_THRESHOLD_DELTA_CHI2),
                         "status":"PREDECLARED_HEURISTIC_NOT_CALIBRATED_P_VALUE"},
            "cases":rows,"gates":{"FIXED_MODEL_TOF_CALCULATION":"PASS",
                                    "LABEL_BLIND_SCORING":"PASS_PROCESS_SEPARATION_NOT_CRYPTOGRAPHIC_BLINDING",
                                    "SOURCE_LAG_FREE_FIT":"NONIDENTIFIABLE_BY_ONE_PHOTON_PAIR",
                                    "ACTUAL_ASTRONOMICAL_EVENT_DATA":"NOT_INGESTED",
                                    "QUADRATIC_LIV_OBSERVED":"NOT_ESTABLISHED",
                                    "SIMULATION_HYPOTHESIS_OBSERVED":"NOT_ESTABLISHED",
                                    "BHB_QUANTUM_FIREWALL":"OPEN_UNCONFIRMED",
                                    "BHB_EXPLOSION":"NOT_ESTABLISHED",
                                    "OSIRIS_PROOF_AUTHORITY":"NO_CHANGE"},
            "claim_ceiling":"SYNTHETIC_IDENTIFIABILITY_AND_NEGATIVE_CONTROLS_ONLY_NO_NEW_OBSERVATIONAL_LIMIT"}
    result["result_sha256"]=digest(result)
    return result


def reveal(protocol:dict[str,Any], observations:dict[str,Any], scoring:dict[str,Any],
           truth:dict[str,Any], protocol_hash:str, truth_hash:str) -> dict[str,Any]:
    if observations.get("sealed_truth_file_sha256") != truth_hash:
        raise AuditError("SEALED_SYNTHETIC_TRUTH_SHA256_MISMATCH")
    if scoring.get("schema")!=SCORE_SCHEMA or scoring.get("protocol_file_sha256")!=protocol_hash or scoring.get("observations_internal_sha256")!=observations.get("result_sha256"):
        raise AuditError("score must match exactly frozen observations")
    if scoring.get("result_sha256") != digest({k:v for k,v in scoring.items() if k!="result_sha256"}):
        raise AuditError("scoring receipt tampered")
    regenerated=prepare(protocol,truth,protocol_hash,truth_hash)
    if regenerated!=observations:
        raise AuditError("observations not exactly reproduced by truth injection")
    names={x["case_id"]:x for x in truth["injections"]}
    rows=[]
    counters={"correct_preferences":0,"incorrect_preferences":0,"inconclusive":0,"unidentifiable":0,"known_source_lag_mismatch_false_positive_controls":0}
    for r in scoring["cases"]:
        t=names[r["case_id"]]
        d=r["model_preference"]
        true=t["injected_model"]
        preferred={"LATTICE_PREFERRED_SYNTHETIC":"LATTICE_DISPERSION",
                   "LORENTZ_PREFERRED_SYNTHETIC":"LORENTZ_INVARIANT"}.get(d)
        if preferred is None:
            c="unidentifiable" if d.startswith("UNIDENTIFIABLE") else "inconclusive"
        elif preferred==true:
            c="correct_preferences"
        else:
            c="incorrect_preferences"
        counters[c]+=1
        if t.get("intrinsic_source_lag_s")=="MIMIC_LATTICE" and preferred=="LATTICE_DISPERSION" and true=="LORENTZ_INVARIANT":
            counters["known_source_lag_mismatch_false_positive_controls"]+=1
        rows.append({"case_id":r["case_id"],"injected_model":true,
                     "intrinsic_source_lag_generation":t["intrinsic_source_lag_s"],
                     "synthetic_noise_s":t["noise_s"],"blind_model_preference":d,
                     "recovery_class":c})
    result={"schema":REVEAL_SCHEMA,"origin":"DETERMINISTIC_SYNTHETIC_ONLY",
            "protocol_file_sha256":protocol_hash,"truth_file_sha256":truth_hash,
            "observations_internal_sha256":observations["result_sha256"],
            "blind_score_internal_sha256":scoring["result_sha256"],
            "cases":rows,"counters":counters,
            "gates":{"SEALED_TRUTH_REPLAY":"PASS_BYTE_EXACT",
                    "SOURCE_INTRINSIC_LAG_CAN_MIMIC_LIV":"DEMONSTRATED_BY_NEGATIVE_CONTROL",
                    "SYNTHETIC_INJECTION_RECOVERY":"PASS_TEST_ONLY",
                    "ACTUAL_LIV_OBSERVED":"NOT_ESTABLISHED",
                    "QUANTUM_FIREWALL":"OPEN_UNCONFIRMED",
                    "SIMULATED_UNIVERSE":"NOT_ESTABLISHED"}}
    result["result_sha256"]=digest(result)
    return result


def main(argv:list[str]|None=None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    flags=ap.add_mutually_exclusive_group(required=True)
    flags.add_argument("--prepare",action="store_true")
    flags.add_argument("--score",action="store_true")
    flags.add_argument("--reveal",action="store_true")
    ap.add_argument("--protocol",type=Path,required=True)
    ap.add_argument("--truth",type=Path)
    ap.add_argument("--observations",type=Path)
    ap.add_argument("--blind-score",type=Path)
    ap.add_argument("--output",type=Path)
    ap.add_argument("--v1",type=Path,default=V1_FILE)
    args=ap.parse_args(argv)
    try:
        proto=json.loads(args.protocol.read_text(encoding="utf-8"))
        validate_protocol(proto)
        proto_hash=file_sha256(args.protocol)
        if args.prepare:
            if args.truth is None:
                raise AuditError("--prepare requires --truth")
            result=prepare(proto,json.loads(args.truth.read_text(encoding="utf-8")),proto_hash,file_sha256(args.truth))
        else:
            if args.observations is None:
                raise AuditError("--observations required")
            obs=json.loads(args.observations.read_text(encoding="utf-8"))
            if args.score:
                result=score(proto,obs,proto_hash,args.v1)
            else:
                if args.truth is None or args.blind_score is None:
                    raise AuditError("--reveal requires --truth and --blind-score")
                s=json.loads(args.blind_score.read_text(encoding="utf-8"))
                result=reveal(proto,obs,s,json.loads(args.truth.read_text(encoding="utf-8")),proto_hash,file_sha256(args.truth))
        content=json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n"
        if args.output:
            args.output.write_text(content,encoding="utf-8")
        else:
            sys.stdout.write(content)
        return 0
    except (AuditError,ValueError,OSError,KeyError) as exc:
        print(f"JANUS-TOF-V2-ERROR: {exc}",file=sys.stderr)
        return 2


if __name__=="__main__":
    raise SystemExit(main())
