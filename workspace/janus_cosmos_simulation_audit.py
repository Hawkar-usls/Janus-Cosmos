#!/usr/bin/env python3
"""JANUS Cosmos V1: falsifiable physics controls for an ALI video claim audit.

This is *not* a detector of a simulated universe. Each mode is a defined toy
or established-theory benchmark. No philosophical inference is a physics gate.
The original video's full transcript was unavailable; timestamp attributions
refer to an identified third-party summary, not verified verbatim quotes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "janus.cosmos.video-physics-audit.v1"
PRECISION = 70
PI = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286")
C_SI = 299792458  # exact m/s in the SI definition, not an inferred CPU tick
FROZEN_R3A_PATH = ROOT / "data/bhb/BHB-R3A-QED-COLLISION-RESULT-2026-10-18.json"
FROZEN_R3A_SHA256 = "f5bd22dc7864b54775011b8cfb9cce5918b58574b121f3e0f8f1a7c5ad092502"
FROZEN_R3A_INTERNAL = "1d772de42de04ac8d6079e1c0d72d0012f42ae2c9efe1c954aeb7019dab18e1b"
ALLOWED_ASSESSMENTS = frozenset({
    "ESTABLISHED_THEORY_NOT_SIMULATION_EVIDENCE",
    "MATHEMATICAL_COUNTEREXAMPLE_NOT_COSMOLOGY_OBSERVATION",
    "OPEN_THEORETICAL_PROPOSAL",
    "CONDITIONAL_PHILOSOPHICAL_ARGUMENT",
    "INTERPRETATION_NOT_EMPIRICAL_DETECTION",
    "MISLEADING_IF_PRESENTED_AS_PROOF",
    "MODEL_SPECIFIC_TEST_POSSIBLE_NOT_DETECTION",
    "UNESTABLISHED_EMPIRICAL_CLAIM",
})
MODELS = frozenset({"lorentz", "schwarzschild", "bh_area_entropy", "flat_torus", "lattice_wave", "conditional_observer_count"})


class AuditError(ValueError):
    """Explicit invalid premise, missing provenance, or scientific promotion."""


def digest(value: Any) -> str:
    blob = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def num(raw: Any, label: str) -> Decimal:
    if isinstance(raw, bool) or not isinstance(raw, (str, int)):
        raise AuditError(f"{label}: use decimal string or integer")
    try:
        value = Decimal(raw)
    except InvalidOperation as exc:
        raise AuditError(f"{label}: invalid decimal") from exc
    if not value.is_finite():
        raise AuditError(f"{label}: finite value required")
    return value


def positive_int(raw: Any, label: str, zero_ok: bool = False) -> int:
    if isinstance(raw, bool) or not isinstance(raw, (str, int)):
        raise AuditError(f"{label}: integer string required")
    s = str(raw)
    if not s.isascii() or not s.isdecimal():
        raise AuditError(f"{label}: non-negative integer required")
    n = int(s)
    if n < 0 or (n == 0 and not zero_ok) or n > 10**30:
        raise AuditError(f"{label}: out of permitted range")
    return n


def show(d: Decimal) -> str:
    return str(+d)


def sine_cosine(x: Decimal) -> tuple[Decimal, Decimal]:
    """Deterministic Decimal Taylor series for 0 <= x <= pi/2."""
    if not Decimal(0) <= x <= PI / 2:
        raise AuditError("series argument outside [0,pi/2]")
    sin_term = x
    cos_term = Decimal(1)
    s, c = sin_term, cos_term
    tolerance = Decimal("1e-66")
    for n in range(1, 90):
        sin_term *= -x * x / (Decimal(2*n) * Decimal(2*n+1))
        cos_term *= -x * x / (Decimal(2*n-1) * Decimal(2*n))
        s += sin_term
        c += cos_term
        if abs(sin_term) < tolerance and abs(cos_term) < tolerance:
            return s, c
    raise AuditError("trigonometric series did not converge")


def verify_frozen_r3a(path: Path = FROZEN_R3A_PATH) -> dict[str, str]:
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != FROZEN_R3A_SHA256:
        raise AuditError("BHB_R3A_FROZEN_SHA256_MISMATCH")
    record = json.loads(raw)
    if record.get("result_sha256") != FROZEN_R3A_INTERNAL:
        raise AuditError("BHB_R3A_INTERNAL_SHA256_MISMATCH")
    gates = record.get("gates", {})
    if gates.get("AMPS_FIREWALL_MICROPHYSICS") != "NOT_MODELED" or gates.get("OBSERVED_BHB_RADIATION") != "NOT_ESTABLISHED":
        raise AuditError("BHB_R3A_SCIENTIFIC_GATE_PROMOTED")
    return {"file_sha256": FROZEN_R3A_SHA256, "internal_sha256": FROZEN_R3A_INTERNAL,
            "observation_gate": "NOT_ESTABLISHED", "firewall_gate": "NOT_MODELED"}


def evaluate_case(case: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(case, dict) or case.get("model") not in MODELS:
        raise AuditError("case object with supported model required")
    name = case.get("case_id")
    if not isinstance(name, str) or not name or len(name) > 120:
        raise AuditError("case_id required")
    model = case["model"]
    with localcontext() as ctx:
        ctx.prec = PRECISION
        p: dict[str, Any] = {}
        if model == "lorentz":
            beta = num(case.get("speed_over_c"), "speed_over_c")
            u = num(case.get("signal_speed_over_c", "1"), "signal_speed_over_c")
            if not Decimal(-1) < beta < Decimal(1) or abs(u) > 1:
                raise AuditError("Lorentz case requires |beta|<1, |signal/c|<=1")
            factor = (1 - beta * beta).sqrt()
            gamma = 1 / factor
            added = (u + beta) / (1 + u * beta)
            if abs(added) > 1:
                raise AuditError("EINSTEIN_VELOCITY_ADDITION_SUPERLUMINAL")
            p = {"beta": show(beta), "signal_in_rest_frame_over_c": show(u),
                 "gamma": show(gamma), "proper_time_over_coordinate_time": show(factor),
                 "contracted_length_over_rest_length": show(factor),
                 "transformed_signal_speed_over_c": show(added),
                 "exact_si_light_speed_m_per_s": C_SI,
                 "interpretation": "SPECIAL_RELATIVITY_NOT_CPU_FRAME_RATE"}
        elif model == "schwarzschild":
            x = num(case.get("r_over_rs"), "r_over_rs")
            e = num(case.get("specific_energy", "1"), "specific_energy")
            if x <= 1 or e < 1:
                raise AuditError("static exterior r/r_s>1 and infall-from-infinity e>=1 required")
            alpha = ((x - 1) / x).sqrt()
            beta = (1 - alpha*alpha/(e*e)).sqrt()
            p = {"r_over_rs":show(x), "static_lapse":show(alpha),
                 "infall_speed_relative_to_static_observer_over_c":show(beta),
                 "one_minus_infall_speed_over_c":show(alpha*alpha/(e*(e+(e*e-alpha*alpha).sqrt()))),
                 "relative_lorentz_gamma":show(e/alpha),
                 "outgoing_radial_photon_local_speed_over_c":"1",
                 "outgoing_radial_photon_E_infinity_over_E_local":show(alpha),
                 "horizon_crossing_kinematics":"STATIC_TETRAD_DEFINED_ONLY_OUTSIDE_HORIZON",
                 "firewall_existence":"OPEN_UNCONFIRMED"}
        elif model == "bh_area_entropy":
            radius = num(case.get("schwarzschild_radius_over_planck_length"), "radius/l_planck")
            if radius <= 0:
                raise AuditError("positive radius required")
            area = 4*PI*radius*radius
            entropy = area/4
            p = {"radius_over_planck_length":show(radius),
                 "area_over_planck_length_squared":show(area),
                 "bekenstein_hawking_entropy_over_kB":show(entropy),
                 "dimensionless_entropy_in_bits":show(entropy/Decimal(2).ln()),
                 "area_law":"S/kB=A/(4*l_P^2)",
                 "interpretation":"AREA_SCALING_NOT_COMPUTER_SCREEN_EVIDENCE"}
        elif model == "flat_torus":
            lengths = [num(case.get(k), k) for k in ("Lx", "Ly", "Lz")]
            if min(lengths) <= 0:
                raise AuditError("positive side lengths required")
            p = {"identified_side_lengths":[show(v) for v in lengths],
                 "volume":show(lengths[0]*lengths[1]*lengths[2]),
                 "shortest_noncontractible_loop":show(min(lengths)),
                 "intrinsic_spatial_curvature":"0",
                 "compact":True, "has_boundary":False,
                 "interpretation":"FLAT_AND_FINITE_MATHEMATICALLY_COMPATIBLE_NO_MEASURED_UNIVERSE_TOPOLOGY"}
        elif model == "lattice_wave":
            theta = num(case.get("wave_number_times_spacing"), "ka")
            if not 0 < theta <= PI:
                raise AuditError("lattice nearest-neighbor model requires 0<ka<=pi")
            sine, cosine = sine_cosine(theta/2)
            p = {"wave_number_times_spacing":show(theta),
                 "angular_frequency_times_spacing_over_c":show(2*sine),
                 "phase_velocity_over_c":show(2*sine/theta),
                 "group_velocity_over_c":show(cosine),
                 "low_k_continuum_limit":"1",
                 "dispersion_relation":"omega=(2*c/a)*sin(k*a/2)",
                 "interpretation":"CHOSEN_NEAREST_NEIGHBOR_LATTICE_TOY_NOT_PLANCK_PIXELS_OR_UNIVERSE_FPS"}
        elif model == "conditional_observer_count":
            original = positive_int(case.get("non_simulated_observers"), "non_simulated_observers")
            simulated = positive_int(case.get("simulated_observers"), "simulated_observers", zero_ok=True)
            p = {"non_simulated_observers":str(original),"simulated_observers":str(simulated),
                 "illustrative_simulation_fraction":show(Decimal(simulated)/Decimal(original+simulated)),
                 "assumption":"EQUALLY_WEIGHTED_COUNTING_WITH_STIPULATED_POPULATIONS_ONLY",
                 "empirical_likelihood_or_prior":None,
                 "interpretation":"CONDITIONAL_COUNT_NOT_ESTIMATED_PROBABILITY_OF_OUR_UNIVERSE"}
        result = {"case_id":name,"model":model,"outputs":p,
                  "numerical_gate":"PASS_DEFINED_MODEL_ONLY", "simulation_detection_gate":"NOT_ESTABLISHED"}
        result["case_sha256"] = digest(result)
        return result


def evaluate_batch(doc: dict[str, Any], r3a_path: Path = FROZEN_R3A_PATH) -> dict[str, Any]:
    if not isinstance(doc, dict) or doc.get("schema") != "janus.cosmos.video-physics-cases.v1":
        raise AuditError("versioned input schema required")
    video = doc.get("video", {})
    if video.get("youtube_id") != "QBpuorUxuSM" or video.get("transcript_status") != "UNAVAILABLE_THIRD_PARTY_TIMESTAMPED_SUMMARY":
        raise AuditError("video source identity and uncertainty required")
    sources = doc.get("sources")
    if not isinstance(sources, dict) or not {"video_summary","nist_si","amps","bostrom","lattice_paper","topology_paper","black_hole_entropy"} <= sources.keys():
        raise AuditError("source provenance incomplete")
    raw_cases, claims = doc.get("cases"), doc.get("claims")
    if not isinstance(raw_cases, list) or not 6 <= len(raw_cases) <= 100:
        raise AuditError("need 6..100 cases")
    if not isinstance(claims, list) or not 12 <= len(claims) <= 100:
        raise AuditError("need 12..100 video-claim records")
    ids = [c.get("case_id") for c in raw_cases if isinstance(c, dict)]
    if len(set(ids)) != len(raw_cases):
        raise AuditError("case IDs must be unique")
    rows = [evaluate_case(c) for c in raw_cases]
    if {r["model"] for r in rows} != MODELS:
        raise AuditError("all six independent mathematical control modes required")
    claim_ids: set[str] = set()
    for c in claims:
        if not isinstance(c, dict) or not isinstance(c.get("claim_id"), str) or not c["claim_id"] or c["claim_id"] in claim_ids:
            raise AuditError("distinct nonempty claim IDs required")
        claim_ids.add(c["claim_id"])
        if c.get("assessment") not in ALLOWED_ASSESSMENTS or c.get("evidence_status") != "SECONDARY_SUMMARY_NOT_VERBATIM":
            raise AuditError("cannot promote unverified summary or philosophical claim to established observation")
        if not isinstance(c.get("timecode"), str) or not isinstance(c.get("paraphrase_ru"), str) or not c.get("assessment_detail_ru") or not c.get("falsification_or_control_ru"):
            raise AuditError("detailed timestamp/assessment/control fields required")
        if not isinstance(c.get("reference_ids"), list) or not c["reference_ids"] or any(s not in sources for s in c["reference_ids"]):
            raise AuditError("unknown or empty supporting sources")
        if any(cid not in ids for cid in c.get("linked_case_ids", [])):
            raise AuditError("claim points to unknown case")
    frozen = verify_frozen_r3a(r3a_path)
    # Pinned continuity check: same local-frame Schwarzschild lapse/gamma
    # for R3A's seven-case receipt. This does NOT grant any firewall physics.
    r3a = json.loads(r3a_path.read_text(encoding="utf-8"))
    comparison = []
    with localcontext() as ctx:
        ctx.prec = PRECISION
        for row in rows:
            if row["model"] != "schwarzschild":
                continue
            prior = next((p for p in r3a["cases"] if p["r_over_rs"] == row["outputs"]["r_over_rs"] and p["electron_specific_energy_from_infinity"] == "1"), None)
            if prior is None:
                raise AuditError("Schwarzschild video audit case missing in frozen R3A benchmark")
            a = Decimal(prior["static_lapse"])
            b = Decimal(row["outputs"]["static_lapse"])
            # Existing frozen R3A receipt intentionally prints only ~25
            # significant decimal digits; compare to its stated precision.
            if abs(a-b) > Decimal("1e-23") * abs(b):
                raise AuditError("video-audit lapse contradicts frozen R3A")
            comparison.append({"case_id":row["case_id"],"r3a_scenario_id":prior["scenario_id"],
                               "lapse_absolute_difference":show(abs(a-b)),"gate":"PASS"})
    counts = {s:sum(c["assessment"] == s for c in claims) for s in sorted(ALLOWED_ASSESSMENTS) if any(c["assessment"] == s for c in claims)}
    result = {"schema":SCHEMA,"experiment_id":doc.get("experiment_id","JANUS-ALI-VIDEO-V1"),
              "video":{"youtube_id":video["youtube_id"],"title":video.get("title"),
                       "transcript_status":video["transcript_status"],
                       "attribution_scope":"SECONDARY_TIMESTAMPED_SUMMARY_NOT_VIDEO_VERBATIM"},
              "models":sorted(MODELS),"cases":rows,"case_count":len(rows),
              "claims_count":len(claims),"claims_sha256":digest(claims),
              "assessment_counts":counts,"r3a_pin":frozen,"r3a_cross_checks":comparison,
              "gates":{"SIX_MODEL_NUMERICAL_CONTROLS":"PASS",
                       "FROZEN_BHB_R3A_CONTINUITY":"PASS",
                       "VIDEO_FULL_TRANSCRIPT_VERIFIED":"NO_SECONDARY_SUMMARY_ONLY",
                       "AMPS_FIREWALL_PHYSICAL_EXISTENCE":"OPEN_UNCONFIRMED",
                       "HOLOGRAPHY_IMPLIED_COMPUTER_SIMULATION":"NOT_ESTABLISHED",
                       "LATTICE_SPACING_OBSERVED":"NOT_ESTABLISHED",
                       "SIMULATION_HYPOTHESIS_OBSERVED":"NOT_ESTABLISHED",
                       "BHB_EXPLOSION_OBSERVED":"NOT_ESTABLISHED",
                       "OSIRIS_PROOF_AUTHORITY":"NO_CHANGE"},
              "scientific_claim_ceiling":"CHECKED_SELECTED_MODELS_AND_SECONDARY_VIDEO_SUMMARY_NOT_PROOF_OF_SIMULATION_FIREWALL_OR_BHB_BOOM"}
    result["result_sha256"] = digest(result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--r3a", type=Path, default=FROZEN_R3A_PATH)
    args = parser.parse_args(argv)
    try:
        output = evaluate_batch(json.loads(args.batch.read_text(encoding="utf-8")), args.r3a)
        content = json.dumps(output, indent=2, ensure_ascii=False, sort_keys=True) + "\n"
        if args.output:
            args.output.write_text(content, encoding="utf-8")
        else:
            sys.stdout.write(content)
    except (AuditError, ValueError, KeyError, OSError) as exc:
        print(f"JANUS-VIDEO-PHYSICS-ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
