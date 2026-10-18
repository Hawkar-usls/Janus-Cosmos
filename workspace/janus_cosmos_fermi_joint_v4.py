#!/usr/bin/env python3
"""JANUS Cosmos V4: multi-GRB Gaussian timing likelihood, SYNTHETIC ONLY.

This is a numerical identifiability and coverage harness, *not* a calibrated
likelihood for raw Fermi LAT counts. Real photons need the LAT IRF (Aeff, PSF,
energy dispersion, exposure, background), GBM prompt source evolution and a
justified source-emission-lag prior before any astrophysical likelihood ratio.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from functools import lru_cache
from pathlib import Path
from typing import Any

from janus_cosmos_lattice_tof import lattice_delay_seconds
from janus_cosmos_simulation_audit import AuditError, digest

SCHEMA="janus.cosmos.fermi-multiburst-identifiability-synthetic.v4"
COSMOLOGY={"H0_km_s_Mpc":"70","Omega_m":"0.3","Omega_lambda":"0.7"}
A_REF="4e-27"  # A test normalization, not measured or fitted astrophysically.
E_REF_GEV=1.
ETA_RANGE=[0.,2.]
INPUTS=[("090510", "0.903", 0.4), ("090902B","1.822",1.1), ("130427A","0.34",2.3)]


@lru_cache(maxsize=32)
def _k_scale(z:str) -> float:
    q=lattice_delay_seconds({"redshift":z,"lattice_spacing_m":A_REF,
        "energy_low_eV":"1e9","energy_high_eV":"1e11"},COSMOLOGY)
    return float(q)/(100**2-E_REF_GEV**2)


def _k_at_energy(z:str,energy_gev:float) -> float:
    """V2 small-q leading-order delay; use its E^2 response as toy coefficient."""
    if energy_gev<1 or not math.isfinite(energy_gev):raise AuditError("toy energy invalid")
    return _k_scale(z)*(energy_gev**2-E_REF_GEV**2)


def toy_bursts(seed:int=20261018,eta_true:float=0.8,source_beta:float=0.,sigma_s:float=0.012) -> list[dict[str,Any]]:
    rng=random.Random(seed)
    result=[]
    for name,z,offset in INPUTS:
        events=[]
        for i in range(42):
            # Broad energy lever arm is an intentionally favorable laboratory toy.
            energy=2+118*(i+0.5)/42
            lag=_k_at_energy(z,energy)
            events.append({"energy_GeV":energy,"k_s":lag,
                           "t_s":offset+(eta_true+source_beta)*lag+rng.gauss(0,sigma_s)})
        rng.shuffle(events)
        result.append({"grb":name,"redshift":z,"sigma_s":sigma_s,"events":events})
    return result


def loss_at(bursts:list[dict[str,Any]],eta:float,source_mode:str,
            beta_prior_sigma:float|None=None) -> float:
    """Profile intercept and burst-specific source lag: beta_j*k(E,z).

    A completely free intrinsic beta_j is algebraically degenerate with eta,
    because the observed coefficient is beta_j+eta for each burst. The
    Gaussian beta prior is a modeling assumption, not a measured constraint.
    """
    if not math.isfinite(eta):raise AuditError("nonfinite toy eta")
    if source_mode not in {"known_zero","gaussian_prior","unbounded"}:
        raise AuditError("invalid source-lag mode")
    if source_mode=="gaussian_prior" and (beta_prior_sigma is None or beta_prior_sigma<=0):
        raise AuditError("missing positive source prior")
    total=0.
    for b in bursts:
        ev=b["events"]
        if len(ev)<3: raise AuditError("too few synthetic photons per burst")
        sig=float(b["sigma_s"])
        if sig<=0 or not math.isfinite(sig):raise AuditError("invalid toy measurement sigma")
        xs=[float(r["k_s"]) for r in ev]
        ys=[float(r["t_s"]) for r in ev]
        if any(not math.isfinite(n) for n in xs+ys):raise AuditError("nonfinite toy photon")
        mx=sum(xs)/len(xs); my=sum(ys)/len(ys)
        xc=[v-mx for v in xs];yc=[v-my for v in ys]
        info=sum(x*x for x in xc)/(sig*sig)
        cross=sum(x*y for x,y in zip(xc,yc))/(sig*sig)
        if info<=0:raise AuditError("zero energy lever arm")
        if source_mode=="known_zero":beta=0.
        elif source_mode=="unbounded":beta=cross/info-eta
        else:beta=(cross-eta*info)/(info+1/(beta_prior_sigma**2))
        resid=[(y-(eta+beta)*x)/sig for x,y in zip(xc,yc)]
        total+=sum(x*x for x in resid)
        if source_mode=="gaussian_prior": total+=(beta/beta_prior_sigma)**2
    return total


def fit_grid(bursts:list[dict[str,Any]],source_mode:str,
             beta_prior_sigma:float|None=None) -> dict[str,Any]:
    # Profiling intercept and linear intrinsic slope yields an exact quadratic
    # in eta. Three evaluations recover its optimum with no grid-resolution
    # artifacts in the boundary-null calibration.
    c0=loss_at(bursts,0.,source_mode,beta_prior_sigma)
    c1=loss_at(bursts,1.,source_mode,beta_prior_sigma)
    c2=loss_at(bursts,2.,source_mode,beta_prior_sigma)
    if source_mode=="unbounded":
        # Exact structural non-identifiability, irrespective of roundoff.
        if max(c0,c1,c2)-min(c0,c1,c2)>1e-7*max(1,min(c0,c1,c2)):
            raise AuditError("unbounded lag unexpectedly constrained eta")
        return {"eta_mle":None,"twice_delta_log_likelihood":0.,
                "identifiability":"EXACT_DEGENERACY_BETA_PER_BURST_PLUS_ETA"}
    a=(c2-2*c1+c0)/2
    if a<=0:raise AuditError("non-convex conditional profile")
    b=c1-c0-a
    best=max(ETA_RANGE[0],min(ETA_RANGE[1],-b/(2*a)))
    cmin=loss_at(bursts,best,source_mode,beta_prior_sigma)
    return {"eta_mle":round(best,8),"twice_delta_log_likelihood":round(max(0.,c0-cmin),9),
            "identifiability":"CONDITIONAL_ON_ASSUMED_SOURCE_PRIOR"}


def run(seed:int=20261018,trials:int=160) -> dict[str,Any]:
    if trials<50 or trials>2000:raise AuditError("null toy trial count outside 50..2000")
    positive=toy_bursts(seed=seed,eta_true=0.8)
    conditional=fit_grid(positive,"known_zero")
    free=fit_grid(positive,"unbounded")
    prior=fit_grid(positive,"gaussian_prior",beta_prior_sigma=0.5)
    confounded=toy_bursts(seed=seed,eta_true=0.,source_beta=0.8)
    confounded_wrong=fit_grid(confounded,"known_zero")
    confounded_free=fit_grid(confounded,"unbounded")
    null_toys=[]
    for i in range(trials):
        case=toy_bursts(seed=seed+100+i,eta_true=0.)
        null_toys.append(fit_grid(case,"known_zero")["twice_delta_log_likelihood"])
    exceed=sum(ts>=conditional["twice_delta_log_likelihood"] for ts in null_toys)
    false_flags=sum(ts>=3.841458820694124 for ts in null_toys)
    result={"schema":SCHEMA,"experiment_id":"JANUS-COSMOS-ALI-V4-FERMI-LAT-RAW-JOINT-2026-10-18",
        "data_kind":"DETERMINISTIC_GAUSSIAN_SYNTHETIC_TIMING_ONLY_NOT_OBSERVED_FERMI_EVENTS",
        "seed":seed,"n_bursts":len(positive),"n_synthetic_photons":sum(len(b["events"]) for b in positive),
        "synthetic_injected_eta":0.8,"eta_definition":"(a/4e-27_m)^2_leading_order_V2_relative_to_1GeV",
        "lattice_spacing_reference_m":A_REF,"cosmology":COSMOLOGY,
        "known_zero_intrinsic_lag":{"fit":conditional,"interpretation":"Only identifiable if source delays are independently known to vanish"},
        "gaussian_intrinsic_lag_prior":{"source_beta_sigma":0.5,"fit":prior,
                "interpretation":"Bayesian/ridge assumption drives an optional conditional result; prior is not empirically calibrated"},
        "unbounded_intrinsic_lag":{"fit":free,"interpretation":"Exact eta + beta_j degeneracy; no observed speed bound"},
        "intrinsic_mimic_negative_control":{"injected_propagation_eta":0.,"injected_source_beta_per_burst":0.8,
             "incorrect_source_fixed_zero_fit":confounded_wrong,
             "correct_source_free_fit":confounded_free},
        "synthetic_null_calibration":{"trials":trials,"seed_base":seed+100,
            "test":"TOY_ONLY_BOUNDARY_PROFILE_TS_ETA_GE_ZERO", "conditional_injection_empirical_p_upper":(exceed+1)/(trials+1),
            "null_fraction_TS_ge_3p841":round(false_flags/trials,6),
            "not_an_observed_p_value":True},
        "real_lat_likelihood_result":None,
        "gates":{"SYNTHETIC_JOINT_PROFILE":"PASS_CONDITIONAL_ONLY", "SYNTHETIC_NULL_CALIBRATION":"PASS_TOY_ONLY",
            "SOURCE_LAG_UNBOUNDED":"EXACT_NONIDENTIFIABILITY_DEMONSTRATED",
            "ASTROPHYSICAL_MULTIBURST_LRT":"BLOCKED_NO_FORWARD_FOLDED_INSTRUMENT_AND_SOURCE_MODEL",
            "LIV_OBSERVED":"NOT_ESTABLISHED","SPACETIME_DISCRETENESS":"NOT_ESTABLISHED",
            "AMPS_FIREWALL":"OPEN_UNCONFIRMED","BHB_BOOM":"NOT_ESTABLISHED"},
        "claim_ceiling":"A synthetic likelihood engine and identifiability stress test; no measured limit, discovery statistic, astrophysical p-value or black-hole-firewall inference."}
    result["result_sha256"]=digest(result)
    return result


def main(argv:list[str]|None=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed",type=int,default=20261018)
    parser.add_argument("--trials",type=int,default=160)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args(argv)
    try:
        body=json.dumps(run(args.seed,args.trials),indent=2,sort_keys=True,ensure_ascii=False)+"\n"
        if args.output:args.output.write_text(body,encoding="utf-8")
        else:print(body,end="")
    except (AuditError,ValueError,OSError) as exc:
        print(f"JANUS-FERMI-V4-JOINT-ERROR: {exc}",file=__import__("sys").stderr)
        return 2
    return 0


if __name__=="__main__":
    raise SystemExit(main())
