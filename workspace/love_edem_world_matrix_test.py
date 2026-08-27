#!/usr/bin/env python3
"""
LOVE -- EDEM robustness test for JANUS COSMOS WORLD MATRIX.

This test does NOT search for a new object. It asks a narrower preregisterable
question: do the already-frozen catalog nearest-neighbour distances collapse
to an exact positional counterpart under an explicitly assumed sub-arcsecond
center-jitter model?

The positional sigma used here is a stress-test prior unless an upstream source
provides a measured uncertainty. The result therefore has a narrow claim ceiling.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path
from typing import Dict, List, Sequence

from janus_world_matrix import DimensionSpec, JanusWorldMatrix, _quantile, _sha256_obj

SCHEMA = "janus.cosmos.love-edem.world-matrix-robustness.v1"
DEFAULT_SIGMA_ARCSEC = 0.7
DEFAULT_SAMPLES = 8192
DEFAULT_SEED = 16016
GATES_ARCSEC = (1.5, 3.0)


def angular_separation_arcsec(ra1_deg: float, dec1_deg: float, ra2_deg: float, dec2_deg: float) -> float:
    ra1,dec1,ra2,dec2 = map(math.radians,(ra1_deg,dec1_deg,ra2_deg,dec2_deg))
    sd = math.sin((dec2-dec1)/2.0)
    sr = math.sin((ra2-ra1)/2.0)
    a = sd*sd + math.cos(dec1)*math.cos(dec2)*sr*sr
    a = min(1.0,max(0.0,a))
    return math.degrees(2.0*math.asin(math.sqrt(a)))*3600.0


def _stats(xs: Sequence[float]) -> Dict[str,float]:
    return {
        "mean": statistics.fmean(xs),
        "std": statistics.pstdev(xs),
        "q05": _quantile(xs,.05),
        "q50": _quantile(xs,.50),
        "q95": _quantile(xs,.95),
        "min": min(xs),
        "max": max(xs),
    }


def run_target(name: str, target: Dict[str,object], samples: int, seed: int, sigma_arcsec: float) -> Dict[str,object]:
    center = target["center"]
    nearest = target["nearest_catalog_detection"]
    if nearest is None:
        return {"target": name, "status": "I_DO_NOT_KNOW", "reason": "NO_NEAREST_CATALOG_DETECTION_IN_INPUT"}

    dims = (
        DimensionSpec("ra_offset_arcsec",0.0,sigma_arcsec),
        DimensionSpec("dec_offset_arcsec",0.0,sigma_arcsec),
    )
    cloud = JanusWorldMatrix(dims, samples=samples, seed=seed,
        namespace=f"LOVE_EDEM__{name}__POSITIONAL_ROBUSTNESS", sampling="latin_hypercube").generate()

    ra0=float(center["ra_deg"]); dec0=float(center["dec_deg"])
    src_ra=float(nearest["ra_deg"]); src_dec=float(nearest["dec_deg"])
    frozen_sep=float(nearest["separation_arcsec"])
    recomputed=angular_separation_arcsec(ra0,dec0,src_ra,src_dec)

    seps: List[float]=[]
    center_shifts: List[float]=[]
    for dra_as,ddec_as in cloud["world_matrix"]:
        pra=ra0+float(dra_as)/3600.0
        pdec=dec0+float(ddec_as)/3600.0
        seps.append(angular_separation_arcsec(pra,pdec,src_ra,src_dec))
        center_shifts.append(angular_separation_arcsec(ra0,dec0,pra,pdec))

    gate_counts={str(g):sum(1 for x in seps if x<=g) for g in GATES_ARCSEC}
    max_shift=max(center_shifts)
    triangle_lower_bound=max(0.0,frozen_sep-max_shift)
    return {
        "target": name,
        "status": "ROBUSTNESS_EVALUATED",
        "frozen_center": center,
        "nearest_catalog_detection": nearest,
        "frozen_nearest_separation_arcsec": frozen_sep,
        "recomputed_nearest_separation_arcsec": recomputed,
        "assumed_center_jitter": {
            "model": "INDEPENDENT_NORMAL_RA_DEC_OFFSETS",
            "sigma_arcsec_each_coordinate": sigma_arcsec,
            "sampling": "DETERMINISTIC_LATIN_HYPERCUBE",
            "seed": seed,
            "samples": samples,
            "measured_uncertainty": False,
            "interpretation": "STRESS_TEST_PRIOR_ONLY",
        },
        "sampled_nearest_separation_arcsec": _stats(seps),
        "sampled_center_displacement_arcsec": _stats(center_shifts),
        "max_sampled_center_displacement_arcsec": max_shift,
        "triangle_inequality_lower_bound_for_any_catalog_source_arcsec": triangle_lower_bound,
        "gate_counts": gate_counts,
        "gate_fractions": {k:v/samples for k,v in gate_counts.items()},
        "all_catalog_sources_excluded_from_gate_by_triangle_bound": {str(g): bool(triangle_lower_bound>g) for g in GATES_ARCSEC},
        "exact_counterpart_promoted": False,
    }


def run(summary: Dict[str,object], samples: int=DEFAULT_SAMPLES, seed: int=DEFAULT_SEED,
        sigma_arcsec: float=DEFAULT_SIGMA_ARCSEC) -> Dict[str,object]:
    targets=summary["targets"]
    result_targets={}
    for offset,name in enumerate(("LOVE","EDEM_SEARCH_CENTER_ZP")):
        result_targets[name]=run_target(name,targets[name],samples,seed+offset,sigma_arcsec)
    core_input={
        "input_experiment_id":summary.get("experiment_id"),
        "input_claim_ceiling":summary.get("claim_ceiling"),
        "targets":{name:{"center":targets[name]["center"],"nearest_catalog_detection":targets[name]["nearest_catalog_detection"]}
                   for name in ("LOVE","EDEM_SEARCH_CENTER_ZP")},
    }
    out={
        "schema":SCHEMA,
        "experiment":"LOVE_EDEM_WORLD_MATRIX_ROBUSTNESS",
        "formula":"RESPICIENS_ET_PROSPICIENS",
        "input_experiment_id":summary.get("experiment_id"),
        "input_claim_ceiling":summary.get("claim_ceiling"),
        "core_input_sha256":_sha256_obj(core_input),
        "samples_per_target":samples,
        "sigma_arcsec_each_coordinate":sigma_arcsec,
        "targets":result_targets,
        "epistemic_firewall":{
            "simulation_is_evidence":False,
            "simulation_is_discovery":False,
            "stress_test_prior_is_measured_uncertainty":False,
            "negative_result_is_valid":True,
            "i_do_not_know_is_valid":True,
        },
        "claim_ceiling":"ROBUSTNESS_TO_ASSUMED_SUB_ARCSECOND_CENTER_JITTER_ONLY__NOT_ANOMALY_EVIDENCE",
    }
    out["input_sha256"]=_sha256_obj(summary)
    out["result_sha256"]=_sha256_obj(out)
    return out


def self_test() -> None:
    sample={
        "experiment_id":"LOVE-EDEM-CATALOG-CROSSCHECK-v1",
        "claim_ceiling":"CATALOG_SOURCE_INVENTORY_AROUND_FROZEN_DIRECTIONS_ONLY",
        "targets":{
            "LOVE":{"center":{"ra_deg":204.30267916666668,"dec_deg":-36.78240527777778,"frame":"ICRS"},
                "nearest_catalog_detection":{"ra_deg":204.2960668,"dec_deg":-36.7804774,"separation_arcsec":20.289587649207107,"catalog":"ALLWISE","source_id":"J133711.05-364649.7"}},
            "EDEM_SEARCH_CENTER_ZP":{"center":{"ra_deg":139.22409686590188,"dec_deg":30.26038779947318,"frame":"ICRS"},
                "nearest_catalog_detection":{"ra_deg":139.2249571,"dec_deg":30.262781,"separation_arcsec":9.021198661597598,"catalog":"ALLWISE","source_id":"J091653.98+301546.0"}}
        }
    }
    a=run(sample,samples=512,seed=16016,sigma_arcsec=.7)
    b=run(sample,samples=512,seed=16016,sigma_arcsec=.7)
    assert a["result_sha256"]==b["result_sha256"]
    for name in ("LOVE","EDEM_SEARCH_CENTER_ZP"):
        t=a["targets"][name]
        assert t["gate_counts"]["3.0"]==0
        assert t["all_catalog_sources_excluded_from_gate_by_triangle_bound"]["3.0"] is True
    print("LOVE_EDEM_WORLD_MATRIX_SELF_TEST=PASS")


def main() -> int:
    p=argparse.ArgumentParser()
    p.add_argument("--input",default="data/love/LOVE-EDEM-CATALOG-CROSSCHECK-v1-LATEST-SUMMARY.json")
    p.add_argument("--output",default=None)
    p.add_argument("--samples",type=int,default=DEFAULT_SAMPLES)
    p.add_argument("--seed",type=int,default=DEFAULT_SEED)
    p.add_argument("--sigma-arcsec",type=float,default=DEFAULT_SIGMA_ARCSEC)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        self_test(); return 0
    summary=json.loads(Path(a.input).read_text(encoding="utf-8"))
    out=run(summary,a.samples,a.seed,a.sigma_arcsec)
    text=json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)
    if a.output: Path(a.output).write_text(text+"\n",encoding="utf-8")
    else: print(text)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
