#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, random
from pathlib import Path

BINS=[
    ("INNER",0.0,3000.0),
    ("MID",3000.0,4500.0),
    ("OUTER",4500.0,5500.0),
    ("EXTREME",5500.0,float("inf")),
]

def load_library(path:Path):
    d=json.loads(path.read_text(encoding="utf-8"))
    assert d["claim_ceiling"]=="CD169_EMPIRICAL_MULTIBEAM_RESIDUAL_CALIBRATION_LIBRARY_ONLY__NO_UNSEEN_PREDICTIVE_SKILL"
    return d

def bin_for_across(abs_across_m:float)->str:
    x=abs(float(abs_across_m))
    for name,lo,hi in BINS:
        if lo<=x<hi:return name
    raise ValueError("unreachable across-track bin")

def channel(lib,lineage,reference,bin_name):
    try:c=lib["library"][lineage][reference][bin_name]
    except KeyError as e:raise ValueError(f"unknown calibration channel: {e}") from e
    if int(c.get("support_cell_count",0))<=0:
        raise ValueError("selected empirical calibration channel has no support")
    return c

def draw_one(rng,c,regime):
    signed=[float(x) for x in c["signed_residual_m"]]
    if regime=="CD169_REPLAY":
        return float(rng.choice(signed)),"SIGNED_EMPIRICAL_BOOTSTRAP"
    mags=[abs(float(x)) for x in signed]
    if regime=="GENERIC_UNSEEN_CALIBRATED":
        m=float(rng.choice(mags))
        sign=-1.0 if rng.random()<0.5 else 1.0
        return sign*m,"EMPIRICAL_MAGNITUDE_BOOTSTRAP_RANDOMIZED_SIGN"
    if regime=="ADVERSARIAL_REAL_TAILS":
        tail=[m for m in mags if m>=100.0]
        if not tail:
            raise ValueError("no empirically observed >=100 m residuals in selected channel/bin; choose another empirical channel rather than inventing a tail")
        m=float(rng.choice(tail))
        sign=-1.0 if rng.random()<0.5 else 1.0
        return sign*m,"EMPIRICAL_GE100_TAIL_STRESS_RANDOMIZED_SIGN"
    if regime=="IDEAL":
        return 0.0,"IDEAL_ZERO"
    raise ValueError("unknown regime")

def sample(lib,regime,lineage,reference,across_values,seed):
    rng=random.Random(seed);out=[]
    for i,x in enumerate(across_values):
        b=bin_for_across(x)
        if regime=="IDEAL":
            residual=0.0;src="IDEAL_ZERO";support=None
        else:
            c=channel(lib,lineage,reference,b)
            residual,src=draw_one(rng,c,regime)
            support=int(c["support_cell_count"])
        out.append({
            "index":i,
            "abs_across_track_m":abs(float(x)),
            "bin":b,
            "depth_error_m":residual,
            "sampling_source":src,
            "calibration_support_cells":support
        })
    return out

def parse_across(args):
    vals=[]
    if args.across_m:
        vals.extend(float(x) for x in args.across_m)
    if args.across_file:
        txt=args.across_file.read_text(encoding="utf-8")
        for token in txt.replace(","," ").split():vals.append(float(token))
    if not vals:raise ValueError("provide --across-m and/or --across-file")
    return vals

def main():
    ap=argparse.ArgumentParser(description="JANUS KUSTO empirical multibeam nuisance sampler")
    ap.add_argument("--library",required=True,type=Path)
    ap.add_argument("--regime",required=True,choices=["IDEAL","CD169_REPLAY","GENERIC_UNSEEN_CALIBRATED","ADVERSARIAL_REAL_TAILS"])
    ap.add_argument("--lineage",default="LINE32_ERA",choices=["LINE18_ERA","LINE32_ERA"])
    ap.add_argument("--reference",default="KN192-07_2008",choices=["KN192-07_2008","KNOX15RR_2008"])
    ap.add_argument("--across-m",type=float,nargs="*")
    ap.add_argument("--across-file",type=Path)
    ap.add_argument("--seed",type=int,default=20260923)
    ap.add_argument("--output",type=Path)
    a=ap.parse_args()
    lib=load_library(a.library)
    across=parse_across(a)
    samples=sample(lib,a.regime,a.lineage,a.reference,across,a.seed)
    payload={
      "artifact_type":"JANUS_KUSTO_SYNTHETIC_MULTIBEAM_NUISANCE_SAMPLES",
      "regime":a.regime,
      "calibration_lineage":a.lineage,
      "independent_reference":a.reference,
      "seed":a.seed,
      "n":len(samples),
      "samples":samples,
      "firewall":{
        "historical_runs_modified":False,
        "generic_unseen_signed_direction_inherited_from_cd169":False,
        "adversarial_tail_probability_forecast":False,
        "unseen_predictive_validation_claim":False
      },
      "claim_ceiling":"CALIBRATED_SYNTHETIC_MULTIBEAM_NUISANCE_ONLY__NOT_UNSEEN_VALIDATED"
    }
    s=json.dumps(payload,indent=2)
    if a.output:a.output.write_text(s+"\n",encoding="utf-8")
    else:print(s)

if __name__=="__main__":
    main()
