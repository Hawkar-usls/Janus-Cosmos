#!/usr/bin/env python3
"""JANUS R0M: deterministic mirrored numerical variants, sky directions and controls."""
from pathlib import Path
import concurrent.futures
import gzip
import hashlib
import json
import math
import random
import re
import sys
import aggregate_verify as audit

P=json.loads(Path("research/numeric_sky/mirror_protocol.json").read_text())
SRC=Path("runtime/numeric_sky/collected")
OUT=Path("runtime/numeric_sky/mirror")
OUT.mkdir(parents=True,exist_ok=True)
COUNT=8

def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,allow_nan=False,indent=2)+"\n")
def float_ok(v):
    return isinstance(v,(float,int)) and not isinstance(v,bool) and math.isfinite(v)
def simple_decimal(v):
    if not float_ok(v) or abs(v)>1e12:return None
    if float(v).is_integer():return str(int(v))
    s=repr(v)
    return s if "e" not in s.lower() and re.fullmatch(r"-?\d+\.\d+",s) else None
def digit_reversal(v):
    s=simple_decimal(v)
    if s is None:return None
    neg=s.startswith("-")
    d=s[1:] if neg else s
    # Reverse the decimal-character sequence; decimal point is kept but reversed in position.
    t=d[::-1]
    if not re.fullmatch(r"\d+(?:\.\d+)?",t):return None
    try:
        val=float(("-" if neg else "")+t)
    except (ValueError,OverflowError):return None
    if not math.isfinite(val):return None
    return {"decimal":("-" if neg else "")+t,"value":val}
def wrap_deg(ra):return (ra%360.0+360.0)%360.0
def transform_sky(ra,dec,operator):
    if operator=="original":return [wrap_deg(ra),dec]
    if operator=="ra_meridian_reflection":return [wrap_deg(-ra),dec]
    if operator=="equatorial_reflection":return [wrap_deg(ra),-dec]
    if operator=="antipode":return [wrap_deg(ra+180),-dec]
    raise ValueError(operator)
SKY_OPS=["original","ra_meridian_reflection","equatorial_reflection","antipode"]
NUM_OPS=["original","minus_a","minus_b","reverse_both_base10","swap_ab"]
def project(a,b):
    # Keep the exact fixed, arbitrary sky projection from R0.
    return [wrap_deg(a), math.degrees(math.asin(max(-1,min(1,2*((abs(b)/180)%1)-1))))]
def num_variant(a,b,method):
    if method=="original":return a,b
    if method=="minus_a":return -a,b
    if method=="minus_b":return a,-b
    if method=="reverse_both_base10":
        r1=digit_reversal(a);r2=digit_reversal(b)
        return (r1["value"],r2["value"]) if r1 and r2 else None
    if method=="swap_ab":return b,a
    raise ValueError(method)

def stream_numeric_mirror():
    counts={"originals":0,"numeric_negatives":0,"reversed_digits":0,"unrepresentable_reversals":0,
            "palindromic_or_equal_digit_reversals":0,"source_shards":8}
    with gzip.open(OUT/"mirrored_numeric_all.jsonl.gz","wt",encoding="utf-8",compresslevel=6) as sink:
        for i in range(COUNT):
            source=SRC/f"numsky-shard-{i}"/"numbers.jsonl.gz"
            assert source.is_file(),str(source)
            with gzip.open(source,"rt",encoding="utf-8") as sourcefile:
                for line in sourcefile:
                    x=json.loads(line);v=x["number"]
                    if not float_ok(v):continue
                    m=digit_reversal(v)
                    datum={"origin":x,"mirrors":{"negated":{"value":-v,"operator":"NEGATE"},
                       "digit_reversal":m},"physical_coordinate_claim":False}
                    sink.write(json.dumps(datum,ensure_ascii=False,separators=(",",":"),allow_nan=False)+"\n")
                    counts["originals"]+=1;counts["numeric_negatives"]+=1
                    if m:
                        counts["reversed_digits"]+=1
                        if m["value"]==v:counts["palindromic_or_equal_digit_reversals"]+=1
                    else:counts["unrepresentable_reversals"]+=1
    save("numeric_mirror_counts.json",counts)
    print("NUMSKY_MIRROR_NUMERIC",json.dumps(counts),flush=True)
    return counts
def sky_mirror(cands):
    points=[]
    for c in cands:
        for op in SKY_OPS:
            ra,dec=transform_sky(c["ra_deg"],c["dec_deg"],op)
            points.append({"ra_deg":ra,"dec_deg":dec,"source_ra_deg":c["ra_deg"],"source_dec_deg":c["dec_deg"],
                "operator":op,"repo":c["repo"],"path":c["path"],"revision":c["revision"],"origin":c["origin"],
                "source_sha256":c["source_sha256"],"original_kind":c["kind"],
                "scientific_qualified_source":c["scientific_qualified"],
                "mirrored_direction_is_new_observation":False})
    save("mirrored_explicit_icrs_directions.json",{"source_coordinates":len(cands),
       "generated_directions":len(points),"mirror_operators":SKY_OPS,"directions":points})
    return points
def projected_variants(projections):
    allpoints=[]
    for p in projections:
        for op in NUM_OPS:
            z=num_variant(p["a"],p["b"],op)
            if z is None:
                allpoints.append({"source":p,"operator":op,"valid":False,"reason":"DECIMAL_REVERSE_UNREPRESENTABLE"})
                continue
            a,b=z;ra,dec=project(a,b)
            allpoints.append({"ra_deg":ra,"dec_deg":dec,"source":p,"operator":op,
                              "valid":True,"a":a,"b":b,"scientific_qualified":False,
                              "physical_direction":False,"note":"Arbitrary synthetic projection, not an observation"})
    save("mirrored_synthetic_math_directions.json",{"sources":len(projections),
           "operators":NUM_OPS,"generated":len(allpoints),"directions":allpoints})
    return allpoints
def run_catalog(candidates,projections):
    explicit=[c for c in candidates if c["scientific_qualified"]][:P["catalog_gate"]["sample_explicit"]]
    synthetic=projections[:P["catalog_gate"]["sample_synthetic"]]
    rnd=random.Random(20261121)
    base_random=[(rnd.random()*360, math.degrees(math.asin(2*rnd.random()-1))) for _ in range(len(explicit))]
    synth_random=[(rnd.random()*360, math.degrees(math.asin(2*rnd.random()-1))) for _ in range(len(synthetic))]
    queries=[]
    for p,q in zip(explicit,base_random):
        for op in SKY_OPS:
            ra,dec=transform_sky(p["ra_deg"],p["dec_deg"],op)
            queries.append({"arm":"explicit_"+op,"ra_deg":ra,"dec_deg":dec,"source":p["path"],"operator":op})
            ra,dec=transform_sky(*q,op)
            queries.append({"arm":"control_explicit_"+op,"ra_deg":ra,"dec_deg":dec,"source":None,"operator":op})
    for p,q in zip(synthetic,synth_random):
        for op in NUM_OPS:
            pair=num_variant(p["a"],p["b"],op)
            # Still count missing decimal mirror explicitly rather than silently omit.
            if pair is not None:
                ra,dec=project(*pair)
                queries.append({"arm":"synthetic_"+op,"ra_deg":ra,"dec_deg":dec,"source":p["path"],"operator":op})
                qr,qd=q
                # Controls by fixed sky design; extra mirrored projection is NOT an a priori physical model.
                qr,qd=transform_sky(qr,qd,SKY_OPS[NUM_OPS.index(op)%4])
                queries.append({"arm":"control_synthetic_"+op,"ra_deg":qr,"dec_deg":qd,"source":None,"operator":op})
    def fetch_one(q):
        result=audit.simbad_probe(q)
        return {**q,**result}
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results=list(pool.map(fetch_one,queries))
    arms={}
    for x in results:
        arm=x["arm"]
        a=arms.setdefault(arm,{"queries":0,"available":0,"nearby":0,"unavailable":0})
        a["queries"]+=1
        if x["status"]=="OK":
            a["available"]+=1;a["nearby"]+=bool(x["object_count"])
        else:a["unavailable"]+=1
    result={"source":"SIMBAD TAP Harvard/Strasbourg mirrors","radius_arcsec":30,
            "catalog_presence_only":True,"physical_signal_detected":False,
            "preregistered_operator_comparison":True,"numerical_trials_counted":len(queries),
            "no_significance_claim_due_to_posthoc_selection_and_sky_density":True,
            "arms":arms,"records":results}
    save("mirror_catalog_check.json",result)
    print("NUMSKY_MIRROR_CATALOG",json.dumps({"queries":len(queries),"arms":arms}),flush=True)
    return result

def main():
    sources,candidates,shards=audit.collect()
    uniq=audit.unique_candidates(candidates)
    synth,meta=audit.math_projection()
    counts=stream_numeric_mirror()
    sky=sky_mirror(uniq)
    vari=projected_variants(synth)
    checks=run_catalog(uniq,synth)
    summary={"stage":"JANUS-NUMERIC-SKY-AUDIT-R0M","status":"COMPLETE_MIRROR_GENERATION__CATALOG_CHECK_BOUNDED",
             "source_run":P["source_run"],"source_shards":COUNT,
             "source_repositories":len(sources),"source_code_numbers":counts["originals"],
             "numeric_mirror_rows":counts["originals"],"unique_original_coordinate_pairs":len(uniq),
             "explicit_original_plus_3_mirrors":len(sky),
             "math_projected_original_plus_4_variants":len(vari),
             "mirrored_sky_catalog_queries":checks["numerical_trials_counted"],
             "catalog_arm_counts":checks["arms"],"production_unchanged":True,
             "artificial_object_detected":False,
             "interpretation":"No evidence for probes or extraterrestrial signals from mirroring. Null and multiplicity controls mandatory."}
    save("mirror_summary.json",summary)
    (OUT/"README.md").write_text(
      "# JANUS Numeric Sky Mirror R0M\n\n"
      "This is a deterministic exploratory transform of the frozen R0 public source corpus. "
      "Every original number has an accompanying negated value and, when representable, a "
      "reversed base-10 decimal representation. Original rows retain source paths and SHA-256.\n\n"
      "Every source-declared celestial direction has the original plus RA-meridian, "
      "equatorial, and antipodal reflections. These are mathematical images, NOT observed "
      "sky positions or candidate extraterrestrial signals.\n\n"
      "Selected arbitrary number-pair projections also receive sign/reversal/swap variants. "
      "All are hypothesis-free and must not be treated as decoded messages.\n\n"
      "SIMBAD lookup with uniform-sky controls is only a catalog proximity test, not a discovery "
      "or a matched-sky statistical test. All potential catalog failures and multiplicity are disclosed.\n")
    print("NUMSKY_MIRROR_FINAL",json.dumps({k:v for k,v in summary.items() if k!="catalog_arm_counts"}),flush=True)
if __name__=="__main__":main()
