#!/usr/bin/env python3
"""Aggregate 8 independent numeric scans; audit ICRS sky positions against external SIMBAD TAP."""
import concurrent.futures
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import re
import requests
import statistics

PROTO=json.loads(Path("research/numeric_sky/audit_protocol.json").read_text())
ROOT=Path("runtime/numeric_sky/collected")
OUT=Path("runtime/numeric_sky/final")
OUT.mkdir(parents=True,exist_ok=True)
NEAR=30.0
SEED=20261121
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def matchscore(c):
    path=c["path"].lower();place=c["origin"].lower();repo=c["repo"].lower()
    s=4 if c.get("scientific_qualified") else -6
    if c["kind"] in ("python_named_constants","python_skycoord"):s+=5
    if c["kind"]=="icrs_unit_vector":s+=4
    if "janus-cosmos" in repo:s+=4
    if "meta-registry" in repo:s+=2
    if any(v in path+" "+place for v in ("osiris","abydos","orion","love","edem","sirius","beacon")):s+=9
    if any(v in path+" "+place for v in ("frozen","target","center","position")):s+=5
    if any(v in path+" "+place for v in ("catalog_results","detections","nearest_20","catalog_crosscheck","candidate_predictions")):s-=15
    if any(v in path for v in ("test_","tests/","archive/","outbox/","receipt-old")):s-=4
    return s
def coordkey(c):return (round(c["ra_deg"],6),round(c["dec_deg"],6))
def coord_valid(c):
    return isinstance(c.get("ra_deg"),(int,float)) and isinstance(c.get("dec_deg"),(int,float)) and math.isfinite(c["ra_deg"]) and math.isfinite(c["dec_deg"]) and 0<=c["ra_deg"]<360 and -90<=c["dec_deg"]<=90
def collect():
    entries=[];cands=[];shards=[]
    for i in range(8):
        p=ROOT/f"numsky-shard-{i}"
        inv=json.loads((p/"inventory.json").read_text())
        co=json.loads((p/"explicit_icrs_candidates.json").read_text())
        shards.append(json.loads((p/"audit_summary.json").read_text()))
        assert len(inv["repos"])==len([r for j,r in enumerate(PROTO["provenance"]["public_repositories"]) if j%8==i])
        entries.extend(inv["repos"])
        cands.extend(co["candidates"])
    names={r["repo"] for r in entries}
    assert len(entries)==len(names)==len(PROTO["provenance"]["public_repositories"])
    for r in entries:
        assert r["repo"] in PROTO["provenance"]["public_repositories"]
    return entries,[x for x in cands if coord_valid(x)],shards
def unique_candidates(cands):
    grouped={}
    for c in cands:
        key=coordkey(c)
        grouped.setdefault(key,[]).append(c)
    result=[]
    for key,src in grouped.items():
        ss=sorted(src,key=lambda c:(-matchscore(c),c["repo"],c["path"],c["origin"]))
        best=ss[0].copy()
        best["references"]=len(src)
        best["duplicate_sources"]=[{"repo":s["repo"],"path":s["path"],"origin":s["origin"]} for s in ss[:10]]
        best["priority_score"]=matchscore(best)
        result.append(best)
    result.sort(key=lambda c:(-c["priority_score"],c["repo"],c["path"],c["origin"]))
    return result
def math_projection():
    """A fixed, non-physical projection from adjacent math literals to ICRS.
    This deliberately requires null controls. No physical association is inferred.
    """
    arr=[]
    pat=re.compile(r"(?i)\b(phi|golden|fibonacc?i|prime|euler|pi|tau|sin|cos|sqrt|taylor|harmonic|period|factorial|newton|modular|mersenne|zeta|math)\b")
    for i in range(8):
        file=ROOT/f"numsky-shard-{i}"/"numbers.jsonl.gz"
        previous={}
        with gzip.open(file,"rt",encoding="utf-8") as f:
            for ln in f:
                try:x=json.loads(ln)
                except json.JSONDecodeError:continue
                if x["kind"] not in ("python_number","json_number","lexical_number"):continue
                n=x["number"]
                if not isinstance(n,(int,float)) or not math.isfinite(n) or abs(n)<1e-9 or abs(n)>1e8:continue
                if not pat.search((x.get("context") or "")+" "+x.get("field","")):continue
                ky=(x["repo"],x["path"],x["source_sha256"])
                p=previous.get(ky)
                previous[ky]=x
                if not p or p["number"]==x["number"]:continue
                a,b=p["number"],n
                ra=a%360
                frac=(abs(b)/180)%1
                dec=math.degrees(math.asin(max(-1,min(1,2*frac-1))))
                arr.append({"ra_deg":ra,"dec_deg":dec,"repo":x["repo"],"path":x["path"],"revision":x["revision"],
                    "a":a,"b":b,"origin":"consecutive_math_literals","method":"ra=a mod 360; dec=asin(2*frac(|b|/180)-1) deg",
                    "physical_direction":False,"eligible_as_discovery":False})
                if len(arr)>100000:break
    arr.sort(key=lambda x:hashlib.sha256((str(SEED)+"|"+x["repo"]+"|"+x["path"]+"|"+str(x["a"])+"|"+str(x["b"])).encode()).hexdigest())
    return arr[:64],{"total_generated_before_selection":len(arr),"selected":min(len(arr),64),"selection":"SHA256 sort seeded by 20261121, before any catalog query"}
def adql_query(ra,dec):
    rad=NEAR/3600
    return ("SELECT TOP 3 main_id, ra, dec, otype, "
            "DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',"+str(ra)+","+str(dec)+")) AS sep_deg "
            "FROM basic WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
            "CIRCLE('ICRS',"+str(ra)+","+str(dec)+","+str(rad)+")) ORDER BY sep_deg ASC")
def simbad_probe(t):
    base=["https://simbad.cfa.harvard.edu/simbad/sim-tap/sync","https://simbad.cds.unistra.fr/simbad/sim-tap/sync"]
    errors=[]
    for url in base:
        try:
            response=requests.get(url,params={"request":"doQuery","lang":"adql","format":"json","query":adql_query(t["ra_deg"],t["dec_deg"])},
                                  headers={"User-Agent":"JANUS-Research-Numeric-Sky-Audit/0.1 (read-only)"},timeout=(5,9))
            response.raise_for_status()
            dat=response.json()
            rows=dat.get("data")
            if rows is None:raise ValueError("No data field in TAP JSON")
            # TAP JSON generally encodes each row as an array. Keep the service's own metadata too.
            meta=dat.get("metadata",[])
            cols=[m.get("name","") for m in meta if isinstance(m,dict)]
            normal=[]
            for row in rows[:3]:
                d={k:v for k,v in zip(cols,row)} if isinstance(row,list) else row if isinstance(row,dict) else {"raw":row}
                normal.append(d)
            return {"status":"OK","service":url,"matches":normal,"object_count":len(rows),"errors":errors}
        except Exception as e:errors.append(str(e)[:160])
    return {"status":"UNAVAILABLE","service":None,"matches":[],"object_count":None,"errors":errors}
def catalog_scan(explicit,synthetic):
    anchors=[x for x in explicit if x["scientific_qualified"]][:24]
    samples=synthetic[:8]
    rng=random.Random(SEED)
    controls=[{"ra_deg":rng.random()*360,"dec_deg":math.degrees(math.asin(2*rng.random()-1)),
               "kind":"uniform_random_control","physical_direction":False} for _ in range(len(anchors)+len(samples))]
    todo=[("explicit",x) for x in anchors]+[("synthetic",x) for x in samples]+[("control",x) for x in controls]
    def work(item):
        group,c=item
        v=simbad_probe(c)
        return {"kind":group,"ra_deg":c["ra_deg"],"dec_deg":c["dec_deg"],"origin":c.get("origin"),
                "path":c.get("path"),"repo":c.get("repo"),"provenance_label":c.get("kind"),
                "query_radius_arcsec":NEAR,**v}
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results=list(pool.map(work,todo))
    available=[x for x in results if x["status"]=="OK"]
    counts={k:{"queried":sum(x["kind"]==k for x in results),"available":sum(x["kind"]==k for x in available),
               "nearby_catalog_objects":sum(x["kind"]==k and x["object_count"]>0 for x in available)} for k in ("explicit","synthetic","control")}
    return {"catalog":"SIMBAD TAP (CDS/Harvard mirror, service-dependent object inventory)",
            "radius_arcsec":NEAR,"null_model":"uniform ICRS sphere with seed 20261121; not matched to Galactic latitude or local catalog completeness",
            "selection":"Top 24 pre-ranked frozen semantically explicit directions + 8 deterministic speculative numerical projections; controls matched in count and radius.",
            "counts":counts,"records":results,
            "scientific_verdict":"CATALOG_CROSSMATCH_ONLY_NOT_DISCOVERY__NO_TECHNOSIGNATURE_OR_PROBE_PROVEN",
            "warning":"If SIMBAD unavailable, no crossmatch conclusion. A catalog neighbor does not identify a signal, anomaly, or artificial source. 30 arcsec is a broad illustrative cone, not astrometric identity."}
def main():
    inv,cands,shards=collect()
    uniq=unique_candidates(cands)
    syn,stat=math_projection()
    check=catalog_scan(uniq,syn)
    totals={"owner":"Hawkar-usls","public_repositories_preregistered":len(PROTO["provenance"]["public_repositories"]),
        "public_repositories_scanned":len(inv),"repos_complete":sum(r["status"]=="COMPLETE_SELECTED_CODE" for r in inv),
        "repos_partial":sum(r["status"]=="PARTIAL" for r in inv),"repos_failed":sum(r["status"]=="FAILED" for r in inv),
        "source_files_selected":sum(r["selected_files"] for r in inv),"numbers_extracted":sum(r["numbers"] for r in inv),
        "labelled_coordinate_occurrences":len(cands),"unique_coordinate_pairs_rounded_microdegree":len(uniq),
        "qualified_explicit_directions":sum(x["scientific_qualified"] for x in uniq),
        "private_unscanned":PROTO["provenance"]["not_scannable_private"],
        "skipped_counts":{reason:sum(r["skipped"].get(reason,0) for r in inv) for reason in inv[0]["skipped"]},
        "revisions":{r["repo"]:r["head_sha"] for r in inv},
        "methodical_limit":"Single default-branch snapshot per repo, source extension whitelist, 1 MiB/file, no vendored/generated, no historical branches or private repos"}
    status="PARTIAL_WITH_REPO_ERRORS" if totals["repos_failed"] or totals["repos_partial"] else "PASS_LITERAL_EXTRACTION_WITHOUT_DISCOVERY"
    result={"stage":"JANUS-NUMERIC-SKY-AUDIT-R0","status":status,"totals":totals,"shards":shards,
            "external_catalog":check["counts"],"projection":stat,
            "science_status":"NO_INDEPENDENT_PROOF_OF_TECHNOSIGNATURE_OR_PROBE",
            "production_change":False}
    save("audit_summary.json",result)
    save("repo_inventory.json",{"repositories":inv})
    save("explicit_icrs_candidates.json",{"total_unique":len(uniq),"candidates":uniq})
    save("synthetic_projections.json",{"hypothesis":"ARBITRARY_DETERMINISTIC_PROJECTION_WITH_NULL_CONTROLS","metadata":stat,"points":syn})
    save("external_catalog_checks.json",check)
    readme=["# JANUS Numeric Sky Audit — research snapshot",
            "Scope: "+str(totals["public_repositories_scanned"])+" public repo default branches; 3 private repos are unscanned.",
            "Numerical literals: "+str(totals["numbers_extracted"])+"; source files: "+str(totals["source_files_selected"])+".",
            "Explicit coordinate pairs (deduplicated): "+str(len(uniq))+".",
            "The 64 exploratory sky projections are ARBITRARY mappings from numeric constants, not astronomical measurements.",
            "The external SIMBAD crosscheck (if available) establishes nearby catalog entries, not unusual physics.",
            "No complete all-git-history guarantee; no instrumental evidence for artificial probes.",
            "All raw numerical sources are saved in compressed numbers.jsonl.gz within each shard artifact."]
    (OUT/"README.md").write_text("\n\n".join(readme)+"\n")
    print("NUMSKY_FINAL",json.dumps({k:v for k,v in result.items() if k!="shards"}),flush=True)
if __name__=="__main__":main()
