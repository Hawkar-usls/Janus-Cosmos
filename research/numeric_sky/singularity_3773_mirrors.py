#!/usr/bin/env python3
"""Зеркало Януса: mirror the frozen SHA-derived 3773rd point into 3774–3776.
Real catalog lookups are read-only. Mirroring alone cannot predict physical objects.
"""
import concurrent.futures
import hashlib
import json
import math
import os
from pathlib import Path
import requests
import collapse_3773 as core
import direct_gaia_irsa_3773 as multi

SPEC=json.loads(Path("research/numeric_sky/singularity_3773_mirrors_protocol.json").read_text())
O=Path("runtime/numeric_sky/singularity_mirrors")
O.mkdir(parents=True,exist_ok=True)
M=SPEC["source"]
H={"User-Agent":"JANUS-COSMOS-Zerkalo-3773-mirrors/1.0 (read-only)"}
def save(name,x):
    (O/name).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
def close(a,b):return abs(a-b)<1e-10
def mirror(ra,dec,op):
    if op=="original":return ra,dec
    if op=="ra_meridian_reflection":return (360-ra)%360,dec
    if op=="equatorial_reflection":return ra,-dec
    if op=="antipode":return (ra+180)%360,-dec
    raise AssertionError("Unknown mirror")
def load_pinned_meta():
    u="https://raw.githubusercontent.com/"+M["meta_repository"]+"/"+M["meta_commit"]+"/"+M["meta_path"]
    r=requests.get(u,headers=H,timeout=(8,60));r.raise_for_status()
    p=json.loads(r.content)
    assert p["source"]["fingerprint_sha256"]==M["base_fingerprint_sha256"]
    assert len(p["points"])==4 and [x["ordinal"] for x in p["points"]]==[3773,3774,3775,3776]
    return p
def verify_source():
    root=core.compute()  # recompute all 3772 source atlas directions; fail on source SHA mismatch.
    assert root["fingerprint_sha256"]==M["base_fingerprint_sha256"]
    assert root["sources"]==3772 and root["direction"]["ordinal"]==3773
    base=root["direction"]
    raw=load_pinned_meta()
    out=[]
    for p in raw["points"]:
        ra,dec=mirror(base["ra_deg"],base["dec_deg"],p["operator"])
        assert close(ra,p["ra_deg"]) and close(dec,p["dec_deg"])
        assert 0<=ra<360 and -90<=dec<=90
        out.append({**p,"frame":"ICRS","source_fingerprint_sha256":root["fingerprint_sha256"],
            "is_new_observation":False,"is_discovery_claim":False})
    payload={"schema":"janus.cosmos.zerkalo_janusa.four_directions.v1",
        "state":"FOUR_REPRODUCIBLE_MATH_DIRECTIONS", "source_fingerprint":root["fingerprint_sha256"],
        "source_original_declared":943,"source_mirror_records":3772,
        "points":out,"geometry":"Three reflections of SHA-derived pseudorandom point 3773",
        "celestial_catalog_checks_separate":True,"not_a_physical_source_detection":True}
    save("3773_3776_directions.json",payload)
    print("JANUS_MIRROR_FOUR",json.dumps({"points":out},ensure_ascii=False),flush=True)
    return out
def gaia_adql(p):
    ra=p["ra_deg"];dec=p["dec_deg"];rad=60/3600
    return ("SELECT TOP 100 source_id, ra, dec, phot_g_mean_mag, parallax, pmra, pmdec "
         "FROM gaiadr3.gaia_source "
         f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{rad}))")
def wise_adql(p):
    ra=p["ra_deg"];dec=p["dec_deg"];rad=60/3600
    return ("SELECT TOP 100 designation, ra, dec, w1mpro, w2mpro, w3mpro, w4mpro "
         "FROM allwise_p3as_psd "
         f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{rad}))")
def catalog_for_target(p):
    # SIMBAD same immutable query implementation as original ordinal 3773 experiment.
    sim=core.simbad(p)
    gaia=multi.check("Gaia DR3 ESA",p,
        ["https://gea.esac.esa.int/tap-server/tap/sync","https://gaia.ari.uni-heidelberg.de/tap/sync"],
        gaia_adql(p),"source_id",["phot_g_mean_mag","parallax","pmra","pmdec"])
    wise=multi.check("AllWISE NASA/IPAC",p,
        ["https://irsa.ipac.caltech.edu/TAP/sync"],wise_adql(p),"designation",
        ["w1mpro","w2mpro","w3mpro","w4mpro"])
    x={"ordinal":p["ordinal"],"operator":p["operator"],"ra_deg":p["ra_deg"],"dec_deg":p["dec_deg"],
       "simbad":sim,"gaia":gaia,"allwise":wise,
       "identification_status":"CATALOG_NEIGHBOR_ONLY_NO_UNIQUE_ASTROPHYSICAL_ASSOCIATION"}
    print("JANUS_MIRROR_TARGET",json.dumps({
        "ordinal":p["ordinal"],"operator":p["operator"],
        "simbad":{"status":sim["status"],"radii":sim.get("num_within_arcsec"),
                  "nearest":sim.get("catalog_entries_top40",[])[:3]},
        "gaia":{"status":gaia["status"],"nearest":gaia.get("nearest",[])[:2]},
        "allwise":{"status":wise["status"],"nearest":wise.get("nearest",[])[:2]}
    },ensure_ascii=False),flush=True)
    return x
def control_audit(points):
    # Controls fixed from previously frozen SHA256 (first 64 out of 256), transformed identically.
    base={"fingerprint_sha256":M["base_fingerprint_sha256"]}
    controls=core.simulate_random_controls(base,n=256)
    selected=controls[:SPEC["catalog"]["random_controls_per_new_mirror"]]
    todo=[]
    for p in points[1:]:
        op=p["operator"]
        for q in selected:
            ra,dec=mirror(q["ra_deg"],q["dec_deg"],op)
            todo.append({"target":p["ordinal"],"operator":op,"control_index":q["control_index"],"ra_deg":ra,"dec_deg":dec})
    assert len(todo)==3*64
    def check(c):
        sim=core.simbad(c)
        return {"target":c["target"],"control_index":c["control_index"],
          "ra_deg":c["ra_deg"],"dec_deg":c["dec_deg"],
          "status":sim["status"],"within_arcsec":sim.get("num_within_arcsec"),
          "service":sim.get("service"),"errors":sim.get("errors",[])}
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        res=list(ex.map(check,todo))
    summary={}
    for target in [3774,3775,3776]:
        subset=[v for v in res if v["target"]==target]
        good=[v for v in subset if v["status"]=="OK"]
        part={"requested":len(subset),"available":len(good),"unavailable":len(subset)-len(good)}
        for arcsec in (2,10,30,60):
            part["with_source_within_"+str(arcsec)+"arcsec"]=sum(
                int(v["within_arcsec"][str(arcsec)]>0) for v in good)
        summary[str(target)]=part
    save("3774_3776_random_control_results.json",{
        "source_seed":"original source SHA256 fingerprint; first 64 of 256 uniform ICRS random controls, then same spatial mirror operator",
        "number_of_controls":len(res),"results":res,"grouped":summary,
        "limitations":"Isotropic random controls, not matched by Galactic latitude or instrument footprint. Mirror groups share correlated base controls."})
    print("JANUS_MIRROR_CONTROLS",json.dumps(summary),flush=True)
    return summary
def main():
    points=verify_source()
    results=[catalog_for_target(p) for p in points[1:]]
    save("3774_3776_external_catalog.json",{
        "schema":"janus.cosmos.zerkalo_janusa.3774_3776.catalog.v1",
        "coordinate_source_sha256":M["base_fingerprint_sha256"],
        "mirror_points":results,
        "individual_catalog_errors_are_not_absences":True,
        "technosignature_or_probe_confirmed":False,
        "telescope_pointing_or_original_astro_measurement":False})
    controls=control_audit(points)
    summary={"schema":"janus.cosmos.zerkalo_janusa.3773_3776.summary.v1",
       "point_count":4,"new_mirrored_points":3,"points":points,
       "catalog_results":[{"ordinal":v["ordinal"],"simbad_status":v["simbad"]["status"],
          "nearest_simbad":v["simbad"].get("catalog_entries_top40",[])[:3],
          "within_arcsec":v["simbad"].get("num_within_arcsec"),
          "gaia_status":v["gaia"]["status"],"gaia_nearest":v["gaia"].get("nearest",[])[:3],
          "allwise_status":v["allwise"]["status"],"allwise_nearest":v["allwise"].get("nearest",[])[:3]}
          for v in results],
       "controls":controls,
       "status":"CATALOG_SCREEN_COMPLETE_WITH_EXPLICIT_UNAVAILABLE_STATUSES",
       "scientific_claim":"NO_NOVEL_TECHNOSIGNATURE_OR_PROBE_FROM_COORDINATE_MIRRORING",
       "production_changed":False}
    save("3773_3776_result.json",summary)
    print("JANUS_MIRROR_SUMMARY",json.dumps(summary,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
