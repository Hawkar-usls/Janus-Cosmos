#!/usr/bin/env python3
"""JANUS COSMOS / Зеркало Януса — preregistered 3773rd sky coordinate and controls.
The SHA-256 fingerprint direction is not a physical inference or a technosignature claim.
"""
import concurrent.futures
import hashlib
import json
import math
import os
from pathlib import Path
import requests
import statistics

P=json.loads(Path("research/numeric_sky/singularity_3773_protocol.json").read_text())
O=Path("runtime/numeric_sky/singularity")
O.mkdir(parents=True,exist_ok=True)
S=P["source"]
H={"User-Agent":"JANUS-COSMOS-ZerkaloJanusa-SkyAudit/1.0 (catalog read-only)",
   "Accept":"application/json"}

def sha(x):return hashlib.sha256(x).hexdigest()
def save(name,obj):(O/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def get_source():
    src="https://raw.githubusercontent.com/"+S["repository"]+"/"+os.environ.get("META_COMMIT","main")+"/"+S["file"]
    r=requests.get(src,timeout=(10,120),headers={"User-Agent":H["User-Agent"]})
    r.raise_for_status()
    raw=r.content
    assert sha(raw)==S["source_sha256"],"FAIL_CLOSED_SOURCE_ATLAS_SHA256"
    data=json.loads(raw)
    xs=data["directions"]
    assert len(xs)==data["generated_directions"]==S["all_mirrored_count"]==3772
    assert data["source_coordinates"]==S["original_count"]==943
    ops=P["source"].get("expected_operators",["original","ra_meridian_reflection","equatorial_reflection","antipode"])
    from collections import Counter
    counts=Counter(v["operator"] for v in xs)
    assert len(counts)==4 and all(counts[q]==943 for q in ops)
    for x in xs:
        assert 0<=x["ra_deg"]<360 and -90<=x["dec_deg"]<=90 and x["operator"] in ops
    return xs,counts

def compute():
    xs,counts=get_source()
    rule=P["fingerprint_collapse"]
    lines=sorted(f'{x["ra_deg"]:.12f},{x["dec_deg"]:.12f},{x["operator"]}' for x in xs)
    digest=hashlib.sha256(rule["domain_separator"].encode("ascii")+("\n".join(lines)+"\n").encode("utf-8")).digest()
    ra=360*int.from_bytes(digest[:8],"big")/2**64
    dec=math.degrees(math.asin(2*int.from_bytes(digest[8:16],"big")/2**64-1))
    xx=math.fsum(math.cos(math.radians(x["dec_deg"]))*math.cos(math.radians(x["ra_deg"])) for x in xs)
    yy=math.fsum(math.cos(math.radians(x["dec_deg"]))*math.sin(math.radians(x["ra_deg"])) for x in xs)
    zz=math.fsum(math.sin(math.radians(x["dec_deg"])) for x in xs)
    magnitude=math.sqrt(xx*xx+yy*yy+zz*zz)
    geometry={"ra_deg":math.degrees(math.atan2(yy,xx))%360,"dec_deg":math.degrees(math.atan2(zz,math.hypot(xx,yy))),
      "sum_vector":[xx,yy,zz],"mean_resultant_fraction":magnitude/3772,
      "interpretation":"Mirror rule produces an equatorial meridian bias; the geometric centroid cannot be interpreted as an unknown physical sky position."}
    result={"schema":"janus.cosmos.zerkalo_janusa.singularity.v1","name_ru":"Зеркало Януса — Сингулярность №3773",
        "source_sha256":S["source_sha256"],"sources":3772,"originals":943,"operators":dict(counts),
        "canonical_line_count":len(lines),"fingerprint_sha256":digest.hex(),
        "direction":{"ra_deg":ra,"dec_deg":dec,"frame":"ICRS","ordinal":3773,
             "construction":"SHA256 coordinate fingerprint; pseudorandom uniform sky position, no astronomical source identity"},
        "geometric_result":geometry,
        "scientific_status":"DETERMINISTIC_POSITION_ONLY__NO_ASTRONOMICAL_DETECTION",
        "production_adoption":False}
    # Independent agreed fingerprint from an earlier numerical calculation.
    assert digest.hex()=="1c42dc53a9591f477cfba9678121c7ab331cc2afa2286b7d5187ef376fced080"
    assert abs(ra-39.742276528418486)<1e-12 and abs(dec+1.3505804808702793)<1e-12
    save("3773_coordinate.json",result)
    print("JANUS_3773_COORDINATE",json.dumps(result,ensure_ascii=False),flush=True)
    return result

def simulate_random_controls(result,n=256):
    # Uniform in sin(Dec), not uniform Dec. Test controls are deterministic once the source fingerprint is locked.
    key=P["fingerprint_collapse"]["domain_separator"].encode("ascii")+bytes.fromhex(result["fingerprint_sha256"])
    out=[]
    for i in range(n):
        h=hashlib.sha256(key+b":control:"+str(i).encode("ascii")).digest()
        ra=360*int.from_bytes(h[:8],"big")/2**64
        dec=math.degrees(math.asin(2*int.from_bytes(h[8:16],"big")/2**64-1))
        out.append({"control_index":i,"ra_deg":ra,"dec_deg":dec})
    save("preregistered_sky_controls.json",{"seed":"SHA256 of frozen coordinate fingerprint and index","n":n,"controls":out})
    return out

def adql_query(t,radius_arcsec=60.):
    ra=t["ra_deg"];dec=t["dec_deg"]
    deg=radius_arcsec/3600
    return ("SELECT TOP 40 main_id, ra, dec, otype, "
        "DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',"+str(ra)+","+str(dec)+")) AS sep_deg "
        "FROM basic WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
        "CIRCLE('ICRS',"+str(ra)+","+str(dec)+","+str(deg)+")) ORDER BY sep_deg")
def simbad(t):
    urls=["https://simbad.cfa.harvard.edu/simbad/sim-tap/sync",
      "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"]
    errors=[]
    for url in urls:
        try:
            r=requests.get(url,params={"request":"doQuery","lang":"adql","format":"json","query":adql_query(t)},
                           headers=H,timeout=(8,30))
            r.raise_for_status()
            data=r.json();rows=data.get("data")
            if rows is None:raise ValueError("SIMBAD TAP missing data")
            fields=[a.get("name","") for a in data.get("metadata",[])]
            parsed=[]
            for row in rows:
                x=dict(zip(fields,row)) if isinstance(row,list) else row
                s=float(x["sep_deg"])*3600
                if not 0<=s<=60.001:continue
                parsed.append({"name":str(x["main_id"]),"type":str(x.get("otype")),"ra_deg":float(x["ra"]),
                    "dec_deg":float(x["dec"]),"separation_arcsec":s})
            parsed.sort(key=lambda a:a["separation_arcsec"])
            # This is only a bounded catalogue inventory: TOP 40 can truncate dense fields.
            n={str(k):sum(x["separation_arcsec"]<=k for x in parsed) for k in (2,10,30,60)}
            return {"status":"OK","service":url,"catalog_entries_top40":parsed[:12],"num_within_arcsec":n,
                    "truncated_possible":len(rows)==40,"errors":errors}
        except Exception as e:errors.append(str(e)[:160])
    return {"status":"UNAVAILABLE","errors":errors,"num_within_arcsec":None,"catalog_entries_top40":[]}

def independent_multiband(ra,dec):
    # Separate catalog ecosystems: Gaia astrometry and AllWISE IR.
    result={}
    try:
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from astroquery.vizier import Vizier
        anchor=SkyCoord(ra=ra*u.deg,dec=dec*u.deg,frame="icrs")
        for label,catalog in (("GAIA_DR3","I/355/gaiadr3"),("ALLWISE","II/328/allwise")):
            try:
                v=Vizier(columns=["*"],row_limit=100)
                v.TIMEOUT=90
                tables=v.query_region(anchor,radius=60*u.arcsec,catalog=catalog)
                pts=[]
                for tb in tables:
                    cols=tb.colnames
                    ra_col=next((x for x in ("RA_ICRS","RAJ2000","RAdeg","RA_ICRS_","RA") if x in cols),None)
                    de_col=next((x for x in ("DE_ICRS","DEJ2000","DEdeg","DE_ICRS_","DEC","DE") if x in cols),None)
                    if ra_col is None or de_col is None:raise ValueError("missing catalog coordinates")
                    for row in tb:
                        try:
                            x=float(row[ra_col]);y=float(row[de_col])
                            sep=float(anchor.separation(SkyCoord(x*u.deg,y*u.deg,frame="icrs")).arcsecond)
                            ident=next((str(row[k]) for k in ("Source","AllWISE","WISE","_2MASS") if k in cols),"UNKNOWN")
                            pts.append({"ra_deg":x,"dec_deg":y,"separation_arcsec":sep,"catalog_id":ident})
                        except (ValueError,TypeError):continue
                pts.sort(key=lambda x:x["separation_arcsec"])
                result[label]={"status":"OK","catalog":catalog,"objects_within_60_arcsec":len(pts),
                    "nearest":pts[:8],"within_2":sum(x["separation_arcsec"]<=2 for x in pts),
                    "within_10":sum(x["separation_arcsec"]<=10 for x in pts),
                    "within_30":sum(x["separation_arcsec"]<=30 for x in pts)}
            except Exception as e:result[label]={"status":"UNAVAILABLE","error":str(e)[:240]}
    except Exception as e:
        result={k:{"status":"UNAVAILABLE","error":str(e)[:240]} for k in ("GAIA_DR3","ALLWISE")}
    return result

def catalog_verify(result,controls):
    target=result["direction"]
    geo=result["geometric_result"]
    events=[("singularity",target),("geometric_mirror_artifact",geo)]+[("random_control",x) for x in controls]
    def one(item):
        group,coords=item
        r=simbad(coords)
        return {"group":group,"coordinates":coords,**r}
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex: obs=list(ex.map(one,events))
    # Target is always first. An unavailable control does not become a missing catalog object.
    census={}
    for group in ("singularity","geometric_mirror_artifact","random_control"):
        relevant=[x for x in obs if x["group"]==group]
        good=[x for x in relevant if x["status"]=="OK"]
        census[group]={"requested":len(relevant),"available":len(good),"unavailable":len(relevant)-len(good)}
        for rad in (2,10,30,60):
            census[group]["near_"+str(rad)+"arcsec"]=sum(x["num_within_arcsec"][str(rad)]>0 for x in good)
    multi=independent_multiband(target["ra_deg"],target["dec_deg"])
    record={"schema":"janus.cosmos.zerkalo_janusa.3773.catalog.v1",
      "date":"2026-11-21","source_fingerprint":result["fingerprint_sha256"],"target":target,
      "observations":obs,"counts":census,"multi_catalog":multi,
      "match_type":"independent_catalog_query_not_direct_telescope_observation",
      "caveats":["This coordinate is a SHA-256 pseudorandom output; proximity to catalog entries is expected in some random directions.",
        "SIMBAD TAP 40-row cap can truncate dense fields and status UNAVAILABLE is not evidence of absence.",
        "Uniform random ICRS controls are not matched on Galactic latitude or survey coverage.",
        "Gaia and WISE share sky source populations and are not automatically independent confirmations of unusual physics.",
        "No repeated imaging or validated telescope pointing at the target has been carried out.",
        "No exclusive novelty, probe or technosignature claims are supported by sky coincidences."],
      "discovery_claim":False,"production_change":False}
    save("3773_catalog_checks.json",record)
    print("JANUS_3773_CATALOG",json.dumps({"target":target,"counts":census,"multi_catalog":multi},ensure_ascii=False),flush=True)
    return record

def main():
    result=compute()
    controls=simulate_random_controls(result)
    checks=catalog_verify(result,controls)
    sig=checks["counts"]["singularity"]
    report={"name":"Зеркало Януса — Сингулярность №3773","stage":"JANUS-COSMOS-ZERKALO-3773",
       "run":"R2", "coordinate":result["direction"],
       "status":"CATALOG_READ_ONLY_COMPLETE" if sig["available"]==1 else "CATALOG_UNAVAILABLE__HOLD",
       "simbad_60arcsec_sources_near_target":checks["observations"][0]["catalog_entries_top40"],
       "simbad_controls":checks["counts"]["random_control"],
       "gaia":checks["multi_catalog"]["GAIA_DR3"],"allwise":checks["multi_catalog"]["ALLWISE"],
       "novel_technosignature_found":False,"probe_found":False,
       "scientific_conclusion":"No unknown astronomical source is inferred from a deterministic hash of source coordinates. Catalogue neighbors require new independent evidence."}
    save("3773_result.json",report)
    print("JANUS_3773_RESULT",json.dumps(report,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
