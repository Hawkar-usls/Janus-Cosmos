#!/usr/bin/env python3
"""JANUS COSMOS R3775: archival multi-epoch observations and matched field/source controls.
All original sky coordinates, source-anchor choice, thresholds and catalog endpoints
are frozen before inspecting new epoch photometry. Research only, no production actions.
"""
import csv
import hashlib
import io
import json
import math
import pathlib
import statistics
import time
import urllib.parse

import requests

P=json.loads(pathlib.Path("research/numeric_sky/3775_independent_epochs_protocol.json").read_text())
O=pathlib.Path("runtime/numeric_sky/3775_epochs");O.mkdir(parents=True,exist_ok=True)
HDR={"User-Agent":"JANUS-COSMOS-R3775-Archival-Audit/1.0 (public research read-only)"}
META_COMMIT="10ce2787fde76bbe8660369175f0195d8b29208f"
IRSA="https://irsa.ipac.caltech.edu/TAP/sync"
MAST="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2"
def save(n,obj):
    (O/n).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
def val(row,key):
    if not row:return None
    try:
        v=row.get(key)
        if v is None or str(v).strip().lower() in ("","null","nan","none","--"):return None
        n=float(v)
        return n if math.isfinite(n) else None
    except (TypeError,ValueError):return None
def distance(ra,dec,x,y):
    a,b,c,d=[math.radians(v) for v in (ra,dec,x,y)]
    x=math.sin((b-d)/2)**2+math.cos(b)*math.cos(d)*math.sin((a-c)/2)**2
    return 2*math.degrees(math.asin(min(1,math.sqrt(max(0,x)))))*3600
def get(url,params,timeout=(8,60),post=False):
    t0=time.monotonic()
    if post:r=requests.post(url,headers=HDR,data=params,timeout=timeout)
    else:r=requests.get(url,headers=HDR,params=params,timeout=timeout)
    r.raise_for_status()
    return r,time.monotonic()-t0
def tap(query,columns,limit=5000):
    errors=[]
    for method in ("GET","POST"):
        try:
            p={"REQUEST":"doQuery","LANG":"ADQL","QUERY":query,"FORMAT":"csv"}
            r,sec=get(IRSA,p,timeout=(9,110),post=(method=="POST"))
            s=r.text
            if s.lstrip().startswith(("<html","<?xml","<!DOCTYPE","ERROR")):
                raise ValueError("Unexpected TAP response, not CSV")
            stream=io.StringIO("\n".join(x for x in s.splitlines() if not x.startswith("#")))
            reader=csv.DictReader(stream)
            if not reader.fieldnames or not set(columns).issubset(set(reader.fieldnames)):
                raise ValueError("Unexpected CSV fields "+str(reader.fieldnames)[:200])
            rows=list(reader)
            if len(rows)>limit:raise ValueError("TAP response over safety limit")
            return {"status":"OK","rows":rows,"count":len(rows),"columns":reader.fieldnames,"request_seconds":round(sec,3),
                    "query_sha256":hashlib.sha256(query.encode()).hexdigest(),"method":method,"errors":errors}
        except Exception as e:errors.append(method+": "+str(e)[:230])
    return {"status":"UNAVAILABLE","rows":[],"count":0,"errors":errors,"query_sha256":hashlib.sha256(query.encode()).hexdigest()}
def fixed_source():
    p=P["origin"];t=P["preselected_allwise_neighbor"]
    url="https://raw.githubusercontent.com/Hawkar-usls/janus-meta-registry/"+META_COMMIT+"/data/JANUS-ZERKALO-JANUSA-SINGULARITY-3773-MIRRORS-2026-11-21.json"
    r,_=get(url,None,timeout=(8,40))
    j=r.json()
    q=next(x for x in j["points"] if x["ordinal"]==3775)
    assert abs(q["ra_deg"]-p["ra_deg"])<1e-12 and abs(q["dec_deg"]-p["dec_deg"])<1e-12
    assert j["source"]["fingerprint_sha256"]==p["source_fingerprint_sha256"]
    sep=distance(q["ra_deg"],q["dec_deg"],t["ra_deg"],t["dec_deg"])
    assert abs(sep-t["separation_arcsec"])<.005
    output={"status":"SOURCE_FROZEN","meta_commit":META_COMMIT,"ordinal":3775,
        "coordinate":{"ra_deg":q["ra_deg"],"dec_deg":q["dec_deg"]},
        "preselected_neighbor":t,"angular_separation_arcsec":sep,
        "science_warning":"A 6.66 arcsec source is not a positional match to the mathematical center; this is a post-hoc catalog neighbor."}
    save("preflight.json",output);print("R3775_FROZEN",json.dumps(output),flush=True);return output
def query_allwise():
    t=P["preselected_allwise_neighbor"];ra=t["ra_deg"];dec=t["dec_deg"]
    c="designation,ra,dec,w1mpro,w2mpro,w1sigmpro,w2sigmpro,ph_qual,cc_flags,ext_flg"
    query=("SELECT TOP 100 "+c+" FROM allwise_p3as_psd "
           f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{15/3600}))")
    result=tap(query,["designation","ra","dec","w1mpro","w2mpro"],limit=120)
    result["rows"]=sorted(result["rows"],key=lambda z:distance(ra,dec,val(z,"ra") or -900,val(z,"dec") or -900))
    match=next((x for x in result["rows"] if x.get("designation")==t["id"]),None)
    if match is None or distance(ra,dec,val(match,"ra"),val(match,"dec"))>1.0:
        result["anchor_status"]="HOLD_ORIGINAL_ALLWISE_SOURCE_NOT_REPRODUCED"
    else:
        result["anchor_status"]="PASS_FROZEN_NEIGHBOR_MATCH";result["anchor"]=match
    result["within_10arcsec"]=sum(distance(ra,dec,val(x,"ra"),val(x,"dec"))<=10 for x in result["rows"])
    save("frozen_neighbor_verification.json",result)
    print("R3775_ALLWISE_ANCHOR",json.dumps({k:v for k,v in result.items() if k!="rows"}),flush=True)
    return result
def field_controls():
    ra=P["origin"]["ra_deg"];dec=P["origin"]["dec_deg"];q=[]
    for radius in P["field_controls"]["offset_degrees"]:
        for angle in P["field_controls"]["position_angle_degrees"]:
            a=math.radians(angle)
            x=(ra+radius*math.sin(a)/math.cos(math.radians(dec)))%360
            y=dec+radius*math.cos(a)
            sep=distance(ra,dec,x,y)/3600
            assert abs(sep-radius)<.001
            q.append({"control_id":len(q)+1,"radius_deg":radius,"position_angle_deg":angle,
                "ra_deg":x,"dec_deg":y})
    assert len(q)==16
    save("matched_field_controls.json",{"status":"FROZEN_LOCAL_SPATIAL_CONTROLS","controls":q,
        "limitation":"Similar sky location, not rigorously matched survey exposure density or Galactic latitude."})
    return q
def select_source_controls():
    t=P["preselected_allwise_neighbor"]
    ra=P["origin"]["ra_deg"];dec=P["origin"]["dec_deg"]
    cfg=P["source_controls"];r=cfg["annulus_deg"][1]
    query=("SELECT TOP 5000 designation,ra,dec,w1mpro,w2mpro,w1sigmpro,w2sigmpro,ph_qual,cc_flags,ext_flg "
        "FROM allwise_p3as_psd "
        f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{r})) "
        f"AND w1mpro BETWEEN {cfg['magnitude_match']['w1_range'][0]} AND {cfg['magnitude_match']['w1_range'][1]} "
        f"AND w2mpro BETWEEN {cfg['magnitude_match']['w2_range'][0]} AND {cfg['magnitude_match']['w2_range'][1]}")
    query_result=tap(query,["designation","ra","dec","w1mpro","w2mpro"],limit=5200)
    rows=query_result["rows"];candidates=[]
    for x in rows:
        try:
            a=val(x,"ra");b=val(x,"dec");u=val(x,"w1mpro");v=val(x,"w2mpro")
            if None in (a,b,u,v):continue
            sep=distance(ra,dec,a,b)
            if sep<cfg["annulus_deg"][0]*3600 or sep>r*3600:continue
            if x.get("designation")==t["id"]:continue
            if x.get("cc_flags") not in (None,"","0000"):continue
            if x.get("ext_flg") not in (None,"","0","0.0"):continue
            if x.get("ph_qual") and ("U" in x["ph_qual"][:2].upper()):continue
            w1err=val(x,"w1sigmpro");w2err=val(x,"w2sigmpro")
            color=u-v
            score=(abs(u-t["w1_mag"])/.6)**2+(abs(v-t["w2_mag"])/.6)**2+(abs(color-(t["w1_mag"]-t["w2_mag"]))/.6)**2
            candidates.append({"id":x["designation"],"ra_deg":a,"dec_deg":b,"w1_mag":u,"w2_mag":v,
                "w1_sigma":w1err,"w2_sigma":w2err,"cc_flags":x.get("cc_flags"),
                "ph_qual":x.get("ph_qual"),"field_sep_arcsec":sep,"match_cost":score})
        except (ValueError,TypeError):continue
    candidates.sort(key=lambda x:(x["match_cost"],x["id"]))
    chosen=[]
    for x in candidates:
        if all(distance(x["ra_deg"],x["dec_deg"],z["ra_deg"],z["dec_deg"])>30 for z in chosen):
            chosen.append(x)
        if len(chosen)>=cfg["count"]:break
    result={"status":"PASS_CONTROL_SELECTION" if len(chosen)>=8 else "HOLD_INSUFFICIENT_MATCHED_CONTROLS",
        "catalog_query_status":query_result["status"],"selected":chosen,"count":len(chosen),
        "matching_candidates_available":len(candidates),"source_query_rows":len(rows),
        "catalog_result_truncated_possible":len(rows)>=5000,"query_sha256":query_result.get("query_sha256"),
        "catalog_errors":query_result.get("errors",[]),
        "matching_method":"Frozen chi-square scaled W1/W2/color distance within 0.35 deg and no repeated source within 30 arcsec; no epoch results used"}
    save("selected_source_controls.json",result)
    print("R3775_CONTROLS",json.dumps({k:v for k,v in result.items() if k!="selected"}),flush=True)
    return chosen
def density_controls(fields):
    out=[]
    for f in fields:
        ra=f["ra_deg"];dec=f["dec_deg"];rad=60/3600
        query=("SELECT COUNT(*) AS n FROM allwise_p3as_psd "
               f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{rad}))")
        r=tap(query,["n"],limit=3)
        out.append({"control_id":f["control_id"],"status":r["status"],"counts_allwise_within_60arcsec":
            int(r["rows"][0]["n"]) if r["status"]=="OK" and r["rows"] else None,
            "errors":r.get("errors",[])})
        print("R3775_FIELD",json.dumps(out[-1]),flush=True)
    save("field_catalog_densities.json",{"results":out,"no_instrumental_completeness_guarantee":True})
    return out
def neowise_one(source):
    ra=source["ra_deg"];dec=source["dec_deg"]
    rad=P["archived_epoch_observations"]["match_radius_arcsec"]/3600
    columns="mjd,ra,dec,w1mpro,w1sigmpro,w2mpro,w2sigmpro,ph_qual,cc_flags,w1rchi2,w2rchi2"
    query=("SELECT TOP 900 "+columns+" FROM neowiser_p1bs_psd "
          f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{rad})) ORDER BY mjd")
    r=tap(query,["mjd","ra","dec","w1mpro","w1sigmpro","w2mpro","w2sigmpro"],limit=1000)
    out={"id":source["id"],"ra_deg":ra,"dec_deg":dec,"status":r["status"],
         "query_sha256":r.get("query_sha256"),"query_errors":r.get("errors",[]),
         "raw_exposures":r["count"],"query_truncated_possible":r["count"]>=900,
         "observations":[],"bands":{}}
    for rec in r["rows"]:
        mjd=val(rec,"mjd");a=val(rec,"ra");b=val(rec,"dec")
        if None in (mjd,a,b):continue
        sep=distance(ra,dec,a,b)
        if sep>P["archived_epoch_observations"]["match_radius_arcsec"]+.01:continue
        entry={"mjd":mjd,"ra_deg":a,"dec_deg":b,"sep_arcsec":sep}
        for band,index in (("w1",0),("w2",1)):
            mag=val(rec,band+"mpro");err=val(rec,band+"sigmpro")
            chi=val(rec,band+"rchi2")
            quality=rec.get("ph_qual")
            cc=rec.get("cc_flags")
            good=(mag is not None and err is not None and 0<err<=.35)
            if chi is not None and chi>3:good=False
            if quality and len(quality)>index and quality[index].upper() not in ("A","B","C"):good=False
            if cc and len(cc)>index and cc[index]!="0":good=False
            entry[band]={"mag":mag,"err_mag":err,"chi2":chi,"good":good,"ph_qual":quality,
                "cc_flags":cc}
        out["observations"].append(entry)
    out["observations"].sort(key=lambda x:x["mjd"])
    for band in ("w1","w2"):
        obs=[q for q in out["observations"] if q[band]["good"]]
        seasons=[]
        for q in obs:
            if not seasons or q["mjd"]-seasons[-1]["last_mjd"]>90:
                seasons.append({"first_mjd":q["mjd"],"last_mjd":q["mjd"],"mags":[],"errs":[]})
            block=seasons[-1];block["last_mjd"]=q["mjd"];block["mags"].append(q[band]["mag"]);block["errs"].append(q[band]["err_mag"])
        season_rows=[]
        for s in seasons:
            n=len(s["mags"]);med=statistics.median(s["mags"])
            err=math.sqrt((statistics.median(s["errs"])**2)/n+.05**2)
            season_rows.append({"first_mjd":s["first_mjd"],"last_mjd":s["last_mjd"],
                "median_mag":med,"floor_stderr_mag":err,"n":n})
        qualified=[s for s in season_rows if s["n"]>=3]
        if qualified:
            low=min(qualified,key=lambda s:s["median_mag"]);high=max(qualified,key=lambda s:s["median_mag"])
            delta=high["median_mag"]-low["median_mag"]
            z=delta/math.sqrt(high["floor_stderr_mag"]**2+low["floor_stderr_mag"]**2)
        else:delta=None;z=None
        flag=(len(qualified)>=3 and delta is not None and delta>=.4 and z>=5)
        out["bands"][band]={"good_exposures":len(obs),"season_rows":season_rows,
            "qualified_seasons":len(qualified),"delta_median_mag":delta,"estimated_z":z,
            "variability_screen_positive":flag}
    out["repeatability_status"]="REPEATABLE_SEASON_COVERAGE" if any(
        out["bands"][b]["qualified_seasons"]>=3 for b in ("w1","w2")) else "INSUFFICIENT_RELIABLE_SEASONS"
    save("epoch_"+source["id"].replace("/","_")+".json",out)
    print("R3775_EPOCH",json.dumps({k:v for k,v in out.items() if k not in ("observations",)},ensure_ascii=False),flush=True)
    return out
def panstarrs_one(source):
    params={"ra":source["ra_deg"],"dec":source["dec_deg"],"radius":2.5/3600,"pagesize":1000}
    out={"id":source["id"],"source_coordinates":{k:source[k] for k in ("ra_deg","dec_deg")}}
    for kind in ("mean","detection"):
        url=MAST+"/"+kind+".csv"
        try:
            r,elapsed=get(url,params,timeout=(8,50))
            reader=csv.DictReader(io.StringIO(r.text))
            fields=reader.fieldnames or []
            if not any(k.lower()=="objid" for k in fields):raise ValueError("Not a standard PanSTARRS CSV: "+str(fields)[:180])
            raw=list(reader)
            if len(raw)>1000:raise ValueError("unbounded PanSTARRS result")
            cols={k.lower():k for k in fields}
            selected=[]
            for x in raw:
                d={k.lower():v for k,v in x.items()}
                row={k:d.get(k) for k in ("objid","detectid","obstime","filterid","psfmag","psfmagerr","psfflux","psffluxerr","ndetections","ramean","decmean") if d.get(k) not in (None,"")}
                selected.append(row)
            out[kind]={"status":"OK","rows":len(raw),"columns":fields[:50],
                "sample":selected[:80],"request_seconds":round(elapsed,3)}
        except Exception as e:out[kind]={"status":"UNAVAILABLE","error":str(e)[:270]}
    save("panstarrs_"+source["id"].replace("/","_")+".json",out)
    print("R3775_PANSTARRS",json.dumps({"source":source["id"],"mean":out["mean"]["status"],"detections":out["detection"]["status"],
        "n_detections":out["detection"].get("rows")}),flush=True)
    return out
def main():
    source=fixed_source()
    anchor=query_allwise()
    field=field_controls()
    controls=select_source_controls()
    field_rows=density_controls(field)
    selected=[{"id":P["preselected_allwise_neighbor"]["id"],"ra_deg":P["preselected_allwise_neighbor"]["ra_deg"],
      "dec_deg":P["preselected_allwise_neighbor"]["dec_deg"]}]+controls
    if anchor["anchor_status"]!="PASS_FROZEN_NEIGHBOR_MATCH":
        result={"status":"HOLD_ANCHOR_CATALOG_NOT_REPRODUCED","new_observation":False,"raw_data_changed":False}
        save("audit_result.json",result);return
    epoch=[neowise_one(s) for s in selected]
    optical=[panstarrs_one(s) for s in selected[:5]]
    good_controls=[x for x in epoch[1:] if x["status"]=="OK"]
    target=epoch[0]
    def comparable(band):
        viable=[x for x in good_controls if x["bands"][band]["qualified_seasons"]>=3]
        flag=sum(x["bands"][band]["variability_screen_positive"] for x in viable)
        return {"eligible_controls":len(viable),"flagged_controls":flag,
            "target_flagged":target["bands"][band]["variability_screen_positive"],
            "target_qualified_seasons":target["bands"][band]["qualified_seasons"],
            "target_delta_mag":target["bands"][band]["delta_median_mag"],
            "target_estimated_z":target["bands"][band]["estimated_z"]}
    stats={b:comparable(b) for b in ("w1","w2")}
    target_repeat=target["repeatability_status"]=="REPEATABLE_SEASON_COVERAGE"
    stats_adequate=len(controls)>=8 and all(stats[b]["eligible_controls"]>=8 for b in ("w1","w2"))
    out={"schema":"janus.cosmos.zerkalo.3775.archival_repeat.v1",
       "status":"QUALIFIED_REPEATABILITY_SCREEN_WITH_POSTHOC_CATALOG_NEIGHBOR" if target_repeat and stats_adequate else "HOLD_REPEATABILITY_OR_MATCHED_CONTROLS_INSUFFICIENT",
       "date":"2026-11-21","coordinate":{"ra_deg":P["origin"]["ra_deg"],"dec_deg":P["origin"]["dec_deg"]},
       "source_association":"NOT_ESTABLISHED_6.66_ARCSEC_OFFSET",
       "anchor_status":anchor["anchor_status"],
       "source_control_count":len(controls),"field_control_count":len(field_rows),
       "matched_field_valid_catalog_counts":sum(r["status"]=="OK" for r in field_rows),
       "neowise_target":{k:v for k,v in target.items() if k!="observations"},
       "empirical_control_screens":stats,
       "panstarrs_target":optical[0],
       "allwise_close_neighbours":anchor["within_10arcsec"],
       "preselected_posthoc":True,"no_blind_holdout":True,
       "no_new_observations_requested":True,"no_new_object_or_probe_discovery":True,
       "scientific_limitations":["NEOWISE-R single-exposure only reports detected sources; absent detections are not image-level upper limits.",
       "Catalog crossmatch radius 2.5 arcsec; crowded fields cause blending and false variability.",
       "Control groups are same-neighborhood and brightness-matched but not rigorously matched exposure/depth.",
       "WISE and NEOWISE reuse telescope instrumentation; Pan-STARRS is an independent optical instrument but may not see a red source.",
       "These are archived preexisting observations, not newly commissioned telescope measurements.",
       "A 6.66 arcsec catalog neighbor cannot be promoted to the arbitrary SHA-derived sky address without independent astrometry."],
       "production_change":False,"radio_R42B_still_HOLD":True}
    save("audit_result.json",out)
    save("source_epoch_photometry.json",{"sources":[{k:v for k,v in s.items() if k!="observations"} for s in epoch]})
    save("panstarrs_optical.json",{"sources":optical})
    print("R3775_FINAL",json.dumps(out,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
