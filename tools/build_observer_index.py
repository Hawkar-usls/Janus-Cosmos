#!/usr/bin/env python3
from __future__ import annotations

import argparse, csv, hashlib, json, re
from datetime import datetime, timezone
from pathlib import Path

LAT_KEYS = {"lat","latitude","lat_deg","latitude_deg"}
LON_KEYS = {"lon","lng","long","longitude","lon_deg","lng_deg","longitude_deg"}
RA_KEYS = {"ra","ra_deg","right_ascension","right_ascension_deg"}
DEC_KEYS = {"dec","dec_deg","declination","declination_deg"}
LABEL_KEYS = ("site","location_name","place","name","target","label","title","artifact_id","row_id","id")
DATE_KEYS = ("date","event_date","created_date","created_date_local","coverage_start","timestamp","datetime")
SKIP_PARTS = {".git","node_modules","archive","vendor"}
SKIP_TOKENS = ("cousteau","kusto")
# These artifacts are valid scientific metadata representations, but their nested
# coordinates duplicate entities rendered by a dedicated canonical layer.
# They remain available to Bridge/metadata logic and are suppressed only from the
# generic rendered-coordinate scan.
SKY_RENDER_QUOTIENT_FILES = {"PALOMAR_OBSERVATION_DAY_LEDGER_v0.1.json"}
COORD_RE = re.compile(r'(?P<lat>\d{1,2}(?:\.\d+)?)\s*°?\s*(?P<ns>[NS])\s*[,;/ ]+\s*(?P<lon>\d{1,3}(?:\.\d+)?)\s*°?\s*(?P<ew>[EW])', re.I)

def norm_key(k):
    return str(k).strip().lower().replace("-", "_").replace(" ", "_")

def as_float(v):
    try:
        if isinstance(v, bool): return None
        return float(v)
    except (TypeError, ValueError):
        return None

def label_of(d, fallback):
    for k in LABEL_KEYS:
        v=d.get(k)
        if isinstance(v,(str,int,float)) and 0 < len(str(v)) <= 180:
            return str(v)
    return fallback

def date_of(d):
    for k in DATE_KEYS:
        v=d.get(k)
        if isinstance(v,(str,int,float)) and str(v).strip():
            return str(v)
    return None

def source_url_of(d):
    for k in ("source_url","url","primary_url","display_url"):
        v=d.get(k)
        if isinstance(v,str) and v.startswith(("http://","https://")):
            return v
    return None

def valid_earth(lat, lon):
    return lat is not None and lon is not None and -90 <= lat <= 90 and -180 <= lon <= 180

def valid_sky(ra, dec):
    return ra is not None and dec is not None and -360 <= ra <= 360 and -90 <= dec <= 90

def parse_registry(path):
    if not path.exists(): return {}, []
    d=json.loads(path.read_text(encoding="utf-8"))
    aliases={}; sites=[]
    for s in d.get("sites",[]):
        for a in set([s.get("canonical_name","")] + list(s.get("aliases") or [])):
            if a: aliases[a.strip().casefold()]=s
        sites.append(s)
    return aliases, sites

def file_sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20), b""): h.update(chunk)
    return h.hexdigest()

def skip_path(path):
    parts={p.lower() for p in path.parts}
    if parts & SKIP_PARTS: return True
    low=str(path).lower()
    return any(t in low for t in SKIP_TOKENS)

def scan_json_file(path, source_root, source_family, aliases):
    earth=[]; sky=[]; unmapped=[]
    try:
        if path.stat().st_size > 3_000_000: return earth,sky,unmapped,{"state":"SKIPPED_SIZE"}
        data=json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        return earth,sky,unmapped,{"state":"PARSE_FAIL","error":str(e)[:160]}
    rel=str(path.relative_to(source_root)).replace("\\","/"); fallback=path.stem

    def emit_earth(d, lat, lon, method, trail):
        if not valid_earth(lat,lon): return
        earth.append({
            "point_id": f"{source_family}:{rel}:{trail}", "site": label_of(d,fallback),
            "lat": round(lat,7), "lng": round(lon,7),
            "coordinate_precision": d.get("coordinate_precision") or d.get("precision") or method,
            "row_kind": d.get("row_kind") or d.get("kind") or "SOURCE_COORDINATE", "date": date_of(d),
            "region": d.get("region") or d.get("country") or d.get("area"),
            "classification_state": d.get("classification_state") or d.get("status") or "SOURCE_BOUND_COORDINATE",
            "source_grade": d.get("source_grade") or ("JANUS_META_PUBLIC_JSON" if source_family=="JANUS_META" else "JANUS_COSMOS_JSON"),
            "source_url": source_url_of(d), "source_family": source_family, "source_path": rel, "coordinate_method": method
        })

    def emit_sky(d, ra, dec, method, trail):
        if not valid_sky(ra,dec): return
        sky.append({
            "target_id": f"{source_family}:{rel}:{trail}", "label": label_of(d,fallback),
            "ra": round(ra%360.0,7), "dec": round(dec,7),
            "row_kind": d.get("row_kind") or d.get("kind") or "JANUS_SKY_COORDINATE", "date": date_of(d),
            "classification_state": d.get("classification_state") or d.get("status") or "SOURCE_BOUND_SKY_COORDINATE",
            "source_family": source_family, "source_path": rel, "source_url": source_url_of(d), "coordinate_method": method
        })

    def walk(x, trail="$"):
        if isinstance(x,dict):
            nk={norm_key(k): k for k in x}
            lat=lon=None
            for k in LAT_KEYS:
                if k in nk: lat=as_float(x[nk[k]]); break
            for k in LON_KEYS:
                if k in nk: lon=as_float(x[nk[k]]); break
            if valid_earth(lat,lon): emit_earth(x,lat,lon,"EXPLICIT_LAT_LON",trail)
            ra=dec=None
            for k in RA_KEYS:
                if k in nk: ra=as_float(x[nk[k]]); break
            for k in DEC_KEYS:
                if k in nk: dec=as_float(x[nk[k]]); break
            if valid_sky(ra,dec): emit_sky(x,ra,dec,"EXPLICIT_RA_DEC",trail)
            site = x.get("site") or x.get("location_name")
            if isinstance(site,str) and site.strip():
                reg=aliases.get(site.strip().casefold())
                if reg:
                    emit_earth(x,float(reg["lat"]),float(reg["lon"]),"SITE_REGISTRY_ALIAS",trail)
                elif not any(p["point_id"].endswith(trail) for p in earth):
                    unmapped.append({"site":site,"row_id":x.get("row_id"),"date":date_of(x),"source_family":source_family,"source_path":rel,"json_path":trail,"coordinate_state":"UNRESOLVED"})
            for k,v in x.items():
                if isinstance(v,str):
                    m=COORD_RE.search(v)
                    if m:
                        la=float(m.group("lat")) * (-1 if m.group("ns").upper()=="S" else 1)
                        lo=float(m.group("lon")) * (-1 if m.group("ew").upper()=="W" else 1)
                        emit_earth(x,la,lo,"DIRECTIONAL_COORDINATE_STRING",trail+"."+str(k))
                walk(v,trail+"."+str(k))
        elif isinstance(x,list):
            for i,v in enumerate(x): walk(v,f"{trail}[{i}]")
    walk(data)
    return earth,sky,unmapped,{"state":"OK","sha256":file_sha256(path)}

def dedupe(points, kind):
    out=[]; seen=set()
    for p in points:
        if kind=="earth": key=(round(float(p["lat"]),6),round(float(p["lng"]),6),p.get("source_path"),p.get("site"),p.get("date"))
        else: key=(round(float(p["ra"]),6),round(float(p["dec"]),6),p.get("source_path"),p.get("label"),p.get("date"))
        if key in seen: continue
        seen.add(key); out.append(p)
    return out

def parse_plates(path):
    rows=[]
    if not path or not path.exists(): return rows
    with path.open(encoding="utf-8",newline="") as f:
        for r in csv.DictReader(f):
            ra=as_float(r.get("ra_deg")); dec=as_float(r.get("dec_deg"))
            if not valid_sky(ra,dec): continue
            rows.append({"target_id":"POSS1:"+str(r.get("plate_id")),"plate_id":r.get("plate_id"),"label":r.get("plate_id"),"ra":round(ra%360,7),"dec":round(dec,7),"row_kind":"POSS1_PLATE_CENTER","classification_state":"SURVEY_POSITION_NOT_ANOMALY","source_family":"POSS1_PUBLIC_MANIFEST","source_path":"jannefi/poss1-plate-slice:data/plate_manifest.csv"})
    return rows

def load_day_ledger(path):
    if not path or not path.exists(): return {}
    try: data=json.loads(path.read_text(encoding="utf-8"))
    except Exception: return {}
    return {r.get("date"):r for r in data.get("rows",[]) if isinstance(r.get("date"),str)}

def enrich_plates_from_days(plates, day_by_date):
    by_plate={}
    for ds,row in day_by_date.items():
        for pid in row.get("plate_ids") or []:
            by_plate[str(pid)]={
                "observation_date":ds,
                "date":ds,
                "day_candidate_count":row.get("candidate_count"),
                "day_plate_count":row.get("plate_count"),
                "day_tile_count":row.get("tile_count"),
                "day_artifact_state":row.get("artifact_state"),
                "day_detector_quality":row.get("detector_quality"),
                "day_exposure":row.get("exposure")
            }
    for p in plates:
        p.update(by_plate.get(str(p.get("plate_id")),{}))
    return len(by_plate)

def match_state(p, day_by_date):
    d=str(p.get("date") or "")
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}",d):
        if d < "1949-11-19" or d > "1957-04-28": return "NULL_TEMPORAL_NONOVERLAP"
        if day_by_date:
            row=day_by_date.get(d)
            if row is None: return "PALOMAR_NO_OBSERVATION_ROW"
            if row.get("observed_flag") is False: return "PALOMAR_NOT_OBSERVED"
            return "EXACT_DAY_JOIN_AVAILABLE"
        return "PENDING_PALOMAR_312_DAY_TABLE"
    if re.fullmatch(r"\d{4}",d):
        y=int(d); return "BLOCKED_DATE_PRECISION" if 1949 <= y <= 1957 else "NULL_TEMPORAL_NONOVERLAP"
    return "NOT_EXACT_DATE_JOINABLE"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo-root",type=Path,default=Path(".")); ap.add_argument("--meta-root",type=Path)
    ap.add_argument("--plate-manifest",type=Path); ap.add_argument("--day-ledger",type=Path); ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args(); repo=args.repo_root.resolve()
    aliases, bootstrap_sites=parse_registry(repo/"domains/earth/SITE_COORDINATE_REGISTRY_v0.2.json")
    earth=[]; sky=[]; unmapped=[]; scan_receipts=[]; suppressed_sky=[]
    roots=[(repo/"domains","JANUS_COSMOS")]
    if args.meta_root and args.meta_root.exists(): roots.append((args.meta_root.resolve(),"JANUS_META"))
    for root,family in roots:
        if not root.exists(): continue
        for path in sorted(root.rglob("*.json")):
            if skip_path(path): continue
            e,s,u,receipt=scan_json_file(path,root,family,aliases)
            earth.extend(e); unmapped.extend(u)
            if family=="JANUS_COSMOS" and path.name in SKY_RENDER_QUOTIENT_FILES:
                if s:
                    suppressed_sky.append({"source_family":family,"path":str(path.relative_to(root)).replace("\\","/"),"suppressed_coordinate_rows":len(s),"reason":"OBSERVATION_DAY_METADATA_DUPLICATES_CANONICAL_POSS1_PLATE_CENTER_LAYER"})
            else:
                sky.extend(s)
            if receipt.get("state")!="OK": scan_receipts.append({"source_family":family,"path":str(path.relative_to(root)),"receipt":receipt})
    earth=dedupe(earth,"earth"); sky=dedupe(sky,"sky"); plates=parse_plates(args.plate_manifest); day=load_day_ledger(args.day_ledger)
    enriched_plate_count=enrich_plates_from_days(plates,day)
    for p in earth: p["palomar_match_state"]=match_state(p,day)
    states={}
    for p in earth: states[p["palomar_match_state"]]=states.get(p["palomar_match_state"],0)+1
    out={
        "schema":"janus.cosmos.observer.autosync_index.v0.4","generated_at_utc":datetime.now(timezone.utc).isoformat(),
        "authority_boundary":"GENERATED_PRESENTATION_INDEX__SOURCE_JSON_REMAINS_CANONICAL",
        "autosync":{"janus_cosmos_scan":True,"janus_meta_scan":bool(args.meta_root and args.meta_root.exists()),"protected_exclusions":["cousteau","kusto"],"plate_manifest_rows":len(plates),"palomar_day_rows":len(day)},
        "representation_quotient":{
            "law":"SAME_SKY_ENTITY_VIA_MULTIPLE_REPRESENTATIONS -> ONE_RENDERED_ENTITY",
            "suppressed_generic_sky_layers":suppressed_sky,
            "suppressed_sky_coordinate_rows":sum(x["suppressed_coordinate_rows"] for x in suppressed_sky),
            "canonical_palomar_render_layer":"POSS1_PUBLIC_MANIFEST plate centers",
            "palomar_day_ledger_role":"TEMPORAL_AND_DETECTOR_METADATA_ENRICHMENT__NOT_SECOND_SKY_POINT_LAYER",
            "plate_centers_enriched_with_observation_day":enriched_plate_count
        },
        "counts":{"earth_points":len(earth),"earth_unmapped_source_rows":len(unmapped),"janus_sky_targets":len(sky),"poss1_plate_centers":len(plates),"all_space_points":len(sky)+len(plates)},
        "earth_points":earth,"earth_unmapped":unmapped[:2000],"janus_sky_targets":sky,"poss1_plate_centers":plates,
        "bridge":{"palomar_window":["1949-11-19","1957-04-28"],"palomar_unique_observation_days_expected":312,"palomar_day_rows_materialized":len(day),"match_state_counts":states,"exact_day_join_gate":"OPEN" if day and len(day)>=312 else "BLOCKED_312_DAY_TABLE","laws":["NO_OBSERVATION != ZERO","TEMPORAL_NONOVERLAP != PHYSICAL_NONCORRELATION","YEAR_PRECISION != DAY_PRECISION","REPORT_COUNT != PHYSICAL_RATE"]},
        "bootstrap_site_registry_count":len(bootstrap_sites),"scan_receipts":scan_receipts,
        "canonical_seal":"THE OBSERVER DISCOVERS SOURCE COORDINATES AUTOMATICALLY, BUT RENDERS EACH DECISION-RELEVANT ENTITY ONCE. OBSERVATION-DAY METADATA ENRICHES POSS-I PLATES; IT DOES NOT CREATE A SECOND SKY."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True); args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("OBSERVER_INDEX",json.dumps(out["counts"],sort_keys=True)); print("BRIDGE_STATES",json.dumps(states,sort_keys=True)); print("REPRESENTATION_QUOTIENT",json.dumps(out["representation_quotient"],sort_keys=True))

if __name__=="__main__": main()
