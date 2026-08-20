from __future__ import annotations

import json, math
from datetime import datetime, timezone
from pathlib import Path
import requests
from astropy.coordinates import SkyCoord, FK4
import astropy.units as u
from astroquery.vizier import Vizier

OUT=Path("data/love/LOVE-EDEM-STARGATE-PALOMAR-PLATE-PROVENANCE-v1-LATEST-RECEIPT.json")
TARGETS={
 "LOVE":(204.30267916666668,-36.78240527777778),
 "EDEM":(139.22409686590188,30.26038779947318),
 "STARGATE_ABYDOS_GEOMETRY":(223.415064157,33.979315670),
}
XE={"EDEM":"XE313","STARGATE_ABYDOS_GEOMETRY":"XE272"}
BASE="https://irsa.ipac.caltech.edu/data/DSS/images/dss1red"

def parse_fits_cards(blob:bytes):
    out={}
    for i in range(0,len(blob)-79,80):
        card=blob[i:i+80].decode('ascii','replace')
        key=card[:8].strip()
        if key=='END': break
        if len(card)>=10 and card[8]=='=':
            raw=card[10:80].split('/',1)[0].strip()
            if raw.startswith("'") and "'" in raw[1:]:
                val=raw[1:raw[1:].find("'")+1]
            else: val=raw
            out[key]=val
    return out

def probe_xe(pid):
    url=f"{BASE}/dss1red_{pid}.fits"
    item={"plate_id":pid,"url":url}
    try:
        r=requests.get(url,headers={"Range":"bytes=0-65535"},stream=True,timeout=60)
        item["http_status"]=r.status_code
        item["content_range"]=r.headers.get("Content-Range")
        buf=b''
        for chunk in r.iter_content(8192):
            buf+=chunk
            if len(buf)>=65536: break
        r.close()
        hdr=parse_fits_cards(buf)
        keep={k:v for k,v in hdr.items() if any(t in k.upper() for t in ["DATE","PLATE","REGION","EPOCH","EXPOS","RA","DEC","OBJECT","SURVEY"]) or k in ["CRVAL1","CRVAL2"]}
        item["header_subset"]=keep
        item["header_card_count"]=len(hdr)
        item["status"]="OK" if hdr else "NO_HEADER_PARSED"
    except Exception as e:
        item["status"]="QUERY_ERROR"; item["error"]=f"{type(e).__name__}: {e}"
    return item

def val(row,name):
    try:
        x=row[name]
        if getattr(x,'mask',False): return None
        if hasattr(x,'item'): x=x.item()
        if isinstance(x,bytes): x=x.decode()
        return x
    except Exception: return None

def full_poss_catalog():
    out={"status":"UNKNOWN","rows":[]}
    try:
        v=Vizier(columns=['*'],row_limit=-1)
        tabs=v.get_catalogs('VI/25')
        allrows=[]
        for tab in tabs:
            for row in tab:
                ra=val(row,'RArad'); de=val(row,'DErad')
                if ra is None or de is None: continue
                try:
                    c=SkyCoord(float(ra)*u.rad,float(de)*u.rad,frame=FK4(equinox='B1950')).icrs
                except Exception: continue
                d={
                  "POSS":str(val(row,'POSS')).strip(),"MLP":int(val(row,'MLP')) if val(row,'MLP') is not None else None,
                  "center_icrs_ra_deg":float(c.ra.deg),"center_icrs_dec_deg":float(c.dec.deg),
                  "Obs.Y":val(row,'Obs.Y'),"Obs.M":val(row,'Obs.M'),"Obs.D":val(row,'Obs.D'),
                  "fObs.M":val(row,'fObs.M'),"fObs.D":val(row,'fObs.D'),
                  "ObsE.h":val(row,'ObsE.h'),"ObsE.m":val(row,'ObsE.m'),"Eexp":val(row,'Eexp'),
                }
                allrows.append(d)
        out["status"]="OK"; out["row_count"]=len(allrows); out["rows"]=allrows
    except Exception as e:
        out["status"]="QUERY_ERROR"; out["error"]=f"{type(e).__name__}: {e}"
    return out

def angsep(ra1,de1,ra2,de2):
    return float(SkyCoord(ra1*u.deg,de1*u.deg).separation(SkyCoord(ra2*u.deg,de2*u.deg)).deg)

def main():
    xe={name:probe_xe(pid) for name,pid in XE.items()}
    cat=full_poss_catalog()
    target_matches={}
    if cat.get('status')=='OK':
        rows=cat['rows']
        for name,(ra,de) in TARGETS.items():
            ranked=[]
            for r in rows:
                rr=dict(r); rr['center_separation_deg']=angsep(ra,de,r['center_icrs_ra_deg'],r['center_icrs_dec_deg'])
                rr['whiteoak_extension']=bool((r.get('MLP') or 0)>=938)
                ranked.append(rr)
            ranked.sort(key=lambda x:x['center_separation_deg'])
            target_matches[name]={"nearest_10_plate_centers":ranked[:10],"nearest":ranked[0]}
    payload={
      "schema":"janus.cosmos.love_edem_stargate.palomar_plate_provenance.receipt.v1",
      "experiment_id":"LOVE-EDEM-STARGATE-PALOMAR-PLATE-PROVENANCE-v1","status":"COMPLETE","run_time_utc":datetime.now(timezone.utc).isoformat(),
      "xe_header_probes":xe,"cds_vi25":{"status":cat.get('status'),"row_count":cat.get('row_count'),"error":cat.get('error')},
      "target_plate_center_matches":target_matches,
      "firewall":{"nearest_plate_center_is_exact_footprint_proof":False,"plate_date_is_reflection_cause":False,"whiteoak_coverage_is_artificial_mirror_evidence":False,"claim_ceiling":"PLATE_PROVENANCE_ONLY"}
    }
    OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(payload,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"xe":xe,"targets":{k:v.get('nearest') for k,v in target_matches.items()}},ensure_ascii=False,indent=2,default=str))
if __name__=='__main__': main()
