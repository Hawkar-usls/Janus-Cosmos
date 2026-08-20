from __future__ import annotations
import json
from datetime import datetime, timezone
from pathlib import Path
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

OUT=Path('data/love/LOVE-EDEM-STARGATE-PALOMAR-VIZIER-DATES-v1-LATEST-RECEIPT.json')
TARGETS={
 'LOVE':(204.30267916666668,-36.78240527777778),
 'EDEM':(139.22409686590188,30.26038779947318),
 'STARGATE_ABYDOS_GEOMETRY':(223.415064157,33.979315670),
}

def clean(v):
    try:
        if getattr(v,'mask',False): return None
    except Exception: pass
    try:
        if hasattr(v,'item'): v=v.item()
    except Exception: pass
    if isinstance(v,bytes): v=v.decode('utf-8','replace')
    if isinstance(v,(str,int,float,bool)) or v is None: return v
    return str(v)

def query(name,ra,dec):
    out={'status':'UNKNOWN','target':{'ra_deg':ra,'dec_deg':dec}}
    try:
        viz=Vizier(columns=['**','+_r'],row_limit=100)
        tabs=viz.query_region(SkyCoord(ra*u.deg,dec*u.deg,frame='icrs'),radius=8*u.deg,catalog='VI/25/nposs')
        rows=[]
        for tab in tabs:
            cols=list(tab.colnames)
            out['columns']=cols
            for row in tab:
                d={c:clean(row[c]) for c in cols}
                rows.append(d)
        out['rows']=rows[:30]
        out['row_count_returned']=len(rows)
        out['status']='OK'
    except Exception as e:
        out['status']='QUERY_ERROR'; out['error']=f'{type(e).__name__}: {e}'
    return out

def main():
    results={k:query(k,*v) for k,v in TARGETS.items()}
    payload={
      'schema':'janus.cosmos.love_edem_stargate.palomar_vizier_dates.receipt.v1',
      'experiment_id':'LOVE-EDEM-STARGATE-PALOMAR-VIZIER-DATES-v1',
      'run_time_utc':datetime.now(timezone.utc).isoformat(),'status':'COMPLETE',
      'catalog':'VI/25/nposs','results':results,
      'firewall':{'nearest_plate_is_not_reflection_evidence':True,'whiteoak_is_not_artificial_mirror_evidence':True,'plate_date_is_not_causation':True,'claim_ceiling':'PALOMAR_PLATE_METADATA_ONLY'}
    }
    OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(payload,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:{'status':v['status'],'columns':v.get('columns'),'first_rows':v.get('rows',[])[:3],'error':v.get('error')} for k,v in results.items()},ensure_ascii=False,indent=2))
if __name__=='__main__': main()
