from __future__ import annotations
import json, urllib.parse, urllib.request
from datetime import datetime, timezone
from pathlib import Path

OUT=Path('data/love/LOVE-EDEM-STARGATE-PALOMAR-VIZIER-HTTP-v1-LATEST-RECEIPT.json')
TARGETS={
 'LOVE':(204.30267916666668,-36.78240527777778),
 'EDEM':(139.22409686590188,30.26038779947318),
 'STARGATE_ABYDOS_GEOMETRY':(223.415064157,33.979315670),
}
BASE='https://vizier.cds.unistra.fr/viz-bin/asu-tsv'

def fetch_target(name,ra,dec):
    params={
      '-source':'VI/25/nposs','-out.all':'','-out.max':'50','-c':f'{ra} {dec}','-c.rm':'600','-sort':'_r'
    }
    url=BASE+'?'+urllib.parse.urlencode(params)
    out={'url':url,'status':'UNKNOWN'}
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'Janus-Cosmos/1.0'})
        with urllib.request.urlopen(req,timeout=45) as r:
            txt=r.read().decode('utf-8','replace')
        out['status']='OK'; out['raw_tsv']=txt[:100000]; out['bytes']=len(txt.encode())
    except Exception as e:
        out['status']='QUERY_ERROR'; out['error']=f'{type(e).__name__}: {e}'
    return out

def main():
    results={k:fetch_target(k,*v) for k,v in TARGETS.items()}
    payload={'schema':'janus.cosmos.palomar.vizier_http.receipt.v1','status':'COMPLETE','run_time_utc':datetime.now(timezone.utc).isoformat(),'targets':TARGETS,'results':results,'firewall':{'plate_metadata_only':True,'no_reflection_claim':True}}
    OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(payload,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:{'status':v['status'],'bytes':v.get('bytes'),'error':v.get('error'),'preview':v.get('raw_tsv','')[:2000]} for k,v in results.items()},ensure_ascii=False,indent=2))
if __name__=='__main__': main()
