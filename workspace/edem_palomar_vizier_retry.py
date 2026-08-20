from __future__ import annotations
import json, urllib.parse, urllib.request
from datetime import datetime, timezone
from pathlib import Path
OUT=Path('data/love/EDEM-PALOMAR-VIZIER-RETRY-v1-LATEST-RECEIPT.json')
ra,dec=139.22409686590188,30.26038779947318
base='https://vizier.cds.unistra.fr/viz-bin/asu-tsv'
params={'-source':'VI/25/nposs','-out.all':'','-out.max':'20','-c':f'{ra} {dec}','-c.rm':'360','-sort':'_r'}
url=base+'?'+urllib.parse.urlencode(params)
p={'schema':'janus.cosmos.edem.palomar_vizier_retry.v1','run_time_utc':datetime.now(timezone.utc).isoformat(),'url':url,'status':'UNKNOWN'}
try:
    req=urllib.request.Request(url,headers={'User-Agent':'Janus-Cosmos/1.0'})
    with urllib.request.urlopen(req,timeout=90) as r: txt=r.read().decode('utf-8','replace')
    p['status']='OK'; p['raw_tsv']=txt; p['bytes']=len(txt.encode())
except Exception as e:
    p['status']='QUERY_ERROR'; p['error']=f'{type(e).__name__}: {e}'
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(p,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':p['status'],'bytes':p.get('bytes'),'error':p.get('error'),'preview':p.get('raw_tsv','')[:5000]},ensure_ascii=False,indent=2))
