#!/usr/bin/env python3
from __future__ import annotations
import json, re, hashlib
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
import requests

OUT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-EA-TPHASE-MISSING-900-PROVENANCE-PROBE-2026-08-21-v1.0.json')
URLS={
 'mgds_landing':'https://www.marine-geo.org/tools/files/30497',
 'mgds_uid':'https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid=30497',
 'datacite':'https://api.datacite.org/dois/10.26022/IEDA/330497',
 'crossref':'https://api.crossref.org/works/10.1029/2022JB024008',
}
WAYBACK_TARGETS=[
 'https://www.marine-geo.org/tools/files/30497',
 'https://www.marine-geo.org/tools/search/entry.php?id=EA_Hydroacoustics',
 'http://www.marine-geo.org/tools/search/entry.php?id=EA_Hydroacoustics',
]

def sha(b): return hashlib.sha256(b).hexdigest()
def meta(r):
 return {'requested_url':r.request.url,'final_url':r.url,'status':r.status_code,'content_type':r.headers.get('content-type'),'bytes':len(r.content),'sha256':sha(r.content)}

def main():
 s=requests.Session(); s.headers['User-Agent']='Mozilla/5.0 Janus-Echo-Cousteau/1.0 missing-900 provenance audit'
 rep={'artifact_id':'JANUS-ECHO-COUSTEAU-EA-TPHASE-MISSING-900-PROVENANCE-PROBE-2026-08-21-v1.0','created_utc':datetime.now(timezone.utc).isoformat(),'known_contract':{'paper_reported_events':6843,'deposited_rows':5943,'difference':900,'deposited_hydrophone_counts':[4,5,6,7,8]},'current_sources':{},'wayback':[],'candidate_explanations':[]}
 for k,u in URLS.items():
  try:
   r=s.get(u,timeout=60,allow_redirects=True); item=meta(r)
   if 'json' in (r.headers.get('content-type') or '').lower():
    try:item['json']=r.json()
    except Exception:item['text_prefix']=r.text[:20000]
   else:
    text=r.text; item['text_prefix']=text[:30000]
    item['numeric_hits']={n:[m.start() for m in re.finditer(str(n),text)][:20] for n in [30497,2504732,5943,6843,900]}
    item['file_uid_like']=list(dict.fromkeys(re.findall(r'(?:data[_-]?uid|file[_-]?uid)[^0-9]{0,20}([0-9]{5,})',text,re.I)))[:50]
   rep['current_sources'][k]=item
  except Exception as e: rep['current_sources'][k]={'error':f'{type(e).__name__}: {e}'}

 for target in WAYBACK_TARGETS:
  cdx='https://web.archive.org/cdx/search/cdx?url='+quote(target,safe='')+'&output=json&fl=timestamp,original,statuscode,digest,mimetype&filter=statuscode:200&filter=collapse:digest&from=2021&to=2026'
  ent={'target':target,'cdx_url':cdx}
  try:
   r=s.get(cdx,timeout=90); ent['cdx']=meta(r); rows=r.json() if r.status_code==200 else []
   ent['snapshots']=rows[1:] if isinstance(rows,list) and rows and isinstance(rows[0],list) else []
   probes=[]; snaps=ent['snapshots']; picks=[]
   if snaps:
    picks=[snaps[0],snaps[-1]]
    if len(snaps)>2:picks.append(snaps[len(snaps)//2])
   seen=set()
   for row in picks:
    ts=row[0]
    if ts in seen: continue
    seen.add(ts); wu=f'https://web.archive.org/web/{ts}id_/{target}'
    try:
     rr=s.get(wu,timeout=90,allow_redirects=True); txt=rr.text
     probes.append({'timestamp':ts,'url':wu,**meta(rr),'numeric_hits':{n:[m.start() for m in re.finditer(str(n),txt)][:20] for n in [30497,2504732,5943,6843,900]},'file_uid_like':list(dict.fromkeys(re.findall(r'(?:data[_-]?uid|file[_-]?uid)[^0-9]{0,20}([0-9]{5,})',txt,re.I)))[:50],'text_prefix':txt[:12000]})
    except Exception as e: probes.append({'timestamp':ts,'url':wu,'error':f'{type(e).__name__}: {e}'})
   ent['snapshot_probes']=probes
  except Exception as e: ent['error']=f'{type(e).__name__}: {e}'
  rep['wayback'].append(ent)
 rep['candidate_explanations']=[
  {'id':'H3_STATION_ROWS_OMITTED_FROM_DEPOSIT','status':'STRONGLY_COMPATIBLE_NOT_VERIFIED','reason':'paper catalog allows origins from >=3 arrival picks while current deposited rows have hydrophone counts 4-8 and difference is exactly 900'},
  {'id':'POST_PUBLICATION_DEPOSIT_VERSION_DRIFT','status':'OPEN','reason':'requires historical file UID or archive snapshot evidence'},
  {'id':'PAPER_COUNT_OR_DEPOSIT_DOCUMENTATION_ERROR','status':'OPEN','reason':'cannot be excluded until historical provenance or author clarification is recovered'}
 ]
 rep['status']='PROBE_COMPLETE__PROVENANCE_GATE_REMAINS_OPEN'
 OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(rep,indent=2,ensure_ascii=False),encoding='utf-8')
 print(json.dumps({'status':rep['status'],'wayback_targets':len(rep['wayback'])},indent=2)); return 0

if __name__=='__main__': raise SystemExit(main())
