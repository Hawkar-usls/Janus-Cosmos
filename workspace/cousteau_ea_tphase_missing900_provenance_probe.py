#!/usr/bin/env python3
from __future__ import annotations
import json, re, hashlib
from datetime import datetime, timezone
from pathlib import Path
import requests

OUT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-EA-TPHASE-MISSING-900-PROVENANCE-PROBE-2026-08-21-v1.1.json')
URLS={
 'mgds_landing':'https://www.marine-geo.org/tools/files/30497',
 'mgds_uid':'https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid=30497',
 'mgds_filespage_js':'https://www.marine-geo.org/tools/search/js/filespage.js?a=20200721',
 'mgds_filesmodel_js':'https://www.marine-geo.org/tools/search/js/filesmodel.js?a=20200721',
 'datacite':'https://api.datacite.org/dois/10.26022/IEDA/330497',
 'crossref':'https://api.crossref.org/works/10.1029/2022JB024008',
}
WAYBACK_TARGETS=[
 'www.marine-geo.org/tools/files/30497',
 'www.marine-geo.org/tools/search/entry.php?id=EA_Hydroacoustics',
 'marine-geo.org/tools/search/entry.php?id=EA_Hydroacoustics',
]

def sha(b): return hashlib.sha256(b).hexdigest()
def meta(r): return {'requested_url':r.request.url,'final_url':r.url,'status':r.status_code,'content_type':r.headers.get('content-type'),'bytes':len(r.content),'sha256':sha(r.content)}
def text_extract(txt):
 return {
  'numeric_hits':{str(n):len(re.findall(str(n),txt)) for n in [30497,2504732,5943,6843,900]},
  'file_uid_like':list(dict.fromkeys(re.findall(r'(?:data[_-]?uid|file[_-]?uid)[^0-9]{0,30}([0-9]{5,})',txt,re.I)))[:100],
  'endpoint_like':list(dict.fromkeys(re.findall(r'["\']([^"\']*(?:file_uids|download|fileinfo|file_info|data_set_uid|data_uid)[^"\']*)["\']',txt,re.I)))[:100],
  'text_prefix':txt[:20000]
 }

def main():
 s=requests.Session(); s.headers['User-Agent']='Mozilla/5.0 Janus-Echo-Cousteau/1.1 missing-900 provenance audit'
 rep={'artifact_id':'JANUS-ECHO-COUSTEAU-EA-TPHASE-MISSING-900-PROVENANCE-PROBE-2026-08-21-v1.1','created_utc':datetime.now(timezone.utc).isoformat(),'supersedes':'...v1.0','known_contract':{'paper_reported_events':6843,'deposited_rows':5943,'difference':900,'deposited_hydrophone_counts':[4,5,6,7,8]},'current_sources':{},'wayback':[]}
 for k,u in URLS.items():
  try:
   r=s.get(u,timeout=60,allow_redirects=True); item=meta(r)
   ctype=(r.headers.get('content-type') or '').lower()
   if 'json' in ctype:
    try:item['json']=r.json()
    except Exception:item.update(text_extract(r.text))
   else:item.update(text_extract(r.text))
   rep['current_sources'][k]=item
  except Exception as e: rep['current_sources'][k]={'error':f'{type(e).__name__}: {e}'}

 for target in WAYBACK_TARGETS:
  ent={'target':target}
  params={'url':target,'output':'json','fl':'timestamp,original,statuscode,digest,mimetype','filter':['statuscode:200','collapse:digest'],'from':'2021','to':'2026'}
  try:
   r=s.get('https://web.archive.org/cdx/search/cdx',params=params,timeout=90); ent['cdx']=meta(r); ent['cdx_body_prefix']=r.text[:2000]
   rows=r.json() if r.status_code==200 else []
   snaps=rows[1:] if isinstance(rows,list) and rows and isinstance(rows[0],list) else []
   ent['snapshots']=snaps
   picks=[]
   if snaps:
    picks=[snaps[0],snaps[-1]]
    if len(snaps)>2:picks.append(snaps[len(snaps)//2])
   probes=[]; seen=set()
   for row in picks:
    ts,orig=row[0],row[1]
    if ts in seen: continue
    seen.add(ts); wu=f'https://web.archive.org/web/{ts}id_/{orig}'
    try:
     rr=s.get(wu,timeout=90,allow_redirects=True); probes.append({'timestamp':ts,'original':orig,**meta(rr),**text_extract(rr.text)})
    except Exception as e: probes.append({'timestamp':ts,'original':orig,'error':f'{type(e).__name__}: {e}'})
   ent['snapshot_probes']=probes
  except Exception as e: ent['error']=f'{type(e).__name__}: {e}'
  rep['wayback'].append(ent)

 dc=rep['current_sources'].get('datacite',{}).get('json',{})
 attrs=((dc.get('data') or {}).get('attributes') or {}) if isinstance(dc,dict) else {}
 rep['datacite_version_signal']={
  'created':attrs.get('created'),'registered':attrs.get('registered'),'updated':attrs.get('updated'),'published':attrs.get('published'),
  'version':attrs.get('version'),'related_identifiers':attrs.get('relatedIdentifiers'),
  'interpretation':'UPDATED_TIMESTAMP_IS_NOT_BY_ITSELF_FILE_VERSION_PROOF'
 }
 rep['candidate_explanations']=[
  {'id':'H3_STATION_ROWS_OMITTED_FROM_DEPOSIT','status':'STRONGLY_COMPATIBLE_NOT_VERIFIED','reason':'paper admits origins with >=3 arrival picks; current file contains only 4-8 recording-hydrophone rows; arithmetic gap is exactly 900'},
  {'id':'POST_PUBLICATION_DEPOSIT_VERSION_DRIFT','status':'TESTABLE_OPEN','reason':'DataCite metadata has later update timestamp; requires historical snapshot/file UID or checksum evidence'},
  {'id':'PAPER_COUNT_OR_DEPOSIT_DOCUMENTATION_ERROR','status':'OPEN','reason':'cannot be excluded without historical file provenance or explicit author/repository clarification'}
 ]
 rep['status']='PROBE_COMPLETE__WAYBACK_QUERY_REPAIRED__PROVENANCE_GATE_REMAINS_OPEN'
 OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(rep,indent=2,ensure_ascii=False),encoding='utf-8')
 print(json.dumps({'status':rep['status'],'wayback':[(x.get('target'),len(x.get('snapshots',[])),x.get('cdx',{}).get('status')) for x in rep['wayback']]},indent=2)); return 0

if __name__=='__main__': raise SystemExit(main())
