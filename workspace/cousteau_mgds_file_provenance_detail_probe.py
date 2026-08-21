#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,re
from datetime import datetime,timezone
from pathlib import Path
import requests
from bs4 import BeautifulSoup

OUT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-MGDS-FILE-2504732-PROVENANCE-DETAIL-PROBE-2026-08-21-v1.0.json')
URLS={
 'datafileservice':'https://www.marine-geo.org/services/xml/datafileservice.php?data_set_uid=30497',
 'datasetiso':'https://www.marine-geo.org/services/xml/datasetisoservice.php?data_set_uid=30497',
 'download_dataset':'https://api.marine-geo.org/services/download/download.php?data_set_uid=30497',
 'file_uids':'https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid=30497',
 'filespage_js':'https://www.marine-geo.org/tools/search/js/filespage.js?a=20200721',
 'filesmodel_js':'https://www.marine-geo.org/tools/search/js/filesmodel.js?a=20200721'
}

def sha(b):return hashlib.sha256(b).hexdigest()
def meta(r):return {'status':r.status_code,'requested_url':r.request.url,'final_url':r.url,'content_type':r.headers.get('content-type'),'content_disposition':r.headers.get('content-disposition'),'bytes':len(r.content),'sha256':sha(r.content)}

def extract_text_fields(text):
 fields={}
 # retain small contexts around exact UID and likely provenance/date terms
 for pat in ['2504732','EA_CTBTO_catalog_all','date','created','modified','updated','uploaded','submit','size','file_name','filename','data_uid','data_set_uid']:
  hits=[]
  for m in re.finditer(pat,text,re.I):
   hits.append(text[max(0,m.start()-350):min(len(text),m.end()+700)])
   if len(hits)>=12:break
  if hits:fields[pat]=hits
 return fields

def main():
 s=requests.Session();s.headers['User-Agent']='Mozilla/5.0 Janus-Echo-Cousteau/1.0 file provenance audit'
 rep={'artifact_id':'JANUS-ECHO-COUSTEAU-MGDS-FILE-2504732-PROVENANCE-DETAIL-PROBE-2026-08-21-v1.0','created_utc':datetime.now(timezone.utc).isoformat(),'dataset_uid':30497,'file_uid':'2504732','responses':{}}
 discovered=[]
 for key,u in URLS.items():
  try:
   r=s.get(u,timeout=90,allow_redirects=True);item=meta(r);ctype=(r.headers.get('content-type') or '').lower()
   if 'json' in ctype:
    try:item['json']=r.json()
    except Exception:item['text_prefix']=r.text[:5000]
   else:
    txt=r.text;item['contexts']=extract_text_fields(txt);item['text_prefix']=txt[:12000]
    # surface endpoint-like strings from JS/XML/HTML
    discovered += re.findall(r'''["']([^"']*(?:file|data|download|xml|service)[^"']*)["']''',txt,re.I)
   rep['responses'][key]=item
  except Exception as e:rep['responses'][key]={'error':f'{type(e).__name__}: {e}'}

 candidates=[]
 for raw in discovered:
  if len(raw)>300:continue
  low=raw.lower()
  if any(x in low for x in ['uid','file','download','service']):candidates.append(raw)
 rep['discovered_endpoint_strings']=list(dict.fromkeys(candidates))[:200]

 # Parse verbose MGDS XML semantically when available.
 x=rep['responses'].get('datafileservice',{})
 if x.get('status')==200:
  try:
   text=requests.get(URLS['datafileservice'],timeout=90).text
   soup=BeautifulSoup(text,'xml')
   uid_nodes=[n for n in soup.find_all(string=re.compile('2504732'))]
   ancestors=[]
   for node in uid_nodes[:20]:
    p=node.parent
    for _ in range(4):
     if p and p.parent:p=p.parent
    ancestors.append(str(p)[:30000] if p else str(node)[:3000])
   rep['file_uid_xml_ancestor_contexts']=ancestors
  except Exception as e:rep['file_uid_xml_parse_error']=f'{type(e).__name__}: {e}'

 rep['interpretation']={'historical_version_provenance_closed':False,'current_file_identity_verified':rep['responses'].get('file_uids',{}).get('json')==['2504732'],'rule':'CURRENT_METADATA_CANNOT_BY_ITSELF_PROVE_WHETHER_900_ROWS_WERE_REMOVED_AFTER_PUBLICATION'}
 rep['status']='PROBE_COMPLETE'
 OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(rep,indent=2,ensure_ascii=False),encoding='utf-8')
 print(json.dumps({'status':rep['status'],'current_file_identity_verified':rep['interpretation']['current_file_identity_verified'],'xml_status':rep['responses'].get('datafileservice',{}).get('status')},indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
