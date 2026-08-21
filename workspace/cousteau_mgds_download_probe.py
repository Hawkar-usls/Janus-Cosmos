#!/usr/bin/env python3
import json, re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup

OUT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-MGDS-EA-HYDROACOUSTICS-DOWNLOAD-PROBE-2026-08-21-v1.0.json')
URLS=[
 'https://www.marine-geo.org/tools/files/30497',
 'https://www.marine-geo.org/doi/10.26022/IEDA/330497',
 'https://www.marine-geo.org/inc/swagger/listings/searchinfo.json',
 'https://www.marine-geo.org/services/search/datasets?query=EA_Hydroacoustics',
 'https://www.marine-geo.org/services/search/datasets?q=EA_Hydroacoustics',
 'https://www.marine-geo.org/services/search/datasets?id=EA_Hydroacoustics',
]
s=requests.Session(); s.headers['User-Agent']='Janus-Echo-Cousteau/1.0 protocol probe'
report={'artifact_id':'JANUS-ECHO-COUSTEAU-MGDS-EA-HYDROACOUSTICS-DOWNLOAD-PROBE-2026-08-21-v1.0','created_at_utc':datetime.now(timezone.utc).isoformat(),'responses':[],'script_probes':[]}
script_urls=[]
for u in URLS:
 try:
  r=s.get(u,timeout=60,allow_redirects=True)
  item={'requested':u,'final':r.url,'status':r.status_code,'content_type':r.headers.get('content-type'),'bytes':len(r.content),'headers':{k:v for k,v in r.headers.items() if k.lower() in ['content-type','content-disposition','location']},'prefix':r.text[:5000]}
  if 'html' in (r.headers.get('content-type') or '').lower():
   soup=BeautifulSoup(r.text,'html.parser')
   item['forms']=[]
   for f in soup.find_all('form'):
    inputs=[]
    for x in f.find_all(['input','button','select']):
     inputs.append({'tag':x.name,'type':x.get('type'),'name':x.get('name'),'value':x.get('value'),'id':x.get('id'),'class':x.get('class')})
    item['forms'].append({'action':urljoin(r.url,f.get('action') or ''),'method':f.get('method'),'id':f.get('id'),'class':f.get('class'),'inputs':inputs})
   item['links']=[urljoin(r.url,a.get('href')) for a in soup.find_all('a') if a.get('href')][:300]
   for sc in soup.find_all('script'):
    if sc.get('src'): script_urls.append(urljoin(r.url,sc.get('src')))
    txt=sc.string or sc.get_text() or ''
    if any(k.lower() in txt.lower() for k in ['30497','download','filedownload','archiveserver','data_set_uid','file_uid']):
     item.setdefault('inline_script_hits',[]).append(txt[:12000])
  report['responses'].append(item)
 except Exception as e:
  report['responses'].append({'requested':u,'error':type(e).__name__+': '+str(e)})

for u in list(dict.fromkeys(script_urls))[:80]:
 try:
  r=s.get(u,timeout=30)
  txt=r.text
  hits=[]
  for pat in ['FileDownloadServer','ArchiveDownloadServer','FileServer','downloadSelected','data_set_uid','file_uid','data_uid','30497']:
   for m in re.finditer(pat,txt,re.I):
    lo=max(0,m.start()-700); hi=min(len(txt),m.end()+1200); hits.append(txt[lo:hi])
    if len(hits)>=20: break
   if len(hits)>=20: break
  if hits: report['script_probes'].append({'url':u,'status':r.status_code,'bytes':len(r.content),'hits':hits})
 except Exception as e:
  report['script_probes'].append({'url':u,'error':str(e)})
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print('wrote',OUT,'responses',len(report['responses']),'script hits',len(report['script_probes']))
