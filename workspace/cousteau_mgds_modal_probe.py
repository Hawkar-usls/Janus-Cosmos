#!/usr/bin/env python3
import json, re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup

OUT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-MGDS-DOWNLOAD-MODAL-PROBE-2026-08-21-v1.0.json')
LANDING='https://www.marine-geo.org/tools/files/30497'
UIDS='https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid=30497'
MODAL='https://www.marine-geo.org/services/download/download_modal.php'
s=requests.Session(); s.headers.update({'User-Agent':'Janus-Echo-Cousteau/1.0 modal probe','Referer':LANDING})
rep={'artifact_id':'JANUS-ECHO-COUSTEAU-MGDS-DOWNLOAD-MODAL-PROBE-2026-08-21-v1.0','created_at_utc':datetime.now(timezone.utc).isoformat()}

r=s.get(UIDS,timeout=45); rep['uid_lookup']={'status':r.status_code,'url':r.url,'text':r.text,'headers':dict(r.headers)}; uids=r.json(); uid=str(uids[0])
m=s.post(MODAL,data={'FileDownload':uid,'data_set_uid':'30497'},timeout=60,allow_redirects=True)
rep['modal']={'status':m.status_code,'url':m.url,'headers':dict(m.headers),'body':m.text}
soup=BeautifulSoup(m.text,'html.parser')
rep['forms']=[]
for i,f in enumerate(soup.find_all('form')):
 vals=[]
 for x in f.find_all(['input','button','select','textarea']):
  vals.append({'tag':x.name,'type':x.get('type'),'name':x.get('name'),'value':x.get('value'),'id':x.get('id'),'class':x.get('class'),'text':x.get_text(' ',strip=True)[:500]})
 rep['forms'].append({'index':i,'action':urljoin(m.url,f.get('action') or ''),'method':f.get('method'),'id':f.get('id'),'class':f.get('class'),'fields':vals,'html':str(f)[:20000]})
rep['links']=[{'href':urljoin(m.url,a.get('href')) if a.get('href') else None,'id':a.get('id'),'class':a.get('class'),'text':a.get_text(' ',strip=True)[:500]} for a in soup.find_all('a')]
rep['buttons']=[{'id':b.get('id'),'class':b.get('class'),'name':b.get('name'),'value':b.get('value'),'type':b.get('type'),'onclick':b.get('onclick'),'text':b.get_text(' ',strip=True)[:500]} for b in soup.find_all('button')]
rep['scripts']=[]
for sc in soup.find_all('script'):
 txt=sc.string or sc.get_text() or ''
 if txt.strip(): rep['scripts'].append(txt[:30000])
# Surface endpoint-like strings anywhere in body.
rep['endpoint_strings']=list(dict.fromkeys(re.findall(r'''["']([^"']*(?:download|file|archive)[^"']*)["']''',m.text,re.I)))[:300]
OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(rep,indent=2,ensure_ascii=False),encoding='utf-8')
print('uid',uid,'modal status',m.status_code,'bytes',len(m.content),'forms',len(rep['forms']),'links',len(rep['links']))
