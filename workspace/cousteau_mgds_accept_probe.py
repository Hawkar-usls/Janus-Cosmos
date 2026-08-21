#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup

OUT = Path('data/cousteau/JANUS-ECHO-COUSTEAU-MGDS-DOWNLOAD-ACCEPT-PROBE-2026-08-21-v1.0.json')
LANDING = 'https://www.marine-geo.org/tools/files/30497'
UIDS = 'https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid=30497'
MODAL = 'https://www.marine-geo.org/services/download/download_modal.php'


def snap(r, body_limit=30000):
    ctype=(r.headers.get('content-type') or '').lower()
    textual='text' in ctype or 'html' in ctype or 'json' in ctype or len(r.content) < 50000
    return {
        'requested_url': r.request.url if r.request else None,
        'method': r.request.method if r.request else None,
        'final_url': r.url,
        'status': r.status_code,
        'history': [{'status':x.status_code,'url':x.url,'location':x.headers.get('location')} for x in r.history],
        'headers': {k:v for k,v in r.headers.items() if k.lower() in ['content-type','content-disposition','location','content-length']},
        'bytes': len(r.content),
        'magic_hex': r.content[:32].hex(),
        'body_prefix': r.text[:body_limit] if textual else None,
    }

s=requests.Session()
s.headers.update({'User-Agent':'Janus-Echo-Cousteau/1.0 download-accept diagnostic','Referer':LANDING})
rep={'artifact_id':'JANUS-ECHO-COUSTEAU-MGDS-DOWNLOAD-ACCEPT-PROBE-2026-08-21-v1.0','created_at_utc':datetime.now(timezone.utc).isoformat(),'steps':[]}

r=s.get(UIDS,timeout=45); r.raise_for_status(); rep['steps'].append({'stage':'uid_lookup',**snap(r)}); uids=[str(x) for x in r.json()]
if '2504732' not in uids: raise RuntimeError(f'UID drift: {uids}')

m=s.post(MODAL,data={'FileDownload':'2504732','data_set_uid':'30497'},timeout=60,allow_redirects=True); m.raise_for_status(); rep['steps'].append({'stage':'modal',**snap(m)})
soup=BeautifulSoup(m.text,'html.parser'); f=soup.find('form',id='data_link') or soup.find('form')
if not f: raise RuntimeError('download form missing')
action=urljoin(m.url,f.get('action') or '')
payload={'purpose':'Research','client':'DataLink','force_download':'1','data_uids':'2504732'}
rep['accept_contract']={'action':action,'payload':payload}

# First capture the raw server decision without auto-following redirects.
a=s.post(action,data=payload,timeout=120,allow_redirects=False); rep['steps'].append({'stage':'accept_no_redirect',**snap(a)})

# Then reproduce browser-like redirect following using the same session.
if a.is_redirect or a.is_permanent_redirect:
    location=urljoin(a.url,a.headers.get('location',''))
    rr=s.get(location,timeout=180,allow_redirects=True); rep['steps'].append({'stage':'accept_redirect_follow',**snap(rr)})
else:
    # Repeat once with redirects enabled to surface any method-specific behavior.
    rr=s.post(action,data=payload,timeout=180,allow_redirects=True); rep['steps'].append({'stage':'accept_follow_enabled',**snap(rr)})

# Parse returned HTML/JSON for endpoint hints, but do not download new content in this diagnostic.
for st in rep['steps'][-2:]:
    body=st.get('body_prefix') or ''
    if body:
        ss=BeautifulSoup(body,'html.parser')
        st['links']=[urljoin(st['final_url'],x.get('href')) for x in ss.find_all('a') if x.get('href')][:100]
        st['forms']=[{'action':urljoin(st['final_url'],x.get('action') or ''),'method':x.get('method'),'html':str(x)[:12000]} for x in ss.find_all('form')][:20]

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(rep,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(json.dumps({'uid':uids,'action':action,'last_steps':[{'stage':x['stage'],'status':x['status'],'bytes':x['bytes'],'content_type':x['headers'].get('Content-Type') or x['headers'].get('content-type'),'location':x['headers'].get('Location') or x['headers'].get('location'),'magic':x['magic_hex']} for x in rep['steps'][-2:]]},indent=2))
