#!/usr/bin/env python3
from __future__ import annotations
import hashlib, io, json, zipfile
from datetime import datetime, timezone
from pathlib import Path
import requests

OUT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-AUTHOR-MAR-DATASET-S1-ACQUISITION-PROBE-2026-08-21-v1.0.json')
DOI='10.1029/2022JB024008'
FILENAME='2022JB024008-sup-0003-Data Set SI-S01.zip'
URLS=[
 f'https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2022JB024008&file=2022JB024008-sup-0003-Data+Set+SI-S01.zip',
 f'https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2022JB024008&file=2022JB024008-sup-0003-Data+Set+SI-S01.zip',
 f'https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2022JB024008&file=2022JB024008-sup-0003-Data%20Set%20SI-S01.zip',
]

def sha(b:bytes): return hashlib.sha256(b).hexdigest()

def main():
 s=requests.Session(); s.headers['User-Agent']='Mozilla/5.0 Janus-Echo-Cousteau/1.0 independent MAR Data Set S1 replication'
 s.headers['Referer']='https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022JB024008'
 rep={'artifact_id':'JANUS-ECHO-COUSTEAU-AUTHOR-MAR-DATASET-S1-ACQUISITION-PROBE-2026-08-21-v1.0','created_utc':datetime.now(timezone.utc).isoformat(),'paper_doi':DOI,'expected_filename':FILENAME,'responses':[],'status':'BLOCKED'}
 for u in URLS:
  try:
   r=s.get(u,timeout=90,allow_redirects=True)
   item={'requested_url':u,'final_url':r.url,'status':r.status_code,'content_type':r.headers.get('content-type'),'content_disposition':r.headers.get('content-disposition'),'bytes':len(r.content),'sha256':sha(r.content)}
   rep['responses'].append(item)
   if r.status_code!=200 or len(r.content)<1000: continue
   try:
    z=zipfile.ZipFile(io.BytesIO(r.content))
   except Exception as e:
    item['zip_error']=f'{type(e).__name__}: {e}'; item['prefix_hex']=r.content[:32].hex(); continue
   members=[]
   for n in z.namelist():
    raw=z.read(n)
    m={'name':n,'bytes':len(raw),'sha256':sha(raw)}
    if len(raw)<2_000_000:
     try:
      txt=raw.decode('utf-8',errors='replace')
      lines=txt.splitlines()
      m['line_count']=len(lines)
      m['first_30_lines']=lines[:30]
      m['last_5_lines']=lines[-5:]
     except Exception: pass
    members.append(m)
   item['zip_members']=members
   rep['status']='SUCCESS'
   rep['selected_download']={'url':r.url,'bytes':len(r.content),'sha256':sha(r.content),'members':[m['name'] for m in members]}
   break
 OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(rep,indent=2,ensure_ascii=False),encoding='utf-8')
 print(json.dumps({'status':rep['status'],'selected':rep.get('selected_download')},indent=2))
 return 0 if rep['status']=='SUCCESS' else 2

if __name__=='__main__': raise SystemExit(main())
