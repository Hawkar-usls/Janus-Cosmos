#!/usr/bin/env python3
from __future__ import annotations

import hashlib, io, json, zipfile
from datetime import datetime, timezone
from pathlib import Path
import requests

OUT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-TURN5-MAR-SUPPORT-DATA-PROBE-2026-08-21-v1.0.json')
URL='https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2022JB024008&file=2022JB024008-sup-0003-Data+Set+SI-S01.zip'

def sha(b:bytes)->str: return hashlib.sha256(b).hexdigest()
def now(): return datetime.now(timezone.utc).isoformat()

def main():
    OUT.parent.mkdir(parents=True,exist_ok=True)
    s=requests.Session(); s.headers.update({'User-Agent':'Janus-Echo-Cousteau/Turn5 scientific data probe'})
    result={'artifact_id':'JANUS-ECHO-COUSTEAU-TURN5-MAR-SUPPORT-DATA-PROBE-2026-08-21-v1.0','created_utc':now(),'source_url':URL,'source_article_doi':'10.1029/2022JB024008','expected_filename':'2022JB024008-sup-0003-Data Set SI-S01.zip','status':'STARTED'}
    try:
        r=s.get(URL,timeout=120,allow_redirects=True)
        result['http']={'status':r.status_code,'final_url':r.url,'content_type':r.headers.get('content-type'),'content_disposition':r.headers.get('content-disposition'),'bytes':len(r.content),'sha256':sha(r.content)}
        r.raise_for_status()
        if not zipfile.is_zipfile(io.BytesIO(r.content)):
            result['status']='BLOCKED_NOT_ZIP'
            result['body_prefix']=r.content[:500].decode('utf-8',errors='replace')
        else:
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                members=[]
                for info in z.infolist():
                    item={'name':info.filename,'bytes':info.file_size,'compressed_bytes':info.compress_size}
                    if not info.is_dir():
                        b=z.read(info)
                        item['sha256']=sha(b)
                        if len(b)<=200000:
                            item['preview']=b[:1600].decode('utf-8',errors='replace')
                    members.append(item)
                result['zip_members']=members
            result['status']='SUCCESS_DATASET_MATERIALIZED'
    except Exception as e:
        result['status']='BLOCKED_DOWNLOAD_OR_PARSE'
        result['error_type']=type(e).__name__; result['error']=str(e)
    OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'http':result.get('http'),'members':[m['name'] for m in result.get('zip_members',[])]},indent=2,ensure_ascii=False))
    return 0

if __name__=='__main__': raise SystemExit(main())
