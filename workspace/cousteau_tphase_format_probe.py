#!/usr/bin/env python3
from __future__ import annotations

import gzip, io, json, tarfile
from datetime import datetime, timezone
from pathlib import Path
import requests
from bs4 import BeautifulSoup

OUT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-EA-TPHASE-INNER-CATALOG-FORMAT-PROBE-2026-08-21-v1.0.json')
UID='2504732'; DS='30497'
LANDING=f'https://www.marine-geo.org/tools/files/{DS}'
UIDS=f'https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid={DS}'
MODAL='https://www.marine-geo.org/services/download/download_modal.php'

s=requests.Session(); s.headers.update({'User-Agent':'Janus-Echo-Cousteau/1.0 format probe','Referer':LANDING})
rep={'artifact_id':'JANUS-ECHO-COUSTEAU-EA-TPHASE-INNER-CATALOG-FORMAT-PROBE-2026-08-21-v1.0','created_at_utc':datetime.now(timezone.utc).isoformat()}

u=s.get(UIDS,timeout=45); u.raise_for_status(); rep['uids']=u.json()
m=s.post(MODAL,data={'FileDownload':UID,'data_set_uid':DS},timeout=60); m.raise_for_status(); soup=BeautifulSoup(m.text,'html.parser'); f=soup.find('form',id='data_link') or soup.find('form'); action=f.get('action')
r=s.post(action,data={'purpose':'Research','client':'DataLink','force_download':'1','data_uids':UID},timeout=180); r.raise_for_status()
rep['download']={'status':r.status_code,'content_type':r.headers.get('content-type'),'content_disposition':r.headers.get('content-disposition'),'bytes':len(r.content),'magic_hex':r.content[:32].hex()}
rep['members']=[]
with tarfile.open(fileobj=io.BytesIO(r.content),mode='r:*') as tf:
    for member in tf.getmembers():
        if not member.isfile(): continue
        fobj=tf.extractfile(member)
        if fobj is None: continue
        b=fobj.read()
        row={'name':member.name,'bytes':len(b),'magic_hex':b[:24].hex()}
        if b.startswith(b'\x1f\x8b'):
            raw=gzip.decompress(b)
            text=raw.decode('utf-8',errors='replace')
            lines=text.splitlines()
            row.update({'nested':'gzip','uncompressed_bytes':len(raw),'line_count':len(lines),'first_40_lines':lines[:40],
                        'last_3_lines':lines[-3:] if len(lines)>=3 else lines,
                        'line_lengths_first_40':[len(x) for x in lines[:40]]})
        elif len(b)<2_000_000:
            text=b.decode('utf-8',errors='replace'); lines=text.splitlines()
            row.update({'nested':'plain','line_count':len(lines),'first_40_lines':lines[:40],'line_lengths_first_40':[len(x) for x in lines[:40]]})
        rep['members'].append(row)

OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(rep,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(json.dumps({'download':rep['download'],'members':[{'name':x['name'],'bytes':x['bytes'],'nested':x.get('nested'),'uncompressed_bytes':x.get('uncompressed_bytes'),'line_count':x.get('line_count'),'first_5_lines':x.get('first_40_lines',[])[:5]} for x in rep['members']]},indent=2,ensure_ascii=False))
