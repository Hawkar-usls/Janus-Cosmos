#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, hashlib, json, re
from datetime import datetime, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk'
BASE='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI'
README=BASE+'/sd11282/README.doc'

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True,type=Path); a=ap.parse_args()
    r={'schema':'janus.cosmos.cousteau.hannah_bodc.tobi_format_probe.v1','created_at_utc':datetime.now(timezone.utc).isoformat(),'status':'STARTED','source_readme':README,'raw_scientific_data_downloaded':False}
    ftp=ftplib.FTP(timeout=25)
    try:
        ftp.connect(HOST,21); ftp.login('anonymous','janus-probe@example.invalid')
        chunks=[]; ftp.retrbinary('RETR '+README,chunks.append); raw=b''.join(chunks)
        text=raw.decode('utf-8',errors='replace'); lines=text.splitlines()
        # Format-only excerpt: C struct declarations, offsets and explanatory lines.
        selected=[]
        for i,line in enumerate(lines,1):
            if 15 <= i <= 145 or 155 <= i <= 190:
                selected.append({'line':i,'text':line.rstrip()[:500]})
        r['readme_sha256']=hashlib.sha256(raw).hexdigest(); r['readme_size_bytes']=len(raw); r['format_lines']=selected
        r['detected_block_size_candidates']=[40960]
        r['status']='REAL_FORMAT_METADATA_READY'
    except Exception as e:
        r['status']='FORMAT_PROBE_FAILED'; r['error_type']=type(e).__name__; r['error']=str(e)
    finally:
        try: ftp.quit()
        except Exception: pass
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':r['status'],'readme_sha256':r.get('readme_sha256'),'line_count':len(r.get('format_lines',[]))},indent=2))
    return 0 if r['status']=='REAL_FORMAT_METADATA_READY' else 2
if __name__=='__main__': raise SystemExit(main())
