#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, json, posixpath, re
from datetime import datetime, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk'
ROOT='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI'
PAT=re.compile(r'(jd0?59|(^|[^0-9])059([^0-9]|$)|line[_ -]?19|tobi[_ -]?0?2|2005[-_]?0?2[-_]?28|readme|\.log$|\.txt$)', re.I)
MAX=12000

def listdir(ftp,path):
    try:
        return [(n,f) for n,f in ftp.mlsd(path, facts=['type','size','modify']) if n not in {'.','..'}]
    except Exception:
        out=[]
        for raw in ftp.nlst(path):
            n=posixpath.basename(raw.rstrip('/')); p=raw if raw.startswith('/') else posixpath.join(path,n)
            cur=ftp.pwd(); typ='file'; size=None
            try: ftp.cwd(p); typ='dir'
            except Exception:
                try: size=ftp.size(p)
                except Exception: pass
            finally:
                try: ftp.cwd(cur)
                except Exception: pass
            out.append((n,{'type':typ,'size':str(size) if size is not None else None}))
        return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True,type=Path); a=ap.parse_args()
    r={'schema':'janus.cosmos.cousteau.hannah_bodc.cd169_target_scout.v1','created_at_utc':datetime.now(timezone.utc).isoformat(),'root':ROOT,'matches':[],'directories_seen':0,'files_seen':0,'status':'STARTED','scientific_claim':False}
    ftp=ftplib.FTP(timeout=25)
    try:
        ftp.connect(HOST,21,timeout=25); ftp.login('anonymous','janus-probe@example.invalid'); ftp.cwd(ROOT); base=ftp.pwd()
        stack=[base]; seen=set(); count=0
        while stack and count<MAX:
            d=stack.pop()
            if d in seen: continue
            seen.add(d); r['directories_seen']+=1
            try: children=listdir(ftp,d)
            except Exception as e:
                r.setdefault('errors',[]).append({'path':d,'error':str(e)}); continue
            # Prioritize likely directories but traverse all.
            children.sort(key=lambda x:(0 if PAT.search(x[0]) else 1,x[0].lower()))
            for n,f in children:
                p=posixpath.join(d.rstrip('/'),n); rel=posixpath.relpath(p,base); typ=f.get('type','unknown')
                count+=1
                if typ=='dir': stack.append(p)
                elif typ=='file': r['files_seen']+=1
                if PAT.search(rel):
                    try: sz=int(f.get('size')) if f.get('size') not in {None,'None'} else None
                    except Exception: sz=None
                    r['matches'].append({'relative_path':rel,'type':typ,'size_bytes':sz,'modify':f.get('modify')})
        r['truncated']=count>=MAX; r['status']='TARGET_SCOUT_READY'
    except Exception as e:
        r['status']='TARGET_SCOUT_FAILED'; r['error']=str(e); r['error_type']=type(e).__name__
    finally:
        try: ftp.quit()
        except Exception: pass
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':r['status'],'directories_seen':r['directories_seen'],'files_seen':r['files_seen'],'match_count':len(r['matches'])},indent=2))
    return 0 if r['status']=='TARGET_SCOUT_READY' else 2
if __name__=='__main__': raise SystemExit(main())
