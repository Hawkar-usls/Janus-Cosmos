#!/usr/bin/env python3
"""Search BODCREQ-9406 for authentic CD169 PRISM navigation artifacts.

Reads only archive metadata and small non-sonar files. Raw TOBI.DAT, CDF imagery
arrays, EM12 XYZ bodies and large archives are not inspected here.
"""
from __future__ import annotations
import argparse, ftplib, hashlib, json, posixpath, re
from datetime import datetime, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk'; ROOT='/bodc/bodc/data/BODCREQ-9406'; MAX_DEPTH=7; MAX_ENTRIES=10000; MAX_READ=8_000_000
PATTERNS=[
  r'\.veh_nav\b',r'navfile\.veh_nav',r'\bwireout\b',r'\bmrgnav_inertia\b',r'commands\.cfg',
  r'\bPRISM\b',r'\bumbilical\b',r'\binertia\b',r'\bdrag\b',r'\bviscosity\b',r'caten',
  r'\bvehicle navigation\b',r'\btowfish navigation\b'
]
SKIP_BODY_EXT={'.dat','.cdf','.ascii','.zip','.gz','.tgz','.tar','.img','.all','.ix1','.ix2'}

def ftp():
    f=ftplib.FTP(timeout=90);f.connect(HOST,21);f.login('anonymous','janus-probe@example.invalid');f.voidcmd('TYPE I');return f

def nlst(f,path):
    out=[]
    for n in f.nlst(path):
        if n in {'.','..',path,path+'/'}:continue
        if not n.startswith('/'):n=posixpath.join(path,n)
        out.append(posixpath.normpath(n))
    return sorted(set(out))

def isdir(f,path):
    old=f.pwd()
    try:f.cwd(path);return True
    except Exception:return False
    finally:
        try:f.cwd(old)
        except Exception:pass

def size(f,path):
    try:
        x=f.size(path);return int(x) if x is not None else None
    except Exception:return None

def ext(path):
    b=posixpath.basename(path).lower();i=b.rfind('.');return b[i:] if i>=0 else ''

def get(path):
    f=ftp();chunks=[]
    try:f.retrbinary('RETR '+path,chunks.append,262144)
    finally:
        try:f.quit()
        except Exception:
            try:f.close()
            except Exception:pass
    return b''.join(chunks)

def extract_text(raw:bytes)->str:
    # Combine permissive decoding with printable ASCII/UTF-16LE string extraction.
    parts=[raw.decode('latin-1','replace')]
    ascii_runs=re.findall(rb'[\x20-\x7e\t\r\n]{4,}',raw)
    parts.append('\n'.join(x.decode('ascii','replace') for x in ascii_runs))
    try:parts.append(raw.decode('utf-16le','ignore'))
    except Exception:pass
    return '\n'.join(parts)

def hits(text):
    found=[]
    for p in PATTERNS:
        for m in re.finditer(p,text,re.I):
            lo=max(0,m.start()-180);hi=min(len(text),m.end()+300)
            found.append({'pattern':p,'context':re.sub(r'\s+',' ',text[lo:hi]).strip()[:700]})
            if len(found)>=100:return found
    return found

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    out={'schema':'janus.cosmos.cousteau.hannah_cd169.native_veh_nav_archive_scout.v3','created_at_utc':datetime.now(timezone.utc).isoformat(),'root':ROOT,'status':'STARTED','sonar_body_inspected':False,'image_data_inspected':False,'entries':[],'content_hits':[],'name_hits':[]}
    f=ftp()
    try:
        stack=[(ROOT,0)];seen=set()
        while stack and len(out['entries'])<MAX_ENTRIES:
            path,d=stack.pop()
            if path in seen or d>MAX_DEPTH:continue
            seen.add(path)
            try:children=nlst(f,path)
            except Exception as e:
                out['entries'].append({'path':path,'type':'dir','list_error':f'{type(e).__name__}: {e}'});continue
            for ch in children:
                if len(out['entries'])>=MAX_ENTRIES:break
                if isdir(f,ch):
                    row={'relative_path':posixpath.relpath(ch,ROOT),'type':'dir'};out['entries'].append(row);stack.append((ch,d+1))
                else:
                    sz=size(f,ch);rel=posixpath.relpath(ch,ROOT);row={'relative_path':rel,'type':'file','size_bytes':sz};out['entries'].append(row)
                    nh=hits(rel)
                    if nh:out['name_hits'].append({'relative_path':rel,'hits':nh})
                    if (sz is not None and sz<=MAX_READ and ext(ch) not in SKIP_BODY_EXT):
                        try:
                            raw=get(ch);txt=extract_text(raw);hh=hits(txt)
                            if hh:out['content_hits'].append({'relative_path':rel,'size_bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'hits':hh})
                        except Exception as e:
                            row['read_error']=f'{type(e).__name__}: {e}'
        out['summary']={
          'entry_count':len(out['entries']),
          'name_hit_files':len(out['name_hits']),
          'content_hit_files':len(out['content_hits']),
          'native_veh_nav_name_found':any('.veh_nav' in x['relative_path'].lower() for x in out['name_hits']),
          'commands_cfg_name_found':any('commands.cfg' in x['relative_path'].lower() for x in out['name_hits']),
          'wireout_text_found':any(any(h['pattern']=='\\bwireout\\b' for h in x['hits']) for x in out['content_hits']),
          'authority':'ARCHIVE_EXISTENCE_AND_SMALL_FILE_TEXT_SEARCH_ONLY'
        }
        out['status']='ARCHIVE_NATIVE_NAV_SCOUT_READY'
    except Exception as e:
        out['status']='ARCHIVE_NATIVE_NAV_SCOUT_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
    finally:
        try:f.quit()
        except Exception:
            try:f.close()
            except Exception:pass
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'summary':out.get('summary'),'name_hits':[x['relative_path'] for x in out.get('name_hits',[])],'content_hits':[x['relative_path'] for x in out.get('content_hits',[])]},indent=2))
    return 0 if out['status']=='ARCHIVE_NATIVE_NAV_SCOUT_READY' else 2
if __name__=='__main__':raise SystemExit(main())
