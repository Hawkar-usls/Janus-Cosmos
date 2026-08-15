#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
import urllib.parse
from pathlib import Path

import requests

ROOT=Path(__file__).resolve().parent
M=json.loads((ROOT/'SKY_MANIFEST_v2_0.json').read_text(encoding='utf-8'))
DATA=ROOT/'external_data'
PROV=DATA/'download_provenance_v2_0.json'
ERRORS=DATA/'download_errors_v2_0.json'
ENDPOINTS=[
    'https://alasky.cds.unistra.fr/hips-image-services/hips2fits',
    'https://alaskybis.cds.unistra.fr/hips-image-services/hips2fits',
]
SESSION=requests.Session()


def sha256_file(p:Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def fits_ok(p:Path):
    if not p.exists() or p.stat().st_size<2880:return False
    with p.open('rb') as f:head=f.read(80)
    return head.startswith(b'SIMPLE') or head.startswith(b'XTENSION')


def hips_url(endpoint,hips,ra,dec,fov,pixels):
    q={
        'hips':hips,'format':'fits','width':str(pixels),'height':str(pixels),
        'fov':str(fov),'projection':'TAN','coordsys':'icrs','rotation_angle':'0',
        'ra':f'{ra:.12f}','dec':f'{dec:.12f}'
    }
    return endpoint+'?'+urllib.parse.urlencode(q)


def download(urls,dst:Path,retries=3):
    dst.parent.mkdir(parents=True,exist_ok=True)
    if fits_ok(dst):
        return {'status':'cached','bytes':dst.stat().st_size,'sha256':sha256_file(dst),'url':'CACHE'}
    last=None
    for url in urls:
        for attempt in range(1,retries+1):
            tmp=dst.with_suffix(dst.suffix+'.part')
            try:
                with SESSION.get(url,stream=True,timeout=(20,240),headers={'User-Agent':'Janus-Cosmos-v2.0/astronomy-validation'}) as r:
                    r.raise_for_status()
                    with tmp.open('wb') as f:
                        for chunk in r.iter_content(1024*1024):
                            if chunk:f.write(chunk)
                if not fits_ok(tmp):
                    sample=tmp.read_bytes()[:200] if tmp.exists() else b''
                    raise RuntimeError('server response is not FITS: '+repr(sample))
                os.replace(tmp,dst)
                return {'status':'downloaded','bytes':dst.stat().st_size,'sha256':sha256_file(dst),'url':url}
            except Exception as e:
                last=e
                try:
                    if tmp.exists():tmp.unlink()
                except Exception:pass
                if attempt<retries:time.sleep(attempt*2)
    raise RuntimeError(f'download failed: {last}')


def plan():
    rows=[]
    o=M['orion'];c=o['center_j2000']
    for s in o['surveys']:
        rows.append({'kind':'ORION','id':f"ORION_{s['family']}_{s['band']}",'dst':DATA/'orion'/s['filename'],
                     'urls':[hips_url(ep,s['hips'],c['ra_deg'],c['dec_deg'],o['fov_deg'],o['pixels']) for ep in ENDPOINTS]})
    bc=M['blind_controls']
    for center in bc['centers']:
        for s in bc['surveys']:
            name=f"{center['id'].lower()}_{s['family'].lower()}_{s['band'].lower()}.fits".replace('2mass','tmass')
            rows.append({'kind':'CONTROL','id':f"{center['id']}_{s['family']}_{s['band']}",'dst':DATA/'controls'/name,
                         'urls':[hips_url(ep,s['hips'],center['ra_deg'],center['dec_deg'],bc['fov_deg'],bc['pixels']) for ep in ENDPOINTS],
                         'center':center,'survey':s})
    for s in M['ngc1425']['surveys']:
        rows.append({'kind':'NGC1425','id':f"NGC1425_{s['band']}",'dst':DATA/'ngc1425'/s['filename'],'urls':[s['url']]})
    return rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dry-run',action='store_true');ap.add_argument('--skip-ngc',action='store_true');a=ap.parse_args()
    rows=[];errors=[]
    for item in plan():
        if a.skip_ngc and item['kind']=='NGC1425':continue
        if a.dry_run:
            print(f"[DRY] {item['id']} -> {item['dst'].relative_to(ROOT)}")
            for u in item['urls']:print('      '+u)
            continue
        try:
            info=download(item['urls'],item['dst'])
            rec={'kind':item['kind'],'id':item['id'],'file':str(item['dst'].relative_to(ROOT)),**info}
            if 'center' in item:rec['center']=item['center'];rec['survey']=item['survey']
            rows.append(rec);print(f"[OK] {item['id']} {info['status']} {info['bytes']} bytes",flush=True)
        except Exception as e:
            errors.append({'kind':item['kind'],'id':item['id'],'error':f'{type(e).__name__}: {e}'})
            print(f"[ERROR] {item['id']}: {e}",flush=True)
    if a.dry_run:
        print(f'DRY-RUN PASS: {len(plan())} source products planned')
        return 0
    PROV.write_text(json.dumps({'schema':'janus.cosmos.download_provenance.v2.0','records':rows,'errors':errors},indent=2,ensure_ascii=False),encoding='utf-8')
    ERRORS.write_text(json.dumps(errors,indent=2,ensure_ascii=False),encoding='utf-8')
    print('DOWNLOAD PASS' if not errors else f'DOWNLOAD PARTIAL: {len(errors)} error(s)')
    return 0 if not errors else 2

if __name__=='__main__':raise SystemExit(main())
