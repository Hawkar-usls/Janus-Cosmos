#!/usr/bin/env python3
"""Metadata-only probe for CD169 PRISM NetCDF files.

Run only after S2 sample extraction. It downloads each archived cd169p*.cdf and
records dimensions, variable names/attributes, and time/navigation scalar ranges.
It intentionally does not emit imagery arrays or classify morphology.
"""
from __future__ import annotations
import argparse, ftplib, hashlib, json, os, re, tempfile
from pathlib import Path

HOST='livftp.noc.ac.uk'
BASE='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI'
FILES=[f'cd169p{i}.cdf' for i in range(1,12)]
INTEREST=('time','date','lat','lon','nav','head','gyro','alt','depth','press','cable','wire','res','pixel','range','sample','slant','ping')

def ftp():
    f=ftplib.FTP(timeout=90); f.connect(HOST,21); f.login('anonymous','janus-probe@example.invalid'); f.voidcmd('TYPE I'); return f

def download(remote,local):
    h=hashlib.sha256(); n=0; f=ftp()
    try:
        with open(local,'wb') as o:
            def cb(b):
                nonlocal n; o.write(b); h.update(b); n+=len(b)
            f.retrbinary('RETR '+remote,cb,blocksize=1048576)
    finally:
        try:f.quit()
        except Exception:
            try:f.close()
            except Exception:pass
    return n,h.hexdigest()

def clean(v):
    try:
        if hasattr(v,'item'):v=v.item()
    except Exception:pass
    if isinstance(v,(str,int,float,bool)) or v is None:return v
    try:return [clean(x) for x in list(v)[:20]]
    except Exception:return str(v)[:500]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True,type=Path); a=ap.parse_args()
    import netCDF4
    out={'schema':'janus.cosmos.cousteau.cd169_prism_cdf_metadata_probe.v1','status':'STARTED','imagery_arrays_emitted':False,'scientific_claim':False,'files':[]}
    try:
      with tempfile.TemporaryDirectory() as td:
        for fn in FILES:
          lp=os.path.join(td,fn); size,sha=download(BASE+'/'+fn,lp); row={'file':fn,'size_bytes':size,'sha256':sha}
          ds=netCDF4.Dataset(lp,'r')
          try:
            row['data_model']=ds.data_model; row['dimensions']={k:(None if d.isunlimited() else len(d)) for k,d in ds.dimensions.items()}; row['global_attributes']={k:clean(ds.getncattr(k)) for k in ds.ncattrs()}
            vars={}
            for name,v in ds.variables.items():
              q={'dimensions':list(v.dimensions),'shape':list(v.shape),'dtype':str(v.dtype),'attributes':{k:clean(v.getncattr(k)) for k in v.ncattrs()}}
              low=name.lower(); interesting=any(k in low for k in INTEREST) or any(any(k in str(x).lower() for k in INTEREST) for x in q['attributes'].values())
              if interesting and v.size and v.ndim<=2:
                try:
                  flat=v[:].reshape(-1); q['first_values']=[clean(x) for x in flat[:5]]; q['last_values']=[clean(x) for x in flat[-5:]]
                  if v.size<=200000 and getattr(v.dtype,'kind','') in 'iuf':
                    import numpy as np
                    arr=np.asarray(flat,dtype=float); arr=arr[np.isfinite(arr)]
                    if arr.size:q['numeric_min']=float(arr.min()); q['numeric_max']=float(arr.max())
                except Exception as e:q['value_probe_error']=f'{type(e).__name__}: {e}'
              vars[name]=q
            row['variables']=vars
          finally:ds.close()
          out['files'].append(row)
      out['status']='PRISM_CDF_METADATA_READY'
    except Exception as e:out['status']='PRISM_CDF_METADATA_FAILED'; out['error_type']=type(e).__name__; out['error']=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'files':[{ 'file':x['file'],'dimensions':x.get('dimensions'),'vars':list(x.get('variables',{}))[:30]} for x in out.get('files',[])]},indent=2))
    return 0 if out['status']=='PRISM_CDF_METADATA_READY' else 2
if __name__=='__main__':raise SystemExit(main())
