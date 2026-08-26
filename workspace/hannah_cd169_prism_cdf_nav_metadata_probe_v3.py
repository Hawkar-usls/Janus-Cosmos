#!/usr/bin/env python3
"""V3 metadata/navigation-only probe for archived CD169 PRISM CDF files.

Purpose: recover authentic processed navigation provenance before any new model.
This program never ranks or inspects sonar intensity/image variables. It inventories
NetCDF metadata and candidate navigation/header variables only.
"""
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, os, re, tempfile
from datetime import datetime, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk'
BASE='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI'
FILES=[f'cd169p{i}.cdf' for i in range(1,12)]
TARGET='2005-02-28T01:07:25Z'

NAV_KEYS=('lat','lon','nav','position','pos','east','north','ship','veh','vehicle','tow','fish','time','date','year','month','day','hour','minute','second','sec','head','gyro','cable','wire','depth','press','alt')
IMAGE_KEYS=('image','imagery','sidescan','side_scan','port','stbd','starboard','amplitude','backscatter','pixel','sample','sonar','swath_data')

def ftp():
    f=ftplib.FTP(timeout=90); f.connect(HOST,21); f.login('anonymous','janus-probe@example.invalid'); f.voidcmd('TYPE I'); return f

def download(path: str, dest: str) -> dict:
    h=hashlib.sha256(); n=0; f=ftp()
    try:
        with open(dest,'wb') as out:
            def cb(b:bytes):
                nonlocal n; out.write(b); h.update(b); n+=len(b)
            f.retrbinary('RETR '+path,cb,blocksize=1024*1024)
    finally:
        try:f.quit()
        except Exception:
            try:f.close()
            except Exception:pass
    return {'size_bytes':n,'sha256':h.hexdigest()}

def clean_scalar(x):
    try:
        if hasattr(x,'item'): x=x.item()
    except Exception: pass
    if isinstance(x,(bytes,bytearray)):
        return x.decode('utf-8','replace')[:500]
    if isinstance(x,float):
        if not math.isfinite(x): return str(x)
        return float(x)
    if isinstance(x,(int,str,bool)) or x is None:return x
    return str(x)[:500]

def sample_variable(v, max_points=7):
    shape=tuple(int(x) for x in getattr(v,'shape',()))
    size=1
    for x in shape:size*=x
    out={'shape':list(shape),'size':size,'dtype':str(getattr(v,'dtype','unknown'))}
    attrs={}
    for a in getattr(v,'ncattrs',lambda:[])():
        try: attrs[a]=clean_scalar(getattr(v,a))
        except Exception as e: attrs[a]=f'<attr_error:{type(e).__name__}>'
    out['attrs']=attrs
    if size<=0:return out
    vals=[]
    try:
        if len(shape)==0:
            vals=[clean_scalar(v[...])]
        elif len(shape)==1:
            idxs=sorted(set([0,max(0,shape[0]//4),max(0,shape[0]//2),max(0,(3*shape[0])//4),shape[0]-1]))
            vals=[{'index':i,'value':clean_scalar(v[i])} for i in idxs]
        elif len(shape)==2 and min(shape)==1:
            n=max(shape); idxs=sorted(set([0,n//2,n-1]));
            vals=[]
            for i in idxs:
                val=v[0,i] if shape[0]==1 else v[i,0]
                vals.append({'linear_index':i,'value':clean_scalar(val)})
        else:
            out['values_not_read_reason']='multidimensional_or_non_header_candidate'
    except Exception as e:
        out['sample_error']=f'{type(e).__name__}: {e}'
    if vals: out['sample_values']=vals
    # For one-dimensional navigation/header candidates only, safe numeric summary.
    if len(shape)==1 and size<=2_000_000:
        try:
            import numpy as np
            a=np.ma.asarray(v[:])
            if np.issubdtype(a.dtype,np.number):
                c=a.compressed() if np.ma.isMaskedArray(a) else np.asarray(a).ravel()
                c=c[np.isfinite(c)]
                if c.size:
                    out['numeric_summary']={'n_finite':int(c.size),'min':float(np.min(c)),'max':float(np.max(c))}
        except Exception as e:
            out['numeric_summary_error']=f'{type(e).__name__}: {e}'
    return out

def is_nav_candidate(name, v):
    s=' '.join([name]+[str(getattr(v,a,'')) for a in getattr(v,'ncattrs',lambda:[])()]).lower()
    return any(k in s for k in NAV_KEYS) and not any(k in name.lower() for k in IMAGE_KEYS)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True,type=Path); a=ap.parse_args()
    out={
      'schema':'janus.cosmos.cousteau.hannah_cd169.prism_cdf_nav_metadata_probe.v3',
      'target_utc':TARGET,
      'status':'STARTED',
      'contract':'JANUS-HANNAH-CD169-NATIVE-OR-CALIBRATED-VEH-NAV-RECOVERY-V3-CONTRACT-2026-08-26-v1.0',
      'sonar_intensity_inspected':False,
      'image_variables_read':False,
      'navigation_model_fitted':False,
      'files':[]
    }
    try:
        from netCDF4 import Dataset
        with tempfile.TemporaryDirectory() as td:
            for fn in FILES:
                remote=BASE+'/'+fn; local=os.path.join(td,fn)
                row={'file':fn,'remote_path':remote,'download':None,'open_status':'STARTED'}
                try:
                    row['download']=download(remote,local)
                    ds=Dataset(local,'r')
                    try:
                        row['data_model']=getattr(ds,'data_model',None)
                        row['global_attributes']={}
                        for att in ds.ncattrs():
                            try:row['global_attributes'][att]=clean_scalar(getattr(ds,att))
                            except Exception as e:row['global_attributes'][att]=f'<attr_error:{type(e).__name__}>'
                        row['dimensions']={k:{'size':len(v),'unlimited':bool(v.isunlimited())} for k,v in ds.dimensions.items()}
                        row['variable_inventory']=[]; row['navigation_candidates']={}
                        for name,v in ds.variables.items():
                            inv={'name':name,'dtype':str(v.dtype),'dimensions':list(v.dimensions),'shape':list(v.shape)}
                            attrs={}
                            for att in v.ncattrs():
                                try:attrs[att]=clean_scalar(getattr(v,att))
                                except Exception:pass
                            inv['attrs']=attrs
                            inv['navigation_candidate']=is_nav_candidate(name,v)
                            inv['image_like_by_name']=any(k in name.lower() for k in IMAGE_KEYS)
                            row['variable_inventory'].append(inv)
                            if inv['navigation_candidate']:
                                row['navigation_candidates'][name]=sample_variable(v)
                        row['open_status']='CDF_METADATA_READY'
                    finally: ds.close()
                except Exception as e:
                    row['open_status']='FAILED'; row['error_type']=type(e).__name__; row['error']=str(e)
                out['files'].append(row)
        ok=[r for r in out['files'] if r['open_status']=='CDF_METADATA_READY']
        out['summary']={
          'cdf_files_requested':len(FILES),
          'cdf_files_opened':len(ok),
          'candidate_variable_names':sorted(set(n for r in ok for n in r.get('navigation_candidates',{}))),
          'explicit_vehicle_named_variables':sorted(set(n for r in ok for n in r.get('navigation_candidates',{}) if any(k in n.lower() for k in ('veh','vehicle','fish','tow')))),
          'explicit_ship_named_variables':sorted(set(n for r in ok for n in r.get('navigation_candidates',{}) if 'ship' in n.lower())),
          'interpretation_authority':'METADATA_AND_HEADER_VARIABLE_DISCOVERY_ONLY'
        }
        out['status']='PRISM_CDF_NAV_METADATA_READY' if ok else 'NO_CDF_OPENED'
    except Exception as e:
        out['status']='PROBE_FAILED'; out['error_type']=type(e).__name__; out['error']=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'summary':out.get('summary'),'errors':[{'file':r['file'],'error':r.get('error')} for r in out.get('files',[]) if r.get('open_status')=='FAILED']},indent=2))
    return 0 if out['status']=='PRISM_CDF_NAV_METADATA_READY' else 2
if __name__=='__main__':raise SystemExit(main())
