#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, os, statistics, struct, tempfile
from pathlib import Path

HOST='livftp.noc.ac.uk'
RAW='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11285/TOBI.DAT'
CDF='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI/cd169p5.cdf'
BLOCK=40960; TARGET=9812
RAW_HASH='039b839861e515efe339542b8c57e5830bd484def961804f9c1fb61387132c3f'
OFF={'port':0x0240,'stbd':0x2180}

def ftp():
    f=ftplib.FTP(timeout=90); f.connect(HOST,21); f.login('anonymous','janus-probe@example.invalid'); f.voidcmd('TYPE I'); return f

def raw_target():
    f=ftp(); buf=bytearray(); idx=0
    try:
      s=f.transfercmd('RETR '+RAW)
      try:
        while idx<=TARGET:
          b=s.recv(1048576)
          if not b:break
          buf.extend(b)
          while len(buf)>=BLOCK:
            block=bytes(buf[:BLOCK]); del buf[:BLOCK]
            if idx==TARGET:return block
            idx+=1
      finally:
        try:s.close()
        except Exception:pass
    finally:
      try:f.close()
      except Exception:pass
    raise RuntimeError('target raw block not reached')

def download(remote,local):
    f=ftp(); h=hashlib.sha256(); n=0
    try:
      with open(local,'wb') as o:
        def cb(b):
          nonlocal n;o.write(b);h.update(b);n+=len(b)
        f.retrbinary('RETR '+remote,cb,blocksize=1048576)
    finally:
      try:f.quit()
      except Exception:
        try:f.close()
        except Exception:pass
    return n,h.hexdigest()

def pearson(a,b):
    if len(a)!=len(b) or len(a)<3:return None
    ma=sum(a)/len(a);mb=sum(b)/len(b);da=[x-ma for x in a];db=[x-mb for x in b]
    den=math.sqrt(sum(x*x for x in da)*sum(x*x for x in db))
    return None if den==0 else sum(x*y for x,y in zip(da,db))/den

def ranks(x):
    # average ranks for ties
    order=sorted(range(len(x)),key=lambda i:x[i]); r=[0.0]*len(x);k=0
    while k<len(order):
      j=k+1
      while j<len(order) and x[order[j]]==x[order[k]]:j+=1
      avg=(k+j-1)/2
      for q in range(k,j):r[order[q]]=avg
      k=j
    return r

def spearman(a,b):return pearson(ranks(a),ranks(b))
def norm_corr(cdf,cand):
    variants={
      'raw_signed':cand,
      'raw_abs':[abs(x) for x in cand],
      'raw_positive':[max(0,x) for x in cand],
      'raw_negative_abs':[abs(min(0,x)) for x in cand],
    }
    out=[]
    for vn,v in variants.items():
      for rev in (False,True):
        q=list(reversed(v)) if rev else v
        out.append({'transform':vn,'reversed':rev,'pearson':pearson(cdf,q),'spearman':spearman(cdf,q)})
    return out

def candidate_sequences(port,stbd):
    c={}
    # whole side, every 4 => 1000
    for side,x in [('port',port),('stbd',stbd)]:
      for phase in range(4):c[f'{side}_every4_phase{phase}']=x[phase::4][:1000]
    # 500 samples from each side: first/last 2000 decimated by 4, all orientation combinations handled by correlation reverse only globally.
    for pseg,pv in [('first2000',port[:2000]),('last2000',port[2000:])]:
      for sseg,sv in [('first2000',stbd[:2000]),('last2000',stbd[2000:])]:
        for pp in range(4):
          p=pv[pp::4][:500]
          for sp in range(4):
            s=sv[sp::4][:500]
            if len(p)==500 and len(s)==500:
              c[f'port_{pseg}_p{pp}__stbd_{sseg}_p{sp}']=p+s
              c[f'stbd_{sseg}_p{sp}__port_{pseg}_p{pp}']=s+p
    return c

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    out={'schema':'janus.cosmos.cousteau.cd169_raw_to_prism_target_crosscheck.v1','status':'STARTED','target_line':TARGET,'morphology_classification':False,'scientific_target_claim':False}
    try:
      raw=raw_target();rh=hashlib.sha256(raw).hexdigest();out['raw_target_sha256']=rh
      if rh!=RAW_HASH:raise RuntimeError('raw hash mismatch')
      port=list(struct.unpack_from('<4000h',raw,OFF['port']));stbd=list(struct.unpack_from('<4000h',raw,OFF['stbd']))
      with tempfile.TemporaryDirectory() as td:
        lp=os.path.join(td,'p5.cdf');sz,sha=download(CDF,lp);out['cdf_integrity']={'size_bytes':sz,'sha256':sha}
        import netCDF4,numpy as np
        ds=netCDF4.Dataset(lp)
        try:
          if len(ds.dimensions['nl'])!=14536 or len(ds.dimensions['ns'])!=1000:raise RuntimeError('unexpected p5 dimensions')
          img=np.asarray(ds.variables['image'][TARGET,:],dtype=float).reshape(-1).tolist()
          date=np.asarray(ds.variables['date'][TARGET,:]).reshape(-1).tolist();time=np.asarray(ds.variables['time'][TARGET,:]).reshape(-1).tolist();sec=float(np.asarray(ds.variables['seconds'][TARGET]).reshape(-1)[0])
          ssa=np.asarray(ds.variables['ss_attributes'][TARGET,:],dtype=float).reshape(-1).tolist();sst=np.asarray(ds.variables['ss_attributes_tobi'][TARGET,:],dtype=float).reshape(-1).tolist();ship=np.asarray(ds.variables['ship_latlon'][TARGET,:],dtype=float).reshape(-1).tolist();sonar=np.asarray(ds.variables['latlon'][TARGET,:],dtype=float).reshape(-1).tolist();depth=np.asarray(ds.variables['depths'][TARGET,:],dtype=float).reshape(-1).tolist()
          out['cdf_target_header']={'date':date,'time':time,'seconds':sec,'ss_attributes_fish_altitude_heading_pitch_roll_yaw':ssa,'ss_attributes_tobi_water_pressure_cable_length':sst,'ship_latlon':ship,'sonar_latlon':sonar,'depths_corrected_uncorrected':depth,'image_length':len(img),'image_sha256_float64_le':hashlib.sha256(np.asarray(img,dtype='<f8').tobytes()).hexdigest()}
          scores=[]
          for name,cand in candidate_sequences(port,stbd).items():
            if len(cand)!=len(img):continue
            for q in norm_corr(img,cand):scores.append({'candidate':name,**q})
          def key(r):
            s=r['spearman'];p=r['pearson'];return -(max(abs(s) if s is not None else -1,abs(p) if p is not None else -1))
          scores.sort(key=key)
          out['top_raw_to_cdf_correspondence_candidates']=scores[:40]
          out['best_candidate']=scores[0] if scores else None
          out['cdf_global_attributes']={k:str(ds.getncattr(k))[:500] for k in ds.ncattrs()}
        finally:ds.close()
      out['status']='RAW_TO_PRISM_TARGET_CROSSCHECK_READY'
    except Exception as e:out['status']='RAW_TO_PRISM_TARGET_CROSSCHECK_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'header':out.get('cdf_target_header'),'best':out.get('best_candidate')},indent=2))
    return 0 if out['status']=='RAW_TO_PRISM_TARGET_CROSSCHECK_READY' else 2
if __name__=='__main__':raise SystemExit(main())
