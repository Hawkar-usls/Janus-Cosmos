#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, os, statistics, struct, tempfile
from pathlib import Path
HOST='livftp.noc.ac.uk'; RAW='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11285/TOBI.DAT'; CDF='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI/cd169p5.cdf'; BLOCK=40960; TARGET=9812; TARGET_HASH='039b839861e515efe339542b8c57e5830bd484def961804f9c1fb61387132c3f'; PORT=0x0240; STBD=0x2180

def ftp():
 f=ftplib.FTP(timeout=90);f.connect(HOST,21);f.login('anonymous','janus-probe@example.invalid');f.voidcmd('TYPE I');return f

def raw_target():
 f=ftp();buf=bytearray();idx=0
 try:
  s=f.transfercmd('RETR '+RAW)
  try:
   while idx<=TARGET:
    b=s.recv(1048576)
    if not b:break
    buf.extend(b)
    while len(buf)>=BLOCK:
     q=bytes(buf[:BLOCK]);del buf[:BLOCK]
     if idx==TARGET:return q
     idx+=1
  finally:
   try:s.close()
   except:pass
 finally:
  try:f.close()
  except:pass
 raise RuntimeError('raw target not reached')

def download(remote,local):
 f=ftp();h=hashlib.sha256();n=0
 try:
  with open(local,'wb') as o:
   def cb(b):
    nonlocal n;o.write(b);h.update(b);n+=len(b)
   f.retrbinary('RETR '+remote,cb,1048576)
 finally:
  try:f.quit()
  except:
   try:f.close()
   except:pass
 return n,h.hexdigest()

def transform(x,mode):
 if mode=='SIGNED':return list(map(float,x))
 if mode=='ABS':return [float(abs(v)) for v in x]
 if mode=='POSITIVE_ONLY':return [float(max(0,v)) for v in x]
 raise ValueError(mode)

def reduce8(x,kind):
 out=[]
 for i in range(0,4000,8):
  a=x[i:i+8]
  if kind=='MEAN':v=sum(a)/8
  elif kind=='MEDIAN':v=statistics.median(a)
  elif kind=='SUM':v=sum(a)
  elif kind=='RMS':v=math.sqrt(sum(q*q for q in a)/8)
  elif kind=='MAX':v=max(a)
  elif kind=='MIN':v=min(a)
  else:raise ValueError(kind)
  out.append(float(v))
 return out

def pearson(a,b):
 ma=sum(a)/len(a);mb=sum(b)/len(b);da=[x-ma for x in a];db=[x-mb for x in b];den=math.sqrt(sum(x*x for x in da)*sum(y*y for y in db));return None if den==0 else sum(x*y for x,y in zip(da,db))/den

def ranks(x):
 order=sorted(range(len(x)),key=lambda i:x[i]);r=[0.]*len(x);k=0
 while k<len(order):
  j=k+1
  while j<len(order) and x[order[j]]==x[order[k]]:j+=1
  v=(k+j-1)/2
  for z in range(k,j):r[order[z]]=v
  k=j
 return r

def spear(a,b):return pearson(ranks(a),ranks(b))
def orient(x,rev):return list(reversed(x)) if rev else x

def add(scores,cdf,name,p,s,order,pr,sr):
 p=orient(p,pr);s=orient(s,sr);q=(p+s) if order=='PORT_THEN_STBD' else (s+p);scores.append({'candidate':name,'side_order':order,'port_reversed':pr,'stbd_reversed':sr,'spearman':spear(cdf,q),'pearson':pearson(cdf,q)})

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();out={'schema':'janus.cosmos.cousteau.cd169.factor8_prism_mapping_test.v1','status':'STARTED','purpose':'FORMAT_MAPPING_ONLY','scientific_claim':False}
 try:
  raw=raw_target();rh=hashlib.sha256(raw).hexdigest();out['raw_target_sha256']=rh
  if rh!=TARGET_HASH:raise RuntimeError('raw hash mismatch')
  po=list(struct.unpack_from('<4000h',raw,PORT));st=list(struct.unpack_from('<4000h',raw,STBD))
  with tempfile.TemporaryDirectory() as td:
   lp=os.path.join(td,'p5.cdf');sz,sha=download(CDF,lp);out['cdf']={'size_bytes':sz,'sha256':sha}
   import netCDF4,numpy as np
   ds=netCDF4.Dataset(lp)
   try:
    cdf=np.asarray(ds.variables['image'][TARGET,:],dtype=float).reshape(-1).tolist();out['cdf_image_length']=len(cdf)
    if len(cdf)!=1000:raise RuntimeError('expected ns=1000')
   finally:ds.close()
  scores=[]
  for mode in ['SIGNED','ABS','POSITIVE_ONLY']:
   p0=transform(po,mode);s0=transform(st,mode)
   for red in ['MEAN','MEDIAN','SUM','RMS','MAX','MIN']:
    p=reduce8(p0,red);s=reduce8(s0,red)
    for order in ['PORT_THEN_STBD','STBD_THEN_PORT']:
     for pr in [False,True]:
      for sr in [False,True]:add(scores,cdf,f'BIN8_{red}_{mode}',p,s,order,pr,sr)
   for pp in range(8):
    p=p0[pp::8][:500]
    for sp in range(8):
     s=s0[sp::8][:500]
     if len(p)!=500 or len(s)!=500:continue
     for order in ['PORT_THEN_STBD','STBD_THEN_PORT']:
      for pr in [False,True]:
       for sr in [False,True]:add(scores,cdf,f'DECIMATE8_{mode}_PORTPH{pp}_STBDPH{sp}',p,s,order,pr,sr)
  scores.sort(key=lambda r:(-abs(r['spearman'] or 0),-abs(r['pearson'] or 0),r['candidate'],r['side_order'],r['port_reversed'],r['stbd_reversed']))
  out['best']=scores[0];out['top50']=scores[:50];out['tested_candidate_count']=len(scores);out['status']='FACTOR8_FORMAT_MAPPING_TEST_READY'
 except Exception as e:out['status']='FACTOR8_FORMAT_MAPPING_TEST_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps({'status':out['status'],'best':out.get('best'),'tested':out.get('tested_candidate_count')},indent=2));return 0 if out['status']=='FACTOR8_FORMAT_MAPPING_TEST_READY' else 2
if __name__=='__main__':raise SystemExit(main())
