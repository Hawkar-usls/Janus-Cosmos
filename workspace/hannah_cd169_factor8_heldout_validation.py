#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, os, statistics, struct, tempfile
from pathlib import Path
HOST='livftp.noc.ac.uk'; RAW='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11285/TOBI.DAT'; CDF='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI/cd169p5.cdf'; BLOCK=40960; PORT=0x0240; STBD=0x2180; INDICES=[0,2422,4845,7267,9690,12112,14535]

def ftp():
 f=ftplib.FTP(timeout=90);f.connect(HOST,21);f.login('anonymous','janus-probe@example.invalid');f.voidcmd('TYPE I');return f

def capture_raw():
 wanted=set(INDICES); got={}; f=ftp(); buf=bytearray(); idx=0; h=hashlib.sha256(); size=0
 try:
  s=f.transfercmd('RETR '+RAW)
  try:
   while True:
    b=s.recv(1048576)
    if not b:break
    h.update(b);size+=len(b);buf.extend(b)
    while len(buf)>=BLOCK:
     q=bytes(buf[:BLOCK]);del buf[:BLOCK]
     if idx in wanted:got[idx]=q
     idx+=1
  finally:
   try:s.close()
   except:pass
 finally:
  try:f.close()
  except:pass
 return got,idx,size,h.hexdigest()

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

def map_line(block):
 p=list(struct.unpack_from('<4000h',block,PORT));s=list(struct.unpack_from('<4000h',block,STBD))
 def bins(x):return [sum(abs(v) for v in x[i:i+8])/8.0 for i in range(0,4000,8)]
 return list(reversed(bins(s)))+bins(p)

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
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();o={'schema':'janus.cosmos.cousteau.cd169.factor8_heldout_validation.v1','status':'STARTED','indices':INDICES,'mapping_retuned':False,'scientific_claim':False}
 try:
  blocks,count,size,rsha=capture_raw();o['raw']={'record_count':count,'size_bytes':size,'sha256':rsha,'captured':sorted(blocks)}
  if sorted(blocks)!=INDICES:raise RuntimeError('heldout raw blocks missing')
  with tempfile.TemporaryDirectory() as td:
   lp=os.path.join(td,'p5.cdf');cs,csha=download(CDF,lp);o['cdf']={'size_bytes':cs,'sha256':csha}
   import netCDF4,numpy as np
   ds=netCDF4.Dataset(lp)
   try:
    rows=[]
    for idx in INDICES:
     mapped=map_line(blocks[idx]);img=np.asarray(ds.variables['image'][idx,:],dtype=float).reshape(-1).tolist()
     rows.append({'index':idx,'pearson':pearson(img,mapped),'spearman':spear(img,mapped),'raw_block_sha256':hashlib.sha256(blocks[idx]).hexdigest(),'cdf_image_sha256_float64_le':hashlib.sha256(np.asarray(img,dtype='<f8').tobytes()).hexdigest()})
    o['rows']=rows;o['summary']={'pearson_median':statistics.median(r['pearson'] for r in rows),'pearson_min':min(r['pearson'] for r in rows),'spearman_median':statistics.median(r['spearman'] for r in rows),'spearman_min':min(r['spearman'] for r in rows)}
   finally:ds.close()
  o['status']='FACTOR8_HELDOUT_FORMAT_VALIDATION_READY'
 except Exception as e:o['status']='FACTOR8_HELDOUT_FORMAT_VALIDATION_FAILED';o['error_type']=type(e).__name__;o['error']=str(e)
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(o,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps({'status':o['status'],'summary':o.get('summary'),'rows':o.get('rows')},indent=2));return 0 if o['status']=='FACTOR8_HELDOUT_FORMAT_VALIDATION_READY' else 2
if __name__=='__main__':raise SystemExit(main())
