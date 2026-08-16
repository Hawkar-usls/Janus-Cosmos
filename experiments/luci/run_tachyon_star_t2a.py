#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json,time,urllib.parse,urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime
from pathlib import Path
SOURCE_SHA256='94dee39740c5a2c402832e8d6b1c6681521331fe11cee630402b95e873a0adcb'
EXPECTED_SOURCES=42
TRAIN_PRE=232.5827
TRAIN_POST=432.2595
MIN_GAP=30.0
MAX_GAP=900.0
ZTF_WORKERS=4
ZTF_COLS=('ra','dec','field','ccdid','qid','filtercode','pid','nid','expid','imgtypecode','obsdate','obsjd','exptime','filefracday','seeing','airmass','moonillf','maglimit','infobits','ipac_pub_date')
CLAIM='EXTERNAL_ARCHIVE_METADATA_DISCOVERY_ONLY__NO_TRANSIENT_DETECTION__NO_FTL__NO_RETROCAUSALITY__NO_TACHYON_IDENTITY__NO_NUCLEAR_CAUSALITY'
def sha256_file(p:Path):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(1<<20),b''):h.update(c)
 return h.hexdigest()
def _get(url:str,timeout=120,retries=2)->bytes:
 err=None
 for attempt in range(1,retries+1):
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'Janus-Cosmos-TachyonStar-T2A/1.2','Accept':'application/json,text/csv,text/plain,*/*'})
   with urllib.request.urlopen(req,timeout=timeout) as r:return r.read()
  except Exception as e:
   err=f'{type(e).__name__}: {e}'
   if attempt<retries: time.sleep(2*attempt)
 raise RuntimeError(f'GET failed after {retries} attempts: {url}: {err}')
def _dt(s:str)->datetime:
 s=(s or '').strip().replace('Z','+00:00')
 return datetime.fromisoformat(s)
def eligible_ztf_triples(rows:list[dict])->list[dict]:
 groups=defaultdict(list);seen=set()
 for r in rows:
  ident=(r['src_id'],r.get('filefracday',''),r.get('field',''),r.get('ccdid',''),r.get('qid',''),r.get('filtercode',''))
  if ident in seen:continue
  seen.add(ident)
  try:exptime=round(float(r.get('exptime') or 0),6);dt=_dt(r.get('obsdate',''))
  except Exception:continue
  key=(r['src_id'],r.get('field',''),r.get('ccdid',''),r.get('qid',''),r.get('filtercode',''),exptime);groups[key].append((dt,r))
 out=[]
 for key,xs in groups.items():
  xs.sort(key=lambda x:(x[0],x[1].get('filefracday','')))
  for i in range(1,len(xs)-1):
   a,b,c=xs[i-1],xs[i],xs[i+1];pre=(b[0]-a[0]).total_seconds();post=(c[0]-b[0]).total_seconds()
   if not(MIN_GAP<=pre<=MAX_GAP and MIN_GAP<=post<=MAX_GAP):continue
   out.append({'src_id':key[0],'field':key[1],'ccdid':key[2],'qid':key[3],'filtercode':key[4],'exptime_s':key[5],'a_obsdate':a[1]['obsdate'],'a_filefracday':a[1]['filefracday'],'a_pid':a[1].get('pid',''),'a_ipac_pub_date':a[1].get('ipac_pub_date',''),'b_obsdate':b[1]['obsdate'],'b_filefracday':b[1]['filefracday'],'b_pid':b[1].get('pid',''),'b_ipac_pub_date':b[1].get('ipac_pub_date',''),'c_obsdate':c[1]['obsdate'],'c_filefracday':c[1]['filefracday'],'c_pid':c[1].get('pid',''),'c_ipac_pub_date':c[1].get('ipac_pub_date',''),'delta_pre_s':pre,'delta_post_s':post,'timing_distance_s':abs(pre-TRAIN_PRE)+abs(post-TRAIN_POST)})
 out.sort(key=lambda r:(r['src_id'],r['timing_distance_s'],r['b_obsdate'],r['b_filefracday'],r['field'],r['ccdid'],r['qid']));return out
def primary_ztf(triples:list[dict])->list[dict]:
 best={}
 for r in triples:
  if r['src_id'] not in best:best[r['src_id']]=r
 return [best[k] for k in sorted(best)]
def write_csv(path:Path,rows:list[dict],fields:list[str]|None=None):
 path.parent.mkdir(parents=True,exist_ok=True)
 if fields is None:fields=list(rows[0].keys()) if rows else []
 with path.open('w',encoding='utf-8',newline='') as f:
  if not fields:return
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
def _ztf_query(idx:int,s:dict):
 sid=s['src_id'];ra=s['ra_deg'];dec=s['dec_deg'];safe=f'{idx:02d}_{hashlib.sha256(sid.encode()).hexdigest()[:12]}'
 params={'POS':f'{ra},{dec}','COLUMNS':','.join(ZTF_COLS),'ct':'csv'};url='https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?'+urllib.parse.urlencode(params)
 b=_get(url,timeout=120,retries=2)
 rr=list(csv.DictReader(b.decode('utf-8-sig').splitlines()));rows=[]
 for r in rr:
  if not r:continue
  r={k:(v if v is not None else '') for k,v in r.items()};r['src_id']=sid;r['source_ra_deg']=ra;r['source_dec_deg']=dec;rows.append(r)
 return {'idx':idx,'sid':sid,'safe':safe,'url':url,'bytes':b,'rows':rows}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--sources',default='data/tachyon_star/JANUS-TACHYON-STAR-T2-EXTERNAL-SOURCE-UNIVERSE.csv');ap.add_argument('--output-dir',default='results/tachyon_star_t2a');a=ap.parse_args();src=Path(a.sources);out=Path(a.output_dir);raw=out/'raw';raw.mkdir(parents=True,exist_ok=True)
 if sha256_file(src)!=SOURCE_SHA256:raise RuntimeError('frozen source-universe SHA mismatch')
 sources=list(csv.DictReader(src.open(encoding='utf-8')))
 if len(sources)!=EXPECTED_SOURCES or len({r['src_id'] for r in sources})!=EXPECTED_SOURCES:raise RuntimeError('source-universe cardinality changed')
 failures=[];tess_rows=[];ztf_rows=[];raw_manifest=[]
 # TESS has an explicit <=5 req/s service limit. Keep this arm serial and paced.
 for idx,s in enumerate(sources):
  sid=s['src_id'];ra=s['ra_deg'];dec=s['dec_deg'];safe=f'{idx:02d}_{hashlib.sha256(sid.encode()).hexdigest()[:12]}';tu='https://mast.stsci.edu/tesscut/api/v0.1/sector?'+urllib.parse.urlencode({'ra':ra,'dec':dec})
  try:
   b=_get(tu,timeout=60,retries=2);p=raw/f'tess_{safe}.json';p.write_bytes(b);raw_manifest.append({'src_id':sid,'arm':'TESS','path':str(p.relative_to(out)),'sha256':sha256_file(p),'bytes':len(b),'url':tu});j=json.loads(b.decode('utf-8'));results=j.get('results',[]) if isinstance(j,dict) else []
   if not results:tess_rows.append({'src_id':sid,'ra_deg':ra,'dec_deg':dec,'sector':'','sectorName':'','camera':'','ccd':'','covered':False})
   else:
    for q in results:tess_rows.append({'src_id':sid,'ra_deg':ra,'dec_deg':dec,'sector':q.get('sector',''),'sectorName':q.get('sectorName',''),'camera':q.get('camera',''),'ccd':q.get('ccd',''),'covered':True})
  except Exception as e:failures.append({'src_id':sid,'arm':'TESS','error':f'{type(e).__name__}: {e}'})
  time.sleep(0.25)
 # IRSA metadata queries are independent by frozen source; parallelism is transport-only.
 results=[]
 with ThreadPoolExecutor(max_workers=ZTF_WORKERS) as ex:
  futs={ex.submit(_ztf_query,idx,s):(idx,s['src_id']) for idx,s in enumerate(sources)}
  for fut in as_completed(futs):
   idx,sid=futs[fut]
   try:results.append(fut.result())
   except Exception as e:failures.append({'src_id':sid,'arm':'ZTF','error':f'{type(e).__name__}: {e}'})
 for res in sorted(results,key=lambda x:x['idx']):
  p=raw/f"ztf_{res['safe']}.csv";p.write_bytes(res['bytes']);raw_manifest.append({'src_id':res['sid'],'arm':'ZTF','path':str(p.relative_to(out)),'sha256':sha256_file(p),'bytes':len(res['bytes']),'url':res['url']});ztf_rows.extend(res['rows'])
 raw_manifest.sort(key=lambda r:(r['src_id'],r['arm']))
 triples=eligible_ztf_triples(ztf_rows);primary=primary_ztf(triples)
 write_csv(out/'tess_sector_coverage.csv',tess_rows,['src_id','ra_deg','dec_deg','sector','sectorName','camera','ccd','covered']);write_csv(out/'ztf_normalized_metadata.csv',ztf_rows,['src_id','source_ra_deg','source_dec_deg']+list(ZTF_COLS));write_csv(out/'ztf_eligible_triples.csv',triples);write_csv(out/'ztf_primary_metadata_sequences.csv',primary)
 (out/'raw_manifest.json').write_text(json.dumps(raw_manifest,indent=2,sort_keys=True)+'\n',encoding='utf-8');status='PASS_METADATA_COMPLETE' if not failures else 'BLOCKED_METADATA_QUERY_FAILURE';rec={'schema':'janus.cosmos.tachyon_star.t2a.receipt.v1','experiment_id':'JANUS-TACHYON-STAR-T2A-EXTERNAL-METADATA-DISCOVERY','status':status,'transport_revision':'R2_ZTF_MAX4_PARALLEL_IGNORE_UNREQUESTED_SERVER_COLUMNS_SAME_FROZEN_QUERIES','source_universe_sha256':SOURCE_SHA256,'sources':EXPECTED_SOURCES,'failures':sorted(failures,key=lambda x:(x['src_id'],x['arm'])),'tess':{'rows':len(tess_rows),'sources_with_coverage':len({r['src_id'] for r in tess_rows if r['covered']}),'sectors_unique':len({str(r['sector']) for r in tess_rows if r['covered']})},'ztf':{'metadata_rows':len(ztf_rows),'sources_with_metadata':len({r['src_id'] for r in ztf_rows}),'eligible_triples':len(triples),'sources_with_eligible_triples':len({r['src_id'] for r in triples}),'primary_sequences':len(primary)},'raw_manifest_sha256':sha256_file(out/'raw_manifest.json'),'claim_ceiling':CLAIM,'external_target_pixels_opened':False};(out/'receipt.json').write_text(json.dumps(rec,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps(rec,indent=2));return 0 if not failures else 3
if __name__=='__main__':raise SystemExit(main())
