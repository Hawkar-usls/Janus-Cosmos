#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from urllib.parse import urlencode
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from janus_cosmos.pipeline import EventWriter,download_source,sha256_file
from janus_cosmos.luci_psf_r1 import measure_psf_at,psf_relative_injection_recovery_gate
from experiments.luci.run_palomar_2f_e import local_coordinate_injection

MANIFEST_GIT_BLOB_SHA1='902ac74624b8c30c365544555f5071c78ec4cae2'
MANIFEST_SHA256_OBSERVED_PRE_PIXEL='f8242864c1822ee6636c885cbe3cc0bf9aa53d4829213ce1ea13c2ee068c6ecb'
EXPECTED_TRIALS=41
EXPECTED_QUALITY_ADMISSIBLE=28
BAD_INFOBITS_THRESHOLD=33554432
CUTOUT_ARCSEC=180
MASK_RADIUS_PX=3.0
CLAIM='PREPOINTED_ZTF_ISOLATED_TRANSIENT_REPLICATION_ONLY__NO_FTL__NO_RETROCAUSALITY__NO_TACHYON_PARTICLE_IDENTIFICATION__NO_NUCLEAR_CAUSALITY__NO_UAP_OR_ARTIFICIAL_ORIGIN_CLAIM'

def _rows(path): return list(csv.DictReader(Path(path).open(encoding='utf-8',newline='')))
def _git_blob_sha1(path:Path):
 b=path.read_bytes();return hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()
def _quality_ok(row): return all(row[f'{e}_imgtypecode']=='o' and int(row[f'{e}_infobits'])<BAD_INFOBITS_THRESHOLD for e in 'abc')
def _product_url(row,epoch,suffix):
 ff=row[f'{epoch}_filefracday'];year=ff[:4];mmdd=ff[4:8];frac=ff[8:14];field=str(row['field']).zfill(6);ccd=str(row['ccdid']).zfill(2);qid=str(row['qid']);filt=row['filtercode'];itype=row[f'{epoch}_imgtypecode']
 base=f'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/{year}/{mmdd}/{frac}/ztf_{ff}_{field}_{filt}_c{ccd}_{itype}_q{qid}_{suffix}'
 return base+'?'+urlencode({'center':f"{row['ra_deg']},{row['dec_deg']}",'size':f'{CUTOUT_ARCSEC}arcsec','gzip':'false'})
def _read_image_wcs(path):
 with fits.open(path,memmap=False,lazy_load_hdus=True) as hdul:
  cand=[]
  for i,h in enumerate(hdul):
   d=getattr(h,'data',None)
   if d is None or np.ndim(d)!=2: continue
   arr=np.asarray(d,dtype=np.float32)
   if arr.size<256: continue
   try:w=WCS(h.header)
   except Exception:continue
   cand.append((arr.size,i,arr,w))
  if not cand: raise RuntimeError(f'no usable 2-D WCS image in {path}')
  _,idx,arr,w=max(cand,key=lambda x:x[0])
 return arr,w,idx
def _xy(wcs,ra,dec):
 x,y=wcs.world_to_pixel_values(float(ra),float(dec));return float(x),float(y)
def _mask_clean(mask,wcs,ra,dec):
 x,y=_xy(wcs,ra,dec);yy,xx=np.indices(mask.shape);sel=((xx-x)**2+(yy-y)**2<=MASK_RADIUS_PX**2)&np.isfinite(mask)
 if not np.any(sel): return {'passed':False,'reason':'NO_MASK_PIXELS_IN_TARGET_APERTURE','x':x,'y':y}
 vals=np.asarray(mask[sel],dtype=np.int64);return {'passed':bool(np.all(vals==0)),'nonzero_pixels':int(np.sum(vals!=0)),'aperture_pixels':int(vals.size),'x':x,'y':y,'radius_px':MASK_RADIUS_PX}
def _safe_frame_gate(image,seed):
 try:return psf_relative_injection_recovery_gate(image,seed=seed)
 except Exception as exc:return {'passed':False,'reason':'R1_FRAME_GATE_EXCEPTION_FAIL_CLOSED','error':f'{type(exc).__name__}: {exc}'}
def classify_epoch(image,wcs,mask,mask_wcs,ra,dec,seed):
 x,y=_xy(wcs,ra,dec);h,w=image.shape
 if not (12<=x<w-12 and 12<=y<h-12): return {'status':'BLOCKED_TARGET_EDGE','x':x,'y':y,'shape':[h,w]}
 mg=_mask_clean(mask,mask_wcs,ra,dec);gate=_safe_frame_gate(image,seed)
 if not gate.get('passed'): return {'status':'BLOCKED_FRAME_SENSITIVITY','x':x,'y':y,'mask_gate':mg,'frame_gate':gate}
 target=measure_psf_at(image,y,x)
 if target is not None:return {'status':'SOURCE_PRESENT' if mg.get('passed') else 'SOURCE_PRESENT_MASK_BLOCKED','x':x,'y':y,'mask_gate':mg,'frame_gate':gate,'target':target}
 local=local_coordinate_injection(image,x,y,gate['injection_base_fwhm_px'])
 if not mg.get('passed'): return {'status':'BLOCKED_TARGET_MASK','x':x,'y':y,'mask_gate':mg,'frame_gate':gate,'local':local}
 if not local.get('passed'): return {'status':'BLOCKED_LOCAL_SENSITIVITY','x':x,'y':y,'mask_gate':mg,'frame_gate':gate,'local':local}
 return {'status':'QUALIFIED_ABSENCE','x':x,'y':y,'mask_gate':mg,'frame_gate':gate,'local':local}
def _classify_trial(a,b,c):
 if a['status']=='QUALIFIED_ABSENCE' and b['status']=='SOURCE_PRESENT' and c['status']=='QUALIFIED_ABSENCE': return 'ISOLATED_B_L0'
 if a['status']=='QUALIFIED_ABSENCE' and b['status']=='QUALIFIED_ABSENCE' and c['status']=='QUALIFIED_ABSENCE': return 'NO_ISOLATED_B_EVENT'
 if 'SOURCE_PRESENT' in a['status'] or 'SOURCE_PRESENT' in c['status']: return 'NON_ISOLATED_SOURCE_PATTERN'
 return 'UNRESOLVED_TRIAL'
def _blocked_epoch(exc): return {'status':'BLOCKED_PIXEL_REPLAY_ERROR','error':f'{type(exc).__name__}: {exc}'}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest',default='data/tachyon_star/JANUS-TACHYON-STAR-T2B-ZTF-PREPOINTED-MANIFEST.csv');ap.add_argument('--output-dir',default='results/tachyon_star_t2b_ztf');ap.add_argument('--cache-dir',default='.cache/tachyon_star_t2b_ztf');a=ap.parse_args();manifest=Path(a.manifest);out=Path(a.output_dir);cache=Path(a.cache_dir);out.mkdir(parents=True,exist_ok=True)
 manifest_sha256=sha256_file(manifest);manifest_blob=_git_blob_sha1(manifest)
 if manifest_blob!=MANIFEST_GIT_BLOB_SHA1: raise RuntimeError(f'frozen manifest Git blob mismatch: {manifest_blob}')
 if manifest_sha256!=MANIFEST_SHA256_OBSERVED_PRE_PIXEL: raise RuntimeError(f'frozen manifest observed SHA256 mismatch: {manifest_sha256}')
 rows=_rows(manifest)
 if len(rows)!=EXPECTED_TRIALS: raise RuntimeError(f'manifest cardinality changed: {len(rows)}')
 qa=sum(_quality_ok(r) for r in rows)
 if qa!=EXPECTED_QUALITY_ADMISSIBLE: raise RuntimeError(f'quality-admissible cardinality changed: {qa}')
 events=EventWriter(out/'events.jsonl');results=[]
 for ti,row in enumerate(rows):
  if not _quality_ok(row): results.append({'trial_id':row['trial_id'],'src_id':row['src_id'],'classification':'BLOCKED_METADATA_QUALITY','metadata':row});continue
  ep={}
  for ei,e in enumerate('abc'):
   try:
    sci_url=_product_url(row,e,'sciimg.fits');msk_url=_product_url(row,e,'mskimg.fits');sci,_=download_source(sci_url,cache,events,target=row['trial_id'],filter_name=f'ZTF_{e}_SCI');msk,_=download_source(msk_url,cache,events,target=row['trial_id'],filter_name=f'ZTF_{e}_MASK');image,wcs,sci_hdu=_read_image_wcs(sci);mask,mwcs,mask_hdu=_read_image_wcs(msk);er=classify_epoch(image,wcs,mask,mwcs,row['ra_deg'],row['dec_deg'],20263000+ti*10+ei);er.update({'science_url':sci_url,'mask_url':msk_url,'science_hdu':sci_hdu,'mask_hdu':mask_hdu,'science_sha256':sha256_file(sci),'mask_sha256':sha256_file(msk)})
   except Exception as exc: er=_blocked_epoch(exc)
   ep[e]=er
  cls=_classify_trial(ep['a'],ep['b'],ep['c']);diff=None
  if cls=='ISOLATED_B_L0':
   try:
    du=_product_url(row,'b','scimrefdiffimg.fits.fz');dp,_=download_source(du,cache,events,target=row['trial_id'],filter_name='ZTF_B_DIFF');di,dw,dh=_read_image_wcs(dp);x,y=_xy(dw,row['ra_deg'],row['dec_deg']);dt=measure_psf_at(di,y,x);diff={'status':'DIFF_SOURCE_PRESENT' if dt is not None else 'DIFF_SOURCE_NOT_MEASURED','url':du,'sha256':sha256_file(dp),'hdu':dh,'x':x,'y':y,'target':dt};cls='ISOLATED_B_L1_DIFF_CONFIRMED' if dt is not None else cls
   except Exception as exc: diff={'status':'DIFF_UNAVAILABLE_OR_UNREADABLE','error':f'{type(exc).__name__}: {exc}'}
  results.append({'trial_id':row['trial_id'],'src_id':row['src_id'],'classification':cls,'metadata':row,'epochs':ep,'difference_adjudication':diff})
 admissible=[r for r in results if r['classification']!='BLOCKED_METADATA_QUALITY'];candidates=[r for r in admissible if r['classification'].startswith('ISOLATED_B_L')];unresolved=[r for r in admissible if r['classification']=='UNRESOLVED_TRIAL'];nulls=[r for r in admissible if r['classification']=='NO_ISOLATED_B_EVENT'];nonisolated=[r for r in admissible if r['classification']=='NON_ISOLATED_SOURCE_PATTERN'];bclusters=len({(r['metadata']['b_filefracday'],r['metadata']['field'],r['metadata']['ccdid'],r['metadata']['qid'],r['metadata']['filtercode']) for r in admissible})
 status='PASS_ADMISSIBLE_NULL' if not candidates and not unresolved else ('CANDIDATE_REQUIRES_ADJUDICATION' if candidates else 'BLOCKED_ADMISSIBLE_INCOMPLETE')
 rec={'schema':'janus.cosmos.tachyon_star.t2b.ztf_prepointed.receipt.v1','experiment_id':'JANUS-TACHYON-STAR-T2B-ZTF-PREPOINTED-PIXEL-REPLICATION','status':status,'manifest_git_blob_sha1':manifest_blob,'manifest_sha256_observed_pre_pixel':manifest_sha256,'trials_total':len(rows),'metadata_quality_blocked':len(rows)-len(admissible),'quality_admissible':len(admissible),'unique_b_exposure_clusters_admissible':bclusters,'no_isolated_b_event':len(nulls),'nonisolated_source_pattern':len(nonisolated),'unresolved':len(unresolved),'isolated_b_candidates':len(candidates),'candidate_trial_ids':[r['trial_id'] for r in candidates],'unresolved_trial_ids':[r['trial_id'] for r in unresolved],'results':results,'claim_ceiling':CLAIM}
 (out/'receipt.json').write_text(json.dumps(rec,indent=2,sort_keys=True)+'\n',encoding='utf-8');print(json.dumps({k:rec[k] for k in ['status','manifest_sha256_observed_pre_pixel','trials_total','metadata_quality_blocked','quality_admissible','unique_b_exposure_clusters_admissible','no_isolated_b_event','nonisolated_source_pattern','unresolved','isolated_b_candidates','candidate_trial_ids','unresolved_trial_ids']},indent=2));return 0 if status=='PASS_ADMISSIBLE_NULL' else 3
if __name__=='__main__': raise SystemExit(main())
