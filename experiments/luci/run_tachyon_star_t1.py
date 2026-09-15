#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json
from datetime import datetime, timedelta
from pathlib import Path
from janus_cosmos.luci import read_luci_fits_image
from janus_cosmos.pipeline import EventWriter, download_source, sha256_file
from experiments.luci.run_palomar_2f_f import replay_row

MANIFEST_SHA256='334fe0cacd0b056214d991be5256e05ab678fe5edb072d46a26ad07d4ae803b7'
EXPECTED_TRIALS=9
EXPECTED_PAIRED=6
EXPECTED_NATIVE_SHAPE=(2048,2048)
CLAIM_CEILING='PREREGISTERED_LUCI_HOLDOUT_ISOLATED_TRANSIENT_MOTIF_TEST_ONLY__NO_FTL__NO_RETROCAUSALITY__NO_TACHYON_IDENTITY__NO_NUCLEAR_CAUSALITY__NO_UAP_ORIGIN_CLAIM'

def _rows(path: Path):
    return list(csv.DictReader(path.open('r',encoding='utf-8',newline='')))

def _phase_row(r, phase):
    return {'src_id':r['src_id'],'file_name':r[f'{phase}_file'],'file_url':r[f'{phase}_url'],'instrument':r['instrument'],'filters':r['filters'],'exact_hdu':r[f'{phase}_hdu'],'exact_x':r[f'{phase}_x'],'exact_y':r[f'{phase}_y'],'exact_naxis1':str(EXPECTED_NATIVE_SHAPE[1]),'exact_naxis2':str(EXPECTED_NATIVE_SHAPE[0]),'file_sha256':r[f'{phase}_sha256'],'date_obs':r[f'{phase}_date_obs']}

def _paired_row(r):
    if not r.get('paired_b_file'): return None
    return {'src_id':r['src_id'],'file_name':r['paired_b_file'],'file_url':r['paired_b_url'],'instrument':r['paired_b_instrument'],'filters':r['paired_b_filters'],'exact_hdu':r['paired_b_hdu'],'exact_x':r['paired_b_x'],'exact_y':r['paired_b_y'],'exact_naxis1':str(EXPECTED_NATIVE_SHAPE[1]),'exact_naxis2':str(EXPECTED_NATIVE_SHAPE[0]),'file_sha256':r['paired_b_sha256'],'date_obs':r['paired_b_date_obs']}

def _download_replay(row, cache, events, seed, target):
    path,_=download_source(row['file_url'],cache,events,target=target,filter_name=row.get('filters','UNKNOWN'))
    got=sha256_file(path)
    if got != row['file_sha256']:
        return {'status':'BLOCKED_FILE_SHA_MISMATCH','expected_sha256':row['file_sha256'],'actual_sha256':got}
    image,meta=read_luci_fits_image(path,require_imaging=True,expected_instrument=row.get('instrument'))
    result=replay_row(row,image,meta,seed)
    return {'status':result['status'],'image_meta':meta,'replay':result}

def _is_absence(x): return x.get('status','').startswith('QUALIFIED_NO_COUNTERPART')
def _is_candidate(x): return x.get('status')=='COUNTERPART_CANDIDATE'
def _adjudicable(x): return _is_absence(x) or _is_candidate(x)

def _exposure_overlap(primary, paired):
    try:
        a=datetime.fromisoformat(primary['image_meta']['date_obs']); b=datetime.fromisoformat(paired['image_meta']['date_obs'])
        ea=float(primary['image_meta']['exptime']); eb=float(paired['image_meta']['exptime'])
    except Exception as exc:
        return {'adjudicating':False,'reason':f'EXPOSURE_METADATA_UNAVAILABLE:{type(exc).__name__}'}
    end_a=a+timedelta(seconds=ea); end_b=b+timedelta(seconds=eb)
    overlap=max(0.0,(min(end_a,end_b)-max(a,b)).total_seconds()); frac=overlap/min(ea,eb) if min(ea,eb)>0 else 0.0
    return {'primary_start':a.isoformat(),'paired_start':b.isoformat(),'primary_exptime_s':ea,'paired_exptime_s':eb,'overlap_seconds':overlap,'overlap_fraction_of_shorter':frac,'adjudicating':frac>=0.5}

def classify_trial(a,b,c,paired=None):
    base={'a_status':a['status'],'b_status':b['status'],'c_status':c['status']}
    if not (_adjudicable(a) and _adjudicable(b) and _adjudicable(c)): return {**base,'class':'UNRESOLVED_PRIMARY_TRIAL'}
    isolated=_is_absence(a) and _is_candidate(b) and _is_absence(c)
    if not isolated: return {**base,'class':'NO_ISOLATED_B_EVENT'}
    out={**base,'class':'ISOLATED_B_EVENT_CANDIDATE'}
    if paired is None: return {**out,'paired_class':'NO_PAIRED_INSTRUMENT_WITNESS'}
    ov=_exposure_overlap(b,paired); out['paired_overlap']=ov
    if not ov.get('adjudicating'): return {**out,'paired_class':'PAIRED_NONOVERLAPPING_OR_METADATA_BLOCKED'}
    if _is_candidate(paired): return {**out,'paired_class':'DUAL_INSTRUMENT_COINCIDENCE_CANDIDATE'}
    if _is_absence(paired): return {**out,'paired_class':'PRIMARY_ONLY_WITH_PAIRED_SENSITIVE_ABSENCE'}
    return {**out,'paired_class':'PAIRED_BLOCKED'}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--manifest',default='data/luci_palomar/JANUS-TACHYON-STAR-T1-LOOKBACK-HOLDOUT-MANIFEST.csv'); ap.add_argument('--output-dir',default='results/tachyon_star_t1'); ap.add_argument('--cache-dir',default='.cache/tachyon_star_t1'); a=ap.parse_args()
    manifest=Path(a.manifest); out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    if sha256_file(manifest)!=MANIFEST_SHA256: raise RuntimeError('sealed manifest SHA mismatch')
    rows=_rows(manifest)
    if len(rows)!=EXPECTED_TRIALS: raise RuntimeError(f'trial cardinality changed: {len(rows)}')
    if sum(bool(r.get('paired_b_file')) for r in rows)!=EXPECTED_PAIRED: raise RuntimeError('paired-witness cardinality changed')
    events=EventWriter(out/'events.jsonl'); cache=Path(a.cache_dir); trials=[]
    for i,r in enumerate(rows):
        phases={}
        for j,p in enumerate(('a','b','c')):
            phases[p]=_download_replay(_phase_row(r,p),cache,events,20263000+i*10+j,f'TACHYON_STAR_T1_{r["trial_id"]}_{p.upper()}')
        paired=None
        if _is_candidate(phases['b']) and r.get('paired_b_file'):
            paired=_download_replay(_paired_row(r),cache,events,20263000+i*10+3,f'TACHYON_STAR_T1_{r["trial_id"]}_PAIRED_B')
        cls=classify_trial(phases['a'],phases['b'],phases['c'],paired)
        trials.append({'trial_id':r['trial_id'],'src_id':r['src_id'],'ra_deg':float(r['ra_deg']),'dec_deg':float(r['dec_deg']),'instrument':r['instrument'],'filters':r['filters'],'delta_pre_s':float(r['delta_pre_s']),'delta_post_s':float(r['delta_post_s']),'timing_distance_s':float(r['timing_distance_s']),'sealed_files':{'a':r['a_file'],'b':r['b_file'],'c':r['c_file'],'paired_b':r.get('paired_b_file') or None},'phases':phases,'paired_b':paired,'classification':cls})
    isolated=[t for t in trials if t['classification']['class']=='ISOLATED_B_EVENT_CANDIDATE']; unresolved=[t for t in trials if t['classification']['class']=='UNRESOLVED_PRIMARY_TRIAL']; dual=[t for t in isolated if t['classification'].get('paired_class')=='DUAL_INSTRUMENT_COINCIDENCE_CANDIDATE']; veto=[t for t in isolated if t['classification'].get('paired_class')=='PRIMARY_ONLY_WITH_PAIRED_SENSITIVE_ABSENCE']
    if unresolved: status='BLOCKED'; scientific='PRIMARY_HOLDOUT_INCOMPLETE'; rc=3
    elif isolated: status='CANDIDATE'; scientific='ONE_OR_MORE_PREREGISTERED_ISOLATED_B_EVENTS_OBSERVED'; rc=4
    else: status='PASS_NEGATIVE'; scientific='NO_PREREGISTERED_ISOLATED_B_EVENT_IN_ADJUDICABLE_PRIMARY_HOLDOUT'; rc=0
    receipt={'schema':'janus.cosmos.tachyon_star.lookback_holdout.receipt.v1','experiment_id':'JANUS-TACHYON-STAR-T1-LOOKBACK-HOLDOUT','status':status,'scientific_status':scientific,'manifest_sha256':MANIFEST_SHA256,'primary_trials':len(trials),'isolated_event_count':len(isolated),'unresolved_trial_count':len(unresolved),'dual_instrument_candidate_count':len(dual),'primary_only_paired_sensitive_absence_count':len(veto),'trials':trials,'claim_ceiling':CLAIM_CEILING,'physical_tachyon_identity_admitted':False}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({k:receipt[k] for k in ['status','scientific_status','primary_trials','isolated_event_count','unresolved_trial_count','dual_instrument_candidate_count','primary_only_paired_sensitive_absence_count']},indent=2)); return rc

if __name__=='__main__': raise SystemExit(main())
