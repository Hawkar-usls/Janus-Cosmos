#!/usr/bin/env python3
from __future__ import annotations
import json, math, statistics
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'data'/'cousteau'/'JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json'
OUT=ROOT/'data'/'cousteau'/'JANUS-ECHO-COUSTEAU-HA10-CONFIRMATORY-AUDIT-2026-08-22-v1.0.json'

def finite(x): return isinstance(x,(int,float)) and math.isfinite(x)
def stats(xs):
    x=[float(v) for v in xs if finite(v)]
    if not x:return None
    return {'n':len(x),'min':min(x),'median':statistics.median(x),'max':max(x),'mean':statistics.fmean(x)}

def main():
    r=json.loads(SRC.read_text())
    pairs=[]; checks=Counter(); passed_checks=Counter(); fail_combo=Counter(); peak_centers=[]; raw_snr=[]; corr_snr=[]; contrast=[]; margins=[]
    for ev in r['events']:
        for cid,s in ev['stations'].items():
            if s.get('data_status')!='ANALYZED': continue
            c=s.get('per_station_event_checks',{})
            failed=tuple(sorted(k for k,v in c.items() if not v))
            fail_combo[failed]+=1
            for k,v in c.items():
                checks[k]+=1
                if v: passed_checks[k]+=1
            raw_snr.append(s.get('raw_event_vs_noise_snr_db'))
            corr_snr.append(s.get('corrected_event_vs_noise_snr_db'))
            contrast.append(s.get('corrected_target_sideband_contrast_db'))
            margins.append(s.get('target_minus_max_control_corrected_snr_db'))
            peak_centers.append(s.get('target_peak_frequency_hz'))
            pairs.append({'source_time_code':ev['source_time_code'],'station':cid,'raw_snr_db':s.get('raw_event_vs_noise_snr_db'),'corrected_snr_db':s.get('corrected_event_vs_noise_snr_db'),'sideband_contrast_db':s.get('corrected_target_sideband_contrast_db'),'target_minus_max_control_snr_db':s.get('target_minus_max_control_corrected_snr_db'),'peak_hz':s.get('target_peak_frequency_hz'),'checks':c,'pass':s.get('per_station_event_pass',False)})
    pass_rates={k:{'passed':passed_checks[k],'total':checks[k],'fraction':passed_checks[k]/checks[k] if checks[k] else None} for k in sorted(checks)}
    out={
      'artifact_id':'JANUS-ECHO-COUSTEAU-HA10-CONFIRMATORY-AUDIT-2026-08-22-v1.0',
      'created_utc':datetime.now(timezone.utc).isoformat(),
      'source_artifact':r['artifact_id'],
      'source_summary':r['summary'],
      'analysis_scope':'POSTRUN_AGGREGATION_ONLY__NO_REPROCESSING__NO_NEW_BANDS__NO_THRESHOLD_CHANGES',
      'station_event_pairs':len(pairs),
      'per_station_event_pass_count':sum(1 for p in pairs if p['pass']),
      'check_pass_rates':pass_rates,
      'metric_distributions':{
        'raw_event_vs_noise_snr_db':stats(raw_snr),
        'corrected_event_vs_noise_snr_db':stats(corr_snr),
        'corrected_target_sideband_contrast_db':stats(contrast),
        'target_minus_max_control_corrected_snr_db':stats(margins),
        'target_peak_frequency_hz':stats(peak_centers)
      },
      'failure_combinations':[{'failed_checks':list(k),'count':v} for k,v in fail_combo.most_common()],
      'pairs':pairs,
      'scientific_interpretation':'NEGATIVE_CONFIRMATORY_RESULT_REMAINS__AUDIT_ONLY_IDENTIFIES_WHICH_PREDECLARED_GATES_FAILED',
      'target_identity':'UNCONFIRMED',
      'target_evidence':'NOT_INCREASED',
      'hard_rules':['NO_RESCORING','NO_THRESHOLD_RETUNING','NO_NEW_FREQUENCY_BANDS','NEGATIVE_RESULT_IMMUTABLE']
    }
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'verdict':r['summary']['verdict'],'pairs':len(pairs),'pair_passes':out['per_station_event_pass_count'],'check_pass_rates':pass_rates,'metric_distributions':out['metric_distributions']},indent=2,ensure_ascii=False))

if __name__=='__main__': main()
