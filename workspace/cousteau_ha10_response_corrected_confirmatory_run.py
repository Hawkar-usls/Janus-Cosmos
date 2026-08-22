#!/usr/bin/env python3
from __future__ import annotations

import hashlib, io, json, math, time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import requests
from obspy import read, read_inventory, UTCDateTime
from scipy.signal import welch

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'/'cousteau'
PROTOCOL=DATA/'JANUS-ECHO-COUSTEAU-HA10-PREDECLARED-RESPONSE-CORRECTED-CONFIRMATORY-PROTOCOL-2026-08-22-v1.0.json'
WINDOWS=DATA/'JANUS-ECHO-COUSTEAU-HA10-CONFIRMATORY-WINDOW-FREEZE-2026-08-22-v1.0.json'
OUT=DATA/'JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0.json'
EXPECTED_WINDOW_FREEZE_SHA256='d5f3c29c1dc4f7d7862724d1225688f9ee88460266d32fb0ec99fabe52cf2671'
FDSN='https://service.earthscope.org/fdsnws'


def sha(b:bytes)->str: return hashlib.sha256(b).hexdigest()
def file_sha(p:Path)->str: return sha(p.read_bytes())
def now()->str: return datetime.now(timezone.utc).isoformat()

def parse_channel_id(cid:str):
    net,sta,loc,cha=cid.split('.')
    return net,sta,loc,cha

def get_bytes(session,url,params,tries=4,timeout=120):
    last=None
    for i in range(tries):
        try:
            r=session.get(url,params=params,timeout=timeout)
            if r.status_code==200 and r.content: return r
            last=RuntimeError(f'HTTP {r.status_code}: {r.text[:300]}')
        except Exception as e: last=e
        time.sleep(1.5*(i+1))
    raise last or RuntimeError('download failed')

def fetch_inventory(session,cid,start,end):
    net,sta,loc,cha=parse_channel_id(cid)
    params={'net':net,'sta':sta,'loc':'--' if loc=='' else loc,'cha':cha,'starttime':start,'endtime':end,'level':'response','format':'xml','nodata':404}
    r=get_bytes(session,FDSN+'/station/1/query',params,timeout=90)
    return read_inventory(io.BytesIO(r.content)), {'url':r.url,'bytes':len(r.content),'sha256':sha(r.content),'content_type':r.headers.get('content-type')}

def fetch_trace(session,cid,start,end):
    net,sta,loc,cha=parse_channel_id(cid)
    params={'net':net,'sta':sta,'loc':'--' if loc=='' else loc,'cha':cha,'starttime':start,'endtime':end,'nodata':404,'format':'miniseed'}
    r=get_bytes(session,FDSN+'/dataselect/1/query',params,timeout=180)
    st=read(io.BytesIO(r.content))
    st.select(network=net,station=sta,channel=cha)
    if len(st)==0: raise RuntimeError(f'no trace for {cid}')
    st.merge(method=0,fill_value=None)
    if len(st)!=1: raise RuntimeError(f'multiple unmerged traces for {cid}: {len(st)}')
    tr=st[0]
    if np.ma.isMaskedArray(tr.data) and np.any(np.ma.getmaskarray(tr.data)): raise RuntimeError(f'gaps present for {cid}; interpolation forbidden')
    tr=tr.copy(); tr.detrend('demean'); tr.detrend('linear')
    fs=float(tr.stats.sampling_rate)
    expected=(UTCDateTime(end)-UTCDateTime(start))*fs
    if len(tr.data) < expected*0.98: raise RuntimeError(f'incomplete window {cid}: {len(tr.data)} < 98% of {expected}')
    return tr, {'url':r.url,'bytes':len(r.content),'sha256':sha(r.content),'npts':int(tr.stats.npts),'sampling_rate_hz':fs,'starttime':str(tr.stats.starttime),'endtime':str(tr.stats.endtime)}

def psd_with_response(tr,inv,nperseg,noverlap):
    x=np.asarray(tr.data,dtype=np.float64)
    fs=float(tr.stats.sampling_rate)
    if len(x)<nperseg: raise RuntimeError(f'npts {len(x)} < nperseg {nperseg}')
    f,p=welch(x,fs=fs,window='hann',nperseg=nperseg,noverlap=noverlap,detrend=False,scaling='density',return_onesided=True)
    resp=inv.get_response(tr.id,tr.stats.starttime)
    h,fr=resp.get_evalresp_response(t_samp=1.0/fs,nfft=nperseg,output='DEF')
    if len(fr)!=len(f) or np.max(np.abs(fr-f))>1e-9: raise RuntimeError('response and Welch frequency grids mismatch')
    amp=np.abs(h)
    corr=np.full_like(p,np.nan,dtype=float)
    valid=np.isfinite(amp)&(amp>0)&np.isfinite(p)
    corr[valid]=p[valid]/(amp[valid]**2)
    return f,p,corr,amp

def band_median(f,p,band):
    lo,hi=map(float,band); m=(f>=lo)&(f<=hi)&np.isfinite(p)
    if not np.any(m): raise RuntimeError(f'no finite bins in band {band}')
    return float(np.median(p[m])), int(np.sum(m))

def peak_in_band(f,p,band):
    lo,hi=map(float,band); idx=np.where((f>=lo)&(f<=hi)&np.isfinite(p))[0]
    if len(idx)==0: raise RuntimeError(f'no peak bins in band {band}')
    j=idx[int(np.argmax(p[idx]))]
    return float(f[j]), float(p[j])

def db_ratio(a,b):
    if not (a>0 and b>0 and math.isfinite(a) and math.isfinite(b)): return float('nan')
    return 10.0*math.log10(a/b)

def summarize_window(f,raw,corr,contract):
    target=contract['primary_target_band_hz']; controls=contract['predeclared_control_bands_hz']
    traw,ntr=band_median(f,raw,target); tcorr,_=band_median(f,corr,target)
    out={'target':{'raw_power':traw,'corrected_power_pa2_per_hz':tcorr,'bin_count':ntr},'controls':{}}
    for name,band in controls.items():
        rv,n=band_median(f,raw,band); cv,_=band_median(f,corr,band)
        out['controls'][name]={'band_hz':band,'raw_power':rv,'corrected_power_pa2_per_hz':cv,'bin_count':n}
    pf,pp=peak_in_band(f,corr,target); out['target']['corrected_peak_frequency_hz']=pf; out['target']['corrected_peak_power']=pp
    return out

def main():
    protocol=json.loads(PROTOCOL.read_text()); windows=json.loads(WINDOWS.read_text())
    if windows.get('status')!='WINDOWS_FROZEN_READY_FOR_FFT': raise RuntimeError('window freeze is not ready')
    if windows.get('window_freeze_sha256')!=EXPECTED_WINDOW_FREEZE_SHA256: raise RuntimeError('window freeze hash mismatch; FFT forbidden')
    if protocol.get('status')!='PREREGISTERED_BEFORE_TARGET_BAND_WAVEFORM_INSPECTION': raise RuntimeError('protocol preregistration status mismatch')
    spec=protocol['spectral_estimator']; freqc=protocol['frequency_contract']; th=protocol['admission_thresholds']
    nperseg=int(spec['welch_nperseg_samples']); noverlap=int(spec['welch_noverlap_samples'])
    s=requests.Session(); s.headers['User-Agent']='Janus-Echo-Cousteau/1.0 frozen confirmatory spectral run'

    # Freeze exact StationXML responses before waveform target statistics are computed.
    inventories={}; response_meta={}
    for ch in protocol['channels']:
        cid=ch['id']; inv,meta=fetch_inventory(s,cid,protocol['public_analysis_interval_utc']['start'],protocol['public_analysis_interval_utc']['end'])
        inventories[cid]=inv; response_meta[cid]=meta

    results=[]; data_errors=[]
    for ev in [x for x in windows['selected_events'] if x.get('complete_on_both_stations')]:
        er={'selection_rank':ev['selection_rank'],'source_time_code':ev['source_time_code'],'origin_utc':ev['origin_utc'],'source_magnitude_db':ev['source_magnitude_db'],'stations':{}}
        for cid,sw in ev['stations'].items():
            sr={'predicted_arrival_utc':sw['predicted_arrival_utc'],'event_window':sw['event_window'],'noise_windows':sw['noise_windows']}
            try:
                tr,emeta=fetch_trace(s,cid,sw['event_window']['start_utc'],sw['event_window']['end_utc'])
                f,raw,corr,amp=psd_with_response(tr,inventories[cid],nperseg,noverlap); event_sum=summarize_window(f,raw,corr,freqc)
                noise_sums=[]; noise_meta=[]
                for nw in sw['noise_windows']:
                    nt,nmeta=fetch_trace(s,cid,nw['start_utc'],nw['end_utc'])
                    nf,nraw,ncorr,namp=psd_with_response(nt,inventories[cid],nperseg,noverlap)
                    if np.max(np.abs(nf-f))>1e-9: raise RuntimeError('noise frequency grid mismatch')
                    noise_sums.append(summarize_window(nf,nraw,ncorr,freqc)); noise_meta.append(nmeta)
                target_noise_raw=float(np.median([x['target']['raw_power'] for x in noise_sums])); target_noise_corr=float(np.median([x['target']['corrected_power_pa2_per_hz'] for x in noise_sums]))
                raw_snr=db_ratio(event_sum['target']['raw_power'],target_noise_raw); corr_snr=db_ratio(event_sum['target']['corrected_power_pa2_per_hz'],target_noise_corr)
                near_corr=float(np.median([event_sum['controls']['near_left']['corrected_power_pa2_per_hz'],event_sum['controls']['near_right']['corrected_power_pa2_per_hz']]))
                side_contrast=db_ratio(event_sum['target']['corrected_power_pa2_per_hz'],near_corr)
                control_scores={}
                for name in freqc['predeclared_control_bands_hz']:
                    nr=float(np.median([x['controls'][name]['raw_power'] for x in noise_sums])); nc=float(np.median([x['controls'][name]['corrected_power_pa2_per_hz'] for x in noise_sums]))
                    control_scores[name]={'raw_event_vs_noise_snr_db':db_ratio(event_sum['controls'][name]['raw_power'],nr),'corrected_event_vs_noise_snr_db':db_ratio(event_sum['controls'][name]['corrected_power_pa2_per_hz'],nc)}
                max_control=max(v['corrected_event_vs_noise_snr_db'] for v in control_scores.values() if math.isfinite(v['corrected_event_vs_noise_snr_db']))
                margin=corr_snr-max_control
                peak=event_sum['target']['corrected_peak_frequency_hz']
                pt=th['per_station_event']
                checks={
                  'raw_snr_ge_threshold': raw_snr>=float(pt['minimum_raw_event_vs_noise_snr_db']),
                  'corrected_snr_ge_threshold': corr_snr>=float(pt['minimum_corrected_event_vs_noise_snr_db']),
                  'raw_corrected_snr_agree': abs(raw_snr-corr_snr)<=float(pt['maximum_abs_raw_minus_corrected_snr_db']),
                  'target_sideband_contrast_ge_threshold': side_contrast>=float(pt['minimum_corrected_target_sideband_contrast_db']),
                  'target_peak_center_in_tolerance': float(pt['target_peak_frequency_must_be_within_hz'][0])<=peak<=float(pt['target_peak_frequency_must_be_within_hz'][1])
                }
                sr.update({'data_status':'ANALYZED','event_waveform':emeta,'noise_waveforms':noise_meta,'event_band_summary':event_sum,'noise_band_summaries':noise_sums,'raw_event_vs_noise_snr_db':raw_snr,'corrected_event_vs_noise_snr_db':corr_snr,'raw_minus_corrected_snr_db':raw_snr-corr_snr,'corrected_target_sideband_contrast_db':side_contrast,'control_band_scores':control_scores,'target_minus_max_control_corrected_snr_db':margin,'target_peak_frequency_hz':peak,'per_station_event_checks':checks,'per_station_event_pass':all(checks.values())})
            except Exception as e:
                sr.update({'data_status':'BLOCKED_PAIR','error_type':type(e).__name__,'error':str(e),'per_station_event_pass':False}); data_errors.append({'event':ev['source_time_code'],'channel':cid,'error':str(e)})
            er['stations'][cid]=sr
        # Cross-station replication is evaluated only after both station results exist.
        ids=[x['id'] for x in protocol['channels']]; a=er['stations'][ids[0]]; b=er['stations'][ids[1]]
        peaks_ok=False
        if a.get('data_status')=='ANALYZED' and b.get('data_status')=='ANALYZED':
            peaks_ok=abs(a['target_peak_frequency_hz']-b['target_peak_frequency_hz'])<=float(th['cross_station_replication']['maximum_peak_center_difference_between_stations_hz'])
        er['cross_station_peak_difference_hz']=abs(a.get('target_peak_frequency_hz',float('nan'))-b.get('target_peak_frequency_hz',float('nan'))) if a.get('data_status')=='ANALYZED' and b.get('data_status')=='ANALYZED' else None
        er['replicated_both_stations']=bool(a.get('per_station_event_pass') and b.get('per_station_event_pass') and peaks_ok)
        results.append(er)

    replicated=[x for x in results if x['replicated_both_stations']]
    analyzed_pairs=[sr for ev in results for sr in ev['stations'].values() if sr.get('data_status')=='ANALYZED']
    repl_pairs=[sr for ev in replicated for sr in ev['stations'].values() if sr.get('data_status')=='ANALYZED']
    margin_req=float(th['global_candidate']['target_band_must_outperform_every_predeclared_control_band_by_db'])
    outperf=[sr['target_minus_max_control_corrected_snr_db']>=margin_req for sr in repl_pairs]
    outperf_fraction=(sum(outperf)/len(outperf)) if outperf else 0.0
    min_events=int(th['global_candidate']['minimum_independent_events_replicated_on_both_stations'])
    required_fraction=float(th['global_candidate']['outperformance_required_on_fraction_of_replicated_station_event_pairs'])

    if len(analyzed_pairs)<min_events*2:
        verdict='BLOCKED_DATA_ACCESS_OR_RESPONSE__NO_SCIENTIFIC_RESULT'
    elif len(replicated)==0:
        verdict='NEGATIVE_CONFIRMATORY_HA10_PUBLIC_SLICE'
    elif len(replicated)<min_events:
        verdict='AMBIGUOUS_REPEAT_INSUFFICIENT'
    elif outperf_fraction<required_fraction:
        verdict='REJECT_NON_SPECIFIC_OR_INSTRUMENTAL'
    else:
        verdict='REPLICATED_HA10_119HZ_CONFIRMATORY_CANDIDATE__NON_SPECIFIC__IDENTITY_UNCONFIRMED'

    out={'artifact_id':'JANUS-ECHO-COUSTEAU-HA10-RESPONSE-CORRECTED-CONFIRMATORY-RUN-001-2026-08-22-v1.0','started_and_completed_utc':now(),'protocol_path':str(PROTOCOL.relative_to(ROOT)),'protocol_file_sha256':file_sha(PROTOCOL),'window_freeze_path':str(WINDOWS.relative_to(ROOT)),'window_freeze_file_sha256':file_sha(WINDOWS),'window_freeze_sha256':windows['window_freeze_sha256'],'response_metadata':response_meta,'processing_contract':{'welch_nperseg':nperseg,'welch_noverlap':noverlap,'whole_spectrum_peak_search_performed':False,'only_predeclared_bands_summarized':True,'response_correction':'PSD_DIVIDE_BY_ABS_RESPONSE_SQUARED','thresholds_changed_after_inspection':False},'events':results,'summary':{'frozen_events':len(results),'analyzed_station_event_pairs':len(analyzed_pairs),'blocked_station_event_pairs':len(data_errors),'replicated_event_count':len(replicated),'replicated_source_time_codes':[x['source_time_code'] for x in replicated],'replicated_station_event_pairs':len(repl_pairs),'target_outperforms_all_controls_by_required_margin_fraction':outperf_fraction,'required_outperformance_fraction':required_fraction,'verdict':verdict,'target_identity':'UNCONFIRMED','target_evidence':'NOT_PROMOTED_BY_FREQUENCY_ONLY'},'data_errors':data_errors,'hard_rules':['FROZEN_WINDOWS_ONLY','NO_FFT_FISHING','NO_UNDECLARED_BANDS','NO_THRESHOLD_RETUNING','RAW_AND_RESPONSE_CORRECTED_SNR_BOTH_REPORTED','MATCHED_NOISE_REQUIRED','BOTH_STATIONS_REQUIRED','CONTROL_OUTPERFORMANCE_REQUIRED','PASS_DOES_NOT_IDENTIFY_SOURCE_OR_TARGET','NEGATIVE_AND_AMBIGUOUS_RESULTS_PRESERVED'],'status':'RUN_COMPLETE'}
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(out['summary'],indent=2,ensure_ascii=False))
    return 0

if __name__=='__main__': raise SystemExit(main())
