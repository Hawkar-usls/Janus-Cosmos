#!/usr/bin/env python3
from __future__ import annotations
import hashlib, io, json, math
from datetime import datetime, timezone
from pathlib import Path
import requests
from obspy import read_inventory

OUT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-HA10-RESPONSE-EDGE-CALIBRATION-PROBE-2026-08-21-v1.0.json')
URL='https://service.earthscope.org/fdsnws/station/1/query'
CHANNELS=[('IM','H10N1','--','EDH'),('IM','H10S2','--','EDH')]
FREQS=[1.0,5.0,10.0,40.0,70.0,90.0,100.0,105.0,110.0,115.0,117.0,119.0,120.0,121.0,124.0]

def sha(b):return hashlib.sha256(b).hexdigest()

def main():
 s=requests.Session();s.headers['User-Agent']='Janus-Echo-Cousteau/1.0 response-only band-edge calibration'
 rep={'artifact_id':'JANUS-ECHO-COUSTEAU-HA10-RESPONSE-EDGE-CALIBRATION-PROBE-2026-08-21-v1.0','created_utc':datetime.now(timezone.utc).isoformat(),'policy':'INSTRUMENT_RESPONSE_ONLY__NO_WAVEFORM_SPECTRAL_INSPECTION','frequencies_hz':FREQS,'channels':[]}
 for net,sta,loc,cha in CHANNELS:
  p={'net':net,'sta':sta,'loc':loc,'cha':cha,'starttime':'2015-01-01','endtime':'2015-01-02','level':'response','format':'xml','nodata':'404'}
  item={'id':f'{net}.{sta}..{cha}'}
  try:
   r=s.get(URL,params=p,timeout=90,allow_redirects=True)
   item.update({'status':r.status_code,'bytes':len(r.content),'stationxml_sha256':sha(r.content),'content_type':r.headers.get('content-type')})
   if r.status_code!=200: item['body_prefix']=r.text[:1000];rep['channels'].append(item);continue
   inv=read_inventory(io.BytesIO(r.content),format='STATIONXML')
   ch=inv.select(network=net,station=sta,channel=cha)[0][0][0]
   item['metadata']={'sample_rate_hz':float(ch.sample_rate),'sensor_description':ch.sensor.description if ch.sensor else None,'sensitivity_value':float(ch.response.instrument_sensitivity.value) if ch.response and ch.response.instrument_sensitivity else None,'sensitivity_frequency_hz':float(ch.response.instrument_sensitivity.frequency) if ch.response and ch.response.instrument_sensitivity else None,'input_units':ch.response.instrument_sensitivity.input_units if ch.response and ch.response.instrument_sensitivity else None,'output_units':ch.response.instrument_sensitivity.output_units if ch.response and ch.response.instrument_sensitivity else None}
   vals=ch.response.get_evalresp_response_for_frequencies(FREQS,output='DEF')
   amps=[abs(complex(v)) for v in vals]
   ref=amps[FREQS.index(10.0)] if amps[FREQS.index(10.0)] else 1.0
   evals=[]
   for f,a,v in zip(FREQS,amps,vals):
    ratio=a/ref if ref else None
    db=20*math.log10(ratio) if ratio and ratio>0 else None
    evals.append({'frequency_hz':f,'response_amplitude_abs':a,'relative_to_10hz':ratio,'relative_to_10hz_db':db,'complex_real':float(complex(v).real),'complex_imag':float(complex(v).imag)})
   item['response_eval']=evals
  except Exception as e:item['error']=f'{type(e).__name__}: {e}'
  rep['channels'].append(item)

 ok=[x for x in rep['channels'] if x.get('response_eval')]
 def vals_at(f):
  return [next(e for e in x['response_eval'] if e['frequency_hz']==f)['relative_to_10hz_db'] for x in ok]
 rep['summary']={'channels_evaluated':len(ok),'response_only':True,'waveform_fft_performed':False,'db_at_100hz_vs_10hz':vals_at(100.0) if ok else [],'db_at_117hz_vs_10hz':vals_at(117.0) if ok else [],'db_at_119hz_vs_10hz':vals_at(119.0) if ok else [],'db_at_120hz_vs_10hz':vals_at(120.0) if ok else [],'target_evidence':'NOT_INCREASED','target_identity':'UNCONFIRMED'}
 rep['status']='SUCCESS_RESPONSE_EDGE_QUANTIFIED' if len(ok)==2 else 'PARTIAL_OR_BLOCKED_RESPONSE_EVAL'
 rep['hard_rules']=['INSTRUMENT_RESPONSE_IS_NOT_OCEAN_SPECTRUM','NO_FFT','NO_TARGET_PROMOTION','119_HZ_INTERPRETATION_REQUIRES_RESPONSE_CORRECTION_AND_PASSBAND_VALIDITY']
 OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(rep,indent=2,ensure_ascii=False),encoding='utf-8')
 print(json.dumps({'status':rep['status'],**rep['summary']},indent=2));return 0

if __name__=='__main__':raise SystemExit(main())
