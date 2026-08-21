#!/usr/bin/env python3
from __future__ import annotations
import hashlib, io, json, math
from datetime import datetime, timezone
from pathlib import Path
import requests

OUT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-HA10-EARTHSCOPE-FDSN-RESPONSE-AND-TINY-WAVEFORM-PROBE-2026-08-21-v1.0.json')
BASE='https://service.earthscope.org/fdsnws'
START='2014-12-11T00:00:00'
END='2015-01-13T00:00:00'
PROBE_START='2015-01-01T00:00:00'
PROBE_END='2015-01-01T00:01:00'
STATIONS=['H10N1','H10S2']
NETWORKS=['IM','IR']

def sha(b:bytes)->str:return hashlib.sha256(b).hexdigest()
def rmeta(r):return {'requested_url':r.request.url,'final_url':r.url,'status':r.status_code,'content_type':r.headers.get('content-type'),'content_disposition':r.headers.get('content-disposition'),'bytes':len(r.content),'sha256':sha(r.content)}

def parse_channel_text(text:str):
 rows=[]
 for line in text.splitlines():
  if not line or line.startswith('#'):continue
  p=line.split('|')
  if len(p)<16:continue
  try: sr=float(p[14])
  except Exception: sr=None
  rows.append({'network':p[0],'station':p[1],'location':p[2],'channel':p[3],'latitude':p[4],'longitude':p[5],'sensor_description':p[10] if len(p)>10 else None,'scale':p[11] if len(p)>11 else None,'scale_frequency':p[12] if len(p)>12 else None,'scale_units':p[13] if len(p)>13 else None,'sample_rate_hz':sr,'starttime':p[15] if len(p)>15 else None,'endtime':p[16] if len(p)>16 else None})
 return rows

def try_obspy(raw:bytes):
 try:
  from obspy import read
  st=read(io.BytesIO(raw))
  return [{'id':tr.id,'starttime':str(tr.stats.starttime),'endtime':str(tr.stats.endtime),'sampling_rate_hz':float(tr.stats.sampling_rate),'npts':int(tr.stats.npts)} for tr in st]
 except Exception as e:return {'error':f'{type(e).__name__}: {e}'}

def main():
 s=requests.Session();s.headers['User-Agent']='Janus-Echo-Cousteau/1.0 response-first public HA10 audit'
 rep={'artifact_id':'JANUS-ECHO-COUSTEAU-HA10-EARTHSCOPE-FDSN-RESPONSE-AND-TINY-WAVEFORM-PROBE-2026-08-21-v1.0','created_utc':datetime.now(timezone.utc).isoformat(),'policy':'RESPONSE_AND_ACCESS_ONLY__NO_OCEAN_SPECTRAL_INSPECTION','station_queries':[],'resolved_channels':[],'response_queries':[],'waveform_probes':[]}
 resolved=[]
 for sta in STATIONS:
  for net in NETWORKS:
   params={'net':net,'sta':sta,'starttime':START,'endtime':END,'level':'channel','format':'text','nodata':'404'}
   try:
    r=s.get(BASE+'/station/1/query',params=params,timeout=60,allow_redirects=True)
    item={'network_requested':net,'station_requested':sta,**rmeta(r)}
    if r.status_code==200:
     item['channels']=parse_channel_text(r.text); resolved.extend(item['channels'])
    else:item['body_prefix']=r.text[:1000]
    rep['station_queries'].append(item)
   except Exception as e:rep['station_queries'].append({'network_requested':net,'station_requested':sta,'error':f'{type(e).__name__}: {e}'})
  # If exact station code absent, probe the triplet family without using any science result.
  if not any(x.get('station')==sta for x in resolved):
   family='H10N*' if sta.startswith('H10N') else 'H10S*'
   for net in NETWORKS:
    params={'net':net,'sta':family,'starttime':START,'endtime':END,'level':'channel','format':'text','nodata':'404'}
    try:
     r=s.get(BASE+'/station/1/query',params=params,timeout=60,allow_redirects=True)
     item={'network_requested':net,'station_requested':family,'family_fallback':True,**rmeta(r)}
     if r.status_code==200:item['channels']=parse_channel_text(r.text);resolved.extend(item['channels'])
     else:item['body_prefix']=r.text[:1000]
     rep['station_queries'].append(item)
    except Exception as e:rep['station_queries'].append({'network_requested':net,'station_requested':family,'family_fallback':True,'error':f'{type(e).__name__}: {e}'})

 # De-duplicate channel identities and prioritize sample rates near the reported 250 Hz.
 uniq={}
 for x in resolved:
  key=(x['network'],x['station'],x['location'],x['channel'])
  uniq[key]=x
 resolved=list(uniq.values()); resolved.sort(key=lambda x:(0 if (x.get('sample_rate_hz') or 0)>=200 else 1,x['network'],x['station'],x['channel']))
 rep['resolved_channels']=resolved
 candidates=[x for x in resolved if (x.get('sample_rate_hz') or 0)>=200][:12]
 if not candidates:candidates=resolved[:12]

 for c in candidates:
  common={'net':c['network'],'sta':c['station'],'loc':c['location'] or '--','cha':c['channel'],'starttime':START,'endtime':END,'level':'response','format':'xml','nodata':'404'}
  try:
   rr=s.get(BASE+'/station/1/query',params=common,timeout=90,allow_redirects=True)
   rep['response_queries'].append({'channel_id':f"{c['network']}.{c['station']}.{c['location']}.{c['channel']}",**rmeta(rr),'xml_prefix':rr.text[:500] if rr.status_code==200 else rr.text[:1000]})
  except Exception as e:rep['response_queries'].append({'channel_id':f"{c['network']}.{c['station']}.{c['location']}.{c['channel']}",'error':f'{type(e).__name__}: {e}'})

  dparams={'net':c['network'],'sta':c['station'],'loc':c['location'] or '--','cha':c['channel'],'starttime':PROBE_START,'endtime':PROBE_END,'nodata':'404','format':'miniseed'}
  try:
   wr=s.get(BASE+'/dataselect/1/query',params=dparams,timeout=90,allow_redirects=True)
   wi={'channel_id':f"{c['network']}.{c['station']}.{c['location']}.{c['channel']}",**rmeta(wr),'window':{'start':PROBE_START,'end':PROBE_END}}
   if wr.status_code==200 and len(wr.content)>0:wi['trace_metadata']=try_obspy(wr.content)
   else:wi['body_prefix']=wr.text[:1000]
   rep['waveform_probes'].append(wi)
  except Exception as e:rep['waveform_probes'].append({'channel_id':f"{c['network']}.{c['station']}.{c['location']}.{c['channel']}",'error':f'{type(e).__name__}: {e}'})

 public_wave=[x for x in rep['waveform_probes'] if x.get('status')==200 and x.get('bytes',0)>0]
 response_ok=[x for x in rep['response_queries'] if x.get('status')==200 and x.get('bytes',0)>0]
 rep['summary']={'resolved_channel_count':len(resolved),'candidate_count':len(candidates),'response_payloads_recovered':len(response_ok),'tiny_waveform_payloads_recovered':len(public_wave),'public_waveform_bytes_verified':bool(public_wave),'spectral_inspection_performed':False,'target_evidence':'NOT_INCREASED','target_identity':'UNCONFIRMED'}
 rep['status']='SUCCESS_PUBLIC_HA10_BYTES_RECOVERED' if public_wave and response_ok else ('PARTIAL_METADATA_OR_RESPONSE_ONLY' if resolved or response_ok else 'BLOCKED_NO_FDSN_CHANNEL_RESOLUTION')
 rep['hard_rules']=['NO_FFT_IN_THIS_PROBE','NO_FREQUENCY_TARGET_INSPECTION','RESPONSE_FIRST','TINY_WINDOW_IS_ACCESS_PROOF_NOT_EVENT_EVIDENCE','NO_TARGET_PROMOTION']
 OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(rep,indent=2,ensure_ascii=False),encoding='utf-8')
 print(json.dumps({'status':rep['status'],**rep['summary']},indent=2));return 0

if __name__=='__main__':raise SystemExit(main())
