#!/usr/bin/env python3
from __future__ import annotations

import gzip, hashlib, io, json, math, re, tarfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

import requests
from bs4 import BeautifulSoup

ROOT=Path(__file__).resolve().parents[1]
PROTOCOL=ROOT/'data'/'cousteau'/'JANUS-ECHO-COUSTEAU-HA10-PREDECLARED-RESPONSE-CORRECTED-CONFIRMATORY-PROTOCOL-2026-08-22-v1.0.json'
OUT=ROOT/'data'/'cousteau'/'JANUS-ECHO-COUSTEAU-HA10-CONFIRMATORY-WINDOW-FREEZE-2026-08-22-v1.0.json'
DATASET_UID='30497'; FILE_UID='2504732'
LANDING=f'https://www.marine-geo.org/tools/files/{DATASET_UID}'
UID_URL=f'https://www.marine-geo.org/tools/search/file_uids.php?data_set_uid={DATASET_UID}'
MODAL_URL='https://www.marine-geo.org/services/download/download_modal.php'
ROW_RE=re.compile(r'^\s*(\d{14})\s+(\d+)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s*$')


def sha(b:bytes)->str: return hashlib.sha256(b).hexdigest()
def csha(o)->str: return sha(json.dumps(o,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode())
def iso(dt:datetime)->str: return dt.astimezone(timezone.utc).isoformat().replace('+00:00','Z')

def parse_code(code:str)->datetime:
    # YYYY DDD HH MM SSS, where SSS is seconds with one decimal place.
    y=int(code[:4]); doy=int(code[4:7]); hh=int(code[7:9]); mm=int(code[9:11]); sec10=int(code[11:14])
    sec=sec10/10.0
    return datetime(y,1,1,tzinfo=timezone.utc)+timedelta(days=doy-1,hours=hh,minutes=mm,seconds=sec)

def hav_km(lat1,lon1,lat2,lon2,r=6371.0088):
    p1,p2=math.radians(lat1),math.radians(lat2); dp=math.radians(lat2-lat1); dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*r*math.asin(min(1.0,math.sqrt(a)))

def acquire():
    s=requests.Session(); s.headers.update({'User-Agent':'Janus-Echo-Cousteau/1.0 confirmatory window freeze','Referer':LANDING})
    u=s.get(UID_URL,timeout=45); u.raise_for_status(); uids=[str(x) for x in u.json()]
    if FILE_UID not in uids: raise RuntimeError(f'file UID drift: {uids}')
    m=s.post(MODAL_URL,data={'FileDownload':FILE_UID,'data_set_uid':DATASET_UID},timeout=60); m.raise_for_status()
    soup=BeautifulSoup(m.text,'html.parser'); f=soup.find('form',id='data_link') or soup.find('form')
    if not f: raise RuntimeError('MGDS download form missing')
    r=s.post(f.get('action'),data={'purpose':'Research','client':'DataLink','force_download':'1','data_uids':FILE_UID},timeout=180); r.raise_for_status()
    target=None
    with tarfile.open(fileobj=io.BytesIO(r.content),mode='r:*') as tf:
        for mem in tf.getmembers():
            if mem.isfile() and mem.name.endswith('EA_CTBTO_catalog_all.dat.gz'):
                fo=tf.extractfile(mem)
                if fo: target=(mem.name,fo.read())
    if target is None: raise RuntimeError('catalog gzip not found')
    name,gz=target; raw=gzip.decompress(gz)
    return {'archive_sha256':sha(r.content),'archive_bytes':len(r.content),'member':name,'member_gzip_sha256':sha(gz),'catalog_ascii_sha256':sha(raw),'catalog_ascii_bytes':len(raw)},raw

def parse(raw:bytes):
    rows=[]
    for lineno,line in enumerate(raw.decode('utf-8').splitlines(),1):
        m=ROW_RE.match(line)
        if not m: continue
        t,n,lat,lon,laterr,lonerr,terr,mag=m.groups()
        rows.append({'source_time_code':t,'origin_utc':parse_code(t),'n_hydrophones':int(n),'lat':float(lat),'lon':float(lon),'lat_error_deg':float(laterr),'lon_error_deg':float(lonerr),'source_time_error_s':float(terr),'source_magnitude_db':float(mag),'source_line':lineno})
    if len(rows)<1000: raise RuntimeError(f'parser rows too small: {len(rows)}')
    return rows

def main():
    protocol=json.loads(PROTOCOL.read_text())
    if protocol.get('status')!='PREREGISTERED_BEFORE_TARGET_BAND_WAVEFORM_INSPECTION': raise RuntimeError('protocol is not frozen preregistration')
    meta,raw=acquire(); rows=parse(raw)
    start=datetime.fromisoformat(protocol['public_analysis_interval_utc']['start'].replace('Z','+00:00'))
    end=datetime.fromisoformat(protocol['public_analysis_interval_utc']['end'].replace('Z','+00:00'))
    ac=protocol['arrival_model']; speed=float(ac['nominal_sofar_group_speed_km_s']); half=float(ac['event_window_half_width_s'])
    stations=protocol['channels']; elig=protocol['event_selection_contract']['eligibility']

    # Compute predicted arrivals for all catalog events. This stage never downloads waveform samples and never computes spectra.
    for row in rows:
        row['predicted_arrivals']={}
        for st in stations:
            d=hav_km(row['lat'],row['lon'],st['latitude_deg'],st['longitude_deg'])
            arr=row['origin_utc']+timedelta(seconds=d/speed)
            row['predicted_arrivals'][st['id']]={'distance_km':d,'arrival_utc':arr}

    candidates=[]
    for row in rows:
        if not (start <= row['origin_utc'] < end): continue
        if row['n_hydrophones'] < int(elig['minimum_n_hydrophones']): continue
        if row['source_time_error_s'] > float(elig['maximum_source_time_error_s']): continue
        if row['lat_error_deg'] > float(elig['maximum_lat_error_deg']): continue
        if row['lon_error_deg'] > float(elig['maximum_lon_error_deg']): continue
        if not all(start <= row['predicted_arrivals'][st['id']]['arrival_utc'] < end for st in stations): continue
        candidates.append(row)
    candidates.sort(key=lambda r:(-r['source_magnitude_db'],r['source_time_error_s'],r['source_time_code']))
    selected=candidates[:int(protocol['event_selection_contract']['maximum_selected_events'])]

    # Precompute all predicted arrivals for contamination checks by station, restricted to public interval.
    arrivals_by_station={st['id']:[] for st in stations}
    for row in rows:
        for st in stations:
            t=row['predicted_arrivals'][st['id']]['arrival_utc']
            if start <= t < end: arrivals_by_station[st['id']].append((t,row['source_time_code']))

    offsets=list(protocol['matched_noise_contract']['preferred_center_offsets_s'])+list(protocol['matched_noise_contract']['fallback_center_offsets_s'])
    frozen=[]
    for rank,row in enumerate(selected,1):
        ev={'selection_rank':rank,'source_time_code':row['source_time_code'],'origin_utc':iso(row['origin_utc']),'source_line':row['source_line'],'n_hydrophones':row['n_hydrophones'],'lat':row['lat'],'lon':row['lon'],'lat_error_deg':row['lat_error_deg'],'lon_error_deg':row['lon_error_deg'],'source_time_error_s':row['source_time_error_s'],'source_magnitude_db':row['source_magnitude_db'],'stations':{}}
        pair_complete=True
        for st in stations:
            sid=st['id']; arr=row['predicted_arrivals'][sid]['arrival_utc']; event_start=arr-timedelta(seconds=half); event_end=arr+timedelta(seconds=half)
            noise=[]
            for off in offsets:
                c=arr+timedelta(seconds=float(off)); ns=c-timedelta(seconds=half); ne=c+timedelta(seconds=half)
                if ns < start or ne >= end: continue
                contaminated=False
                for other_t,other_code in arrivals_by_station[sid]:
                    if other_code==row['source_time_code']: continue
                    if abs((other_t-c).total_seconds()) <= 600: contaminated=True; break
                if contaminated: continue
                noise.append({'center_offset_s':off,'start_utc':iso(ns),'center_utc':iso(c),'end_utc':iso(ne)})
                if len(noise)==2: break
            if len(noise)<2: pair_complete=False
            ev['stations'][sid]={'station_lat':st['latitude_deg'],'station_lon':st['longitude_deg'],'source_to_station_distance_km':round(row['predicted_arrivals'][sid]['distance_km'],6),'predicted_arrival_utc':iso(arr),'event_window':{'start_utc':iso(event_start),'end_utc':iso(event_end)},'noise_windows':noise,'pair_complete':len(noise)==2}
        ev['complete_on_both_stations']=pair_complete
        frozen.append(ev)

    complete=[x for x in frozen if x['complete_on_both_stations']]
    freeze_core={'protocol_artifact_id':protocol['artifact_id'],'catalog_ascii_sha256':meta['catalog_ascii_sha256'],'selected_complete_events':complete}
    freeze_hash=csha(freeze_core)
    minimum=int(protocol['event_selection_contract']['minimum_events_required_for_confirmatory_execution'])
    status='WINDOWS_FROZEN_READY_FOR_FFT' if len(complete)>=minimum else 'BLOCKED_UNDERPOWERED__DO_NOT_RELAX_SELECTION_CRITERIA'
    out={'artifact_id':'JANUS-ECHO-COUSTEAU-HA10-CONFIRMATORY-WINDOW-FREEZE-2026-08-22-v1.0','created_utc':datetime.now(timezone.utc).isoformat(),'protocol_path':str(PROTOCOL.relative_to(ROOT)),'protocol_blob_sha_not_embedded':'VERIFY_BY_GIT_PROVENANCE','source':meta,'catalog_rows':len(rows),'eligible_event_count_before_topN':len(candidates),'selected_event_count':len(selected),'complete_event_count':len(complete),'minimum_required':minimum,'selected_events':frozen,'fft_performed':False,'waveform_samples_inspected':False,'target_band_inspected':False,'window_freeze_sha256':freeze_hash,'status':status,'hard_rules':['NO_WAVEFORM_SPECTRAL_INSPECTION_IN_WINDOW_FREEZE','NO_TARGET_BAND_INSPECTION','NO_THRESHOLD_RETUNING','WINDOW_FREEZE_HASH_PRECEDES_FFT']}
    OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'status':status,'catalog_rows':len(rows),'eligible':len(candidates),'selected':len(selected),'complete':len(complete),'window_freeze_sha256':freeze_hash},indent=2))
    return 0 if status=='WINDOWS_FROZEN_READY_FOR_FFT' else 2

if __name__=='__main__': raise SystemExit(main())
