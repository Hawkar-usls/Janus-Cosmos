#!/usr/bin/env python3
from __future__ import annotations
import argparse, bisect, ftplib, hashlib, json, math, statistics, struct
from datetime import datetime, timezone
from pathlib import Path
import hannah_cd169_ground_fixed_nav_reconstruction as v1

TARGET=v1.TARGET
FROZEN=v1.FROZEN
BLOCK=v1.BLOCK
RAW=v1.RAW
SAMPLE_M=0.75
DISCOVERY={'01:00:00','01:29:00'}
OFF=v1.OFF
R=6371008.8


def destination(lat,lon,bearing_deg,distance_m):
    p1=math.radians(lat);l1=math.radians(lon);br=math.radians(bearing_deg);d=distance_m/R
    p2=math.asin(math.sin(p1)*math.cos(d)+math.cos(p1)*math.sin(d)*math.cos(br))
    l2=l1+math.atan2(math.sin(br)*math.sin(d)*math.cos(p1),math.cos(d)-math.sin(p1)*math.sin(p2))
    return {'lat':math.degrees(p2),'lon':((math.degrees(l2)+540)%360)-180}

def reconstruct_v2(nav,cable,h):
    c=v1.interp_rows(cable,h['utc'],'cable_m');p=h['pressure_m_equiv'];total=c+200.0
    if total<=abs(p) or h['heading_deg'] is None:return {'valid':False,'cable_m':c,'pressure_m_equiv':p}
    lay=math.sqrt(total*total-p*p);ship=v1.ship_state(nav,h['utc']);bearing=(h['heading_deg']+180.0)%360;veh=destination(ship['lat'],ship['lon'],bearing,lay)
    return {'valid':True,'cable_m':c,'pressure_m_equiv':p,'total_cable_plus_umbilical_m':total,'horizontal_layback_m':lay,'raw_tobi_heading_deg':h['heading_deg'],'layback_bearing_deg':bearing,'ship':ship,'vehicle':veh}

def map_target(h,rec):
    if not rec.get('valid'):return {'valid':False}
    e,n=v1.en_offset((rec['vehicle']['lat'],rec['vehicle']['lon']),FROZEN);hr=math.radians(h['heading_deg'])
    along=e*math.sin(hr)+n*math.cos(hr);cross=e*math.cos(hr)-n*math.sin(hr);slant=math.sqrt(h['altitude_m']**2+cross**2);idx=int(round(slant/SAMPLE_M));side='stbd' if cross>=0 else 'port'
    return {'valid':True,'east_m':e,'north_m':n,'along_track_miss_m':along,'cross_track_m':cross,'ground_range_m':abs(cross),'slant_range_m':slant,'sample_index':idx,'side':side,'in_4000_samples':0<=idx<4000}

def stats(xs):
    xs=sorted(xs)
    if not xs:return {'n':0,'median_m':None,'p90_m':None,'max_m':None}
    return {'n':len(xs),'median_m':statistics.median(xs),'p90_m':v1.percentile(xs,0.9),'max_m':max(xs)}

def values_for_target(block,m):
    if not block or not m.get('in_4000_samples'):return None
    i=m['sample_index'];side=m['side'];out={}
    for name in [side+'_sidescan',side+'_swath']:
        vals=list(struct.unpack_from('<4000h',block,OFF[name]));out[name]={'sample_index':i,'raw_signed':vals[i],'abs_amplitude':abs(vals[i]),'local_z':v1.local_z(vals,i)}
    return out

def iso(x):return x.isoformat().replace('+00:00','Z')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    out={'schema':'janus.cosmos.cousteau.cd169.veh_nav_directional_repair_v2.v1','status':'STARTED','contract':'JANUS-HANNAH-CD169-VEH-NAV-DIRECTIONAL-REPAIR-V2-CONTRACT-2026-08-26-v1.0','model_generated_after_v1_failure':True,'target_period_is_discovery_not_validation':True,'parameter_retuning':False,'intensity_used_for_nav_or_mapping':False,'native_veh_nav_claim':False}
    try:
        nr=v1.get(v1.NAV);cr=v1.get(v1.CABLE);nav=v1.parse_nav(nr);cable=v1.parse_cable(cr);v1.add_arc(nav)
        headers=[];target_block=None;rawhash=hashlib.sha256();idx=0;buf=bytearray();f=v1.ftp();size=0
        try:
            sock=f.transfercmd('RETR '+RAW)
            try:
                while True:
                    chunk=sock.recv(1048576)
                    if not chunk:break
                    rawhash.update(chunk);size+=len(chunk);buf.extend(chunk)
                    while len(buf)>=BLOCK:
                        b=bytes(buf[:BLOCK]);del buf[:BLOCK];h=v1.raw_header(b,idx)
                        if h['utc']:headers.append(h)
                        if idx==9812:target_block=b
                        idx+=1
            finally:
                try:sock.close()
                except Exception:pass
        finally:
            try:f.close()
            except Exception:pass
        times=[h['utc'].timestamp() for h in headers]
        def nearest_header(t):
            j=bisect.bisect_left(times,t.timestamp());cand=[]
            if j<len(headers):cand.append(headers[j])
            if j>0:cand.append(headers[j-1])
            return min(cand,key=lambda h:abs((h['utc']-t).total_seconds()))
        rows=[]
        for tim,lat,lon,logged_layback in v1.ANCHORS:
            hh,mm,ss=map(int,tim.split(':'));t=datetime(2005,2,28,hh,mm,ss,tzinfo=timezone.utc);h=nearest_header(t)
            r1=v1.reconstruct(nav,cable,h);r2=reconstruct_v2(nav,cable,h)
            e={'anchor_utc':iso(t),'anchor_time':tim,'split':'DISCOVERY' if tim in DISCOVERY else 'HELDOUT','raw_ping_index':h['index'],'raw_ping_utc':iso(h['utc']),'logged_tobi':{'lat':lat,'lon':lon,'layback_m':logged_layback},'raw_heading_deg':h['heading_deg'],'v1_vehicle':r1.get('vehicle'),'v2_vehicle':r2.get('vehicle'),'v1_residual_m':v1.hav_m((r1['vehicle']['lat'],r1['vehicle']['lon']),(lat,lon)) if r1.get('valid') else None,'v2_residual_m':v1.hav_m((r2['vehicle']['lat'],r2['vehicle']['lon']),(lat,lon)) if r2.get('valid') else None,'v2_layback_m':r2.get('horizontal_layback_m')}
            rows.append(e)
        held=[r for r in rows if r['split']=='HELDOUT'];disc=[r for r in rows if r['split']=='DISCOVERY']
        v1s=stats([r['v1_residual_m'] for r in held if r['v1_residual_m'] is not None]);v2s=stats([r['v2_residual_m'] for r in held if r['v2_residual_m'] is not None]);dominance=(v2s['median_m']<v1s['median_m'] and v2s['p90_m']<v1s['p90_m'])
        out['integrity']={'raw_sha256':rawhash.hexdigest(),'raw_size_bytes':size,'record_count':idx,'nav_sha256':hashlib.sha256(nr).hexdigest(),'cable_sha256':hashlib.sha256(cr).hexdigest()}
        out['heldout_validation']={'anchors':held,'v1_same_heldout_stats':v1s,'v2_heldout_stats':v2s,'directional_repair_dominance_pass':dominance,'absolute_accuracy_threshold_m':None,'parameters_tuned_after_results':False}
        out['discovery_region_descriptive_only']={'anchors':disc,'excluded_from_validation':True}
        th=nearest_header(TARGET);rec=reconstruct_v2(nav,cable,th);mapping=map_target(th,rec)
        out['target_mapping_after_validation']={'raw_ping_index':th['index'],'raw_ping_utc':iso(th['utc']),'delta_to_frozen_target_s':(th['utc']-TARGET).total_seconds(),'reconstruction':rec,'mapping':mapping,'measurements_after_geometry_freeze':values_for_target(target_block,mapping)}
        if dominance:
            out['status']='V2_HELDOUT_DIRECTIONAL_REPAIR_DOMINANCE_PASS__HISTORICAL_NAV_REPRODUCTION_IMPROVED'
        else:
            out['status']='V2_HELDOUT_DIRECTIONAL_REPAIR_DOMINANCE_FAIL__MODEL_INSUFFICIENT'
        out['interpretation_ceiling']='HISTORICAL_NAVIGATION_MODEL_REPRODUCTION_ONLY__NOT_NATIVE_VEH_NAV_OR_ABSOLUTE_GROUND_TRUTH'
    except Exception as e:
        out['status']='V2_DIRECTIONAL_REPAIR_EXECUTION_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'heldout':out.get('heldout_validation'),'target':out.get('target_mapping_after_validation')},indent=2));return 0 if out['status'].startswith('V2_HELDOUT_') else 2
if __name__=='__main__':raise SystemExit(main())
