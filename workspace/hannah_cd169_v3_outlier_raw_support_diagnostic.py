#!/usr/bin/env python3
from __future__ import annotations
import argparse, math, statistics
from datetime import datetime, timezone
from pathlib import Path
import json
import hannah_bodc_9406_tobi_timestamp_locator_v3 as loc
import hannah_cd169_ground_fixed_nav_reconstruction as v1
import hannah_cd169_veh_nav_directional_repair_v2 as v2

SCIENCE=[
 ('00:00:00',-3.860,-12.118,-3.907,-12.153,6481),
 ('00:15:00',-3.849,-12.133,-3.895,-12.167,6277),
 ('00:29:00',-3.840,-12.118,-3.886,-12.151,6262),
 ('01:00:00',-3.820,-12.121,-3.869,-12.145,6032),
 ('01:29:00',-3.801,-12.124,-3.855,-12.135,6155),
 ('02:00:00',-3.781,-12.126,-3.835,-12.128,6050),
 ('02:30:00',-3.761,-12.127,-3.814,-12.128,6016),
 ('03:00:00',-3.740,-12.128,-3.795,-12.131,6112),
 ('03:30:00',-3.720,-12.130,-3.774,-12.134,6314),
 ('04:00:00',-3.699,-12.132,-3.752,-12.136,6308),
 ('04:30:00',-3.679,-12.136,-3.734,-12.140,6634),
 ('05:00:00',-3.658,-12.141,-3.715,-12.146,6586),
 ('05:30:00',-3.636,-12.149,-3.683,-12.181,6304),
 ('06:00:00',-3.615,-12.152,-3.670,-12.151,6110),
 ('06:31:00',-3.592,-12.155,-3.649,-12.154,6200)
]

def iso(t):return t.isoformat().replace('+00:00','Z')

def bearing_deg(a,b):
    p1=math.radians(a[0]);p2=math.radians(b[0]);dl=math.radians(b[1]-a[1])
    y=math.sin(dl)*math.cos(p2);x=math.cos(p1)*math.sin(p2)-math.sin(p1)*math.cos(p2)*math.cos(dl)
    return (math.degrees(math.atan2(y,x))+360)%360

def adiff(a,b):return abs((a-b+180)%360-180)

def stats(xs):
    return {'n':len(xs),'median':statistics.median(xs) if xs else None,'max':max(xs) if xs else None}

def actual_raw_header(files, starts, target):
    # Search chronologically plausible file(s), then correct index from observed timestamps.
    candidates=[]
    for i,(st,e) in enumerate(starts):
        nxt=starts[i+1][0] if i+1<len(starts) else None
        est_end=st.timestamp()+4*e['record_count']
        if st<=target and ((nxt and target<nxt) or (nxt is None and target.timestamp()<est_end+3600)):
            candidates.append((st,e))
    # Also include nearest starting file on either side as robust support audit.
    near=sorted(starts,key=lambda x:abs((x[0]-target).total_seconds()))[:2]
    by={e['relative_path']:(st,e) for st,e in candidates+near}
    best=None
    attempts=[]
    for st,e in by.values():
        idx=max(0,min(e['record_count']-1,int(round((target-st).total_seconds()/4))))
        for _ in range(2):
            s=max(0,idx-5);run=loc.stream_capture(e['path'],s,11)
            vals=[]
            for h in run['blocks']:
                d=loc.dt_of(h)
                if d:vals.append((abs((d-target).total_seconds()),d,h))
            if not vals:break
            vals.sort(key=lambda x:x[0]);delta,d,h=vals[0]
            attempts.append({'relative_path':e['relative_path'],'estimated_index':idx,'nearest_index':h['record_index'],'nearest_utc':iso(d),'abs_delta_s':delta})
            if best is None or delta<best[0]:best=(delta,d,h,e)
            signed=(d-target).total_seconds()
            if abs(signed)<=8:break
            idx=max(0,min(e['record_count']-1,h['record_index']-int(round(signed/4))))
    return best,attempts

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    out={'schema':'janus.cosmos.cousteau.hannah_cd169.v3_outlier_raw_support_diagnostic.run.v1','contract':'JANUS-HANNAH-CD169-V3-OUTLIER-AND-RAW-SUPPORT-DIAGNOSTIC-CONTRACT-2026-08-26-v1.0','status':'STARTED','role':'DIAGNOSTIC_ONLY_NOT_PROMOTION_GATE','image_or_intensity_read':False,'model_retuned':False,'outliers_dropped':False}
    try:
        files=loc.data_files();starts=[];file_rows=[]
        for e in files:
            h=loc.first_valid(e['path']);d=loc.dt_of(h)
            file_rows.append({'relative_path':e['relative_path'],'size_bytes':e['size_bytes'],'record_count':e['record_count'],'first_valid':loc.compact(h)})
            if d:starts.append((d,e))
        starts.sort(key=lambda x:x[0]);out['raw_file_chronology']=file_rows
        nr=v1.get(v1.NAV);cr=v1.get(v1.CABLE);nav=v1.parse_nav(nr);cable=v1.parse_cable(cr);v1.add_arc(nav)
        rows=[]
        for tim,slat,slon,tlat,tlon,logged_lay in SCIENCE:
            hh,mm,ss=map(int,tim.split(':'));t=datetime(2005,2,28,hh,mm,ss,tzinfo=timezone.utc)
            best,attempts=actual_raw_header(files,starts,t)
            row={'anchor_utc':iso(t),'science_log_ship':{'lat':slat,'lon':slon},'science_log_tobi':{'lat':tlat,'lon':tlon,'layback_m':logged_lay},'raw_search_attempts':attempts}
            ext=v1.ship_state(nav,t);row['external_cd169_nav_ship']={'lat':ext['lat'],'lon':ext['lon']};row['external_ship_vs_science_log_ship_m']=v1.hav_m((ext['lat'],ext['lon']),(slat,slon))
            logged_b=bearing_deg((slat,slon),(tlat,tlon));logged_d=v1.hav_m((slat,slon),(tlat,tlon));row['science_log_ship_to_tobi']={'distance_m':logged_d,'bearing_deg':logged_b}
            tags=[]
            if row['external_ship_vs_science_log_ship_m']>=1000:tags.append('SHIP_NAV_SOURCE_DISAGREEMENT_CANDIDATE')
            if best:
                delta,d,h,e=best;row['raw_support']={'relative_path':e['relative_path'],'record_index':h['record_index'],'utc':iso(d),'delta_s':(d-t).total_seconds(),'abs_delta_s':delta}
                if delta<=10:tags.append('RAW_TIME_SUPPORT_GOOD')
                elif delta>60:tags.append('RAW_TIME_SUPPORT_BOUNDARY')
                rec=v2.reconstruct_v2(nav,cable,h);rb=(h['heading_deg']+180)%360 if h.get('heading_deg') is not None else None
                row['raw_header']={'heading_deg':h.get('heading_deg'),'pressure_m_equiv':h.get('pressure_m_equiv'),'altitude_raw':h.get('altitude_raw')}
                row['reciprocal_raw_heading_deg']=rb;row['bearing_mismatch_deg']=None if rb is None else adiff(rb,logged_b)
                if row['external_ship_vs_science_log_ship_m']<=250 and row['bearing_mismatch_deg'] is not None and row['bearing_mismatch_deg']>=20:tags.append('HEADING_OR_TOW_DYNAMICS_CANDIDATE')
                if rec.get('valid'):
                    row['v2_vehicle']=rec['vehicle'];row['v2_residual_to_science_log_tobi_m']=v1.hav_m((rec['vehicle']['lat'],rec['vehicle']['lon']),(tlat,tlon));row['v2_horizontal_layback_m']=rec['horizontal_layback_m'];row['v2_layback_delta_vs_science_log_m']=rec['horizontal_layback_m']-logged_lay
            else:row['raw_support']=None;tags.append('RAW_SUPPORT_NOT_LOCATED')
            row['diagnostic_tags']=tags;rows.append(row)
        out['anchors']=rows
        key={r['anchor_utc'][11:19]:r for r in rows}
        out['focused_diagnostics']={
          '00:15':key.get('00:15:00'),
          '05:30':key.get('05:30:00'),
          '06:31':key.get('06:31:00')
        }
        good=[r['raw_support']['abs_delta_s'] for r in rows if r.get('raw_support')]
        out['summary']={'raw_files_found':len(files),'anchors':len(rows),'raw_support_abs_delta_s':stats(good),'anchors_with_good_raw_support':sum('RAW_TIME_SUPPORT_GOOD' in r['diagnostic_tags'] for r in rows),'ship_source_disagreement_candidates':[r['anchor_utc'] for r in rows if 'SHIP_NAV_SOURCE_DISAGREEMENT_CANDIDATE' in r['diagnostic_tags']],'heading_or_tow_dynamics_candidates':[r['anchor_utc'] for r in rows if 'HEADING_OR_TOW_DYNAMICS_CANDIDATE' in r['diagnostic_tags']]}
        out['status']='V3_OUTLIER_RAW_SUPPORT_DIAGNOSTIC_READY'
    except Exception as e:
        out['status']='V3_OUTLIER_RAW_SUPPORT_DIAGNOSTIC_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'summary':out.get('summary'),'focused':out.get('focused_diagnostics')},indent=2));return 0 if out['status']=='V3_OUTLIER_RAW_SUPPORT_DIAGNOSTIC_READY' else 2
if __name__=='__main__':raise SystemExit(main())
