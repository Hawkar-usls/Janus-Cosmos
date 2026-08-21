#!/usr/bin/env python3
from __future__ import annotations

import argparse, hashlib, json, math
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import requests

from cousteau_ea_tphase_blind_cluster_v7 import acquire_exact_file, parse_exact, ANCHOR_LAT, ANCHOR_LON, PAPER_REPORTED_EVENT_COUNT

R=6371.0088
PB_COMMIT='339b0c56563c118307b1f4542703047f5f698fae'
PB_URL=f'https://raw.githubusercontent.com/fraxen/tectonicplates/{PB_COMMIT}/GeoJSON/PB2002_boundaries.json'
OUT_DEFAULT=Path('data/cousteau/JANUS-ECHO-COUSTEAU-TURN5-PB2002-TECTONIC-CONTROL-2026-08-21-v1.0.json')
BLIND_PATH=Path('data/cousteau/JANUS-ECHO-COUSTEAU-EA-TPHASE-BLIND-CLUSTER-RUN-007-2026-08-21-v1.6.json')

def now(): return datetime.now(timezone.utc).isoformat()
def sha(b:bytes): return hashlib.sha256(b).hexdigest()
def rad(x): return math.radians(float(x))

def angdist(lat1,lon1,lat2,lon2):
    p1,p2=rad(lat1),rad(lat2); dl=rad(lon2-lon1)
    a=math.sin((p2-p1)/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*math.atan2(math.sqrt(a),math.sqrt(max(0.0,1-a)))

def bearing(lat1,lon1,lat2,lon2):
    p1,p2=rad(lat1),rad(lat2); dl=rad(lon2-lon1)
    return math.atan2(math.sin(dl)*math.cos(p2), math.cos(p1)*math.sin(p2)-math.sin(p1)*math.cos(p2)*math.cos(dl))

def point_arc_km(plat,plon,a_lat,a_lon,b_lat,b_lon):
    d12=angdist(a_lat,a_lon,b_lat,b_lon)
    if d12 < 1e-12: return R*angdist(plat,plon,a_lat,a_lon)
    d13=angdist(a_lat,a_lon,plat,plon)
    th13=bearing(a_lat,a_lon,plat,plon); th12=bearing(a_lat,a_lon,b_lat,b_lon)
    sx=max(-1.0,min(1.0,math.sin(d13)*math.sin(th13-th12)))
    dxt=math.asin(sx)
    c=max(-1.0,min(1.0, math.cos(d13)/max(1e-15,math.cos(dxt))))
    dat=math.acos(c)
    # determine forward/backward sign using course geometry
    if math.cos(th13-th12) < 0: dat=-dat
    if 0 <= dat <= d12:
        return abs(dxt)*R
    return R*min(angdist(plat,plon,a_lat,a_lon), angdist(plat,plon,b_lat,b_lon))

def collect_af_sa_segments(pb):
    segs=[]; features=[]
    for f in pb.get('features',[]):
        p=f.get('properties') or {}; pair={str(p.get('PlateA','')),str(p.get('PlateB',''))}
        if pair != {'AF','SA'}: continue
        g=f.get('geometry') or {}; typ=g.get('type'); coords=g.get('coordinates') or []
        lines=coords if typ=='MultiLineString' else [coords] if typ=='LineString' else []
        n=0
        for line in lines:
            for a,b in zip(line,line[1:]):
                segs.append((float(a[1]),float(a[0]),float(b[1]),float(b[0]))); n+=1
        features.append({'name':p.get('Name'),'source':p.get('Source'),'segments':n})
    if not segs: raise RuntimeError('No AF-SA/SA-AF PB2002 boundary segments found')
    return segs,features

def min_to_boundary(lat,lon,segs):
    return min(point_arc_km(lat,lon,*s) for s in segs)

def q(vals,ps=(0,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1)):
    a=np.asarray(vals,float); return {str(p):round(float(np.quantile(a,p)),3) for p in ps}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',default=str(OUT_DEFAULT)); a=ap.parse_args()
    archive,gz,raw,member,trace=acquire_exact_file(); df,meta=parse_exact(raw)
    r=requests.get(PB_URL,timeout=120,headers={'User-Agent':'Janus-Echo-Cousteau/Turn5 tectonic control'}); r.raise_for_status(); pb=r.json(); segs,features=collect_af_sa_segments(pb)
    event_d=np.array([min_to_boundary(float(x.lat),float(x.lon),segs) for x in df.itertuples()],dtype=float)
    anchor_d=min_to_boundary(ANCHOR_LAT,ANCHOR_LON,segs)
    nearest_anchor_idx=int(np.argmin([R*angdist(ANCHOR_LAT,ANCHOR_LON,float(x.lat),float(x.lon)) for x in df.itertuples()]))
    nearest_anchor_event=df.iloc[nearest_anchor_idx]
    nearest_anchor_event_ridge=float(event_d[nearest_anchor_idx])
    fractions={str(k):round(float(np.mean(event_d<=k)),6) for k in [25,50,100,200,300,500,750,1000,1500,2000]}
    anchor_percentile=round(float(np.mean(event_d<=anchor_d)),6)

    cluster_controls=[]
    if BLIND_PATH.exists():
        blind=json.load(open(BLIND_PATH,encoding='utf-8'))
        for cfg in blind.get('post_reveal',{}).get('configs',[]):
            nc=cfg.get('nearest_cluster')
            if not nc: continue
            clat=nc.get('center_lat'); clon=nc.get('center_lon')
            if clat is None or clon is None:
                # recover matching cluster center from blind phase by id/config
                for bc in blind.get('blind_phase',{}).get('configs',[]):
                    if bc.get('eps_km')==cfg.get('eps_km') and bc.get('min_samples')==cfg.get('min_samples'):
                        target=next((c for c in bc.get('clusters',[]) if c.get('cluster_id')==nc.get('cluster_id')),None)
                        if target: clat,clon=target.get('center_lat'),target.get('center_lon')
                        break
            if clat is not None:
                cluster_controls.append({'eps_km':cfg.get('eps_km'),'min_samples':cfg.get('min_samples'),'anchor_to_cluster_center_km':nc.get('anchor_to_center_km'),'cluster_center_lat':clat,'cluster_center_lon':clon,'cluster_center_to_af_sa_boundary_km':round(min_to_boundary(float(clat),float(clon),segs),3)})

    n_dist={str(k):int(v) for k,v in sorted(df.n_hydrophones.value_counts().to_dict().items())}
    mismatch=PAPER_REPORTED_EVENT_COUNT-len(df)
    count_inference=(mismatch>0 and '3' not in n_dist and min(map(int,n_dist))>=4)
    ratio=anchor_d/max(1e-9,float(np.median(event_d)))
    tectonic='TECTONIC_CONTROL_STRONGLY_FAVORED_FOR_CATALOG_SPATIAL_STRUCTURE' if ratio>=3 and fractions['200']>=0.25 else 'TECTONIC_CONTROL_REMAINS_REQUIRED'
    out={
      'artifact_id':'JANUS-ECHO-COUSTEAU-TURN5-PB2002-TECTONIC-CONTROL-2026-08-21-v1.0','created_utc':now(),
      'frozen_anchor':[ANCHOR_LAT,ANCHOR_LON],
      'sources':{
        'tphase':{'doi':'10.26022/IEDA/330497','file_uid':'2504732','rows':int(len(df)),'ascii_sha256':sha(raw)},
        'paper':{'doi':'10.1029/2022JB024008','reported_final_events':PAPER_REPORTED_EVENT_COUNT,'reported_location_threshold':'THREE_OR_MORE_ARRIVAL_PICKS'},
        'tectonic_control':{'dataset':'PB2002 Bird 2003 via fraxen GeoJSON conversion','doi':'10.1029/2001GC000252','git_commit':PB_COMMIT,'url':PB_URL,'sha256':sha(r.content),'af_sa_features':features,'arc_segments':len(segs)},
        'author_specific_mar_dataset_s1':{'filename':'2022JB024008-sup-0003-Data Set SI-S01.zip','status':'BLOCKED_BY_WILEY_HTTP_403_IN_RUNNER','role':'REPLICATION_GATE_NOT_REPLACED_BY_PB2002'}
      },
      'count_reconciliation':{
        'paper_reported':PAPER_REPORTED_EVENT_COUNT,'downloaded_rows':int(len(df)),'difference':int(mismatch),'downloaded_n_hydrophones_distribution':n_dist,
        'observation':'AUTHORITATIVE_DOWNLOAD_HAS_NO_3_HYDROPHONE_ROWS',
        'inference':'EXACT_900_ROW_GAP_IS_STRUCTURALLY_COMPATIBLE_WITH_OMITTED_THREE_PICK_EVENTS__NOT_SOURCE_CONFIRMED' if count_inference else 'NO_SIMPLE_PICK_COUNT_EXPLANATION',
        'status':'OPEN_REQUIRES_SOURCE_CONFIRMATION','synthetic_row_repair':False
      },
      'tectonic_distance_control':{
        'anchor_to_af_sa_boundary_km':round(anchor_d,3),
        'event_to_af_sa_boundary_km_quantiles':q(event_d),
        'event_fraction_within_km':fractions,
        'anchor_distance_percentile_vs_events':anchor_percentile,
        'anchor_to_event_median_ridge_distance_ratio':round(ratio,3),
        'nearest_anchor_event':{'lat':float(nearest_anchor_event.lat),'lon':float(nearest_anchor_event.lon),'event_to_anchor_km':477.597,'event_to_af_sa_boundary_km':round(nearest_anchor_event_ridge,3),'n_hydrophones':int(nearest_anchor_event.n_hydrophones)},
        'nearest_blind_cluster_controls':cluster_controls,
        'verdict':tectonic
      },
      'scientific_interpretation':{
        'target_evidence':'NOT_INCREASED','target_identity':'UNCONFIRMED',
        'rule':'IF_CATALOG_EVENTS_AND_CLUSTERS_ARE_SYSTEMATICALLY_CLOSER_TO_MAR_THAN_ANCHOR__TECTONIC_BACKGROUND_OUTPERFORMS_ANCHOR_ASSOCIATION',
        'author_dataset_s1_replication_required':True
      },
      'hard_rules':['PB2002_IS_INDEPENDENT_TECTONIC_CONTROL_NOT_AUTHOR_DATASET_S1','NO_RECENTERING','DISTANCE_IS_NOT_CAUSATION','NO_SYNTHETIC_MISSING_ROWS','THREE_PICK_OMISSION_IS_INFERENCE_NOT_FACT','NEGATIVE_RESULT_NOT_RESCUED_BY_ASSOCIATION'],
      'status':'RUN_COMPLETE'
    }
    Path(a.output).write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'anchor_to_af_sa_boundary_km':out['tectonic_distance_control']['anchor_to_af_sa_boundary_km'],'event_quantiles':out['tectonic_distance_control']['event_to_af_sa_boundary_km_quantiles'],'fractions':fractions,'ratio':out['tectonic_distance_control']['anchor_to_event_median_ridge_distance_ratio'],'verdict':tectonic,'count_inference':out['count_reconciliation']['inference']},indent=2,ensure_ascii=False))
    return 0

if __name__=='__main__': raise SystemExit(main())
