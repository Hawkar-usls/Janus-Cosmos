#!/usr/bin/env python3
from __future__ import annotations

import argparse, json, math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from cousteau_ea_tphase_blind_cluster import EARTH_RADIUS_KM, haversine_km
from cousteau_ea_tphase_blind_cluster_v7 import acquire_exact_file, parse_exact, ANCHOR_LAT, ANCHOR_LON

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'/'cousteau'
OUT_DEFAULT=DATA/'JANUS-ECHO-COUSTEAU-EA-TPHASE-EVENT-ERROR-AWARE-DISTANCE-RUN-001-2026-08-21-v1.0.json'


def now(): return datetime.now(timezone.utc).isoformat()


def km_per_lon_degree(lat_deg: float) -> float:
    return (math.pi/180.0)*EARTH_RADIUS_KM*math.cos(math.radians(lat_deg))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',default=str(OUT_DEFAULT)); a=ap.parse_args()
    archive,gz,raw,member,trace=acquire_exact_file(); df,meta=parse_exact(raw)
    d=haversine_km(df.lat.to_numpy(),df.lon.to_numpy(),ANCHOR_LAT,ANCHOR_LON)
    order=np.argsort(d)
    rows=[]
    for rank,idx in enumerate(order[:25],start=1):
        r=df.iloc[int(idx)]
        lat_km=float(r.lat_error_deg)*(math.pi/180.0)*EARTH_RADIUS_KM
        lon_km=float(r.lon_error_deg)*km_per_lon_degree(float(r.lat))
        radial_proxy=math.hypot(lat_km,lon_km)
        axis_sum_proxy=abs(lat_km)+abs(lon_km)
        nominal=float(d[int(idx)])
        rows.append({
          'rank':rank,'catalog_row_zero_based':int(idx),'source_line':int(r.source_line),'source_time_code':str(r.source_time_code),
          'n_hydrophones':int(r.n_hydrophones),'lat':float(r.lat),'lon':float(r.lon),'nominal_anchor_distance_km':round(nominal,3),
          'reported_lat_error_deg':float(r.lat_error_deg),'reported_lon_error_deg':float(r.lon_error_deg),'reported_source_time_error_s':float(r.source_time_error_s),
          'source_magnitude_db':float(r.source_magnitude_db),'lat_error_km_proxy':round(lat_km,3),'lon_error_km_proxy':round(lon_km,3),
          'radial_error_extent_proxy_km':round(radial_proxy,3),'axis_sum_error_extent_proxy_km':round(axis_sum_proxy,3),
          'nominal_minus_radial_proxy_km':round(max(0.0,nominal-radial_proxy),3),
          'nominal_minus_axis_sum_proxy_km':round(max(0.0,nominal-axis_sum_proxy),3)
        })
    nearest=rows[0]
    # These are geometric proxies only; the catalog header does not state a confidence level/covariance for lat/lon errors.
    verdict='NEAREST_EVENT_REMAINS_SEPARATED_BEYOND_CONSERVATIVE_AXIS_SUM_PROXY' if nearest['nominal_minus_axis_sum_proxy_km']>0 else 'REPORTED_ERROR_EXTENT_CAN_REACH_ANCHOR_UNDER_AXIS_SUM_PROXY'
    out={
      'artifact_id':'JANUS-ECHO-COUSTEAU-EA-TPHASE-EVENT-ERROR-AWARE-DISTANCE-RUN-001-2026-08-21-v1.0',
      'created_utc':now(),'frozen_anchor':[ANCHOR_LAT,ANCHOR_LON],
      'source':{'dataset':'EA_Hydroacoustics','doi':'10.26022/IEDA/330497','file_uid':'2504732','member':member,
                'authoritative_rows':int(len(df)),'catalog_ascii_sha256':__import__('hashlib').sha256(raw).hexdigest()},
      'method':{
        'nominal_distance':'great-circle haversine','lat_error_km_proxy':'reported_lat_error_deg * pi/180 * Earth_radius',
        'lon_error_km_proxy':'reported_lon_error_deg * pi/180 * Earth_radius * cos(event_lat)',
        'radial_proxy':'sqrt(lat_km^2 + lon_km^2)','axis_sum_proxy':'abs(lat_km)+abs(lon_km)',
        'confidence_interpretation':'NONE_ASSUMED__CATALOG_HEADER_DOES_NOT_SPECIFY_COVARIANCE_OR_CONFIDENCE_LEVEL',
        'role':'conservative geometry sanity check, not a statistical confidence region'
      },
      'nearest_25':rows,
      'summary':{
        'nearest_nominal_distance_km':nearest['nominal_anchor_distance_km'],'nearest_event_lat':nearest['lat'],'nearest_event_lon':nearest['lon'],
        'nearest_n_hydrophones':nearest['n_hydrophones'],'nearest_lat_error_deg':nearest['reported_lat_error_deg'],'nearest_lon_error_deg':nearest['reported_lon_error_deg'],
        'nearest_radial_error_extent_proxy_km':nearest['radial_error_extent_proxy_km'],'nearest_axis_sum_error_extent_proxy_km':nearest['axis_sum_error_extent_proxy_km'],
        'nearest_nominal_minus_radial_proxy_km':nearest['nominal_minus_radial_proxy_km'],
        'nearest_nominal_minus_axis_sum_proxy_km':nearest['nominal_minus_axis_sum_proxy_km'],'verdict':verdict,
        'target_identity':'UNCONFIRMED'
      },
      'hard_rules':['ERROR_PROXY_IS_NOT_CONFIDENCE_ELLIPSE','NO_COVARIANCE_INFERENCE','DISTANCE_IS_NOT_CAUSATION','NEAREST_EVENT_CANNOT_RESCUE_BLIND_CLUSTER_NEGATIVE','NO_RECENTERING'],
      'status':'RUN_COMPLETE'
    }
    Path(a.output).write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(out['summary'],indent=2,ensure_ascii=False)); return 0

if __name__=='__main__': raise SystemExit(main())
