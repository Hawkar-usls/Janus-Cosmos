#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from datetime import datetime, timezone
from pathlib import Path
import requests

LAT=-47.3
LON=47.9
SURVEY='ANTAC23'
URL='https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/trackline_combined_dynamic/MapServer/1/query'
UA={'User-Agent':'Janus-Echo-Cousteau/1.0 frozen-point nearest-track geometry probe'}
R=6371.0088

def xy(lon,lat):
    return ((lon-LON)*111.195*math.cos(math.radians(LAT)),(lat-LAT)*111.195)

def lonlat(x,y):
    return (LON+x/(111.195*math.cos(math.radians(LAT))),LAT+y/111.195)

def nearest_on_segment(ax,ay,bx,by):
    vx=bx-ax;vy=by-ay;den=vx*vx+vy*vy
    if den<=0:return ax,ay,0.0
    t=max(0.0,min(1.0,-(ax*vx+ay*vy)/den))
    return ax+t*vx,ay+t*vy,t

def hav(lat1,lon1,lat2,lon2):
    p1,p2=map(math.radians,[lat1,lat2]);dp=math.radians(lat2-lat1);dl=math.radians(lon2-lon1)
    a=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*R*math.asin(min(1,math.sqrt(a)))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    p={'where':f"SURVEY_ID='{SURVEY}'",'outFields':'SURVEY_ID,SURVEY_TYPE,INST_SRC,COUNTRY,PLATFORM,PROJECT,CHIEF,START_YR,END_YR,SURVEY_YEAR,DOWNLOAD_URL','returnGeometry':'true','outSR':'4326','f':'json','geometryPrecision':'7','resultRecordCount':'1000'}
    r=requests.get(URL,params=p,headers=UA,timeout=180);r.raise_for_status();js=r.json()
    if 'error' in js:raise RuntimeError(js['error'])
    best=None;segments=0;vertices=0;paths_count=0
    for fi,f in enumerate(js.get('features',[])):
        geom=f.get('geometry') or {}
        for pi,path in enumerate(geom.get('paths',[])):
            paths_count+=1;vertices+=len(path)
            for si in range(len(path)-1):
                a0,b0=path[si],path[si+1]
                ax,ay=xy(float(a0[0]),float(a0[1]));bx,by=xy(float(b0[0]),float(b0[1]))
                qx,qy,t=nearest_on_segment(ax,ay,bx,by);qlon,qlat=lonlat(qx,qy)
                d=hav(LAT,LON,qlat,qlon);segments+=1
                cand=(d,fi,pi,si,t,qlat,qlon,a0,b0)
                if best is None or cand[0]<best[0]:best=cand
    if best is None:raise RuntimeError('No ANTAC23 polyline segment returned')
    d,fi,pi,si,t,qlat,qlon,a0,b0=best
    attrs=js['features'][fi].get('attributes',{})
    tier='TIER_B_VERY_NEAR_SINGLEBEAM_CONTROL_NOT_EXACT_POINT_MEASUREMENT' if d<=1 else ('TIER_B_NEAR_SINGLEBEAM_CONTROL' if d<=10 else 'TIER_C_REGIONAL_TRACKLINE_CONTEXT_ONLY')
    out={
      'artifact_id':'JANUS-COUSTEAU-S47P3000-E47P9000-ANTAC23-NEAREST-GEOMETRY-RUN-001-2026-08-26-v1.0',
      'schema':'janus.cosmos.cousteau.nearest_track_geometry.run.v1',
      'created_at_utc':datetime.now(timezone.utc).isoformat(),
      'status':'COMPLETE',
      'frozen_point':{'latitude':LAT,'longitude':LON,'recentered':False},
      'survey':attrs,
      'geometry_inventory':{'features':len(js.get('features',[])),'paths':paths_count,'vertices':vertices,'segments':segments},
      'nearest_geometry':{
        'distance_km':d,'nearest_latitude':qlat,'nearest_longitude':qlon,
        'feature_index':fi,'path_index':pi,'segment_index':si,'segment_fraction':t,
        'segment_start_lonlat':a0,'segment_end_lonlat':b0
      },
      'tier':tier,
      'bathymetry_values_read':False,
      'morphology_tested':False,
      'point_coverage_claim':False,
      'interpretation':'Nearest single-beam survey centerline geometry only. It does not measure the frozen point unless the line itself intersects the point, and it is not a multibeam swath.',
      'hard_rules':['NO_RECENTERING','NO_Z_VALUES_READ','SINGLEBEAM_TRACK_NE_MULTIBEAM_SWATH','NEARBY_TRACK_NE_EXACT_POINT_MEASUREMENT','NO_MORPHOLOGY_CLAIM']
    }
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'survey':SURVEY,'distance_km':d,'nearest':[qlat,qlon],'tier':tier,'segments':segments},indent=2))
if __name__=='__main__':main()
