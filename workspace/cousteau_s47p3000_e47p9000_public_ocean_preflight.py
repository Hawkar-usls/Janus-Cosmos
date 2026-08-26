#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, tempfile
from datetime import datetime, timezone
from pathlib import Path
import requests

LAT=-47.3
LON=47.9
BUFFERS=[0,1,5,10,25,50]
UA={'User-Agent':'Janus-Echo-Cousteau/1.0 frozen-point provenance preflight'}
TRACK='https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/trackline_combined_dynamic/MapServer'
TRACK_LAYERS={1:'BATHYMETRY',2:'GRAVITY',3:'MAGNETICS',4:'MULTI_CHANNEL_SEISMICS',5:'SEISMIC_REFRACTION',7:'SIDE_SCAN_SONAR',8:'SINGLE_CHANNEL_SEISMICS',9:'SUBBOTTOM_PROFILE'}
MB='https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/multibeam_dynamic/FeatureServer/0'
MBD='https://gis.ngdc.noaa.gov/arcgis/rest/services/multibeam_datasets/MapServer'
GMRT_POINT='https://www.gmrt.org/services/PointServer'
GMRT_GRID='https://www.gmrt.org/services/GridServer'

def compact(attrs):
    keys=['SURVEY_ID','SURVEY_TYPE','INST_SRC','COUNTRY','PLATFORM','PROJECT','CHIEF','START_YR','END_YR','SURVEY_YEAR','DOWNLOAD_URL','DATASET_UID','FILE_UID','NAME','TITLE','CRUISE','FORMAT','FILE_NAME','URL']
    out={k:attrs.get(k) for k in keys if k in attrs and attrs.get(k) not in (None,'')}
    if not out:
        out={k:v for k,v in list(attrs.items())[:25] if k.lower() not in {'shape','objectid'}}
    return out

def arc_query(base,km):
    url=base+'/query'
    p={'where':'1=1','geometry':f'{LON},{LAT}','geometryType':'esriGeometryPoint','inSR':'4326','spatialRel':'esriSpatialRelIntersects','outFields':'*','returnGeometry':'false','f':'json','resultRecordCount':'100'}
    if km>0:
        p['distance']=str(km*1000);p['units']='esriSRUnit_Meter'
    r=requests.get(url,params=p,headers=UA,timeout=90);r.raise_for_status();js=r.json()
    if 'error' in js: raise RuntimeError(js['error'])
    return {'count':len(js.get('features',[])),'features':[compact(x.get('attributes',{})) for x in js.get('features',[])[:100]],'exceededTransferLimit':js.get('exceededTransferLimit',False),'query_url':r.url}

def gmrt_point(lat,lon):
    r=requests.get(GMRT_POINT,params={'latitude':lat,'longitude':lon,'format':'json'},headers=UA,timeout=90);r.raise_for_status()
    try:return {'resolved':True,'json':r.json(),'url':r.url}
    except Exception:return {'resolved':True,'text':r.text[:10000],'url':r.url}

def extract_z(obj):
    if isinstance(obj,(int,float)) and math.isfinite(float(obj)): return float(obj)
    if isinstance(obj,dict):
        for k in ['z','Z','elevation','Elevation','depth','Depth','value','Value']:
            if k in obj:
                z=extract_z(obj[k])
                if z is not None:return z
        for v in obj.values():
            z=extract_z(v)
            if z is not None:return z
    if isinstance(obj,list):
        for v in obj:
            z=extract_z(v)
            if z is not None:return z
    return None

def gmrt_mask_probe():
    # Frozen before results: +/-0.05 degree context, high-resolution topo-mask only.
    p={'north':LAT+0.05,'south':LAT-0.05,'west':LON-0.05,'east':LON+0.05,'layer':'topo-mask','format':'netcdf','resolution':'high'}
    r=requests.get(GMRT_GRID,params=p,headers=UA,timeout=180);r.raise_for_status()
    raw=r.content
    with tempfile.NamedTemporaryFile(suffix='.nc') as f:
        f.write(raw);f.flush()
        import netCDF4, numpy as np
        ds=netCDF4.Dataset(f.name)
        try:
            vars={k:v for k,v in ds.variables.items()}
            zvar=None
            for k,v in vars.items():
                if v.ndim>=2 and getattr(v.dtype,'kind','') in 'fiu': zvar=v;break
            if zvar is None: raise RuntimeError('no numeric 2D grid variable')
            arr=np.ma.asarray(zvar[:],dtype=float)
            finite=np.isfinite(arr.filled(np.nan))
            valid=int(finite.sum());total=int(arr.size)
            lonv=None;latv=None
            for k,v in vars.items():
                lk=k.lower()
                if v.ndim==1 and ('lon' in lk or lk in ('x',)):lonv=np.asarray(v[:],dtype=float)
                if v.ndim==1 and ('lat' in lk or lk in ('y',)):latv=np.asarray(v[:],dtype=float)
            exact=None;nearest_valid=None
            if lonv is not None and latv is not None:
                ix=int(np.argmin(abs(lonv-LON)));iy=int(np.argmin(abs(latv-LAT)))
                aa=arr.filled(np.nan)
                # accommodate [lat,lon] normal ordering
                if aa.shape[-2:]==(len(latv),len(lonv)):
                    zv=float(aa[iy,ix]) if np.isfinite(aa[iy,ix]) else None
                    exact={'grid_lat':float(latv[iy]),'grid_lon':float(lonv[ix]),'value':zv,'finite':zv is not None}
                    ys,xs=np.where(np.isfinite(aa))
                    if len(xs):
                        dy=(latv[ys]-LAT)*111.195
                        dx=(lonv[xs]-LON)*111.195*math.cos(math.radians(LAT))
                        j=int(np.argmin(dx*dx+dy*dy));nearest_valid={'distance_km_approx':float(math.hypot(dx[j],dy[j])),'lat':float(latv[ys[j]]),'lon':float(lonv[xs[j]]),'value':float(aa[ys[j],xs[j]])}
            return {'resolved':True,'http_bytes':len(raw),'content_type':r.headers.get('content-type'),'dimensions':{k:len(v) for k,v in ds.dimensions.items()},'grid_variable':zvar.name,'valid_cells':valid,'total_cells':total,'valid_fraction':valid/total if total else None,'nearest_grid_node_to_frozen_point':exact,'nearest_valid_highres_cell':nearest_valid,'request_url':r.url,'global_attributes':{k:str(ds.getncattr(k))[:500] for k in ds.ncattrs()}}
        finally:ds.close()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    rep={'artifact_id':'JANUS-COUSTEAU-S47P3000-E47P9000-PUBLIC-OCEAN-DATA-PREFLIGHT-RUN-001-2026-08-26-v1.0','schema':'janus.cosmos.cousteau.frozen_ocean_point_public_data_preflight.run.v1','created_at_utc':datetime.now(timezone.utc).isoformat(),'status':'STARTED','frozen_point':{'latitude_deg':LAT,'longitude_deg_east':LON,'no_recenter':True},'queries':{},'scientific_claims':{'morphology':False,'artificial_structure':False,'acoustic_resonance':False}}
    # GMRT exact + 8 fixed neighborhood points (not chosen from results)
    gm=[]
    offsets=[(0,0),(-.05,-.05),(-.05,0),(-.05,.05),(0,-.05),(0,.05),(.05,-.05),(.05,0),(.05,.05)]
    for dy,dx in offsets:
        try:
            q=gmrt_point(LAT+dy,LON+dx);q['offset_deg']=[dy,dx];q['lat']=LAT+dy;q['lon']=LON+dx;q['extracted_numeric_value']=extract_z(q.get('json'))
        except Exception as e:q={'resolved':False,'offset_deg':[dy,dx],'lat':LAT+dy,'lon':LON+dx,'error':f'{type(e).__name__}: {e}'}
        gm.append(q)
    rep['queries']['gmrt_pointserver']={'points':gm,'exact':gm[0]}
    try:rep['queries']['gmrt_highres_topo_mask']=gmrt_mask_probe()
    except Exception as e:rep['queries']['gmrt_highres_topo_mask']={'resolved':False,'error':f'{type(e).__name__}: {e}'}
    # NCEI multi-channel trackline catalog
    track={}
    for lid,name in TRACK_LAYERS.items():
        bb={}
        for km in BUFFERS:
            try:bb[str(km)]={'resolved':True,**arc_query(f'{TRACK}/{lid}',km)}
            except Exception as e:bb[str(km)]={'resolved':False,'error':f'{type(e).__name__}: {e}'}
        track[name]={'layer_id':lid,'buffers':bb}
    rep['queries']['ncei_trackline']=track
    # NCEI multibeam centerline survey catalog
    mb={}
    for km in BUFFERS:
        try:mb[str(km)]={'resolved':True,**arc_query(MB,km)}
        except Exception as e:mb[str(km)]={'resolved':False,'error':f'{type(e).__name__}: {e}'}
    rep['queries']['ncei_multibeam_dynamic']=mb
    # NCEI raw/processed/product dataset geometry layers
    mbd={}
    for lid,name in [(0,'PRODUCTS'),(1,'PROCESSED'),(2,'RAW')]:
        bb={}
        for km in BUFFERS:
            try:bb[str(km)]={'resolved':True,**arc_query(f'{MBD}/{lid}',km)}
            except Exception as e:bb[str(km)]={'resolved':False,'error':f'{type(e).__name__}: {e}'}
        mbd[name]={'layer_id':lid,'buffers':bb}
    rep['queries']['ncei_multibeam_datasets']=mbd
    # Summary is mechanical counts only.
    exact_track=[n for n,v in track.items() if v['buffers']['0'].get('count',0)>0]
    near_track={km:[n for n,v in track.items() if v['buffers'][str(km)].get('count',0)>0] for km in BUFFERS[1:]}
    exact_mb=mb['0'].get('count',0)
    nearest_mb_radius=next((km for km in BUFFERS if mb[str(km)].get('count',0)>0),None)
    exact_mbd=[n for n,v in mbd.items() if v['buffers']['0'].get('count',0)>0]
    nearest_mbd={n:next((km for km in BUFFERS if v['buffers'][str(km)].get('count',0)>0),None) for n,v in mbd.items()}
    mask=rep['queries']['gmrt_highres_topo_mask'];mask_exact=(mask.get('nearest_grid_node_to_frozen_point') or {}).get('finite') if mask.get('resolved') else None
    rep['summary']={'gmrt_exact_numeric_value':gm[0].get('extracted_numeric_value'),'gmrt_highres_mask_exact_finite':mask_exact,'gmrt_highres_mask_nearest_valid':mask.get('nearest_valid_highres_cell') if mask.get('resolved') else None,'ncei_exact_trackline_modalities':exact_track,'ncei_near_trackline_modalities_by_radius_km':near_track,'ncei_multibeam_centerline_exact_count':exact_mb,'ncei_multibeam_centerline_nearest_positive_radius_km':nearest_mb_radius,'ncei_multibeam_dataset_exact_layers':exact_mbd,'ncei_multibeam_dataset_nearest_positive_radius_km':nearest_mbd}
    exact_measured=bool(exact_track or exact_mb or exact_mbd or mask_exact is True)
    any_near=nearest_mb_radius is not None or any(v is not None for v in nearest_mbd.values()) or any(near_track[k] for k in near_track) or (mask.get('nearest_valid_highres_cell') is not None if mask.get('resolved') else False)
    rep['gate_state']={'G0_COORDINATE_FREEZE':'PASS','G1_GLOBAL_BATHYMETRY_CONTEXT':'PASS' if gm[0].get('resolved') else 'UNRESOLVED','G2_NCEI_TRACKLINE_COVERAGE':'EXACT_HIT' if exact_track else ('NEAR_HIT' if any(near_track[k] for k in near_track) else 'NO_HIT_WITHIN_50KM_IN_QUERIED_SERVICE'),'G3_NCEI_MULTIBEAM_COVERAGE':'EXACT_OR_HIGHRES_HIT' if (exact_mb or exact_mbd or mask_exact is True) else ('NEAR_HIT' if any_near else 'NO_HIT_WITHIN_50KM_IN_QUERIED_SERVICES'),'G4_HIGH_RESOLUTION_LOCAL_DATA':'EXACT_HIGHRES_GMRT_MASK' if mask_exact is True else ('NEAR_HIGHRES_GMRT_MASK' if mask.get('nearest_valid_highres_cell') else 'NO_HIGHRES_GMRT_MASK_CELL_IN_WINDOW_OR_UNRESOLVED'),'G5_LOCAL_DATA_FETCH':'OPEN' if exact_measured else 'HOLD','G6_MORPHOLOGY':'BLOCKED_NOT_TESTED','G7_ACOUSTIC_RESONANCE':'BLOCKED_NOT_TESTED'}
    rep['claim_ceiling']='PUBLIC_BATHYMETRY_AND_SURVEY_COVERAGE_PREFLIGHT_ONLY'
    rep['hard_rules']=['NO_RECENTERING','GLOBAL_GRID_NE_LOCAL_MEASURED_SWATH','CENTERLINE_NE_SWATH_FOOTPRINT','NEARBY_NE_EXACT','NO_MORPHOLOGY_CLASSIFICATION','NO_ARTIFICIAL_STRUCTURE_CLAIM','NO_ACOUSTIC_RESONANCE_CLAIM','NEGATIVE_RESULTS_ARE_FIRST_CLASS']
    rep['status']='COMPLETE'
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(rep,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
    print(json.dumps({'status':rep['status'],'summary':rep['summary'],'gate_state':rep['gate_state']},indent=2))
if __name__=='__main__':main()
