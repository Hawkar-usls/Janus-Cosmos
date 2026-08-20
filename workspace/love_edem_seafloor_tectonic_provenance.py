#!/usr/bin/env python3
import json, math, os, re, tempfile, zipfile
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import requests
from netCDF4 import Dataset
import shapefile
from shapely.geometry import shape as shp_shape, Point
from shapely.ops import transform as shp_transform
from pyproj import CRS, Transformer

LAT = -3.8654180644718967
LON = 3.854924373538978
BUFFERS_KM = [0, 1, 5, 10, 25]
OUT = Path('data/love/JANUS-LOVE-EDEM-SEAFLOOR-TECTONIC-SUBMERGENCE-AND-DIRECT-DATA-RUN-001-RECEIPT.json')
UA = {'User-Agent': 'Janus-Cosmos/1.0 scientific provenance audit'}

TID_CODES = {
    10:'singlebeam direct measurement', 11:'multibeam direct measurement',
    12:'seismic direct measurement', 13:'isolated sounding',
    17:'combination of direct methods', 40:'satellite-gravity prediction',
    41:'computer interpolation', 44:'sounding-constrained/satellite-guided interpolation',
    70:'pre-generated mixed-source grid', 71:'unknown source'
}

def get(url, **kwargs):
    r = requests.get(url, headers=UA, timeout=90, **kwargs)
    r.raise_for_status()
    return r

def sample_nc(url, stem):
    p = Path(tempfile.gettempdir()) / stem
    if not p.exists() or p.stat().st_size < 10000:
        r = get(url, stream=True)
        with p.open('wb') as f:
            for chunk in r.iter_content(1024*1024):
                if chunk: f.write(chunk)
    with Dataset(p) as ds:
        vars_ = ds.variables
        def find(names):
            for n in names:
                if n in vars_: return n
            return None
        xn = find(['lon','longitude','x'])
        yn = find(['lat','latitude','y'])
        if xn is None or yn is None:
            raise RuntimeError(f'coordinate variables not found: {list(vars_)}')
        x = np.asarray(vars_[xn][:], dtype=float)
        y = np.asarray(vars_[yn][:], dtype=float)
        qlon = ((LON + 180) % 360) - 180 if np.nanmin(x) < 0 else LON % 360
        ix = int(np.nanargmin(np.abs(x-qlon)))
        iy = int(np.nanargmin(np.abs(y-LAT)))
        candidates=[]
        for name,v in vars_.items():
            if name in [xn,yn] or getattr(v,'ndim',0) < 2: continue
            if v.shape[-2:] == (len(y),len(x)):
                candidates.append(name)
        if not candidates:
            raise RuntimeError(f'2D data variable not found: {list(vars_)}')
        vn = candidates[0]
        z = vars_[vn][iy,ix]
        if np.ma.isMaskedArray(z) and bool(z.mask): val=None
        else: val=float(z)
        return {
            'url':url,'variable':vn,'nearest_grid_lon_deg':float(x[ix]),
            'nearest_grid_lat_deg':float(y[iy]),'value':val,
            'units':getattr(vars_[vn],'units',None),
            'grid_shape':[int(len(y)),int(len(x))]
        }

def gebco_tid_wms():
    base='https://wms.gebco.net/2026/mapserv'
    d=0.01
    # WMS 1.3 + EPSG:4326 axis order is lat,lon.
    common={
      'SERVICE':'WMS','VERSION':'1.3.0','REQUEST':'GetFeatureInfo',
      'LAYERS':'gebco_2026_tid_2','QUERY_LAYERS':'gebco_2026_tid_2',
      'CRS':'EPSG:4326','BBOX':f'{LAT-d},{LON-d},{LAT+d},{LON+d}',
      'WIDTH':'101','HEIGHT':'101','I':'50','J':'50','STYLES':''
    }
    attempts=[]
    for fmt in ['application/json','text/plain','text/html']:
        p=dict(common); p['INFO_FORMAT']=fmt
        try:
            r=get(base, params=p)
            text=r.text[:5000]
            attempts.append({'format':fmt,'status_code':r.status_code,'content_type':r.headers.get('content-type'),'body_head':text[:1000]})
            ints=[int(x) for x in re.findall(r'(?<![\d.])-?\d+(?![\d.])', text)]
            hits=[x for x in ints if x in TID_CODES]
            if hits:
                code=hits[0]
                return {'resolved':True,'tid_code':code,'tid_meaning':TID_CODES[code],'attempts':attempts}
        except Exception as e:
            attempts.append({'format':fmt,'error':repr(e)})
    return {'resolved':False,'tid_code':None,'attempts':attempts}

def gmrt_point():
    url='https://www.gmrt.org/services/PointServer'
    try:
        r=get(url, params={'latitude':LAT,'longitude':LON,'format':'json'})
        try: payload=r.json()
        except Exception: payload={'raw':r.text[:3000]}
        return {'resolved':True,'request_url':r.url,'payload':payload}
    except Exception as e:
        return {'resolved':False,'error':repr(e)}

def local_projection():
    wgs=CRS.from_epsg(4326)
    aeqd=CRS.from_proj4(f'+proj=aeqd +lat_0={LAT} +lon_0={LON} +datum=WGS84 +units=m +no_defs')
    return Transformer.from_crs(wgs,aeqd,always_xy=True).transform

def gmrt_shapefile_hits(url, kind):
    result={'url':url,'kind':kind,'resolved':False,'buffers':{str(k):[] for k in BUFFERS_KM}}
    try:
        zpath=Path(tempfile.gettempdir())/f'gmrt_{kind}.zip'
        if not zpath.exists() or zpath.stat().st_size < 1000:
            zpath.write_bytes(get(url).content)
        dest=Path(tempfile.gettempdir())/f'gmrt_{kind}_shp'; dest.mkdir(exist_ok=True)
        with zipfile.ZipFile(zpath) as z: z.extractall(dest)
        shpfiles=list(dest.rglob('*.shp'))
        if not shpfiles: raise RuntimeError('no .shp in zip')
        rdr=shapefile.Reader(str(shpfiles[0]))
        fields=[f[0] for f in rdr.fields[1:]]
        proj=local_projection(); pt=Point(0,0)
        for sr in rdr.iterShapeRecords():
            try:
                geom=shp_transform(proj, shp_shape(sr.shape.__geo_interface__))
                dist=float(geom.distance(pt))
            except Exception:
                continue
            rec={fields[i]:sr.record[i] for i in range(min(len(fields),len(sr.record)))}
            compact={k:v for k,v in rec.items() if v not in (None,'')}
            compact['distance_m']=dist
            for km in BUFFERS_KM:
                threshold=1.0 if km==0 else km*1000.0
                if dist <= threshold:
                    result['buffers'][str(km)].append(compact)
        for k in result['buffers']:
            result['buffers'][k]=sorted(result['buffers'][k],key=lambda x:x['distance_m'])[:50]
        result.update({'resolved':True,'shapefile':shpfiles[0].name,'fields':fields,'feature_count':len(rdr)})
    except Exception as e:
        result['error']=repr(e)
    return result

def ncei_hits():
    base='https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/multibeam_dynamic/MapServer/0/query'
    out={}
    for km in BUFFERS_KM:
        p={
          'where':'1=1','geometry':f'{LON},{LAT}','geometryType':'esriGeometryPoint',
          'inSR':'4326','spatialRel':'esriSpatialRelIntersects','outFields':'SURVEY_ID,PLATFORM,SURVEY_YEAR,SOURCE,NGDC_ID,INSTRUMENT,TRACK_LENGTH,DOWNLOAD_URL,START_TIME,END_TIME',
          'returnGeometry':'false','f':'json'
        }
        if km>0:
            p['distance']=str(km*1000); p['units']='esriSRUnit_Meter'
        try:
            js=get(base,params=p).json()
            feats=[f.get('attributes',{}) for f in js.get('features',[])]
            out[str(km)]={'resolved':True,'count':len(feats),'features':feats[:100]}
        except Exception as e:
            out[str(km)]={'resolved':False,'error':repr(e)}
    return out

def sph_to_vec(lat,lon):
    la,lo=np.radians([lat,lon]); return np.array([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)])
def vec_to_sph(v):
    v=v/np.linalg.norm(v); return float(np.degrees(np.arcsin(v[2]))), float(np.degrees(np.arctan2(v[1],v[0])))
def quat_from_euler_pole(lat,lon,angle_deg):
    k=sph_to_vec(lat,lon); h=np.radians(angle_deg)/2
    return np.r_[np.cos(h), k*np.sin(h)]
def qslerp(q0,q1,t):
    q0=q0/np.linalg.norm(q0); q1=q1/np.linalg.norm(q1); dot=float(np.dot(q0,q1))
    if dot<0: q1=-q1; dot=-dot
    dot=min(1,max(-1,dot))
    if dot>0.9995:
        q=q0+t*(q1-q0); return q/np.linalg.norm(q)
    th=math.acos(dot); return (math.sin((1-t)*th)/math.sin(th))*q0+(math.sin(t*th)/math.sin(th))*q1

def qrotate(q,v):
    w=q[0]; u=q[1:]
    return 2*np.dot(u,v)*u+(w*w-np.dot(u,u))*v+2*w*np.cross(u,v)
def q_to_euler(q):
    q=q/np.linalg.norm(q); w=float(np.clip(q[0],-1,1)); angle=2*math.acos(w)
    s=math.sin(angle/2)
    if abs(s)<1e-12: return {'pole_lat_deg':90.0,'pole_lon_deg':0.0,'angle_deg':0.0}
    axis=q[1:]/s; plat,plon=vec_to_sph(axis)
    adeg=math.degrees(angle)
    if adeg>180: adeg-=360
    return {'pole_lat_deg':plat,'pole_lon_deg':plon,'angle_deg':adeg}

def parse_rotation(text, moving, fixed):
    rows=[]
    for line in text.splitlines():
        s=line.strip()
        if not s or s.startswith('!'): continue
        parts=s.split()
        try:
            if int(parts[0])==moving and int(parts[5])==fixed:
                rows.append((float(parts[1]),float(parts[2]),float(parts[3]),float(parts[4])))
        except Exception: pass
    return sorted(rows)

def interp_rotation(rows, age):
    exact=[r for r in rows if abs(r[0]-age)<1e-9]
    if exact:
        a,la,lo,ang=exact[0]; q=quat_from_euler_pole(la,lo,ang); return q, {'lower_age_ma':a,'upper_age_ma':a,**q_to_euler(q)}
    lo_r=max((r for r in rows if r[0]<=age), default=None, key=lambda x:x[0])
    hi_r=min((r for r in rows if r[0]>=age), default=None, key=lambda x:x[0])
    if not lo_r or not hi_r: raise RuntimeError('rotation age out of range')
    t=(age-lo_r[0])/(hi_r[0]-lo_r[0])
    q=qslerp(quat_from_euler_pole(*lo_r[1:]),quat_from_euler_pole(*hi_r[1:]),t)
    return q, {'lower_age_ma':lo_r[0],'upper_age_ma':hi_r[0],**q_to_euler(q)}

def great_circle_km(a,b):
    la1,lo1=np.radians(a); la2,lo2=np.radians(b)
    c=np.sin((la2-la1)/2)**2+np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 6371.0088*2*np.arcsin(min(1,math.sqrt(float(c))))

def tectonics(age_result, rate_result, dir_result):
    age=age_result.get('value')
    result={'age_grid':age_result,'formation_full_spreading_rate_grid':rate_result,'formation_spreading_direction_grid':dir_result}
    if age is None:
        result['status']='AGE_CELL_UNRESOLVED'; return result
    rot_url='https://raw.githubusercontent.com/EarthByte/presentday-agegridding/master/AgeGridInput/Global_250-0Ma_Rotations_2019_v3.rot'
    txt=get(rot_url).text
    afr=parse_rotation(txt,701,0)
    sam=parse_rotation(txt,201,701)
    q,meta=interp_rotation(afr,age)
    paleo=vec_to_sph(qrotate(q,sph_to_vec(LAT,LON)))
    sam_meta=None
    try:
        _,sam_meta=interp_rotation(sam,age)
    except Exception as e:
        sam_meta={'error':repr(e)}
    result.update({
      'status':'AGE_AND_FIRST_ORDER_PLATE_RECONSTRUCTION_RESOLVED',
      'crust_formation_age_ma':age,
      'interpretation':'If the Seton 2020 cell is oceanic, the site was created as seafloor at a spreading ridge at approximately this age; this is not a later flooding date of pre-existing continental land.',
      'africa_plate_701_absolute_rotation_at_age':meta,
      'first_order_paleolocation_assuming_africa_plate_701':{'lat_deg':paleo[0],'lon_deg':paleo[1]},
      'great_circle_offset_paleolocation_to_present_km':great_circle_km(paleo,(LAT,LON)),
      'south_america_201_relative_to_africa_701_at_age':sam_meta,
      'rotation_caveat':'The exact oceanic static-polygon plate ID at the point is not independently resolved here. Plate 701 is an African-side first-order assignment; do not use this reconstruction before the crust-formation age or as an exact ridge-axis paleoposition without static-polygon validation.'
    })
    return result

def main():
    earth='https://earthbyte.org/webdav/ftp/earthbyte/agegrid/2020/Grids/'
    age=sample_nc(earth+'age.2020.1.GeeK2007.6m.nc','janus_age.nc')
    rate=sample_nc(earth+'full_rate.2020.1.GeeK2007.6m.nc','janus_rate.nc')
    direc=sample_nc(earth+'dir.2020.1.GeeK2007.6m.nc','janus_dir.nc')
    tid=gebco_tid_wms()
    gp=gmrt_point()
    sw=gmrt_shapefile_hits('https://gmrt.org/shapefiles/gmrt_swath_polygons.zip','swaths')
    tr=gmrt_shapefile_hits('https://gmrt.org/shapefiles/gmrt_cruise_tracks.zip','tracks')
    nc=ncei_hits()
    tect=tectonics(age,rate,direc)
    direct=False
    if tid.get('tid_code') in (10,11,12,13,17): direct=True
    if sw.get('resolved') and sw['buffers'].get('0'): direct=True
    if nc.get('0',{}).get('count',0)>0: direct=True
    receipt={
      'artifact_id':'JANUS-LOVE-EDEM-SEAFLOOR-TECTONIC-SUBMERGENCE-AND-DIRECT-DATA-RUN-001-2026-08-20-v1.0',
      'schema':'janus.cosmos.love_edem.seafloor.tectonic_submergence_direct_data.v1',
      'status':'DIRECT_DATA_PROVENANCE_RESOLVED__MORPHOLOGY_AUTHORIZED' if direct else 'TECTONIC_CLOCK_RESOLVED_OR_ATTEMPTED__DIRECT_DATA_PROVENANCE_NOT_YET_SUFFICIENT__MORPHOLOGY_BLOCKED',
      'frozen_point':{'lat_deg':LAT,'lon_deg_east':LON,'no_recenter':True},
      'frozen_buffers_km':BUFFERS_KM,
      'geological_clock':tect,
      'regional_timeline_context':{
        'south_atlantic_rifting':'Early Cretaceous; regional literature brackets major rifting roughly 145-125 Ma then stronger separation toward ~115-110 Ma.',
        'angola_margin_oceanic_crust_onset':'Published models place oceanic crust generation along Brazil-Angola sector around ~110 Ma, with northward propagation and model uncertainty.',
        'equatorial_gateway':'Shallow/intermediate circulation became possible around ~100 Ma in some reconstructions; deep-water connection developed later over several Myr.',
        'critical_distinction':'REGIONAL_MARINE_TRANSGRESSION != EXACT_POINT_SUBMERGENCE. If the exact point samples oceanic crust, its water-covered history begins with oceanic crust formation, not a later flood of continental terrain.'
      },
      'gebco_2026_tid':tid,
      'gmrt_pointserver':gp,
      'gmrt_curated_swath_polygons':sw,
      'gmrt_cruise_tracks':tr,
      'ncei_multibeam_trackline_query':nc,
      'morphology_gate':{
        'authorized':direct,
        'rule':'Only inspect geometric seabed morphology if at least one independent direct-sounding provenance route reaches the exact point. Nearby-only hits remain context, not a substitute for the point.'
      },
      'scientific_firewall':[
        'OCEANIC_CRUST_AGE_IS_NOT_AN_ARTIFACT_AGE',
        'PLATE_RECONSTRUCTION_IS_NOT_EVIDENCE_OF_ANCIENT_TECHNOLOGY',
        'REGIONAL_TRANSGRESSION_IS_NOT_A_LOCAL_FLOOD_EVENT',
        'NO_RECENTERING_TO_NEARBY_FEATURES',
        'INTERPOLATED_BATHYMETRY_CANNOT_AUTHORIZE_ARTIFICIAL_MORPHOLOGY_CLAIMS'
      ]
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(receipt,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf-8')
    print(json.dumps({'status':receipt['status'],'age_ma':tect.get('crust_formation_age_ma'),'tid':tid.get('tid_code'),'gmrt_swath_exact':len(sw.get('buffers',{}).get('0',[])) if sw.get('resolved') else None,'ncei_exact':nc.get('0',{}).get('count')},indent=2))

if __name__=='__main__': main()
