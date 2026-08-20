#!/usr/bin/env python3
import json
from pathlib import Path
import requests

LAT=-3.8654180644718967
LON=3.854924373538978
BUFFERS=[0,1,5,10,25]
OUT=Path('data/love/JANUS-LOVE-EDEM-NCEI-TRACKLINE-ALL-CHANNELS-RUN-001-RECEIPT.json')
BASE='https://gis.ngdc.noaa.gov/arcgis/rest/services/web_mercator/trackline_combined_dynamic/MapServer'
LAYERS={0:'ALL_SURVEY_TYPES',1:'BATHYMETRY',2:'GRAVITY',3:'MAGNETICS',4:'MULTI_CHANNEL_SEISMICS',5:'SEISMIC_REFRACTION',7:'SIDE_SCAN_SONAR',8:'SINGLE_CHANNEL_SEISMICS',9:'SUBBOTTOM_PROFILE'}
UA={'User-Agent':'Janus-Cosmos/1.0 scientific provenance audit'}

def query(layer,km):
    url=f'{BASE}/{layer}/query'
    p={'where':'1=1','geometry':f'{LON},{LAT}','geometryType':'esriGeometryPoint','inSR':'4326','spatialRel':'esriSpatialRelIntersects','outFields':'*','returnGeometry':'false','f':'json'}
    if km>0:
        p['distance']=str(km*1000);p['units']='esriSRUnit_Meter'
    r=requests.get(url,params=p,headers=UA,timeout=90);r.raise_for_status();js=r.json()
    return [x.get('attributes',{}) for x in js.get('features',[])]

def main():
    results={}
    for lid,name in LAYERS.items():
        b={}
        for km in BUFFERS:
            try:
                f=query(lid,km); b[str(km)]={'resolved':True,'count':len(f),'features':f[:100]}
            except Exception as e:b[str(km)]={'resolved':False,'error':repr(e)}
        results[name]={'layer_id':lid,'buffers':b}
    exact_types=[name for name,r in results.items() if name!='ALL_SURVEY_TYPES' and r['buffers']['0'].get('count',0)>0]
    nearby25=[name for name,r in results.items() if name!='ALL_SURVEY_TYPES' and r['buffers']['25'].get('count',0)>0]
    receipt={'artifact_id':'JANUS-LOVE-EDEM-NCEI-TRACKLINE-ALL-CHANNELS-RUN-001-2026-08-20-v1.0','schema':'janus.cosmos.love_edem.ncei_trackline_all_channels.v1','status':'COMPLETE','frozen_point':{'lat_deg':LAT,'lon_deg_east':LON,'no_recenter':True},'frozen_buffers_km':BUFFERS,'service':BASE,'results':results,'summary':{'exact_point_data_types':exact_types,'within_25km_data_types':nearby25},'interpretation_rule':'Trackline catalog intersection establishes survey provenance only. It does not by itself establish data quality, depth resolution, or artificial morphology. Nearby lines cannot replace exact-point data.','claim_ceiling':'TRACKLINE_PROVENANCE_ONLY__NO_MORPHOLOGY_OR_ARTIFICIAL_STRUCTURE_CLAIM'}
    OUT.parent.mkdir(parents=True,exist_ok=True);OUT.write_text(json.dumps(receipt,ensure_ascii=False,indent=2,default=str)+'\n')
    print(json.dumps(receipt['summary'],indent=2))
if __name__=='__main__':main()
