#!/usr/bin/env python3
import json, math
from pathlib import Path
import pygplates
from gplately import PlateModelManager

LAT=-3.8654180644718967
LON=3.854924373538978
AGE=95.9502944946289
OUT=Path('data/love/JANUS-LOVE-EDEM-STATIC-PLATE-ID-VALIDATION-RUN-001-RECEIPT.json')

def gc_km(lat1,lon1,lat2,lon2):
    r=6371.0088
    a1,a2=map(math.radians,[lat1,lat2]); dl=math.radians(lon2-lon1)
    c=math.sin((a2-a1)/2)**2+math.cos(a1)*math.cos(a2)*math.sin(dl/2)**2
    return r*2*math.asin(min(1,math.sqrt(c)))

def main():
    model=PlateModelManager().get_model('Muller2019', data_dir='plate-model-repo')
    rotations=model.get_rotation_model()
    static=model.get_static_polygons()
    part=pygplates.PlatePartitioner(static, rotations)
    rg=part.partition_point((LAT,LON))
    if rg is None:
        result={'resolved':False,'reason':'NO_STATIC_POLYGON_CONTAINS_POINT'}
    else:
        f=rg.get_feature()
        pid=f.get_reconstruction_plate_id()
        name=f.get_name()
        valid=f.get_valid_time()
        pt=pygplates.PointOnSphere(LAT,LON)
        rot=rotations.get_rotation(AGE,pid)
        paleo=rot*pt
        plat,plon=paleo.to_lat_lon()
        result={
          'resolved':True,'reconstruction_plate_id':int(pid),
          'static_polygon_name':name,'static_polygon_valid_time_ma':[float(valid[0]),float(valid[1])],
          'paleolocation_at_age_ma':AGE,
          'paleolat_deg':float(plat),'paleolon_deg':float(plon),
          'great_circle_offset_to_present_km':gc_km(plat,plon,LAT,LON)
        }
    receipt={
      'artifact_id':'JANUS-LOVE-EDEM-STATIC-PLATE-ID-VALIDATION-RUN-001-2026-08-20-v1.0',
      'schema':'janus.cosmos.love_edem.static_plate_id_validation.v1',
      'frozen_point':{'lat_deg':LAT,'lon_deg_east':LON,'no_recenter':True},
      'age_input_ma_GeeK2007_2arcmin':AGE,
      'model':'Muller2019 via GPlately PlateModelManager / pyGPlates static polygons',
      'result':result,
      'prior_first_order_assumption':{'plate_id':701,'paleolat_deg':-23.74977654071065,'paleolon_deg':-14.548890558166233,'offset_km':2964.6175728514113},
      'comparison_rule':'If static polygon plate ID differs from 701, demote prior first-order reconstruction and use this result for subsequent motion-path calculations. If it equals 701, the prior plate-carrier assumption is independently validated.',
      'claim_ceiling':'STATIC_PLATE_ASSIGNMENT_AND_KINEMATIC_RECONSTRUCTION_ONLY__NO_ARTIFACT_OR_RELAY_CLAIM'
    }
    receipt['status']='STATIC_PLATE_ID_RESOLVED' if result.get('resolved') else 'STATIC_PLATE_ID_UNRESOLVED'
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
