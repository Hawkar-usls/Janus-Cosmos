#!/usr/bin/env python3
import json, math, hashlib, pathlib, requests
import numpy as np

LAT = -3.8654180644718967
LEGACY_LON = 3.854924373538978
PB_URL = 'https://raw.githubusercontent.com/fraxen/tectonicplates/339b0c56563c118307b1f4542703047f5f698fae/GeoJSON/PB2002_boundaries.json'
PREREG='data/cousteau/JANUS-ECHO-COUSTEAU-LOVE-EDEM-LATITUDE-CIRCLE-SIDEREAL-PHASE-REVERSE-SCAN-PREREG-2026-08-21-v1.0.json'
OUT='data/cousteau/JANUS-ECHO-COUSTEAU-LOVE-EDEM-LATITUDE-CIRCLE-SIDEREAL-PHASE-REVERSE-SCAN-RUN-001-2026-08-21-v1.0.json'
R=6371.0088


def hav(lat1,lon1,lat2,lon2):
    p1,p2=np.radians(lat1),np.radians(lat2)
    dp=p2-p1; dl=np.radians(lon2-lon1)
    a=np.sin(dp/2)**2+np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
    return 2*R*np.arcsin(np.sqrt(np.clip(a,0,1)))


def bearing(lat1,lon1,lat2,lon2):
    p1,p2=map(math.radians,[lat1,lat2]); dl=math.radians(lon2-lon1)
    y=math.sin(dl)*math.cos(p2)
    x=math.cos(p1)*math.sin(p2)-math.sin(p1)*math.cos(p2)*math.cos(dl)
    return math.atan2(y,x)


def angdist(lat1,lon1,lat2,lon2):
    return float(hav(lat1,lon1,lat2,lon2))/R


def point_segment_km(lat,lon,a,b):
    lat1,lon1=a; lat2,lon2=b
    d13=angdist(lat1,lon1,lat,lon)
    if d13==0: return 0.0
    th13=bearing(lat1,lon1,lat,lon); th12=bearing(lat1,lon1,lat2,lon2)
    s=math.sin(d13)*math.sin(th13-th12)
    s=max(-1.0,min(1.0,s))
    dxt=math.asin(s)
    dat=math.atan2(math.sin(d13)*math.cos(th13-th12), math.cos(d13))
    d12=angdist(lat1,lon1,lat2,lon2)
    if 0.0 <= dat <= d12:
        return abs(dxt)*R
    return min(float(hav(lat,lon,lat1,lon1)), float(hav(lat,lon,lat2,lon2)))


def min_boundary_km(lat,lon,segs):
    return min(point_segment_km(lat,lon,a,b) for a,b in segs)


def percentile_le(values, x):
    a=np.asarray(values,float)
    return float(np.mean(a <= x))


def main():
    prereg=json.load(open(PREREG,encoding='utf-8'))
    raw=requests.get(PB_URL,timeout=60); raw.raise_for_status()
    pb_sha=hashlib.sha256(raw.content).hexdigest()
    gj=raw.json()
    segs=[]; feats=[]
    for f in gj['features']:
        p=f.get('properties',{})
        name=p.get('Name','')
        if name not in ('AF-SA','SA-AF'): continue
        g=f['geometry']
        lines=[g['coordinates']] if g['type']=='LineString' else g['coordinates']
        n=0
        for line in lines:
            for i in range(len(line)-1):
                lon1,lat1=line[i]; lon2,lat2=line[i+1]
                segs.append(((lat1,lon1),(lat2,lon2))); n+=1
        feats.append({'name':name,'source':p.get('Source'),'segments':n})
    if not segs: raise RuntimeError('No AF-SA/SA-AF segments')

    coarse=np.arange(-180.0,180.0001,0.25)
    coarse_d=np.array([min_boundary_km(LAT,float(lon),segs) for lon in coarse])
    i=int(np.argmin(coarse_d)); seed=float(coarse[i])
    refine=np.arange(seed-1.0,seed+1.0001,0.01)
    refine=((refine+180)%360)-180
    refine_d=np.array([min_boundary_km(LAT,float(lon),segs) for lon in refine])
    j=int(np.argmin(refine_d)); minlon=float(refine[j]); mind=float(refine_d[j])
    legacy_d=float(min_boundary_km(LAT,LEGACY_LON,segs))

    # Search all coarse longitudes for local minima, preserving the full look-elsewhere frame.
    locmins=[]
    for k in range(1,len(coarse)-1):
        if coarse_d[k] <= coarse_d[k-1] and coarse_d[k] <= coarse_d[k+1]:
            locmins.append({'longitude_deg_east':round(float(coarse[k]),6),'distance_km':round(float(coarse_d[k]),3)})
    locmins=sorted(locmins,key=lambda x:x['distance_km'])[:12]

    # Secondary predeclared geometry using already-frozen nearest blind cluster centers from Turn 6 input.
    turn6_path='data/cousteau/JANUS-ECHO-COUSTEAU-TURN5-PB2002-TECTONIC-CONTROL-2026-08-21-v1.0.json'
    turn6=json.load(open(turn6_path,encoding='utf-8'))
    centers=[(x['cluster_center_lat'],x['cluster_center_lon']) for x in turn6['tectonic_distance_control']['nearest_blind_cluster_controls']]
    center_curve=np.array([min(float(hav(LAT,float(lon),la,lo)) for la,lo in centers) for lon in coarse])
    ci=int(np.argmin(center_curve))

    thresholds=[25,50,100,200,500]
    out={
      'artifact_id':'JANUS-ECHO-COUSTEAU-LOVE-EDEM-LATITUDE-CIRCLE-SIDEREAL-PHASE-REVERSE-SCAN-RUN-001-2026-08-21-v1.0',
      'status':'RUN_COMPLETE',
      'preregistration':{'path':PREREG,'status':prereg['status']},
      'frozen_intrinsic_latitude_deg':LAT,
      'legacy_longitude_deg_east':LEGACY_LON,
      'sources':{
        'pb2002':{'url':PB_URL,'sha256':pb_sha,'af_sa_features':feats,'segments':len(segs)},
        'prior_turn_control':turn6_path
      },
      'primary_af_sa_scan':{
        'coarse_step_deg':0.25,
        'refinement_step_deg':0.01,
        'full_longitude_trials_coarse':int(len(coarse)),
        'minimum_distance_km':round(mind,3),
        'minimum_longitude_deg_east':round(minlon,6),
        'legacy_distance_km':round(legacy_d,3),
        'legacy_distance_percentile_low_is_close':round(percentile_le(coarse_d,legacy_d),6),
        'fraction_of_longitudes_within_km':{str(t):round(float(np.mean(coarse_d<=t)),6) for t in thresholds},
        'distance_quantiles_km':{str(q):round(float(np.quantile(coarse_d,q)),3) for q in [0,0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1]},
        'top_local_minima_coarse':locmins
      },
      'secondary_blind_cluster_center_scan':{
        'nearest_possible_distance_km_on_latitude_circle':round(float(center_curve[ci]),3),
        'longitude_deg_east':round(float(coarse[ci]),6),
        'legacy_nearest_center_distance_km':round(float(min(hav(LAT,LEGACY_LON,la,lo) for la,lo in centers)),3)
      },
      'interpretation':{
        'legacy_longitude_is_intrinsic':False,
        'latitude_circle_intersects_or_approaches_ridge_somewhere': bool(mind < 25.0),
        'ridge_intersection_counts_as_target_evidence':False,
        'previous_legacy_longitude_negative_tests_preserved':True,
        'claim_ceiling':'COORDINATE_PROVENANCE_AND_GEOMETRIC_CONTROL_ONLY'
      },
      'hard_rules':prereg['hard_rules']
    }
    pathlib.Path(OUT).parent.mkdir(parents=True,exist_ok=True)
    pathlib.Path(OUT).write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(out['primary_af_sa_scan'],indent=2))

if __name__=='__main__': main()
