#!/usr/bin/env python3
from __future__ import annotations
import json, math, hashlib
from datetime import datetime, timezone, timedelta
from pathlib import Path

OUT=Path('data/cousteau/JANUS-FULL-SYSTEM-POINTING-COUNCIL-RUN-001-2026-08-22-v1.0.json')
LOVE=(204.30267916666668,-36.78240527777778)
EDEM=(139.22409686590188,30.26038779947318)

def jd(dt):
    y,m=dt.year,dt.month
    D=dt.day+(dt.hour+(dt.minute+(dt.second+dt.microsecond/1e6)/60)/60)/24
    if m<=2: y-=1; m+=12
    A=y//100; B=2-A+A//4
    return int(365.25*(y+4716))+int(30.6001*(m+1))+D+B-1524.5

def gmst(dt):
    J=jd(dt); T=(J-2451545.0)/36525
    return (280.46061837+360.98564736629*(J-2451545.0)+0.000387933*T*T-T*T*T/38710000)%360

def altaz(dt,lat,lon,ra,dec):
    lst=(gmst(dt)+lon)%360
    H=(lst-ra+540)%360-180
    ph,de,hr=map(math.radians,(lat,dec,H))
    alt=math.degrees(math.asin(math.sin(ph)*math.sin(de)+math.cos(ph)*math.cos(de)*math.cos(hr)))
    y=-math.sin(hr)*math.cos(de)
    x=math.sin(de)*math.cos(ph)-math.cos(de)*math.sin(ph)*math.cos(hr)
    return {'lst_deg':round(lst,6),'hour_angle_deg':round(H,6),'alt_deg':round(alt,6),'az_deg':round(math.degrees(math.atan2(y,x))%360,6)}

def clamp(x): return max(0.0,min(1.0,x))

def score(c):
    # Authority-bearing lanes
    hrain=.20*c['raw']+.16*c['source']+.16*c['location']+.10*c['time']+.18*c['controls']+.20*c['repeat']
    cosmos=.32*c['time']+.25*c['location']+.23*c['cosmos']+.20*(1-c['posthoc'])
    cousteau=.34*c['raw']+.30*c['acoustic']+.21*c['controls']+.15*c['repeat']
    cat=.30*c['testability']+.30*c['falsifiable']+.20*c['blindable']+.20*c['source']
    fundamentum=.35*c['falsifiable']+.25*c['controls']+.20*c['source']+.20*(1-c['claim_risk'])
    aifc=.35*c['time']+.25*c['location']+.20*c['source']+.20*c['blindable']
    # Association lanes are inspectable, not truth-bearing.
    inaihr=.55*c['novelty']+.45*c['cross_domain']
    aura=.60*c['symbolic']+.40*c['cross_domain']
    weighted=(.20*hrain+.17*cosmos+.18*cousteau+.14*cat+.14*fundamentum+.12*aifc+.05*inaihr)
    penalty=.22*c['posthoc']+.18*c['claim_risk']+.12*c['unverified']
    total=clamp(weighted-penalty)
    return {k:round(v,6) for k,v in {'HRain':hrain,'Cosmos':cosmos,'Cousteau':cousteau,'Fast_CAT':cat,'Fundamentum':fundamentum,'AIFC':aifc,'iNaiHR':inaihr,'Aura_symbolic_zero_authority':aura,'authority_weighted_before_penalty':weighted,'penalty':penalty,'final':total}.items()}

# Inputs are frozen from current JANUS receipts. Historical cases without exact raw acoustic chain are deliberately penalized.
candidates=[
 {'id':'EA_TPHASE_EVENT_2013_076_100526','class':'REAL_HYDROACOUSTIC_EVENT','lat':-2.412,'lon':-0.193,'utc':'2013-03-17T10:05:26Z','raw':.82,'source':1,'location':.96,'time':1,'controls':.96,'repeat':.72,'cosmos':1,'posthoc':0,'acoustic':1,'testability':.95,'falsifiable':1,'blindable':.92,'novelty':.65,'cross_domain':.90,'symbolic':.35,'claim_risk':.05,'unverified':0,
  'notes':'Nearest authoritative EA T-phase catalog event to the legacy Love-Edem point; six hydrophones; nominal distance 477.597 km; location errors do not bridge gap. Likely tectonic control, not unexplained target.'},
 {'id':'HA10_PUBLIC_WAVEFORM_SLICE_2015_01_01','class':'REAL_PUBLIC_ACOUSTIC_CALIBRATION','lat':-7.845673,'lon':-14.48023,'utc':'2015-01-01T00:00:00Z','raw':1,'source':1,'location':1,'time':1,'controls':1,'repeat':1,'cosmos':1,'posthoc':0,'acoustic':1,'testability':1,'falsifiable':1,'blindable':1,'novelty':.55,'cross_domain':.85,'symbolic':.25,'claim_risk':0,'unverified':0,
  'notes':'EarthScope IM.H10N1..EDH public miniSEED + StationXML; 250 Hz; nominal 1-100 Hz channel. Best measurement/calibration target, not an anomaly.'},
 {'id':'LOVE_EDEM_LATITUDE_CIRCLE_MAR_CROSSING','class':'GEOMETRIC_SEARCH_CORRIDOR','lat':-3.8654180645,'lon':-11.79,'utc':None,'raw':.55,'source':.93,'location':.62,'time':0,'controls':1,'repeat':.88,'cosmos':.82,'posthoc':.85,'acoustic':.88,'testability':.82,'falsifiable':.95,'blindable':.45,'novelty':.92,'cross_domain':1,'symbolic':1,'claim_risk':.40,'unverified':0,
  'notes':'Full 360 scan found AF-SA/MAR boundary minimum ~0.52 km near -11.79E at intrinsic latitude. Exact longitude is post-hoc and MUST NOT be retargeted as evidence; useful only as preregistered corridor/control.'},
 {'id':'LEGACY_LOVE_EDEM_OCEAN_POINT','class':'LEGACY_TIMESTAMP_SPECIFIC_TEST','lat':-3.865418,'lon':3.854924,'utc':None,'raw':.30,'source':.90,'location':1,'time':0,'controls':.95,'repeat':.30,'cosmos':.45,'posthoc':0,'acoustic':.55,'testability':.70,'falsifiable':1,'blindable':.85,'novelty':.70,'cross_domain':1,'symbolic':1,'claim_risk':.45,'unverified':0,
  'notes':'Preserved original timestamp-specific test point. Current T-phase spatial result is negative; do not rescue by association.'},
 {'id':'SHAG_HARBOUR_1967','class':'HISTORICAL_TRANSMEDIUM_CONTEXT','lat':43.50,'lon':-65.73,'utc':None,'date':'1967-10-04','raw':.12,'source':.68,'location':.55,'time':.35,'controls':.45,'repeat':.20,'cosmos':.30,'posthoc':0,'acoustic':.28,'testability':.38,'falsifiable':.45,'blindable':.20,'novelty':.90,'cross_domain':.78,'symbolic':.65,'claim_risk':.65,'unverified':.48,
  'notes':'Event/investigation verified; open primary sonar-tracking chain not verified. Exact UTC absent, so no admitted celestial alignment.'},
 {'id':'NIMITZ_2004_11_14','class':'HISTORICAL_MULTI_SENSOR_UAP_CONTEXT','lat':32.5,'lon':-118.0,'utc':None,'date':'2004-11-14','raw':.20,'source':.82,'location':.62,'time':.40,'controls':.55,'repeat':.55,'cosmos':.34,'posthoc':0,'acoustic':.05,'testability':.40,'falsifiable':.55,'blindable':.25,'novelty':.88,'cross_domain':.55,'symbolic':.50,'claim_risk':.60,'unverified':.25,
  'notes':'Radar/FLIR/visual context; no public sonar channel in current Cousteau record. Exact event UTC not frozen here.'},
 {'id':'TREPANG_1971','class':'DISPUTED_MILITARY_EXERCISE_CONTEXT','lat':None,'lon':None,'utc':None,'raw':.08,'source':.45,'location':.20,'time':.20,'controls':.65,'repeat':.20,'cosmos':.10,'posthoc':0,'acoustic':.12,'testability':.28,'falsifiable':.55,'blindable':.20,'novelty':.70,'cross_domain':.45,'symbolic':.45,'claim_risk':.75,'unverified':.70,
  'notes':'Disputed Arctic weapons-test context; no demonstrated public raw sonar-contact chain.'},
 {'id':'GALLAUDET_1980S_RELAYED_SONAR_ACCOUNT','class':'TESTIMONY_ARCHIVE_TARGET','lat':None,'lon':None,'utc':None,'raw':0,'source':.55,'location':0,'time':.10,'controls':.50,'repeat':.15,'cosmos':0,'posthoc':0,'acoustic':.65,'testability':.25,'falsifiable':.50,'blindable':.10,'novelty':.85,'cross_domain':.45,'symbolic':.55,'claim_risk':.65,'unverified':.75,
  'notes':'Sworn testimony relays a second-hand 1980s sonar account; raw track and exact location unavailable.'},
 {'id':'TITANIC_2004_KNOWN_TARGET','class':'KNOWN_TARGET_CALIBRATION','lat':41.7325,'lon':-49.9469444444,'utc':None,'raw':1,'source':1,'location':1,'time':.55,'controls':1,'repeat':.96,'cosmos':.55,'posthoc':0,'acoustic':1,'testability':1,'falsifiable':1,'blindable':.90,'novelty':.20,'cross_domain':.30,'symbolic':.10,'claim_risk':0,'unverified':0,
  'notes':'Known-target real-ocean benchmark with public raw multibeam and SVP; not an anomaly. Complex phase-resolved target-return gate remains open.'},
]

for c in candidates:
    c['scores']=score(c)
    if c.get('utc') and c.get('lat') is not None:
        dt=datetime.fromisoformat(c['utc'].replace('Z','+00:00'))
        c['cosmos_geometry']={'gmst_deg':round(gmst(dt),6),'LOVE':altaz(dt,c['lat'],c['lon'],*LOVE),'EDEM':altaz(dt,c['lat'],c['lon'],*EDEM),
          'rule':'ALT_AZ_IS_GEOMETRY_ONLY__NOT_CAUSAL_ASSOCIATION'}
    else:
        c['cosmos_geometry']={'status':'NOT_ADMITTED_WITHOUT_EXACT_UTC_AND_TRACEABLE_LOCATION'}

rank=sorted(candidates,key=lambda c:c['scores']['final'],reverse=True)
measurement=[c for c in rank if c['class'] in ('REAL_PUBLIC_ACOUSTIC_CALIBRATION','REAL_HYDROACOUSTIC_EVENT','KNOWN_TARGET_CALIBRATION','GEOMETRIC_SEARCH_CORRIDOR')]
historical=[c for c in rank if c['class'].startswith('HISTORICAL') or c['class'] in ('DISPUTED_MILITARY_EXERCISE_CONTEXT','TESTIMONY_ARCHIVE_TARGET')]

result={
 'artifact_id':'JANUS-FULL-SYSTEM-POINTING-COUNCIL-RUN-001-2026-08-22-v1.0',
 'created_utc':datetime.now(timezone.utc).isoformat(),
 'question':'Where should JANUS look next when dates, locations, acoustic anomalies and Cosmos geometry are evaluated together?',
 'council':{
  'HRain':'STRUCTURAL_CONTEXT','iNaiHR':'ASSOCIATIVE_CONTEXT__ASSOCIATION_NOT_EVIDENCE','DemiHead':'PRESERVE_DISAGREEMENT_AND_CLAIM_CEILING','Fast_CAT':'BLINDABLE_REVIEW_AND_NO_POSTHOC_DROPPING','Aura':'SYMBOLIC_REFLECTION_ZERO_EVIDENCE_AUTHORITY','Janus_Cosmos':'CELESTIAL_GEOMETRY_AND_ANTI_PSEUDOREPLICATION','Cousteau':'REAL_ACOUSTIC_DATA_AND_CONTROLS','Fundamentum':'FALSIFICATION_AND_NEGATIVE_RESULT_PRESERVATION','AIFC':'TIME_LOCATION_PROVENANCE_AND_FAIL_CLOSED_ADMISSION'},
 'weights_rule':'AURA_HAS_ZERO_AUTHORITY_WEIGHT__INAIHR_SMALL_ROUTING_WEIGHT__EVIDENCE_LANES_DOMINATE',
 'ranked_candidates':rank,
 'system_pointing':{
  'LOOK_HERE_NEXT':measurement[0]['id'],
  'WHY':measurement[0]['notes'],
  'BEST_REAL_EVENT_FOR_DATE_LOCATION_COSMOS_CROSSCHECK':'EA_TPHASE_EVENT_2013_076_100526',
  'BEST_HISTORICAL_ANOMALY_ARCHIVE_TARGET':historical[0]['id'] if historical else None,
  'SECONDARY_CALIBRATION_TARGET':'TITANIC_2004_KNOWN_TARGET',
  'MAR_CORRIDOR_RULE':'Use the latitude-circle/MAR crossing only as a preregistered corridor/control; never promote the post-hoc -11.79 longitude itself as target evidence.'
 },
 'hard_rules':['EXACT_UTC_REQUIRED_FOR_CELESTIAL_ALIGNMENT','ALT_AZ_ALIGNMENT_IS_NOT_CAUSATION','RAW_ACOUSTIC_DATA_OUTRANKS_HISTORICAL_NARRATIVE','POSTHOC_LONGITUDE_CANNOT_BE_RETARGETED','NEGATIVE_RESULTS_REMAIN_NEGATIVE','AURA_SYMBOLIC_OUTPUT_HAS_ZERO_EVIDENCE_AUTHORITY','NO_UNDERWATER_PYRAMID_DETECTED'],
 'target_identity':'UNCONFIRMED'
}
result['sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'LOOK_HERE_NEXT':result['system_pointing']['LOOK_HERE_NEXT'],'BEST_REAL_EVENT':result['system_pointing']['BEST_REAL_EVENT_FOR_DATE_LOCATION_COSMOS_CROSSCHECK'],'BEST_HISTORICAL':result['system_pointing']['BEST_HISTORICAL_ANOMALY_ARCHIVE_TARGET'],'top5':[(c['id'],c['scores']['final']) for c in rank[:5]],'sha256':result['sha256']},indent=2,ensure_ascii=False))
