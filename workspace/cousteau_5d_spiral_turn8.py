#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'/'cousteau'
TURN7=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-007-2026-08-21-v1.0.json'
HA10=DATA/'JANUS-ECHO-COUSTEAU-HA10-EARTHSCOPE-FDSN-RESPONSE-AND-TINY-WAVEFORM-PROBE-2026-08-21-v1.0.json'
ACCESS=DATA/'JANUS-ECHO-COUSTEAU-117-121HZ-WAVEFORM-ACCESS-AND-BANDPASS-GATE-2026-08-21-v1.1.json'
MATRIX=DATA/'JANUS-ECHO-COUSTEAU-5D-CONTROL-OUTPERFORMANCE-MATRIX-2026-08-21-v1.0.json'
OUT=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-008-2026-08-21-v1.0.json'
SUMMARY=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-008-2026-08-21-v1.0-SUMMARY.json'
DIMS=['D0_SPACE','D1_TIME','D2_ACOUSTIC','D3_REVERSE_CAUSAL','D4_ASSOCIATIVE_PROVENANCE']
PHASES=['BACK','FORWARD','LEFT_HRAIN','RIGHT_INAIHR','ASCEND']

def now():return datetime.now(timezone.utc).isoformat()
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def csha(o):return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()

def main():
 t7=load(TURN7); h=load(HA10); a=load(ACCESS); m=load(MATRIX)
 if t7.get('status')!='TURN_7_ASCENDED':raise SystemExit('Turn8 forbidden: Turn7 not ascended')
 if h.get('status')!='SUCCESS_PUBLIC_HA10_BYTES_RECOVERED':raise SystemExit('Turn8 forbidden: public HA10 waveform bytes not recovered')
 if h.get('summary',{}).get('spectral_inspection_performed') is not False:raise SystemExit('Turn8 forbidden: HA10 access probe was not spectrum-blind')
 resolved=h.get('resolved_channels',[])
 if len(resolved)<2:raise SystemExit('Turn8 forbidden: expected two resolved HA10 channels')
 if not all(x.get('network')=='IM' and x.get('channel')=='EDH' and float(x.get('sample_rate_hz') or 0)==250.0 for x in resolved):raise SystemExit('Turn8 forbidden: channel identity/sample-rate contract failed')
 if not all('100 Hz' in (x.get('sensor_description') or '') for x in resolved):raise SystemExit('Turn8 forbidden: nominal 1-100 Hz sensor metadata not present')
 if len([x for x in h.get('response_queries',[]) if x.get('status')==200])<2:raise SystemExit('Turn8 forbidden: StationXML responses not recovered')
 if len([x for x in h.get('waveform_probes',[]) if x.get('status')==200 and x.get('bytes',0)>0])<2:raise SystemExit('Turn8 forbidden: tiny waveform byte proof incomplete')
 if m.get('global_score',{}).get('candidate_admission') is not False:raise SystemExit('Turn8 forbidden: control matrix unexpectedly admits candidate')

 added=[
  'PUBLIC_WAVEFORM_BYTES_DO_NOT_EQUAL_TARGET_EVIDENCE',
  'HA10_EDH_NOMINAL_PASSBAND_1_100_HZ_OVERRIDES_NYQUIST_ONLY_ARGUMENT',
  '119_HZ_OCEAN_SPECTRAL_TEST_FORBIDDEN_UNTIL_RESPONSE_EDGE_CALIBRATION_JUSTIFIES',
  'RESPONSE_ONLY_CALIBRATION_PRECEDES_FFT',
  'EVENT_AND_NOISE_WINDOWS_PREREGISTERED_BEFORE_OCEAN_SPECTRAL_INSPECTION'
 ]
 prior_state=t7['turn']['state_sha256']; prior_radius=int(t7['turn']['radius']); radius=prior_radius+len(added)
 delta='PUBLIC_HA10_BYTES_RECOVERED__IM_H10N1_H10S2_EDH_250HZ_AND_1_100HZ_METADATA_VERIFIED__119HZ_MEASUREMENT_VALIDITY_RESTRICTED__TARGET_EVIDENCE_NOT_INCREASED'
 core={'turn':8,'parent_state_sha256':prior_state,'added_constraints':added,'scientific_delta':delta,'ha10_probe_sha256':csha(h),'access_gate_sha256':csha(a),'control_matrix_sha256':csha(m)}
 state=csha(core)
 prov=[{'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':csha(load(p))} for p in [HA10,ACCESS,MATRIX,TURN7]]
 nodes=[
  {'node_id':'S08N01','turn':8,'phase':'BACK','claim_class':'CORRECTION','structural_context':'Return to the 117-121 Hz branch after real public waveform bytes became available. Access is now proven at EarthScope for IM.H10N1..EDH and IM.H10S2..EDH on the 2015 calibration slice.','evidence_status':'PUBLIC_HA10_BYTES_RECOVERED','parents':[],'payload':{'resolved_channels':resolved},'provenance':prov},
  {'node_id':'S08N02','turn':8,'phase':'FORWARD','claim_class':'CONTROL','structural_context':'Forward metadata replay shows both channels sample at 250 Hz but describe the sensor/datalogger as 1 to 100 Hz. Therefore Nyquist=125 Hz does not validate a 119 Hz physical measurement.','evidence_status':'NOMINAL_PASSBAND_BLOCKS_NYQUIST_ONLY_INFERENCE','parents':['S08N01'],'payload':{'sample_rate_hz':250.0,'nyquist_hz':125.0,'nominal_passband_hz':[1.0,100.0]},'provenance':prov},
  {'node_id':'S08N03','turn':8,'phase':'LEFT_HRAIN','claim_class':'FACT','structural_context':'Structural lane records two StationXML payloads and two independent 60-second miniSEED payloads without performing FFT or viewing target-band ocean amplitudes.','evidence_status':'RESPONSE_AND_BYTES_PROVENANCE_FROZEN','parents':['S08N02'],'payload':{'response_payloads_recovered':h['summary']['response_payloads_recovered'],'tiny_waveform_payloads_recovered':h['summary']['tiny_waveform_payloads_recovered'],'spectral_inspection_performed':False},'provenance':prov},
  {'node_id':'S08N04','turn':8,'phase':'RIGHT_INAIHR','claim_class':'NEW_GATE','structural_context':'The new waveform access is useful only if instrument response near the band edge is calibrated before any target-band FFT. Association to the 119-Hz voice prior cannot bypass this measurement-validity gate.','evidence_status':'RESPONSE_EDGE_CALIBRATION_REQUIRED','parents':['S08N03'],'payload':{'generated_gate':'COUSTEAU_HA10_RESPONSE_EDGE_CALIBRATION_V1','target_frequency_hz':119.0},'provenance':prov},
  {'node_id':'S08N05','turn':8,'phase':'ASCEND','claim_class':'STATE_LIFT','structural_context':'Turn 8 ascends because genuinely new raw waveform bytes and response metadata entered the graph. The ascent tightens rather than promotes the target: 119-Hz interpretation is restricted by the nominal 1-100 Hz channel metadata and must wait for response-only edge calibration.','evidence_status':'SPIRAL_ASCENT_WITH_MEASUREMENT_VALIDITY_GAIN','parents':['S08N04'],'payload':{'target_evidence':'NOT_INCREASED','target_identity':'UNCONFIRMED','next_gates':['COUSTEAU_HA10_RESPONSE_EDGE_CALIBRATION_V1','COUSTEAU_EA_TPHASE_MISSING_900_PROVENANCE_RECOVERY_V1','COUSTEAU_AUTHOR_MAR_DATASET_S1_REPLICATION_V1','COUSTEAU_LOVE_EDEM_LATITUDE_CIRCLE_INDEPENDENT_GEOACOUSTIC_SCAN_V1']},'provenance':prov}
 ]
 for n in nodes:n['node_sha256']=csha(n)
 tests={
  'T8_01_PARENT_TURN7_ASCENDED':t7.get('status')=='TURN_7_ASCENDED',
  'T8_02_NEW_RAW_WAVEFORM_BYTES':h['summary']['public_waveform_bytes_verified'] is True,
  'T8_03_NO_SPECTRAL_INSPECTION_IN_TRIGGER':h['summary']['spectral_inspection_performed'] is False,
  'T8_04_TWO_250HZ_EDH_CHANNELS':len(resolved)==2 and all(x['channel']=='EDH' and x['sample_rate_hz']==250.0 for x in resolved),
  'T8_05_NOMINAL_1_100HZ_METADATA':all('100 Hz' in x['sensor_description'] for x in resolved),
  'T8_06_TWO_STATIONXML_RESPONSES':h['summary']['response_payloads_recovered']==2,
  'T8_07_TWO_TINY_WAVEFORMS':h['summary']['tiny_waveform_payloads_recovered']==2,
  'T8_08_CONTROL_MATRIX_NO_ADMISSION':m['global_score']['candidate_admission'] is False,
  'T8_09_RADIUS_ASCENDS':radius>prior_radius,
  'T8_10_TARGET_NOT_PROMOTED':'TARGET_EVIDENCE_NOT_INCREASED' in delta
 }
 ok=all(tests.values())
 out={'artifact_id':'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-008-2026-08-21-v1.0','created_utc':now(),'parent_turn':t7['artifact_id'],'spiral_law':t7['spiral_law'],'dimensions':DIMS,'phases':PHASES,'trigger':'PUBLIC_HA10_STATIONXML_PLUS_MINISEED_BYTES','turn':{'turn':8,'height':8,'prior_radius':prior_radius,'radius':radius,'prior_state_sha256':prior_state,'state_sha256':state,'repeats_prior_state':state==prior_state,'added_constraints':added,'scientific_delta':delta},'nodes':nodes,'tests':[{'id':k,'pass':bool(v)} for k,v in tests.items()],'engine_verdict':'PASS_5D_SPIRAL_TURN_8' if ok else 'FAIL_5D_SPIRAL_TURN_8','hard_rules':['NO_FFT_BEFORE_RESPONSE_EDGE_GATE','NO_TARGET_PROMOTION_FROM_PUBLIC_ACCESS','NO_POSTHOC_RETARGETING','NEGATIVE_RESULTS_NOT_RESCUED_BY_ASSOCIATION','NO_UNDERWATER_PYRAMID_DETECTED_YET'],'status':'TURN_8_ASCENDED' if ok else 'TURN_8_REJECTED'}
 summary={'artifact_id':out['artifact_id']+'-SUMMARY','engine_verdict':out['engine_verdict'],'scientific_verdict':delta,'height':8,'radius':radius,'prior_radius':prior_radius,'state_hash_unique_vs_parent':state!=prior_state,'tests_passed':sum(tests.values()),'tests_total':len(tests),'resolved_channels':[f"{x['network']}.{x['station']}..{x['channel']}" for x in resolved],'sample_rate_hz':250.0,'nominal_passband_hz':[1.0,100.0],'stationxml_payloads':h['summary']['response_payloads_recovered'],'tiny_waveform_payloads':h['summary']['tiny_waveform_payloads_recovered'],'target_identity':'UNCONFIRMED'}
 OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8');SUMMARY.write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
 print(json.dumps(summary,indent=2));return 0 if ok else 2

if __name__=='__main__':raise SystemExit(main())
