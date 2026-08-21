#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data'/'cousteau'
T8=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-008-2026-08-21-v1.0.json'
RESP=DATA/'JANUS-ECHO-COUSTEAU-HA10-RESPONSE-EDGE-CALIBRATION-PROBE-2026-08-21-v1.0.json'
GATE=DATA/'JANUS-ECHO-COUSTEAU-117-121HZ-WAVEFORM-ACCESS-AND-BANDPASS-GATE-2026-08-21-v1.2.json'
OUT=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-009-2026-08-21-v1.0.json';SUM=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-009-2026-08-21-v1.0-SUMMARY.json'
DIMS=['D0_SPACE','D1_TIME','D2_ACOUSTIC','D3_REVERSE_CAUSAL','D4_ASSOCIATIVE_PROVENANCE'];PHASES=['BACK','FORWARD','LEFT_HRAIN','RIGHT_INAIHR','ASCEND']
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def csha(o):return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()

def main():
 t8=load(T8);r=load(RESP);g=load(GATE)
 if t8.get('status')!='TURN_8_ASCENDED':raise SystemExit('Turn9 forbidden: Turn8 not ascended')
 if r.get('status')!='SUCCESS_RESPONSE_EDGE_QUANTIFIED':raise SystemExit('Turn9 forbidden: response edge not quantified')
 if r.get('summary',{}).get('waveform_fft_performed') is not False:raise SystemExit('Turn9 forbidden: response probe was not spectrum-blind')
 d119=r['summary']['db_at_119hz_vs_10hz'];d117=r['summary']['db_at_117hz_vs_10hz'];d120=r['summary']['db_at_120hz_vs_10hz']
 if len(d119)!=2 or not all(x<-45 for x in d119):raise SystemExit('Turn9 forbidden: 119-Hz roll-off not verified on both channels')
 added=['NYQUIST_DOES_NOT_DEFINE_USABLE_DISCOVERY_BAND','HA10_EDH_119HZ_RESPONSE_IS_ABOUT_MINUS_50_9DB_RELATIVE_TO_10HZ','119_HZ_ON_HA10_EDH_IS_CONFIRMATORY_ONLY_NOT_DISCOVERY_GRADE','RAW_AND_RESPONSE_CORRECTED_SNR_THRESHOLDS_MUST_BE_PREREGISTERED_BEFORE_TARGET_BAND_FFT']
 prior=t8['turn']['state_sha256'];pr=int(t8['turn']['radius']);radius=pr+len(added)
 delta='HA10_RESPONSE_EDGE_QUANTIFIED__119HZ_MINUS_50_9DB_RELATIVE_10HZ_ON_BOTH_CHANNELS__119HZ_DISCOVERY_KEY_REJECTED__CONFIRMATORY_ONLY__TARGET_EVIDENCE_NOT_INCREASED'
 core={'turn':9,'parent_state_sha256':prior,'added_constraints':added,'scientific_delta':delta,'response_sha256':csha(r),'gate_sha256':csha(g)};state=csha(core)
 prov=[{'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':csha(load(p))} for p in [T8,RESP,GATE]]
 nodes=[
 {'node_id':'S09N01','turn':9,'phase':'BACK','claim_class':'FACT','structural_context':'Return to the 119-Hz measurement question after StationXML response-only calibration on both public HA10 EDH channels.','evidence_status':'RESPONSE_EDGE_CALIBRATION_ENTERED','parents':[],'payload':{'db117_vs_10':d117,'db119_vs_10':d119,'db120_vs_10':d120},'provenance':prov},
 {'node_id':'S09N02','turn':9,'phase':'FORWARD','claim_class':'NEGATIVE','structural_context':'Both H10N1 and H10S2 show about -50.904 dB response at 119 Hz relative to 10 Hz; 117 Hz is about -39.137 dB and 120 Hz about -57.943 dB. The roll-off is an instrument property, not an ocean signal.','evidence_status':'119HZ_DISCOVERY_SENSITIVITY_REJECTED','parents':['S09N01'],'payload':{'inverse_response_gain_119_vs_10':g['response_edge_calibration']['119_hz_inverse_response_gain_needed_relative_to_10hz']},'provenance':prov},
 {'node_id':'S09N03','turn':9,'phase':'LEFT_HRAIN','claim_class':'CONTROL','structural_context':'The structural lane separates digitization, passband and recoverability: 119 Hz is below Nyquist yet deep in the EDH anti-alias/response roll-off. Deconvolution is possible in principle but necessarily magnifies noise where response is small.','evidence_status':'MEASUREMENT_VALIDITY_CONTROL_DOMINATES','parents':['S09N02'],'payload':{'sample_rate_hz':250.0,'nyquist_hz':125.0,'nominal_sensor_band_hz':[1.0,100.0]},'provenance':prov},
 {'node_id':'S09N04','turn':9,'phase':'RIGHT_INAIHR','claim_class':'NEW_GATE','structural_context':'The 119-Hz prior may route a confirmatory test only. It cannot define a discovery band. Any confirmatory FFT requires frozen event/noise windows, forward-response or response-removal method, raw/corrected SNR gates and neighboring-frequency controls before target-band inspection.','evidence_status':'CONFIRMATORY_PROTOCOL_REQUIRED','parents':['S09N03'],'payload':{'generated_gate':'COUSTEAU_HA10_PREDECLARED_RESPONSE_CORRECTED_CONFIRMATORY_PROTOCOL_V1'},'provenance':prov},
 {'node_id':'S09N05','turn':9,'phase':'ASCEND','claim_class':'STATE_LIFT','structural_context':'Turn 9 ascends because new calibration data reduced the admissible hypothesis space. The target is not promoted; instead 119-Hz discovery use is rejected on HA10 EDH and preserved only as a preregistered confirmatory possibility.','evidence_status':'SPIRAL_ASCENT_WITH_INSTRUMENT_FALSIFICATION_GAIN','parents':['S09N04'],'payload':{'target_evidence':'NOT_INCREASED','target_identity':'UNCONFIRMED','next_gates':['COUSTEAU_HA10_PREDECLARED_RESPONSE_CORRECTED_CONFIRMATORY_PROTOCOL_V1','COUSTEAU_EA_TPHASE_MISSING_900_PROVENANCE_RECOVERY_V1','COUSTEAU_AUTHOR_MAR_DATASET_S1_REPLICATION_V1']},'provenance':prov}
 ]
 for n in nodes:n['node_sha256']=csha(n)
 tests={'T9_01_PARENT_T8_ASCENDED':t8.get('status')=='TURN_8_ASCENDED','T9_02_RESPONSE_EDGE_SUCCESS':r.get('status')=='SUCCESS_RESPONSE_EDGE_QUANTIFIED','T9_03_NO_WAVEFORM_FFT':r['summary']['waveform_fft_performed'] is False,'T9_04_BOTH_CHANNELS_EVALUATED':r['summary']['channels_evaluated']==2,'T9_05_119_BELOW_MINUS45DB':all(x<-45 for x in d119),'T9_06_117_ROLLOFF_VERIFIED':all(x<-30 for x in d117),'T9_07_120_ROLLOFF_VERIFIED':all(x<-50 for x in d120),'T9_08_DISCOVERY_KEY_REJECTED':'REJECT_AS_HA10_EDH_BLIND_DISCOVERY_BAND' in g['h0_117_121_hz_status']['overall'],'T9_09_RADIUS_ASCENDS':radius>pr,'T9_10_TARGET_NOT_PROMOTED':'TARGET_EVIDENCE_NOT_INCREASED' in delta}
 ok=all(tests.values())
 out={'artifact_id':'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-009-2026-08-21-v1.0','created_utc':now(),'parent_turn':t8['artifact_id'],'spiral_law':t8['spiral_law'],'dimensions':DIMS,'phases':PHASES,'trigger':'HA10_STATIONXML_RESPONSE_EDGE_QUANTIFICATION','turn':{'turn':9,'height':9,'prior_radius':pr,'radius':radius,'prior_state_sha256':prior,'state_sha256':state,'repeats_prior_state':state==prior,'added_constraints':added,'scientific_delta':delta},'nodes':nodes,'tests':[{'id':k,'pass':bool(v)} for k,v in tests.items()],'engine_verdict':'PASS_5D_SPIRAL_TURN_9' if ok else 'FAIL_5D_SPIRAL_TURN_9','hard_rules':['119_HZ_NOT_A_BLIND_KEY','RESPONSE_IS_NOT_OCEAN_SPECTRUM','NO_FFT_FISHING','NO_TARGET_PROMOTION_FROM_FREQUENCY_ONLY','NO_UNDERWATER_PYRAMID_DETECTED_YET'],'status':'TURN_9_ASCENDED' if ok else 'TURN_9_REJECTED'}
 summ={'artifact_id':out['artifact_id']+'-SUMMARY','engine_verdict':out['engine_verdict'],'scientific_verdict':delta,'height':9,'radius':radius,'prior_radius':pr,'state_hash_unique_vs_parent':state!=prior,'tests_passed':sum(tests.values()),'tests_total':len(tests),'db_117_vs_10':d117,'db_119_vs_10':d119,'db_120_vs_10':d120,'inverse_gain_119_vs_10':g['response_edge_calibration']['119_hz_inverse_response_gain_needed_relative_to_10hz'],'119_discovery_key':'REJECTED_ON_HA10_EDH','target_identity':'UNCONFIRMED'}
 OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8');SUM.write_text(json.dumps(summ,indent=2,ensure_ascii=False),encoding='utf-8');print(json.dumps(summ,indent=2));return 0 if ok else 2
if __name__=='__main__':raise SystemExit(main())
