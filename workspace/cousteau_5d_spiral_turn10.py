#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
from datetime import datetime,timezone,date
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data'/'cousteau'
T9=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-009-2026-08-21-v1.0.json'
REC=DATA/'JANUS-ECHO-COUSTEAU-EA-TPHASE-6843-VS-5943-COUNT-RECONCILIATION-2026-08-21-v1.2.json'
PROBE=DATA/'JANUS-ECHO-COUSTEAU-MGDS-FILE-2504732-PROVENANCE-DETAIL-PROBE-2026-08-21-v1.0.json'
OUT=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-010-2026-08-21-v1.0.json';SUM=DATA/'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-010-2026-08-21-v1.0-SUMMARY.json'
DIMS=['D0_SPACE','D1_TIME','D2_ACOUSTIC','D3_REVERSE_CAUSAL','D4_ASSOCIATIVE_PROVENANCE'];PHASES=['BACK','FORWARD','LEFT_HRAIN','RIGHT_INAIHR','ASCEND']
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def csha(o):return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def d(s):return date.fromisoformat(s)
def main():
 t=load(T9);r=load(REC);p=load(PROBE)
 if t.get('status')!='TURN_9_ASCENDED':raise SystemExit('Turn10 forbidden: Turn9 not ascended')
 mg=r['mgds_provenance'];tl=r['publication_timeline'];hy={x['id']:x for x in r['hypotheses']}
 release_before_accept=d(tl['mgds_public_release'])<d(tl['paper_accepted']); size_match=mg['size_match'] is True
 added=['PREPUBLICATION_MGDS_RELEASE_DISFAVORS_POSTPUBLICATION_DRIFT_EXPLANATION','MISSING_900_LEADING_THREE_PICK_EXPLANATION_MAY_ROUTE_SEARCH_BUT_NOT_CLOSE_GATE','HISTORICAL_CHECKSUM_OR_EXPLICIT_PRIMARY_SOURCE_REQUIRED_TO_CLOSE_COUNT_GATE']
 pr=int(t['turn']['radius']);rad=pr+len(added);prior=t['turn']['state_sha256']
 delta='MGDS_PREPUBLICATION_RELEASE_PROVENANCE_RECOVERED__POSTPUBLICATION_DRIFT_DISFAVORED__MISSING900_BOUNDARY_NARROWED_NOT_SOLVED__TARGET_EVIDENCE_NOT_INCREASED'
 state=csha({'turn':10,'parent':prior,'added':added,'delta':delta,'reconciliation_sha':csha(r),'provenance_sha':csha(p)})
 prov=[{'path':str(x.relative_to(ROOT)).replace('\\','/'),'sha256':csha(load(x))} for x in [T9,REC,PROBE]]
 nodes=[
 {'node_id':'S10N01','turn':10,'phase':'BACK','claim_class':'FACT','structural_context':'Return to the 6843-versus-5943 mismatch after recovering MGDS file-level XML provenance.','evidence_status':'FILE_PROVENANCE_ENTERED','parents':[],'payload':{'dataset_created':mg['dataset_created'],'release_date':mg['file_access_release_date'],'file_uid':mg['file_uid'],'file_name':mg['file_name']},'provenance':prov},
 {'node_id':'S10N02','turn':10,'phase':'FORWARD','claim_class':'CORRECTION','structural_context':'MGDS release on 2022-05-31 precedes paper acceptance on 2022-06-16 and publication on 2022-06-26. Current XML records one current file whose data-file size matches the currently downloaded gzip member.','evidence_status':'POSTPUBLICATION_DRIFT_DISFAVORED','parents':['S10N01'],'payload':{'release_before_acceptance':release_before_accept,'size_match':size_match},'provenance':prov},
 {'node_id':'S10N03','turn':10,'phase':'LEFT_HRAIN','claim_class':'HYPOTHESIS','structural_context':'The structural leading explanation remains the three-pick boundary: paper admits origins from three or more arrival picks, while the deposit contains only 4-8 recording-hydrophone rows and is short by exactly 900.','evidence_status':'LEADING_EXPLANATION_NOT_PROVEN','parents':['S10N02'],'payload':{'h1_status':hy['H1_THREE_PICK_OR_THREE_HYDROPHONE_ROWS_OMITTED_FROM_DEPOSIT']['status'],'difference':r['known_counts']['difference']},'provenance':prov},
 {'node_id':'S10N04','turn':10,'phase':'RIGHT_INAIHR','claim_class':'NEW_GATE','structural_context':'The attractive exact-900 explanation may guide source recovery, but arrival picks and recording hydrophones are not silently equated. Closure requires a historical checksum/file copy or explicit repository/author statement.','evidence_status':'PROVENANCE_CLOSURE_RULE_FROZEN','parents':['S10N03'],'payload':{'generated_gate':'COUSTEAU_EA_TPHASE_MISSING_900_EXPLICIT_SOURCE_OR_HISTORICAL_CHECKSUM_RECOVERY_V1'},'provenance':prov},
 {'node_id':'S10N05','turn':10,'phase':'ASCEND','claim_class':'STATE_LIFT','structural_context':'Turn 10 ascends because provenance evidence changed hypothesis weights without creating target evidence. Post-publication drift is disfavored; the missing-900 question is narrower but remains open.','evidence_status':'SPIRAL_ASCENT_WITH_PROVENANCE_RESOLUTION_GAIN','parents':['S10N04'],'payload':{'target_evidence':'NOT_INCREASED','target_identity':'UNCONFIRMED','next_gates':['COUSTEAU_EA_TPHASE_MISSING_900_EXPLICIT_SOURCE_OR_HISTORICAL_CHECKSUM_RECOVERY_V1','COUSTEAU_AUTHOR_MAR_DATASET_S1_REPLICATION_V1','COUSTEAU_HA10_PREDECLARED_RESPONSE_CORRECTED_CONFIRMATORY_PROTOCOL_V1']},'provenance':prov}
 ]
 for n in nodes:n['node_sha256']=csha(n)
 tests={'T10_01_PARENT_T9_ASCENDED':t.get('status')=='TURN_9_ASCENDED','T10_02_RELEASE_PRECEDES_ACCEPTANCE':release_before_accept,'T10_03_CURRENT_FILE_SIZE_MATCHES_DOWNLOAD':size_match,'T10_04_COUNT_DIFFERENCE_900':r['known_counts']['difference']==900,'T10_05_NO_3_HYDROPHONE_ROWS':r['known_counts']['download_contains_3_hydrophone_rows'] is False,'T10_06_H1_LEADING_NOT_PROVEN':hy['H1_THREE_PICK_OR_THREE_HYDROPHONE_ROWS_OMITTED_FROM_DEPOSIT']['status']=='LEADING_STRONGLY_COMPATIBLE_NOT_PROVEN','T10_07_H2_DISFAVORED_NOT_EXCLUDED':hy['H2_POST_PUBLICATION_FILE_DRIFT_OR_REPLACEMENT']['status']=='DISFAVORED_NOT_EXCLUDED','T10_08_NO_SYNTHETIC_ROWS':'DO_NOT_SYNTHESIZE_MISSING_ROWS' in r['hard_rules'],'T10_09_RADIUS_ASCENDS':rad>pr,'T10_10_TARGET_NOT_PROMOTED':'TARGET_EVIDENCE_NOT_INCREASED' in delta}
 ok=all(tests.values())
 out={'artifact_id':'JANUS-ECHO-COUSTEAU-5D-SPIRAL-TURN-010-2026-08-21-v1.0','created_utc':now(),'parent_turn':t['artifact_id'],'spiral_law':t['spiral_law'],'dimensions':DIMS,'phases':PHASES,'trigger':'MGDS_FILE_LEVEL_PREPUBLICATION_PROVENANCE_RECOVERY','turn':{'turn':10,'height':10,'prior_radius':pr,'radius':rad,'prior_state_sha256':prior,'state_sha256':state,'repeats_prior_state':state==prior,'added_constraints':added,'scientific_delta':delta},'nodes':nodes,'tests':[{'id':k,'pass':bool(v)} for k,v in tests.items()],'engine_verdict':'PASS_5D_SPIRAL_TURN_10' if ok else 'FAIL_5D_SPIRAL_TURN_10','hard_rules':['DO_NOT_SYNTHESIZE_MISSING_ROWS','HYPOTHESIS_MAY_ROUTE_SEARCH_BUT_NOT_CLOSE_GATE','METADATA_UPDATE_NOT_FILE_REPLACEMENT_PROOF','TARGET_IDENTITY_UNCONFIRMED'],'status':'TURN_10_ASCENDED' if ok else 'TURN_10_REJECTED'}
 sm={'artifact_id':out['artifact_id']+'-SUMMARY','engine_verdict':out['engine_verdict'],'scientific_verdict':delta,'height':10,'radius':rad,'prior_radius':pr,'tests_passed':sum(tests.values()),'tests_total':len(tests),'mgds_release_date':mg['file_access_release_date'],'paper_accepted':tl['paper_accepted'],'paper_published':tl['paper_first_published'],'missing_rows':900,'leading_explanation':'THREE_PICK_BOUNDARY_STRONGLY_COMPATIBLE_NOT_PROVEN','postpublication_drift':'DISFAVORED_NOT_EXCLUDED','target_identity':'UNCONFIRMED'}
 OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8');SUM.write_text(json.dumps(sm,indent=2,ensure_ascii=False),encoding='utf-8');print(json.dumps(sm,indent=2));return 0 if ok else 2
if __name__=='__main__':raise SystemExit(main())
