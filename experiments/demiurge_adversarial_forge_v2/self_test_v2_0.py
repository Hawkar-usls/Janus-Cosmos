#!/usr/bin/env python3
from __future__ import annotations
import copy, hashlib, json, math, random, subprocess, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
import janus_cosmos_core_v2 as core
import demiurge_forge_v2 as forge
import download_sky_v2 as dl
import janus_cosmos_v2_0 as sky


def check(cond,msg):
    if not cond:raise RuntimeError(msg)


def main():
    print('JANUS COSMOS v2.0.2 OFFLINE SELF-TEST',flush=True)
    # 1. Core geometry and null generators.
    rng=np.random.default_rng(1234);raw=rng.normal(size=(160,120)).astype(np.float32)
    g=forge.default_genome();x=core.normalize(raw,g);v=core.geometry(x)
    check(x.shape==(128,128),'normalize shape')
    check(v.shape==(14,) and np.all(np.isfinite(v)),'geometry vector')
    for model in ('phase_iaaft','block_shuffle'):
        y=core.surrogate(x,np.random.default_rng(5),model,g)
        check(y.shape==x.shape and np.all(np.isfinite(y)),f'{model} surrogate')
    print('[PASS] core geometry + both null generators')

    # 2. Forge must not import/read sky manifest or target data.
    source=(ROOT/'demiurge_forge_v2.py').read_text(encoding='utf-8')
    for forbidden in ('SKY_MANIFEST','external_data/orion','external_data/ngc1425','parent_v1_6_report'):
        check(forbidden not in source,f'forge source crosses blind wall: {forbidden}')
    print('[PASS] forge/target blind-wall source scan')

    # 3. Re-forge must reproduce the portable detector identity. Raw metrics
    # remain scientific evidence and are checked against a tight tolerance
    # envelope instead of being forced into a cross-platform byte identity.
    expected=json.loads((ROOT/'EXPECTED_FORGE_v2_0.json').read_text(encoding='utf-8'))
    receipt,_,_=forge.forge(write=False,quiet=True)
    deviations=forge.metric_contract_deviations(receipt)
    identity_problems=[]
    if receipt['genome_sha256']!=expected['genome_sha256']:identity_problems.append('genome_sha256')
    if receipt['freeze_sha256']!=expected['freeze_sha256']:identity_problems.append('freeze_sha256')
    if deviations:identity_problems.append('metric_contract')
    if identity_problems:
        print(json.dumps({'portable_verification':'FAIL','problems':identity_problems,
              'actual_genome_sha256':receipt['genome_sha256'],'expected_genome_sha256':expected['genome_sha256'],
              'actual_freeze_sha256':receipt['freeze_sha256'],'expected_freeze_sha256':expected['freeze_sha256'],
              'metric_deviations':deviations,'numeric_environment':receipt['numeric_environment']},indent=2),flush=True)
    check(not identity_problems,'portable forge verification mismatch')
    check(receipt['validation_pass'] is True,'forge validation did not pass')
    check(forge.metric_sha(receipt)==receipt['metrics_sha256'],'metric evidence hash mismatch')

    # Cross-platform numerical evidence may move inside the registered envelope
    # without changing detector identity, but its exact evidence hash changes.
    within=copy.deepcopy(receipt);within['validation_metrics']['fitness']+=5e-7
    within['metrics_sha256']=forge.metric_sha(within)
    check(forge.freeze_identity_sha(within)==receipt['freeze_sha256'],
          'metric evidence contaminated detector identity')
    check(not forge.metric_contract_deviations(within),'registered float tolerance rejected harmless drift')
    check(within['metrics_sha256']!=receipt['metrics_sha256'],'metric evidence hash hid numerical drift')

    material=copy.deepcopy(receipt);material['validation_metrics']['fitness']+=2e-4
    material['metrics_sha256']=forge.metric_sha(material)
    check(forge.metric_contract_deviations(material),'material metric drift escaped tolerance gate')
    print('[PASS] portable detector identity + exact metric receipt + tolerance gate')

    # 4. Frozen payload and fixed backbone verification.
    frozen=sky.verify_forge()
    check(frozen['genome_sha256']==expected['genome_sha256'],'runtime frozen detector mismatch')
    print('[PASS] frozen detector integrity + backbone lock')

    # 5. Monte Carlo resolution check.
    M=json.loads((ROOT/'SKY_MANIFEST_v2_0.json').read_text(encoding='utf-8'))
    N=M['monte_carlo']['test_nulls_per_model'];min_p=1/(N+1)
    alphas=[M['orion']['alpha_corrected'],M['ngc1425']['alpha_corrected'],M['blind_controls']['alpha_corrected_whole'],M['blind_controls']['alpha_corrected_corridor']]
    check(all(min_p<a for a in alphas),f'powered null count cannot resolve alpha: min_p={min_p}, alphas={alphas}')
    print(f'[PASS] Monte Carlo power floor: min p={min_p:.9f}')

    # 6. Blind control coordinates must reproduce from the committed seed string.
    seed=int.from_bytes(hashlib.sha256(b'JANUS_COSMOS_V2_BLIND_CONTROLS_20260814').digest()[:8],'big');rr=random.Random(seed);regen=[]
    for i in range(4):
        ra=rr.random()*360;dec=math.degrees(math.asin(rr.uniform(-.85,.85)));regen.append((ra,dec))
    for got,(ra,dec) in zip(M['blind_controls']['centers'],regen):
        check(abs(got['ra_deg']-ra)<1e-12 and abs(got['dec_deg']-dec)<1e-12,'blind control coordinate drift')
    print('[PASS] blind control coordinate regeneration')

    # 7. Downloader plan: 4 Orion + 16 controls + 2 HST = 22 products.
    plan=dl.plan();check(len(plan)==22,f'downloader plan count={len(plan)}')
    check(all(all(u.startswith('https://') for u in item['urls']) for item in plan),'non-HTTPS data URL')
    print('[PASS] downloader plan = 22 HTTPS FITS products')

    # 8. Parent evidence boundary: v1.6 is prior motivation, never forge training input.
    p=json.loads((ROOT/'prior_evidence/parent_v1_6_report.json').read_text(encoding='utf-8'))
    check(p.get('status')=='PASS','parent v1.6 evidence is not PASS')
    flags=set(p.get('global_flags',[]))
    check('INDEPENDENT_SURVEY_BELT_REPLICATION' in flags and 'SURVEY_ARTIFACT_DISFAVORED' in flags,'unexpected v1.6 parent evidence')
    print('[PASS] parent v1.6 evidence receipt recognized')

    print('SELF-TEST PASS',flush=True)
    return 0

if __name__=='__main__':raise SystemExit(main())
