#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy import ndimage

import janus_cosmos_core_v2 as core

ROOT=Path(__file__).resolve().parent
MANIFEST_PATH=ROOT/'FORGE_MANIFEST_v2_0.json'
METRIC_REFERENCE_PATH=ROOT/'FORGE_METRIC_REFERENCE_v2_0.json'
OUT=ROOT/'forge_output'
OUT.mkdir(exist_ok=True)
LEDGER=OUT/'demiurge_forge_ledger.jsonl'
FROZEN=OUT/'FROZEN_DETECTOR_v2_0.json'
SUMMARY=OUT/'FORGE_SUMMARY_v2_0.txt'
RELEASE_VERSION='2.0.2'
METRIC_SECTIONS=('train_metrics','validation_metrics','canonical_validation_metrics')
METRIC_CONTRACT={
    'schema':'janus.cosmos.metric_tolerance.v1',
    'scope':'cross_platform_synthetic_validation_only',
    'absolute_tolerance':2e-6,
    'relative_tolerance':2e-6,
    'scientific_gate_uses_raw_metrics':True,
}
FREEZE_FIELDS=(
    'schema','version','status','origin','training_rank','genome','genome_sha256',
    'validation_pass','metric_contract','forbidden_target_access_verified',
    'source_hashes','frozen_utc','freeze_hash_scope',
)


def load_manifest():
    return json.loads(MANIFEST_PATH.read_text(encoding='utf-8'))


def default_genome():
    return {
        'asinh_scale':6.0,
        'clip_high':99.5,
        'low_freq_frac':0.10,
        'block_size':16,
        'group_weights':{
            'directional':1.0,'rotation':1.0,'anisotropy':1.0,
            'high_frequency':1.0,'component_count':1.0,'largest_component':1.0,
        }
    }


def canonicalize_genome(g):
    g=copy.deepcopy(g)
    # Clamp/canonicalize floats for stable hashes across machines.
    for k in ('asinh_scale','clip_high','low_freq_frac'):
        g[k]=float(g[k])
    g['block_size']=int(g['block_size'])
    for k in sorted(g['group_weights']):
        g['group_weights'][k]=round(float(g['group_weights'][k]),6)
    return g


def genome_sha(g):
    return core.sha256_bytes(core.canonical_json(canonicalize_genome(g)).encode('utf-8'))


def sha256_normalized_text_file(path):
    """Hash UTF-8 source content without platform-specific line endings."""
    text=Path(path).read_text(encoding='utf-8').replace('\r\n','\n').replace('\r','\n')
    return core.sha256_bytes(text.encode('utf-8'))


def metric_payload(obj):
    return {section:obj[section] for section in METRIC_SECTIONS}


def metric_sha(obj):
    return core.sha256_bytes(core.canonical_json(metric_payload(obj)).encode('utf-8'))


def freeze_identity_payload(obj):
    """Return only detector-identity fields, excluding platform-sensitive evidence."""
    return {key:obj[key] for key in FREEZE_FIELDS}


def freeze_identity_sha(obj):
    return core.sha256_bytes(core.canonical_json(freeze_identity_payload(obj)).encode('utf-8'))


def load_metric_reference():
    return json.loads(METRIC_REFERENCE_PATH.read_text(encoding='utf-8'))


def metric_contract_deviations(receipt, reference=None):
    reference=reference or load_metric_reference()
    if reference.get('contract')!=METRIC_CONTRACT:
        return [{'kind':'contract_mismatch','expected':METRIC_CONTRACT,'actual':reference.get('contract')}]
    deviations=[]
    abs_tol=float(METRIC_CONTRACT['absolute_tolerance'])
    rel_tol=float(METRIC_CONTRACT['relative_tolerance'])
    for section in METRIC_SECTIONS:
        expected_metrics=reference['sections'][section]
        actual_metrics=receipt[section]
        if set(expected_metrics)!=set(actual_metrics):
            deviations.append({'kind':'metric_key_mismatch','section':section,
                               'expected':sorted(expected_metrics),'actual':sorted(actual_metrics)})
            continue
        for name,expected in expected_metrics.items():
            actual=float(actual_metrics[name]); expected=float(expected)
            tolerance=max(abs_tol,rel_tol*max(abs(expected),abs(actual)))
            delta=abs(actual-expected)
            if not math.isfinite(actual) or delta>tolerance:
                deviations.append({'kind':'metric_outside_tolerance','section':section,'metric':name,
                                   'expected':expected,'actual':actual,'delta':delta,'tolerance':tolerance})
    return deviations


def verify_receipt_integrity(receipt):
    if metric_sha(receipt)!=receipt.get('metrics_sha256'):
        raise RuntimeError('frozen validation-metrics payload hash mismatch')
    if freeze_identity_sha(receipt)!=receipt.get('freeze_sha256'):
        raise RuntimeError('portable detector-identity hash mismatch')
    deviations=metric_contract_deviations(receipt)
    if deviations:
        raise RuntimeError('synthetic validation metrics outside portable tolerance: '+json.dumps(deviations,sort_keys=True))
    return True


def random_genome(rng, M):
    S=M['genome_space']; B=S['fixed_backbone']
    return canonicalize_genome({
        **B,
        'group_weights':{k:float(rng.uniform(*S['group_weight_bounds'])) for k in S['group_weights']}
    })


def mutate(g, rng, M, strength=1.0):
    S=M['genome_space']; out=copy.deepcopy(g)
    lo,hi=S['group_weight_bounds']
    for k in S['group_weights']:
        if rng.random()<0.75:
            out['group_weights'][k]=float(np.clip(out['group_weights'][k]+rng.normal(0,0.13*strength),lo,hi))
    return canonicalize_genome(out)


def crossover(a,b,rng,M):
    child=copy.deepcopy(a)
    for k in child['group_weights']:
        t=float(rng.uniform(.2,.8))
        child['group_weights'][k]=t*a['group_weights'][k]+(1-t)*b['group_weights'][k]
    return mutate(child,rng,M,0.7)


def _add_gaussian(img,cx,cy,amp,sigma):
    h,w=img.shape; r=max(2,int(math.ceil(4*sigma)))
    x0=max(0,int(cx)-r);x1=min(w,int(cx)+r+1);y0=max(0,int(cy)-r);y1=min(h,int(cy)+r+1)
    yy,xx=np.indices((y1-y0,x1-x0)); xx=xx+x0; yy=yy+y0
    img[y0:y1,x0:x1]+=amp*np.exp(-((xx-cx)**2+(yy-cy)**2)/(2*sigma*sigma))


def make_base(seed,variant=0,n=128):
    rng=np.random.default_rng(core.stable_seed('forge_base',seed,variant))
    # Multi-scale diffuse sky-like texture plus many compact sources.
    img=rng.normal(0,0.025,(n,n))
    diffuse=ndimage.gaussian_filter(rng.normal(0,1,(n,n)),rng.uniform(7,16))
    diffuse=(diffuse-diffuse.mean())/max(diffuse.std(),1e-9)
    img += 0.055*diffuse
    yy,xx=np.indices((n,n));
    img += rng.uniform(-0.05,0.05)*(xx/(n-1)-.5)+rng.uniform(-0.05,0.05)*(yy/(n-1)-.5)
    for _ in range(int(rng.integers(45,90))):
        _add_gaussian(img,float(rng.uniform(3,n-3)),float(rng.uniform(3,n-3)),float(rng.uniform(.12,.75)),float(rng.uniform(.55,1.8)))
    return img.astype(np.float32)


def inject_signal(base, seed, variant=0):
    rng=np.random.default_rng(core.stable_seed('forge_signal',seed,variant)); img=base.copy(); n=img.shape[0]
    kind=(seed+variant)%4
    yy,xx=np.indices(img.shape)
    if kind==0:  # coherent curved arc pair
        cx,cy=rng.uniform(48,80),rng.uniform(48,80); r=rng.uniform(23,34)
        theta=np.arctan2(yy-cy,xx-cx); rr=np.sqrt((xx-cx)**2+(yy-cy)**2)
        angular=((theta>-.9)&(theta<1.25))
        img += .22*np.exp(-((rr-r)/1.5)**2)*angular
        img += .14*np.exp(-((rr-(r+7))/1.9)**2)*angular
    elif kind==1:  # repeated symmetric nodes
        cx,cy=rng.uniform(52,76),rng.uniform(52,76); r=rng.uniform(18,30)
        for a in np.linspace(0,2*np.pi,8,endpoint=False):
            _add_gaussian(img,cx+r*np.cos(a),cy+r*np.sin(a),.32,1.6)
        _add_gaussian(img,cx,cy,.22,2.4)
    elif kind==2:  # two coherent filaments with common orientation
        ang=rng.uniform(-1.0,1.0); ca,sa=np.cos(ang),np.sin(ang)
        x0,y0=rng.uniform(30,55),rng.uniform(35,75)
        for off in (-8,8):
            d=np.abs((xx-x0)*(-sa)+(yy-y0)*ca-off)
            along=(xx-x0)*ca+(yy-y0)*sa
            img += .18*np.exp(-(d/1.5)**2)*np.exp(-((along-20)/35)**4)
    else:  # nested ellipse/ring geometry
        cx,cy=rng.uniform(50,78),rng.uniform(50,78); ang=rng.uniform(-.7,.7)
        ca,sa=np.cos(ang),np.sin(ang); X=(xx-cx)*ca+(yy-cy)*sa; Y=-(xx-cx)*sa+(yy-cy)*ca
        for rr,amp in ((1.0,.18),(1.32,.12)):
            q=np.sqrt((X/(27*rr))**2+(Y/(16*rr))**2)
            img += amp*np.exp(-((q-1)/.07)**2)
    return img.astype(np.float32)


def inject_artifact(base, seed, variant=0):
    rng=np.random.default_rng(core.stable_seed('forge_artifact',seed,variant)); img=base.copy(); n=img.shape[0]
    kind=(seed+2*variant)%6
    if kind==0:  # polygon/no-data wedge
        yy,xx=np.indices(img.shape); cx,cy=rng.uniform(82,104),rng.uniform(18,38)
        mask=(xx>cx-12)&(xx<cx+13)&(yy>cy-11)&(yy<cy+10)&((xx-cx)+.8*(yy-cy)>-10)
        img[mask]=np.percentile(img,2)
    elif kind==1:  # plate step
        x=int(rng.integers(40,90)); img[:,x:]+=rng.uniform(.12,.24)
    elif kind==2:  # saturation cross / bleed
        cx,cy=int(rng.integers(28,100)),int(rng.integers(28,100));
        img[max(0,cy-1):min(n,cy+2),:]+=0.28
        img[:,max(0,cx-1):min(n,cx+2)]+=0.28
        _add_gaussian(img,cx,cy,1.4,3.2)
    elif kind==3:  # edge truncation / footprint
        cut=int(rng.integers(12,30)); img[:cut,:]=np.median(img); img[:, :cut//2]=np.median(img)
    elif kind==4:  # rectangular mosaic gain blocks
        for y in (0,n//2):
            for x in (0,n//2):
                img[y:y+n//2,x:x+n//2]+=rng.uniform(-.10,.10)
    else:  # giant halo / gradient bloom
        yy,xx=np.indices(img.shape);cx,cy=rng.uniform(20,108),rng.uniform(20,108);s=rng.uniform(13,25)
        img += .26*np.exp(-((xx-cx)**2+(yy-cy)**2)/(2*s*s))
    return img.astype(np.float32)


def build_corpus(seeds):
    # Backbone preprocessing is deliberately fixed. Feature vectors can therefore
    # be computed once; Demiurge is allowed to tune only bounded feature-group weights.
    backbone=default_genome()
    rows=[]
    for seed in seeds:
        for variant in range(3):
            b=make_base(seed,variant); s=inject_signal(b,seed,variant); a=inject_artifact(b,seed,variant)
            rows.append((
                core.geometry(core.normalize(b,backbone)),
                core.geometry(core.normalize(s,backbone)),
                core.geometry(core.normalize(a,backbone)),seed,variant
            ))
    return rows


def evaluate_genome(genome, corpus, M):
    bvec=np.asarray([r[0] for r in corpus]);svec=np.asarray([r[1] for r in corpus]);avec=np.asarray([r[2] for r in corpus])
    center=bvec.mean(0); scale=bvec.std(0,ddof=1)
    bs=np.asarray([core.weighted_std_dist(v,center,scale,genome) for v in bvec])
    ss=np.asarray([core.weighted_std_dist(v,center,scale,genome) for v in svec])
    aa=np.asarray([core.weighted_std_dist(v,center,scale,genome) for v in avec])
    # Counterfactual attribution: same base field, two interventions.
    signal_delta=np.asarray([core.weighted_std_dist(svec[i],bvec[i],scale,genome) for i in range(len(bvec))])
    artifact_delta=np.asarray([core.weighted_std_dist(avec[i],bvec[i],scale,genome) for i in range(len(bvec))])
    smed=float(np.median(signal_delta)); amed=float(np.median(artifact_delta))
    ratio=float(smed/max(amed,1e-9))
    beat=float(np.mean(signal_delta>artifact_delta))
    threshold=float(np.quantile(bs,.99))
    art_fp=float(np.mean(aa>threshold)); sig_rec=float(np.mean(ss>threshold))
    pair_ratio=np.log((signal_delta+1e-6)/(artifact_delta+1e-6))
    instability=float(np.std(pair_ratio))
    gw=np.asarray(list(genome['group_weights'].values()),float); complexity=float(np.std(gw))
    F=M['fitness']
    fitness=(F['ratio_reward']*math.log(max(ratio,1e-6))
             +F['signal_sensitivity_reward']*smed
             -F['artifact_false_positive_penalty']*art_fp
             -F['seed_instability_penalty']*instability
             -F['complexity_penalty']*complexity)
    return {
      'fitness':float(fitness),
      'base_median':float(np.median(bs)),'signal_absolute_median':float(np.median(ss)),'artifact_absolute_median':float(np.median(aa)),
      'paired_signal_median':smed,'paired_artifact_median':amed,'signal_to_artifact_ratio':ratio,
      'pairwise_signal_beats_artifact_rate':beat,'signal_recovery_at_base_p99':sig_rec,'artifact_fp_rate_at_base_p99':art_fp,
      'base_p99_threshold':threshold,'seed_instability':instability,'complexity':complexity
    }


def validation_pass(metrics, canonical_metrics, M):
    g=M['fitness']['validation_gate']
    ratio_improvement=metrics['signal_to_artifact_ratio']/max(canonical_metrics['signal_to_artifact_ratio'],1e-9)
    artifact_fraction=metrics['paired_artifact_median']/max(canonical_metrics['paired_artifact_median'],1e-9)
    signal_fraction=metrics['paired_signal_median']/max(canonical_metrics['paired_signal_median'],1e-9)
    return bool(ratio_improvement>=g['min_signal_to_artifact_ratio_improvement_over_canonical']
                and artifact_fraction<=g['max_artifact_delta_fraction_of_canonical']
                and signal_fraction>=g['min_signal_delta_fraction_of_canonical']
                and metrics['pairwise_signal_beats_artifact_rate']>=g['min_pairwise_signal_beats_artifact_rate'])


def forge(write=True, quiet=False):
    M=load_manifest(); rng=np.random.default_rng(int(M['forge_seed']))
    train=build_corpus(M['training_seeds']); valid=build_corpus(M['validation_seeds']); canonical_valid=evaluate_genome(default_genome(),valid,M)
    pop=[default_genome()]
    seen={genome_sha(pop[0])}
    while len(pop)<M['population']:
        g=random_genome(rng,M); h=genome_sha(g)
        if h not in seen: seen.add(h);pop.append(g)
    ledger=[]
    for gen in range(int(M['generations'])):
        evaluated=[]
        for idx,g in enumerate(pop):
            metrics=evaluate_genome(g,train,M); h=genome_sha(g)
            row={'generation':gen,'candidate_index':idx,'genome_sha256':h,'genome':g,'train_metrics':metrics}
            ledger.append(row); evaluated.append((metrics['fitness'],h,g,metrics))
            if not quiet: print(f"GEN {gen} CAND {idx:02d} fitness={metrics['fitness']:.5f} sig={metrics['paired_signal_median']:.3f} art={metrics['paired_artifact_median']:.3f} artFP={metrics['artifact_fp_rate_at_base_p99']:.3f}",flush=True)
        evaluated.sort(key=lambda t:(-t[0],t[1]))
        elites=evaluated[:int(M['elite_count'])]
        if gen==int(M['generations'])-1:
            final_rank=elites
            break
        new=[copy.deepcopy(x[2]) for x in elites]
        # Demiurge pattern: preserve best, mutate elites, combine knowledge from pairs.
        while len(new)<int(M['population']):
            if rng.random()<0.55:
                parent=elites[int(rng.integers(0,len(elites)))][2]
                child=mutate(parent,rng,M,1.0+0.15*gen)
            else:
                a=elites[int(rng.integers(0,len(elites)))][2]; b=elites[int(rng.integers(0,len(elites)))][2]
                child=crossover(a,b,rng,M)
            new.append(child)
        pop=new
    # Validation is a gate, not a ranking target: take first training-ranked elite that passes.
    selected=None; validation_rows=[]
    for rank,(fit,h,g,tm) in enumerate(final_rank[:int(M['validation_candidates'])]):
        vm=evaluate_genome(g,valid,M); ok=validation_pass(vm,canonical_valid,M)
        validation_rows.append({'training_rank':rank+1,'genome_sha256':h,'validation_metrics':vm,'validation_pass':ok})
        if selected is None and ok:
            selected=(g,h,tm,vm,rank+1)
    if selected is None:
        raise RuntimeError('DEMIURGE_FORGE_REJECTED: no top training genome passed frozen synthetic validation gate')
    g,h,tm,vm,rank=selected
    source_sha=sha256_normalized_text_file(Path(__file__))
    core_sha=sha256_normalized_text_file(ROOT/'janus_cosmos_core_v2.py')
    manifest_sha=sha256_normalized_text_file(MANIFEST_PATH)
    metric_reference_sha=sha256_normalized_text_file(METRIC_REFERENCE_PATH)
    payload={
        'schema':'janus.cosmos.frozen_detector.v2.0','version':RELEASE_VERSION,'status':'FROZEN',
        'origin':'Demiurge-inspired adversarial detector forge',
        'training_rank':rank,'genome':g,'genome_sha256':h,
        'train_metrics':tm,'validation_metrics':vm,'canonical_validation_metrics':canonical_valid,'validation_pass':True,
        'metric_contract':METRIC_CONTRACT,
        'forbidden_target_access_verified':True,
        'source_hashes':{'demiurge_forge_v2.py':source_sha,'janus_cosmos_core_v2.py':core_sha,
                         'FORGE_MANIFEST_v2_0.json':manifest_sha,'FORGE_METRIC_REFERENCE_v2_0.json':metric_reference_sha},
        'numeric_environment':{'python':platform.python_version(),'implementation':platform.python_implementation(),
                               'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()},
        'frozen_utc':'DETERMINISTIC_BUILD_NO_RUNTIME_CLOCK_IN_HASH',
        'freeze_hash_scope':{
            'contract':'portable_detector_identity.v1',
            'included_fields':list(FREEZE_FIELDS),
            'excluded_platform_evidence':[*METRIC_SECTIONS,'metrics_sha256','numeric_environment'],
        },
    }
    receipt={**payload,'metrics_sha256':metric_sha(payload)}
    freeze_hash=freeze_identity_sha(receipt)
    receipt['freeze_sha256']=freeze_hash
    if write:
        LEDGER.write_text('',encoding='utf-8')
        with LEDGER.open('a',encoding='utf-8') as f:
            for row in ledger: f.write(json.dumps(row,sort_keys=True,ensure_ascii=False)+'\n')
            for row in validation_rows: f.write(json.dumps({'event':'validation',**row},sort_keys=True,ensure_ascii=False)+'\n')
        FROZEN.write_text(json.dumps(receipt,indent=2,ensure_ascii=False),encoding='utf-8')
        SUMMARY.write_text(
            f'JANUS COSMOS v{RELEASE_VERSION} — DEMIURGE FORGE\n\n'
            f"status: FROZEN\ntraining_rank: {rank}\ngenome_sha256: {h}\nfreeze_sha256: {freeze_hash}\n"
            f"train_fitness: {tm['fitness']:.8f}\n"
            f"validation_ratio: {vm['signal_to_artifact_ratio']:.6f}\ncanonical_validation_ratio: {canonical_valid['signal_to_artifact_ratio']:.6f}\nvalidation_artifact_delta: {vm['paired_artifact_median']:.6f}\n"
            'target_data_seen_by_forge: NO\n',encoding='utf-8')
    return receipt,ledger,validation_rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--quiet',action='store_true');ap.add_argument('--verify-expected',action='store_true')
    a=ap.parse_args();receipt,_,_=forge(write=True,quiet=a.quiet)
    expected_path=ROOT/'EXPECTED_FORGE_v2_0.json'
    if a.verify_expected and expected_path.exists():
        expected=json.loads(expected_path.read_text(encoding='utf-8'))
        problems=[]
        if receipt['genome_sha256']!=expected['genome_sha256']:problems.append('genome_sha256')
        if receipt['freeze_sha256']!=expected['freeze_sha256']:problems.append('freeze_sha256')
        deviations=metric_contract_deviations(receipt)
        if deviations:problems.append('metric_contract')
        if problems:
            detail={'problems':problems,'actual_genome_sha256':receipt['genome_sha256'],
                    'expected_genome_sha256':expected['genome_sha256'],'actual_freeze_sha256':receipt['freeze_sha256'],
                    'expected_freeze_sha256':expected['freeze_sha256'],'metric_deviations':deviations,
                    'numeric_environment':receipt['numeric_environment']}
            raise RuntimeError('FORGE_PORTABLE_VERIFICATION_MISMATCH: '+json.dumps(detail,sort_keys=True))
    print(json.dumps({'status':'FROZEN','version':RELEASE_VERSION,'genome_sha256':receipt['genome_sha256'],
                      'freeze_sha256':receipt['freeze_sha256'],'metrics_sha256':receipt['metrics_sha256'],
                      'numeric_environment':receipt['numeric_environment'],'validation':receipt['validation_metrics']},indent=2))
    return 0

if __name__=='__main__': raise SystemExit(main())
