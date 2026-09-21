#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, random
from pathlib import Path

def load_pack(path: Path):
    d=json.loads(path.read_text(encoding='utf-8'))
    assert d['firewall']['validated_on_unseen_real_data'] is False
    assert d['firewall']['cd169_can_validate_this_pack'] is False
    return d

def sample(pack, regime, method, n, seed):
    rng=random.Random(seed); out=[]
    if regime=='IDEAL':
        return [{'east_m':0.0,'north_m':0.0,'radius_m':0.0,'source':'IDEAL'} for _ in range(n)]
    if regime=='ADVERSARIAL_REAL_TAILS':
        radii=pack['models']['ADVERSARIAL_REAL_TAILS']['tail_radius_m']
        for _ in range(n):
            r=rng.choice(radii); a=rng.uniform(0,2*math.pi)
            out.append({'east_m':r*math.cos(a),'north_m':r*math.sin(a),'radius_m':r,'source':'ADVERSARIAL_REAL_TAILS'})
        return out
    if regime!='REALISTIC_CALIBRATED': raise ValueError('unknown regime')
    if method=='SPARSE_LOG_GENERIC':
        radii=pack['models']['GENERIC_SPARSE_LOG_TO_NATIVE_RADIAL']['empirical_radius_m']
        for _ in range(n):
            r=rng.choice(radii); a=rng.uniform(0,2*math.pi)
            out.append({'east_m':r*math.cos(a),'north_m':r*math.sin(a),'radius_m':r,'source':'GENERIC_SPARSE_LOG_TO_NATIVE_RADIAL'})
        return out
    if method=='TWO_ENDPOINT_LINEAR_INTERPOLATION':
        m=pack['models']['TWO_ENDPOINT_LINEAR_INTERPOLATION']; aa=m['empirical_along_m']; cc=m['empirical_cross_m_signed_cd169']
        if len(aa)!=len(cc): raise ValueError('paired empirical arrays differ')
        for _ in range(n):
            i=rng.randrange(len(aa)); along=aa[i]; cross=abs(cc[i])*(-1 if rng.random()<0.5 else 1)
            out.append({'along_m':along,'cross_m':cross,'radius_m':math.hypot(along,cross),'source':'TWO_ENDPOINT_LINEAR_INTERPOLATION_BOOTSTRAP_CROSS_SIGN_RANDOMIZED'})
        return out
    raise ValueError('unknown method')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--pack',required=True,type=Path)
    ap.add_argument('--regime',required=True,choices=['IDEAL','REALISTIC_CALIBRATED','ADVERSARIAL_REAL_TAILS'])
    ap.add_argument('--method',default='SPARSE_LOG_GENERIC',choices=['SPARSE_LOG_GENERIC','TWO_ENDPOINT_LINEAR_INTERPOLATION'])
    ap.add_argument('--n',type=int,default=1000)
    ap.add_argument('--seed',type=int,default=20260921)
    ap.add_argument('--output',type=Path)
    a=ap.parse_args()
    data=sample(load_pack(a.pack),a.regime,a.method,a.n,a.seed)
    payload={'regime':a.regime,'method':a.method,'n':a.n,'seed':a.seed,'samples':data,'claim_ceiling':'CALIBRATED_SYNTHETIC_NAVIGATION_NUISANCE_ONLY__NOT_VALIDATED_PREDICTIVE_SKILL'}
    s=json.dumps(payload,indent=2)
    if a.output:a.output.write_text(s+'\n',encoding='utf-8')
    else:print(s)

if __name__=='__main__':main()
