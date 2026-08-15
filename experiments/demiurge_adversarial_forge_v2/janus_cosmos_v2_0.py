#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
import time
import traceback
import zlib
from pathlib import Path

import numpy as np
from scipy import ndimage

import janus_cosmos_core_v2 as core
import demiurge_forge_v2 as forge

VERSION='2.0.2'
ROOT=Path(__file__).resolve().parent
SKY=json.loads((ROOT/'SKY_MANIFEST_v2_0.json').read_text(encoding='utf-8'))
FORGE_EXPECTED_PATH=ROOT/'EXPECTED_FORGE_v2_0.json'
FROZEN_PATH=ROOT/'forge_output'/'FROZEN_DETECTOR_v2_0.json'
DATA=ROOT/'external_data'
OUT=ROOT/'results_v2_0';OUT.mkdir(exist_ok=True)
CHECK=OUT/'checkpoints';CHECK.mkdir(exist_ok=True)
EVENTS=OUT/'janus-cosmos-v2.0-events.jsonl'
REPORT=OUT/'janus-cosmos-v2.0-report.json'
SUMMARY=OUT/'SUMMARY_v2.0.txt'
TERMINAL=OUT/'terminal.log'
LOGFH=None


def log(msg):
    print(msg,flush=True)
    if LOGFH is not None:
        LOGFH.write(msg+'\n');LOGFH.flush()


def emit(event,**kw):
    row={'ts_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'event':event,**kw}
    s=json.dumps(row,sort_keys=True,ensure_ascii=False);log(s)
    with EVENTS.open('a',encoding='utf-8') as f:f.write(s+'\n')


def verify_forge():
    if not FROZEN_PATH.exists():raise RuntimeError('FROZEN_DETECTOR missing: run demiurge_forge_v2.py first')
    frozen=json.loads(FROZEN_PATH.read_text(encoding='utf-8'))
    if frozen.get('status')!='FROZEN' or not frozen.get('validation_pass'):
        raise RuntimeError('forge result is not a validated frozen detector')
    forge.verify_receipt_integrity(frozen)
    if core.sha256_bytes(core.canonical_json(frozen['genome']).encode('utf-8'))!=frozen['genome_sha256']:
        raise RuntimeError('frozen genome hash mismatch')
    if FORGE_EXPECTED_PATH.exists():
        ex=json.loads(FORGE_EXPECTED_PATH.read_text(encoding='utf-8'))
        if frozen['genome_sha256']!=ex['genome_sha256'] or frozen['freeze_sha256']!=ex['freeze_sha256']:
            raise RuntimeError('frozen detector does not match packaged portable identity expectation')
    # Make sure forge did not silently evolve prohibited backbone values.
    fm=json.loads((ROOT/'FORGE_MANIFEST_v2_0.json').read_text(encoding='utf-8'))
    backbone=fm['genome_space']['fixed_backbone']
    for k,v in backbone.items():
        if frozen['genome'][k]!=v:raise RuntimeError(f'forbidden backbone drift: {k}')
    return frozen


def corridor_from_norm_points(x,points,half_width=10,margin=8):
    pts=np.asarray(points,float);c=pts.mean(0);q=pts-c
    _,_,vh=np.linalg.svd(q,full_matrices=False);axis=vh[0]
    if axis[0]<0:axis=-axis
    perp=np.array([-axis[1],axis[0]]);along=q@axis
    amin=float(along.min()-margin);amax=float(along.max()+margin)
    width=max(16,int(math.ceil(amax-amin))+1);height=max(8,int(math.ceil(2*half_width))+1)
    av=np.linspace(amin,amax,width);pv=np.linspace(-half_width,half_width,height)
    A,P=np.meshgrid(av,pv);X=c[0]+A*axis[0]+P*perp[0];Y=c[1]+A*axis[1]+P*perp[1]
    crop=ndimage.map_coordinates(x,[Y,X],order=1,mode='reflect')
    ch,cw=crop.shape;side=max(ch,cw);py=side-ch;px=side-cw
    crop=np.pad(crop,((py//2,py-py//2),(px//2,px-px//2)),mode='reflect')
    if side!=core.IMAGE_SIZE:crop=ndimage.zoom(crop,(core.IMAGE_SIZE/side,core.IMAGE_SIZE/side),order=1)
    return np.clip(crop,0,1).astype(np.float32),{'center_xy':c.tolist(),'axis_xy':axis.tolist(),'sample_shape':[height,width]}


def checkpoint_key(file_sha,genome_sha,label,variant,model,nulls,cal,seeds):
    d={'file_sha':file_sha,'genome_sha':genome_sha,'label':label,'variant':variant,'model':model,'nulls':nulls,'cal':cal,'seeds':seeds,'version':VERSION}
    return hashlib.sha256(core.canonical_json(d).encode()).hexdigest(),d


def run_model(x,genome,file_sha,label,variant,model,nulls,cal,seeds):
    key,settings=checkpoint_key(file_sha,core.sha256_bytes(core.canonical_json(genome).encode()),label,variant,model,nulls,cal,seeds)
    cp=CHECK/f'{label}__{variant}__{model}.json'
    if cp.exists():
        old=json.loads(cp.read_text(encoding='utf-8'))
        if old.get('settings_hash')==key:
            log(f'[CACHE] {label} {variant} {model}');return old['result']
    def progress(seed,i,n):log(f'      {label} {variant} {model} seed={seed}: {i}/{n}')
    res=core.empirical_test(x,genome,model,nulls,cal,seeds,(label,variant),progress=progress)
    cp.write_text(json.dumps({'settings_hash':key,'settings':settings,'result':res},indent=2),encoding='utf-8')
    return res


def render_png(path,a,marks=()):
    a=np.nan_to_num(np.asarray(a,float),nan=0,posinf=0,neginf=0);lo,hi=np.percentile(a,[1,99.7])
    a=np.clip((a-lo)/max(hi-lo,1e-12),0,1);a=np.arcsinh(8*a)/np.arcsinh(8)
    rgb=np.repeat((a*255).astype(np.uint8)[:,:,None],3,axis=2)
    for x,y in marks:
        x=int(round(x));y=int(round(y))
        for d in range(-4,5):
            if 0<=y<rgb.shape[0] and 0<=x+d<rgb.shape[1]:rgb[y,x+d]=255
            if 0<=y+d<rgb.shape[0] and 0<=x<rgb.shape[1]:rgb[y+d,x]=255
    H,W,_=rgb.shape;raw=b''.join(b'\x00'+rgb[y].tobytes() for y in range(H))
    def ch(t,d):return struct.pack('>I',len(d))+t+d+struct.pack('>I',zlib.crc32(t+d)&0xffffffff)
    path.write_bytes(b'\x89PNG\r\n\x1a\n'+ch(b'IHDR',struct.pack('>IIBBBBB',W,H,8,2,0,0,0))+ch(b'IDAT',zlib.compress(raw,9))+ch(b'IEND',b''))


def band_pass(models,alpha):
    return all(models[m]['p_empirical']<alpha for m in ('phase_iaaft','block_shuffle'))


def evaluate_image(path,label,genome,variants,alpha,nulls,cal,seeds):
    raw,h,meta=core.read_primary_fits(path);x=core.normalize(raw,genome);file_sha=core.sha256_file(path)
    out={'file':str(path.relative_to(ROOT)),'sha256':file_sha,'meta':meta,'variants':{}}
    for variant,xv in variants(x,h,raw.shape):
        m={}
        for model in ('phase_iaaft','block_shuffle'):
            emit('model_start',label=label,variant=variant,model=model)
            m[model]=run_model(xv,genome,file_sha,label,variant,model,nulls,cal,seeds)
            emit('model_complete',label=label,variant=variant,model=model,p=m[model]['p_empirical'])
        out['variants'][variant]={'models':m,'robust':band_pass(m,alpha)}
    return out,x,h,raw.shape


def make_orion_variants(genome,stars):
    def fn(x,h,shape):
        cor,cdiag=core.belt_corridor(x,h,shape,stars,half_width=10,margin=8)
        fn.last_corridor=cdiag
        return [('WHOLE',x),('BELT_CORRIDOR',cor)]
    fn.last_corridor=None
    return fn


def make_control_variants(template):
    def fn(x,h,shape):
        cor,cdiag=corridor_from_norm_points(x,template,10,8);fn.last_corridor=cdiag
        return [('WHOLE',x),('PSEUDO_BELT_CORRIDOR',cor)]
    fn.last_corridor=None
    return fn


def make_whole_only():
    return lambda x,h,shape:[('WHOLE',x)]


def control_file(center,survey):
    return DATA/'controls'/f"{center['id'].lower()}_{survey['family'].lower()}_{survey['band'].lower()}.fits".replace('2mass','tmass')


def survey_family_robust(bands,variant,family):
    vals=[b for b in bands.values() if b['family']==family]
    return bool(vals and all(b['analysis']['variants'][variant]['robust'] for b in vals))


def control_field_false_positive(bands,variant):
    dss=any(b['analysis']['variants'][variant]['robust'] for b in bands.values() if b['family']=='DSS2')
    tm=any(b['analysis']['variants'][variant]['robust'] for b in bands.values() if b['family']=='2MASS')
    return bool(dss and tm)


def main():
    global LOGFH
    ap=argparse.ArgumentParser();ap.add_argument('--smoke',action='store_true');ap.add_argument('--self-test',action='store_true');ap.add_argument('--nulls',type=int);ap.add_argument('--cal-nulls',type=int)
    a=ap.parse_args();EVENTS.write_text('',encoding='utf-8');LOGFH=TERMINAL.open('w',encoding='utf-8',buffering=1)
    frozen=verify_forge();genome=frozen['genome'];mc=SKY['monte_carlo']
    nulls=int(a.nulls or (24 if a.smoke else mc['test_nulls_per_model']));cal=int(a.cal_nulls or (12 if a.smoke else mc['calibration_nulls_per_model']));seeds=list(mc['seeds'])
    min_p=1/(nulls+1)
    if not a.smoke:
        required=max(SKY['orion']['alpha_corrected'],SKY['ngc1425']['alpha_corrected'],SKY['blind_controls']['alpha_corrected_whole'],SKY['blind_controls']['alpha_corrected_corridor'])
        # The smallest threshold is the limiting one, not the largest.
        required=min(SKY['orion']['alpha_corrected'],SKY['ngc1425']['alpha_corrected'],SKY['blind_controls']['alpha_corrected_whole'],SKY['blind_controls']['alpha_corrected_corridor'])
        if min_p>=required:raise RuntimeError(f'underpowered Monte Carlo: min_p={min_p} cannot resolve strict alpha={required}')
    report={'schema':'janus.cosmos.v2.0.report','version':VERSION,'status':'RUNNING','smoke_only':bool(a.smoke),
            'frozen_detector':{'genome_sha256':frozen['genome_sha256'],'freeze_sha256':frozen['freeze_sha256'],'validation_metrics':frozen['validation_metrics']},
            'monte_carlo':{'test_nulls_per_model':nulls,'calibration_nulls':cal,'seeds':seeds,'min_resolvable_p':min_p},
            'blind_controls':{},'orion':{},'ngc1425':{},'global_status':{},'errors':[],'claim_ceiling':SKY['claim_ceiling']}
    if a.self_test:
        log('SELF-TEST PASS: frozen detector verified');return 0
    try:
        # ---------------- Blind unrelated sky controls ----------------
        bc=SKY['blind_controls'];template=bc['pseudo_belt_template_norm'];whole_alpha=bc['alpha_corrected_whole'];cor_alpha=bc['alpha_corrected_corridor']
        fp_whole=0;fp_corr=0
        for center in bc['centers']:
            field={'center':center,'bands':{}}
            for s in bc['surveys']:
                p=control_file(center,s)
                if not p.exists():raise RuntimeError(f'missing control file {p}; run download_sky_v2.py')
                label=f"{center['id']}_{s['family']}_{s['band']}"
                # Evaluate same data once, but apply variant-specific alpha after.
                analysis,x,h,shape=evaluate_image(p,label,genome,make_control_variants(template),min(whole_alpha,cor_alpha),nulls,cal,seeds)
                for v in analysis['variants']:
                    alpha=whole_alpha if v=='WHOLE' else cor_alpha
                    analysis['variants'][v]['robust']=band_pass(analysis['variants'][v]['models'],alpha)
                field['bands'][label]={'family':s['family'],'band':s['band'],'analysis':analysis}
            field['false_positive_whole']=control_field_false_positive(field['bands'],'WHOLE')
            field['false_positive_corridor']=control_field_false_positive(field['bands'],'PSEUDO_BELT_CORRIDOR')
            fp_whole+=int(field['false_positive_whole']);fp_corr+=int(field['false_positive_corridor'])
            report['blind_controls'][center['id']]=field
        gate=bc['specificity_gate'];specificity=bool(fp_whole<=gate['max_false_positive_fields_whole'] and fp_corr<=gate['max_false_positive_fields_corridor'])
        report['global_status']['specificity_gate']={'status':'PASS' if specificity else 'FAIL','false_positive_fields_whole':fp_whole,
             'false_positive_fields_corridor':fp_corr,'field_count':len(bc['centers']),'max_whole':gate['max_false_positive_fields_whole'],'max_corridor':gate['max_false_positive_fields_corridor']}

        # ---------------- Orion target, frozen detector ----------------
        o=SKY['orion'];bands={};orion_norm={}
        for s in o['surveys']:
            p=DATA/'orion'/s['filename'];label=f"ORION_{s['family']}_{s['band']}"
            if not p.exists():raise RuntimeError(f'missing Orion file {p}')
            variant_fn=make_orion_variants(genome,o['belt_stars_j2000'])
            analysis,x,h,shape=evaluate_image(p,label,genome,variant_fn,o['alpha_corrected'],nulls,cal,seeds)
            bands[label]={'family':s['family'],'band':s['band'],'analysis':analysis,'corridor':variant_fn.last_corridor}
            orion_norm[label]=x
            marks=[]
            for star in o['belt_stars_j2000']:
                px,py=core.world_to_pixel(h,star['ra_deg'],star['dec_deg']);nx,ny=core.raw_to_norm(px,py,shape,core.IMAGE_SIZE);marks.append((nx,ny))
            render_png(OUT/f"{label.lower()}_preview.png",x,marks)
        family_status={}
        for fam in ('DSS2','2MASS'):
            family_status[fam]={
              'whole_robust':survey_family_robust(bands,'WHOLE',fam),
              'corridor_robust':survey_family_robust(bands,'BELT_CORRIDOR',fam)}
        whole_cross=all(family_status[f]['whole_robust'] for f in ('DSS2','2MASS'))
        corr_cross=all(family_status[f]['corridor_robust'] for f in ('DSS2','2MASS'))
        # Same HiPS projection makes these pixel grids directly comparable.
        corr_pairs={}
        labels=sorted(orion_norm)
        for i in range(len(labels)):
            for j in range(i+1,len(labels)):
                corr_pairs[f'{labels[i]}__{labels[j]}']=core.smooth_correlation(orion_norm[labels[i]],orion_norm[labels[j]],2.0)
        report['orion']={'bands':bands,'family_status':family_status,'whole_cross_survey':whole_cross,'belt_corridor_cross_survey':corr_cross,
                         'cross_band_smooth_correlations':corr_pairs,'status':('SKY_FIXED_MORPHOLOGY_CANDIDATE' if specificity and corr_cross else ('DETECTOR_SPECIFICITY_BLOCKED' if not specificity else 'NOT_REPLICATED'))}

        # ---------------- NGC1425 HST target ----------------
        ng=SKY['ngc1425'];ngbands={}
        for s in ng['surveys']:
            p=DATA/'ngc1425'/s['filename'];label=f"NGC1425_{s['band']}"
            if not p.exists():raise RuntimeError(f'missing NGC1425 file {p}')
            analysis,x,h,shape=evaluate_image(p,label,genome,make_whole_only(),ng['alpha_corrected'],nulls,cal,seeds)
            ngbands[label]={'family':s['family'],'band':s['band'],'analysis':analysis}
            render_png(OUT/f"{label.lower()}_preview.png",x)
        ngrobust=all(v['analysis']['variants']['WHOLE']['robust'] for v in ngbands.values())
        report['ngc1425']={'bands':ngbands,'cross_filter_robust':ngrobust,
                           'status':('HST_CROSS_FILTER_MORPHOLOGY_CANDIDATE' if specificity and ngrobust else ('DETECTOR_SPECIFICITY_BLOCKED' if not specificity else 'NOT_REPLICATED'))}

        report['global_status']['orion_candidate_admitted']=bool(specificity and corr_cross)
        report['global_status']['ngc1425_candidate_admitted']=bool(specificity and ngrobust)
        report['global_status']['claim_ceiling_enforced']=True
        report['status']='PASS'
    except Exception as e:
        report['status']='FAIL';report['errors'].append({'error':f'{type(e).__name__}: {e}','traceback':traceback.format_exc(limit=12)});emit('fatal_error',error=str(e))
    REPORT.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    lines=[f'JANUS COSMOS v{VERSION} — DEMIURGE-ADVERSARIAL COSMOS CHECK','',f"status: {report['status']}",f"smoke_only: {report['smoke_only']}",
           f"frozen_genome: {frozen['genome_sha256']}"]
    if report.get('global_status'):
        sp=report['global_status'].get('specificity_gate',{});lines+=['',f"SPECIFICITY_GATE: {sp.get('status','N/A')}",f"control_false_positive_whole: {sp.get('false_positive_fields_whole','N/A')}",f"control_false_positive_corridor: {sp.get('false_positive_fields_corridor','N/A')}"]
    if report.get('orion'):lines+=['',f"ORION: {report['orion'].get('status')}",f"orion_whole_cross_survey: {report['orion'].get('whole_cross_survey')}",f"orion_belt_corridor_cross_survey: {report['orion'].get('belt_corridor_cross_survey')}"]
    if report.get('ngc1425'):lines+=['',f"NGC1425: {report['ngc1425'].get('status')}",f"ngc1425_cross_filter_robust: {report['ngc1425'].get('cross_filter_robust')}"]
    if report['errors']:lines+=['','ERRORS:']+[f"- {e['error']}" for e in report['errors']]
    lines+=['','Claim ceiling: '+SKY['claim_ceiling']]
    SUMMARY.write_text('\n'.join(lines)+'\n',encoding='utf-8');log('\n'.join(lines))
    if LOGFH:LOGFH.close()
    return 0 if report['status']=='PASS' else 2

if __name__=='__main__':raise SystemExit(main())
