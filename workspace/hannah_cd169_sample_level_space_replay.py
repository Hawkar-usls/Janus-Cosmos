#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, gzip, hashlib, json, math, statistics, struct
from datetime import datetime, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk'
FILE='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11285/TOBI.DAT'
BLOCK=40960
START,END=9767,9857
PRIMARY=set(range(9804,9821)); TARGET=9812
TARGET_HASH='039b839861e515efe339542b8c57e5830bd484def961804f9c1fb61387132c3f'
FROZEN=(-3.8654180644718967,-12.142441475)
PIXEL_M=0.75
ARRAYS={'port_sidescan':0x0240,'stbd_sidescan':0x2180,'profiler':0x40C0,'port_swath':0x6000,'stbd_swath':0x7F40}


def ftp():
    f=ftplib.FTP(timeout=90); f.connect(HOST,21); f.login('anonymous','janus-probe@example.invalid'); f.voidcmd('TYPE I'); return f

def dos_dt(d,t):
    try:return datetime(1980+((d>>9)&127),(d>>5)&15,d&31,(t>>11)&31,(t>>5)&63,(t&31)*2,tzinfo=timezone.utc)
    except Exception:return None

def circ_mean_deg(vals):
    if not vals:return None
    s=sum(math.sin(math.radians(v)) for v in vals); c=sum(math.cos(math.radians(v)) for v in vals)
    if abs(s)<1e-15 and abs(c)<1e-15:return None
    return (math.degrees(math.atan2(s,c))+360)%360

def telemetry(b):
    version,tim,dat,alt=struct.unpack_from('<HHHH',b,0x30)
    gyro=list(struct.unpack_from('<8h',b,0xD4)); pressure=list(struct.unpack_from('<8H',b,0xE4))
    water_path,wire=struct.unpack_from('<hh',b,0x114)
    dt=dos_dt(dat,tim)
    return {'version':version,'datetime_utc':dt.isoformat().replace('+00:00','Z') if dt else None,'altitude_m':float(alt),'gyro_heading_deg':circ_mean_deg([(x/10.0)-10.1 for x in gyro]),'pressure_dbar_median':statistics.median([(x/10.0)-5.0 for x in pressure]),'water_path_ms':float(water_path),'wire_out_m':float(wire)}

def vals(b,off):return list(struct.unpack_from('<4000h',b,off))
def medmad(xs):
    m=statistics.median(xs); mad=statistics.median(abs(x-m) for x in xs); return float(m),float(mad)
def robust_zs(xs):
    m,mad=medmad(xs); scale=1.4826*mad
    return ([None]*len(xs),m,mad) if scale==0 else ([(x-m)/scale for x in xs],m,mad)
def local_z(xs,i,half=20):
    lo=max(0,i-half); hi=min(len(xs),i+half+1); n=xs[lo:i]+xs[i+1:hi]
    if not n:return None
    m,mad=medmad(n); sc=1.4826*mad
    return None if sc==0 else (xs[i]-m)/sc

def ecdf(values):
    finite=sorted(v for v in values if v is not None and math.isfinite(v))
    if not finite:return [None]*len(values)
    import bisect
    n=len(finite); out=[]
    for v in values:
        out.append(None if v is None or not math.isfinite(v) else bisect.bisect_right(finite,v)/n)
    return out

def destination(lat,lon,bearing_deg,dist_m):
    r=6371008.8; d=dist_m/r; br=math.radians(bearing_deg); p1=math.radians(lat); l1=math.radians(lon)
    p2=math.asin(math.sin(p1)*math.cos(d)+math.cos(p1)*math.sin(d)*math.cos(br))
    l2=l1+math.atan2(math.sin(br)*math.sin(d)*math.cos(p1),math.cos(d)-math.sin(p1)*math.sin(p2))
    return {'latitude':math.degrees(p2),'longitude':((math.degrees(l2)+540)%360)-180}

def approx_projection(name,i,alt,heading):
    if name=='profiler':return {'state':'BLOCKED_PROFILER_DIFFERENT_TIMING_GEOMETRY'}
    sl=i*PIXEL_M
    if sl<alt:return {'state':'WATER_COLUMN_BY_CALIBRATION_HYPOTHESIS','slant_range_m':sl,'altitude_m':alt}
    gr=math.sqrt(max(0.0,sl*sl-alt*alt))
    side='port' if name.startswith('port_') else 'starboard'; bearing=(heading-90 if side=='port' else heading+90)%360
    return {'state':'SCIENCE_LOG_APPROX_PROJECTION','side':side,'sample_index':i,'slant_range_m':sl,'ground_range_m':gr,'bearing_deg':bearing,'coordinate':destination(FROZEN[0],FROZEN[1],bearing,gr),'ground_fixed_proof':False}

def stream_blocks():
    f=ftp(); got={}; idx=0; buf=bytearray(); bytes_seen=0
    try:
        s=f.transfercmd('RETR '+FILE)
        try:
            while idx<=END:
                chunk=s.recv(1048576)
                if not chunk:break
                bytes_seen+=len(chunk); buf.extend(chunk)
                while len(buf)>=BLOCK and idx<=END:
                    b=bytes(buf[:BLOCK]); del buf[:BLOCK]
                    if START<=idx<=END: got[idx]=b
                    idx+=1
        finally:
            try:s.close()
            except Exception:pass
    finally:
        try:f.close()
        except Exception:pass
    return got,bytes_seen,idx

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--summary',required=True,type=Path); ap.add_argument('--metrics',required=True,type=Path); a=ap.parse_args()
    out={'schema':'janus.cosmos.cousteau.cd169_sample_level_space_replay.v1','status':'STARTED','source_file':'sd11285/TOBI.DAT','frozen_target_block':TARGET,'sample_values_redistributed':False,'scientific_anomaly_threshold':None,'geometry_threshold':None,'ground_fixed_proof':False}
    try:
        blocks,seen,nextidx=stream_blocks(); out['stream']={'captured_blocks':len(blocks),'bytes_streamed_from_start':seen,'next_block_index':nextidx,'expected_start':START,'expected_end':END}
        if len(blocks)!=(END-START+1): raise RuntimeError(f'expected {END-START+1} blocks, got {len(blocks)}')
        th=hashlib.sha256(blocks[TARGET]).hexdigest(); out['target_block_sha256']=th
        if th!=TARGET_HASH: raise RuntimeError('target block hash mismatch')
        tele={i:telemetry(blocks[i]) for i in sorted(blocks)}; out['target_telemetry']=tele[TARGET]
        out['timestamp_sequence']={'first':tele[START]['datetime_utc'],'target':tele[TARGET]['datetime_utc'],'last':tele[END]['datetime_utc']}
        arr_by_block={}; z_by_block={}; block_stats={}
        for i,b in blocks.items():
            arr_by_block[i]={n:vals(b,o) for n,o in ARRAYS.items()}; z_by_block[i]={}; block_stats[i]={}
            for n,xs in arr_by_block[i].items():
                z,m,mad=robust_zs(xs); z_by_block[i][n]=z; block_stats[i][n]={'median_raw':m,'mad_raw':mad,'finite_robust_z':sum(v is not None for v in z)}
        target_metrics={}; rankings={}; nadir={}
        alt=tele[TARGET]['altitude_m']; heading=tele[TARGET]['gyro_heading_deg']; continuous_nadir=alt/PIXEL_M; discrete_nadir=math.ceil(continuous_nadir)
        out['nadir_prediction']={'altitude_m':alt,'pixel_m_calibration_hypothesis':PIXEL_M,'continuous_sample_index':continuous_nadir,'nearest_discrete_valid_bottom_sample':discrete_nadir,'selection_used_intensity':False,'classification':'EXTERNAL_PRISM_CALIBRATION_PENDING_CD169_CDF_CROSSCHECK'}
        for n in ARRAYS:
            xs=arr_by_block[TARGET][n]; tz=z_by_block[TARGET][n]
            lz=[local_z(xs,i) for i in range(4000)]
            pers=[]; signc=[]
            for j in range(4000):
                vs=[z_by_block[i][n][j] for i in sorted(PRIMARY) if z_by_block[i][n][j] is not None]
                if vs:
                    pers.append(float(statistics.median(abs(v) for v in vs)))
                    signc.append(abs(sum(1 if v>0 else -1 if v<0 else 0 for v in vs)/len(vs)))
                else: pers.append(None); signc.append(None)
            p1=ecdf([abs(v) if v is not None else None for v in tz]); p2=ecdf([abs(v) if v is not None else None for v in lz]); p3=ecdf(pers)
            scores=[]; rows=[]
            for j in range(4000):
                score=None if any(v is None for v in (p1[j],p2[j],p3[j])) else (p1[j]+p2[j]+p3[j])/3
                scores.append(score)
                rows.append({'sample_index':j,'target_robust_z':None if tz[j] is None else round(tz[j],6),'local_contrast_z':None if lz[j] is None else round(lz[j],6),'primary_persistence':None if pers[j] is None else round(pers[j],6),'primary_sign_consistency':None if signc[j] is None else round(signc[j],6),'review_priority_score':None if score is None else round(score,9)})
            ranked=sorted(((-s,j) for j,s in enumerate(scores) if s is not None))[:20]
            top=[]
            for neg,j in ranked:
                r=dict(rows[j]); r['approx_projection']=approx_projection(n,j,alt,heading); top.append(r)
            rankings[n]=top; target_metrics[n]=rows
            if n!='profiler' and 0<=discrete_nadir<4000:
                r=dict(rows[discrete_nadir]); r['approx_projection']=approx_projection(n,discrete_nadir,alt,heading); nadir[n]=r
        out['target_block_array_stats']=block_stats[TARGET]
        out['top20_review_priority_per_array']=rankings
        out['predicted_nadir_sample_metrics']=nadir
        out['profiler_mapping_state']='BLOCKED_PENDING_PROFILER_SAMPLE_RATE_AND_SOUND_SPEED_MODEL'
        out['native_navigation_state']='PRISM_SHIP_NAV_AND_CABLE_INPUTS_RECOVERED__VEH_NAV_OUTPUT_NOT_ARCHIVED_OR_NOT_FOUND'
        out['space_state']='SAMPLE_AXIS_REPLAY_READY__EARTH_FIXED_REPLICATION_BLOCKED_PENDING_VEH_NAV_RECONSTRUCTION_AND_CD169_CDF_RANGE_CROSSCHECK'
        out['hard_rules']=['TOP20_IS_REVIEW_PRIORITY_NOT_ANOMALY','NADIR_INDEX_SELECTED_FROM_ALTITUDE_AND_FROZEN_0.75M_CALIBRATION_NOT_INTENSITY','SCIENCE_LOG_APPROX_PROJECTION_IS_NOT_GROUND_FIXED_PROOF','ADJACENT_BLOCK_SAME_INDEX_PERSISTENCE_IS_NOT_GROUND_FIXED_REPLICATION','NO_GEOMETRY_CLASS']
        out['status']='REAL_SAMPLE_LEVEL_REPLAY_READY__SPACE_GATE_PARTIAL'
        metrics={'schema':'janus.cosmos.cousteau.cd169_sample_metrics.v1','target_block':TARGET,'primary_blocks':sorted(PRIMARY),'array_length':4000,'metrics_by_array':target_metrics,'raw_int16_values_included':False,'authority':'REVIEW_PRIORITY_ONLY'}
        a.metrics.parent.mkdir(parents=True,exist_ok=True)
        with gzip.open(a.metrics,'wt',encoding='utf-8') as g: json.dump(metrics,g,ensure_ascii=False,separators=(',',':'))
        out['metrics_gzip_sha256']=hashlib.sha256(a.metrics.read_bytes()).hexdigest(); out['metrics_gzip_size_bytes']=a.metrics.stat().st_size
    except Exception as e:
        out['status']='SAMPLE_LEVEL_REPLAY_FAILED'; out['error_type']=type(e).__name__; out['error']=str(e)
    a.summary.parent.mkdir(parents=True,exist_ok=True); a.summary.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':out['status'],'nadir':out.get('nadir_prediction'),'space_state':out.get('space_state'),'top_indices':{k:[x['sample_index'] for x in v[:5]] for k,v in out.get('top20_review_priority_per_array',{}).items()}},indent=2))
    return 0 if out['status']=='REAL_SAMPLE_LEVEL_REPLAY_READY__SPACE_GATE_PARTIAL' else 2
if __name__=='__main__':raise SystemExit(main())
