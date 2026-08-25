#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, posixpath, struct
from datetime import datetime, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk'
ROOT='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI'
BLOCK=40960
TARGET=datetime(2005,2,28,1,7,25,tzinfo=timezone.utc)

def conn():
    f=ftplib.FTP(timeout=30); f.connect(HOST,21); f.login('anonymous','janus-probe@example.invalid'); return f

def read_slice(path,offset,n):
    f=conn(); data=b''
    try:
        s=f.transfercmd('RETR '+path, rest=offset)
        try:
            while len(data)<n:
                b=s.recv(min(65536,n-len(data)))
                if not b: break
                data+=b
        finally:
            try: s.close()
            except Exception: pass
        try: f.voidresp()
        except Exception: pass
    finally:
        try: f.close()
        except Exception: pass
    return data

def dos_dt(d,t):
    try:
        year=1980+((d>>9)&0x7f); month=(d>>5)&0x0f; day=d&0x1f
        hour=(t>>11)&0x1f; minute=(t>>5)&0x3f; second=(t&0x1f)*2
        return datetime(year,month,day,hour,minute,second,tzinfo=timezone.utc)
    except Exception: return None

def dd(deg,mins):
    if not math.isfinite(mins) or abs(mins)>=60.5: return None
    return deg-mins/60 if deg<0 else deg+mins/60

def header(block):
    if len(block)<0x138: return {'valid':False,'reason':'short'}
    try:
        heading=block[:48].split(b'\x00',1)[0].decode('ascii',errors='replace').strip()
        version,tim,dat,alt=struct.unpack_from('<HHHH',block,0x30)
        lon_deg,lat_deg,lon_min,lat_min=struct.unpack_from('<hhff',block,0x38)
        magx=list(struct.unpack_from('<8i',block,0x44)); magy=list(struct.unpack_from('<8i',block,0x64)); magz=list(struct.unpack_from('<8i',block,0x84))
        roll=list(struct.unpack_from('<8h',block,0xA4)); pitch=list(struct.unpack_from('<8h',block,0xB4)); gyro=list(struct.unpack_from('<8h',block,0xD4))
        press=list(struct.unpack_from('<8H',block,0xE4)); temp=list(struct.unpack_from('<8H',block,0xF4)); cond=list(struct.unpack_from('<8H',block,0x104))
        water_path,wire_out=struct.unpack_from('<hh',block,0x114)
        lss=list(struct.unpack_from('<8i',block,0x118))
        dt=dos_dt(dat,tim)
        return {
          'valid':dt is not None,
          'heading_text':heading,'version':version,'datetime_utc':dt.isoformat() if dt else None,
          'packed_time_hex':f'{tim:04x}','packed_date_hex':f'{dat:04x}','altitude_raw':alt,
          'ship_position_raw':{'lon_degs':lon_deg,'lat_degs':lat_deg,'lon_mins':lon_min,'lat_mins':lat_min},
          'ship_position_dd':{'longitude':dd(lon_deg,lon_min),'latitude':dd(lat_deg,lat_min)},
          'magx_raw':magx,'magy_raw':magy,'magz_raw':magz,'roll_raw':roll,'pitch_raw':pitch,'gyro_raw':gyro,
          'pressure_raw':press,'temperature_raw':temp,'conductivity_raw':cond,'water_path_raw':water_path,'wire_out_m':wire_out,'lss_raw':lss,
        }
    except Exception as e: return {'valid':False,'reason':type(e).__name__+': '+str(e)}

def valid_dt(h):
    if not h.get('valid') or not h.get('datetime_utc'): return None
    try:
        d=datetime.fromisoformat(h['datetime_utc']); return d if 2004<=d.year<=2006 else None
    except Exception: return None

def probe_record(path,idx):
    b=read_slice(path,idx*BLOCK,BLOCK); h=header(b); h['record_index']=idx; h['block_sha256']=hashlib.sha256(b).hexdigest(); h['bytes_read']=len(b); return h

def file_entries():
    f=conn(); out=[]
    try:
        for sd,facts in f.mlsd(ROOT,facts=['type']):
            if facts.get('type')!='dir' or not sd.startswith('sd'): continue
            d=posixpath.join(ROOT,sd)
            for n,ff in f.mlsd(d,facts=['type','size']):
                if ff.get('type')=='file' and n.upper() in {'TOBI.DAT','TOBIA.DAT'}:
                    sz=int(ff['size']); out.append({'path':posixpath.join(d,n),'relative_path':sd+'/'+n,'size_bytes':sz,'record_count':sz//BLOCK,'size_mod_block':sz%BLOCK})
    finally:
        try:f.quit()
        except Exception:pass
    return sorted(out,key=lambda x:x['relative_path'])

def nearest_valid_boundary(path,n,side):
    inds=range(0,min(8,n)) if side=='first' else range(n-1,max(-1,n-9),-1)
    samples=[]
    for i in inds:
        h=probe_record(path,i); samples.append(h)
        if valid_dt(h): return h,samples
    return None,samples

def locate(path,n,target):
    lo,hi=0,n-1; trace=[]; best=None
    while lo<=hi and len(trace)<32:
        mid=(lo+hi)//2; h=probe_record(path,mid); dt=valid_dt(h); trace.append({'index':mid,'datetime_utc':h.get('datetime_utc'),'valid':bool(dt)})
        if dt is None:
            # probe a tiny neighborhood for a valid timestamp
            found=None
            for j in (mid-1,mid+1,mid-2,mid+2):
                if 0<=j<n:
                    hh=probe_record(path,j); d=valid_dt(hh)
                    if d is not None: found=(j,hh,d); break
            if found is None: break
            mid,h,dt=found
        delta=abs((dt-target).total_seconds())
        if best is None or delta<best[0]: best=(delta,mid,h)
        if dt<target: lo=mid+1
        elif dt>target: hi=mid-1
        else: best=(0,mid,h); break
    candidates=[]
    center=best[1] if best else max(0,min(n-1,lo))
    for i in range(max(0,center-3),min(n,center+4)):
        h=probe_record(path,i); d=valid_dt(h)
        if d: candidates.append({'delta_seconds':(d-target).total_seconds(),'header':h})
    candidates.sort(key=lambda x:abs(x['delta_seconds']))
    return {'trace':trace,'nearest_candidates':candidates}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True,type=Path); a=ap.parse_args()
    r={'schema':'janus.cosmos.cousteau.hannah_bodc.tobi_timestamp_locator.v1','target_utc':TARGET.isoformat(),'block_size_bytes':BLOCK,'status':'STARTED','scientific_claim':False,'files':[]}
    try:
        entries=file_entries(); candidates=[]
        for e in entries:
            first,fs=nearest_valid_boundary(e['path'],e['record_count'],'first'); last,ls=nearest_valid_boundary(e['path'],e['record_count'],'last')
            row={**e,'first_valid_header':first,'last_valid_header':last,'first_probe_count':len(fs),'last_probe_count':len(ls)}
            fd=valid_dt(first or {}); ld=valid_dt(last or {})
            if fd and ld and min(fd,ld)<=TARGET<=max(fd,ld):
                row['contains_target_by_boundary']=True; candidates.append(row)
            else: row['contains_target_by_boundary']=False
            r['files'].append(row)
        r['target_file_candidates']=[x['relative_path'] for x in candidates]
        r['location_runs']=[]
        for c in candidates:
            loc=locate(c['path'],c['record_count'],TARGET); r['location_runs'].append({'relative_path':c['relative_path'],**loc})
        r['status']='REAL_TOBI_TARGET_TIMESTAMP_LOCATED' if r['location_runs'] else 'TARGET_NOT_BRACKETED_BY_RAW_FILES'
    except Exception as e:
        r['status']='TOBI_TIMESTAMP_LOCATOR_FAILED'; r['error_type']=type(e).__name__; r['error']=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':r['status'],'target_file_candidates':r.get('target_file_candidates'),'location_run_count':len(r.get('location_runs',[]))},indent=2))
    return 0 if r['status'] in {'REAL_TOBI_TARGET_TIMESTAMP_LOCATED','TARGET_NOT_BRACKETED_BY_RAW_FILES'} else 2
if __name__=='__main__': raise SystemExit(main())
