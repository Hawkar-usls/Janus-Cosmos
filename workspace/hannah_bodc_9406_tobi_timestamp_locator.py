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
    f=ftplib.FTP(timeout=45); f.connect(HOST,21); f.login('anonymous','janus-probe@example.invalid'); return f


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
    except Exception:return None


def transfer_prefix(path,nbytes):
    f=conn(); out=bytearray()
    try:
        s=f.transfercmd('RETR '+path)
        try:
            while len(out)<nbytes:
                b=s.recv(min(262144,nbytes-len(out)))
                if not b: break
                out.extend(b)
        finally:
            try:s.close()
            except Exception:pass
    finally:
        try:f.close()
        except Exception:pass
    return bytes(out)


def head_probe(path,count=8):
    raw=transfer_prefix(path,count*BLOCK); samples=[]; first=None
    for i in range(min(count,len(raw)//BLOCK)):
        b=raw[i*BLOCK:(i+1)*BLOCK]; h=header(b); h['record_index']=i; h['block_sha256']=hashlib.sha256(b).hexdigest(); samples.append(h)
        if first is None and valid_dt(h): first=h
    return first,samples


def stream_window(path,start_index,count):
    """Sequentially stream from byte zero, discard preceding records, capture window."""
    f=conn(); captured=[]; skipped=0; buf=bytearray(); want_start=start_index*BLOCK; want_end=(start_index+count)*BLOCK; pos=0
    try:
        s=f.transfercmd('RETR '+path)
        try:
            while pos<want_end:
                b=s.recv(262144)
                if not b: break
                next_pos=pos+len(b)
                if next_pos>want_start and pos<want_end:
                    lo=max(0,want_start-pos); hi=min(len(b),want_end-pos); buf.extend(b[lo:hi])
                pos=next_pos
        finally:
            try:s.close()
            except Exception:pass
    finally:
        try:f.close()
        except Exception:pass
    for j in range(len(buf)//BLOCK):
        idx=start_index+j; block=bytes(buf[j*BLOCK:(j+1)*BLOCK]); h=header(block); h['record_index']=idx; h['block_sha256']=hashlib.sha256(block).hexdigest(); captured.append(h)
    return {'bytes_streamed_from_start':pos,'captured':captured,'captured_bytes':len(buf)}


def file_entries():
    f=conn(); out=[]
    try:
        for sd,facts in f.mlsd(ROOT,facts=['type']):
            if facts.get('type')!='dir' or not sd.startswith('sd'):continue
            d=posixpath.join(ROOT,sd)
            for n,ff in f.mlsd(d,facts=['type','size']):
                if ff.get('type')=='file' and n.upper() in {'TOBI.DAT','TOBIA.DAT'}:
                    sz=int(ff['size']); out.append({'path':posixpath.join(d,n),'relative_path':sd+'/'+n,'size_bytes':sz,'record_count':sz//BLOCK,'size_mod_block':sz%BLOCK})
    finally:
        try:f.quit()
        except Exception:pass
    return sorted(out,key=lambda x:x['relative_path'])


def compact_header(h):
    if not h:return None
    return {k:v for k,v in h.items() if k not in {'magx_raw','magy_raw','magz_raw','roll_raw','pitch_raw','gyro_raw','pressure_raw','temperature_raw','conductivity_raw','lss_raw'}}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True,type=Path); a=ap.parse_args()
    r={'schema':'janus.cosmos.cousteau.hannah_bodc.tobi_timestamp_locator.v2','target_utc':TARGET.isoformat(),'block_size_bytes':BLOCK,'cadence_seconds_from_readme':4,'transport_rule':'SEQUENTIAL_STREAM_NO_FTP_REST','status':'STARTED','scientific_claim':False,'files':[]}
    try:
        entries=file_entries(); starts=[]
        for e in entries:
            first,samples=head_probe(e['path'],8); fd=valid_dt(first or {})
            row={**e,'first_valid_header':compact_header(first),'first_probe_count':len(samples),'estimated_end_utc_from_4s_cadence':(fd.timestamp()+4*(e['record_count']-1) if fd else None)}
            if row['estimated_end_utc_from_4s_cadence'] is not None:
                row['estimated_end_utc_from_4s_cadence']=datetime.fromtimestamp(row['estimated_end_utc_from_4s_cadence'],timezone.utc).isoformat()
            r['files'].append(row)
            if fd: starts.append((fd,e))
        starts.sort(key=lambda x:x[0])
        r['chronology']=[{'relative_path':e['relative_path'],'start_utc':dt.isoformat()} for dt,e in starts]
        candidates=[]
        for i,(dt,e) in enumerate(starts):
            next_dt=starts[i+1][0] if i+1<len(starts) else None
            approx_end=dt.timestamp()+4*e['record_count']
            if dt<=TARGET and ((next_dt is not None and TARGET<next_dt) or (next_dt is None and TARGET.timestamp()<approx_end)):
                candidates.append((dt,e))
        r['target_file_candidates']=[e['relative_path'] for _,e in candidates]
        r['location_runs']=[]
        for start,e in candidates:
            estimated=int(round((TARGET-start).total_seconds()/4.0)); estimated=max(0,min(e['record_count']-1,estimated))
            win_start=max(0,estimated-5); run=stream_window(e['path'],win_start,11)
            valid=[]
            for h in run['captured']:
                d=valid_dt(h)
                if d: valid.append({'delta_seconds':(d-TARGET).total_seconds(),'header':h})
            valid.sort(key=lambda x:abs(x['delta_seconds']))
            # One correction pass if cadence estimate misses due to recording gaps.
            correction=None
            if valid and abs(valid[0]['delta_seconds'])>12:
                adjust=int(round(valid[0]['delta_seconds']/4.0)); new_est=max(0,min(e['record_count']-1,valid[0]['header']['record_index']-adjust)); new_start=max(0,new_est-5)
                run2=stream_window(e['path'],new_start,11); valid2=[]
                for h in run2['captured']:
                    d=valid_dt(h)
                    if d:valid2.append({'delta_seconds':(d-TARGET).total_seconds(),'header':h})
                valid2.sort(key=lambda x:abs(x['delta_seconds']))
                correction={'new_estimated_index':new_est,'bytes_streamed_from_start':run2['bytes_streamed_from_start'],'nearest_candidates':valid2}
                if valid2 and (not valid or abs(valid2[0]['delta_seconds'])<abs(valid[0]['delta_seconds'])): valid=valid2
            r['location_runs'].append({'relative_path':e['relative_path'],'file_start_utc':start.isoformat(),'initial_estimated_index':estimated,'bytes_streamed_from_start':run['bytes_streamed_from_start'],'nearest_candidates':valid,'correction_pass':correction})
        if r['location_runs'] and any(x['nearest_candidates'] for x in r['location_runs']):
            r['status']='REAL_TOBI_TARGET_TIMESTAMP_LOCATED'
        else:r['status']='TARGET_NOT_LOCATED_AFTER_SEQUENTIAL_STREAM'
    except Exception as e:
        r['status']='TOBI_TIMESTAMP_LOCATOR_FAILED'; r['error_type']=type(e).__name__; r['error']=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':r['status'],'target_file_candidates':r.get('target_file_candidates'),'location_run_count':len(r.get('location_runs',[])),'error':r.get('error')},indent=2))
    return 0 if r['status'] in {'REAL_TOBI_TARGET_TIMESTAMP_LOCATED','TARGET_NOT_LOCATED_AFTER_SEQUENTIAL_STREAM'} else 2
if __name__=='__main__':raise SystemExit(main())
