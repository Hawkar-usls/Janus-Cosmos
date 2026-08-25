#!/usr/bin/env python3
from __future__ import annotations
import argparse, ftplib, hashlib, json, math, posixpath, struct
from datetime import datetime, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk'; ROOT='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI'; BLOCK=40960
TARGET=datetime(2005,2,28,1,7,25,tzinfo=timezone.utc)

def ftp():
    f=ftplib.FTP(timeout=60); f.connect(HOST,21); f.login('anonymous','janus-probe@example.invalid'); f.voidcmd('TYPE I'); return f

def dos_dt(d,t):
    try:return datetime(1980+((d>>9)&127),(d>>5)&15,d&31,(t>>11)&31,(t>>5)&63,(t&31)*2,tzinfo=timezone.utc)
    except:return None

def dd(deg,mins):
    if not math.isfinite(mins) or abs(mins)>60.5:return None
    return deg-mins/60 if deg<0 else deg+mins/60

def parse(b):
    if len(b)!=BLOCK:return {'valid':False,'reason':'block_length','bytes':len(b)}
    try:
        ver,tm,dt,alt=struct.unpack_from('<HHHH',b,0x30); stamp=dos_dt(dt,tm)
        lond,latd,lonm,latm=struct.unpack_from('<hhff',b,0x38)
        gx=list(struct.unpack_from('<8h',b,0xD4)); roll=list(struct.unpack_from('<8h',b,0xA4)); pitch=list(struct.unpack_from('<8h',b,0xB4))
        press=list(struct.unpack_from('<8H',b,0xE4)); temp=list(struct.unpack_from('<8H',b,0xF4)); cond=list(struct.unpack_from('<8H',b,0x104))
        magx=list(struct.unpack_from('<8i',b,0x44)); magy=list(struct.unpack_from('<8i',b,0x64)); magz=list(struct.unpack_from('<8i',b,0x84)); lss=list(struct.unpack_from('<8i',b,0x118))
        water,wire=struct.unpack_from('<hh',b,0x114)
        return {'valid':bool(stamp and 2004<=stamp.year<=2006),'datetime_utc':stamp.isoformat() if stamp else None,
          'heading_text':b[:48].split(b'\0',1)[0].decode('ascii','replace').strip(),'version':ver,'altitude_raw':alt,
          'ship_position_dd':{'longitude':dd(lond,lonm),'latitude':dd(latd,latm)},'ship_position_raw':{'lon_deg':lond,'lat_deg':latd,'lon_min':lonm,'lat_min':latm},
          'gyro_raw':gx,'roll_raw':roll,'pitch_raw':pitch,'pressure_raw':press,'temperature_raw':temp,'conductivity_raw':cond,
          'magx_raw':magx,'magy_raw':magy,'magz_raw':magz,'water_path_raw':water,'wire_out_m':wire,'lss_raw':lss,
          'packed_time_hex':f'{tm:04x}','packed_date_hex':f'{dt:04x}'}
    except Exception as e:return {'valid':False,'reason':type(e).__name__+': '+str(e)}

def data_files():
    f=ftp(); out=[]
    try:
        for rawd in f.nlst(ROOT):
            sd=posixpath.basename(rawd.rstrip('/'))
            if not sd.startswith('sd'):continue
            d=rawd if rawd.startswith('/') else posixpath.join(ROOT,sd)
            try:names=f.nlst(d)
            except:continue
            for rawp in names:
                n=posixpath.basename(rawp.rstrip('/'))
                if n.upper() not in {'TOBI.DAT','TOBIA.DAT'}:continue
                p=rawp if rawp.startswith('/') else posixpath.join(d,n)
                try:sz=f.size(p)
                except:sz=None
                if sz:out.append({'path':p,'relative_path':sd+'/'+n,'size_bytes':int(sz),'record_count':int(sz)//BLOCK,'size_mod_block':int(sz)%BLOCK})
    finally:
        try:f.quit()
        except:pass
    return sorted(out,key=lambda x:x['relative_path'])

def stream_capture(path,start_index,count):
    want0=start_index*BLOCK; want1=(start_index+count)*BLOCK; pos=0; buf=bytearray(); f=ftp()
    try:
        s=f.transfercmd('RETR '+path)
        try:
            while pos<want1:
                chunk=s.recv(1024*1024)
                if not chunk:break
                nxt=pos+len(chunk)
                if nxt>want0 and pos<want1:
                    buf.extend(chunk[max(0,want0-pos):min(len(chunk),want1-pos)])
                pos=nxt
        finally:
            try:s.close()
            except:pass
    finally:
        try:f.close()
        except:pass
    blocks=[]
    for j in range(len(buf)//BLOCK):
        bb=bytes(buf[j*BLOCK:(j+1)*BLOCK]); h=parse(bb); h['record_index']=start_index+j; h['block_sha256']=hashlib.sha256(bb).hexdigest(); blocks.append(h)
    return {'bytes_streamed':pos,'captured_bytes':len(buf),'blocks':blocks}

def first_valid(path):
    run=stream_capture(path,0,8)
    for h in run['blocks']:
        if h.get('valid'):return h
    return None

def compact(h):
    if not h:return None
    keep=['datetime_utc','heading_text','version','altitude_raw','ship_position_dd','wire_out_m','record_index','block_sha256']
    return {k:h.get(k) for k in keep}

def dt_of(h):
    try:return datetime.fromisoformat(h['datetime_utc']) if h and h.get('valid') else None
    except:return None

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True,type=Path); a=ap.parse_args()
    r={'schema':'janus.cosmos.cousteau.hannah_bodc.tobi_timestamp_locator.v3','target_utc':TARGET.isoformat(),'block_size_bytes':BLOCK,'cadence_seconds_from_readme':4,'transport':'anonymous_FTP_NLST_SIZE_sequential_RETR','status':'STARTED','scientific_claim':False,'files':[]}
    try:
        fs=data_files(); starts=[]
        for e in fs:
            h=first_valid(e['path']); d=dt_of(h); r['files'].append({**e,'first_valid':compact(h)})
            if d:starts.append((d,e))
        starts.sort(key=lambda x:x[0]); r['chronology']=[{'relative_path':e['relative_path'],'start_utc':d.isoformat()} for d,e in starts]
        cand=[]
        for i,(d,e) in enumerate(starts):
            nd=starts[i+1][0] if i+1<len(starts) else None
            estimated_end=d.timestamp()+4*e['record_count']
            if d<=TARGET and ((nd and TARGET<nd) or (nd is None and TARGET.timestamp()<estimated_end)):cand.append((d,e))
        r['target_file_candidates']=[e['relative_path'] for _,e in cand]; r['location_runs']=[]
        for start,e in cand:
            idx=max(0,min(e['record_count']-1,int(round((TARGET-start).total_seconds()/4))))
            def make_run(center):
                st=max(0,center-6); rr=stream_capture(e['path'],st,13); vals=[]
                for h in rr['blocks']:
                    d=dt_of(h)
                    if d:vals.append({'delta_seconds':(d-TARGET).total_seconds(),'header':h})
                vals.sort(key=lambda x:abs(x['delta_seconds'])); return st,rr,vals
            st,rr,vals=make_run(idx); correction=None
            if vals and abs(vals[0]['delta_seconds'])>12:
                corrected=max(0,min(e['record_count']-1,vals[0]['header']['record_index']-int(round(vals[0]['delta_seconds']/4))))
                st2,rr2,vals2=make_run(corrected); correction={'estimated_index':corrected,'window_start_index':st2,'bytes_streamed':rr2['bytes_streamed'],'nearest':vals2}
                if vals2 and abs(vals2[0]['delta_seconds'])<abs(vals[0]['delta_seconds']):vals=vals2
            r['location_runs'].append({'relative_path':e['relative_path'],'file_start_utc':start.isoformat(),'initial_estimated_index':idx,'window_start_index':st,'bytes_streamed':rr['bytes_streamed'],'nearest':vals,'correction':correction})
        r['status']='REAL_TOBI_TARGET_TIMESTAMP_LOCATED' if any(x['nearest'] for x in r['location_runs']) else 'TARGET_NOT_LOCATED'
    except Exception as e:r.update(status='TOBI_TIMESTAMP_LOCATOR_FAILED',error_type=type(e).__name__,error=str(e))
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(r,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':r['status'],'files':len(r['files']),'target_file_candidates':r.get('target_file_candidates'),'error':r.get('error')},indent=2)); return 0 if r['status'] in {'REAL_TOBI_TARGET_TIMESTAMP_LOCATED','TARGET_NOT_LOCATED'} else 2
if __name__=='__main__':raise SystemExit(main())
