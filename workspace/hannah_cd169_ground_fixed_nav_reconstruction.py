#!/usr/bin/env python3
from __future__ import annotations
import argparse, bisect, ftplib, hashlib, json, math, statistics, struct
from datetime import datetime, timedelta, timezone
from pathlib import Path

HOST='livftp.noc.ac.uk'
RAW='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11285/TOBI.DAT'
NAV='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI/cd169.nav'
CABLE='/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/TOBI/cd169.cable'
BLOCK=40960
TARGET=datetime(2005,2,28,1,7,25,tzinfo=timezone.utc)
FROZEN=(-3.8654180644718967,-12.142441475)
WINDOW=timedelta(minutes=5)
UMBILICAL_M=200.0
SAMPLE_M=0.75
OFF={'port_sidescan':0x0240,'stbd_sidescan':0x2180,'port_swath':0x6000,'stbd_swath':0x7F40}
ANCHORS=[
 ('00:00:00',-3.907,-12.153,6481),('00:15:00',-3.895,-12.167,6277),('00:29:00',-3.886,-12.151,6262),
 ('01:00:00',-3.869,-12.145,6032),('01:29:00',-3.855,-12.135,6155),('02:00:00',-3.835,-12.128,6050),
 ('02:30:00',-3.814,-12.128,6016),('03:00:00',-3.795,-12.131,6112),('03:30:00',-3.774,-12.134,6314),
 ('04:00:00',-3.752,-12.136,6308),('04:30:00',-3.734,-12.140,6634),('05:00:00',-3.715,-12.146,6586),
 ('05:30:00',-3.683,-12.181,6304),('06:00:00',-3.670,-12.151,6110),('06:31:00',-3.649,-12.154,6200)
]

def ftp():
    f=ftplib.FTP(timeout=90); f.connect(HOST,21); f.login('anonymous','janus-probe@example.invalid'); f.voidcmd('TYPE I'); return f

def get(path):
    f=ftp(); out=[]
    try:f.retrbinary('RETR '+path,out.append,1048576)
    finally:
        try:f.quit()
        except Exception:
            try:f.close()
            except Exception:pass
    return b''.join(out)

def dos_dt(d,t):
    try:return datetime(1980+((d>>9)&127),(d>>5)&15,d&31,(t>>11)&31,(t>>5)&63,(t&31)*2,tzinfo=timezone.utc)
    except Exception:return None

def circ_mean(vals):
    if not vals:return None
    s=sum(math.sin(math.radians(v)) for v in vals); c=sum(math.cos(math.radians(v)) for v in vals)
    return (math.degrees(math.atan2(s,c))+360)%360 if abs(s)+abs(c)>1e-12 else None

def raw_header(b,idx):
    tm,dt,alt=struct.unpack_from('<HHH',b,0x32); stamp=dos_dt(dt,tm)
    gyro=struct.unpack_from('<8h',b,0xD4); pressure=struct.unpack_from('<8H',b,0xE4)
    return {'index':idx,'utc':stamp,'altitude_m':float(alt),'heading_deg':circ_mean([(x/10.0)-10.1 for x in gyro]),'pressure_m_equiv':float(statistics.median([(x/10.0)-5.0 for x in pressure]))}

def dt_fields(date,hm):
    if len(date)!=6 or len(hm)!=4:return None
    try:return datetime(2000+int(date[:2]),int(date[2:4]),int(date[4:6]),int(hm[:2]),int(hm[2:4]),tzinfo=timezone.utc)
    except:return None

def parse_nav(raw):
    out=[]
    for line in raw.decode('utf-8','replace').splitlines():
        p=line.split()
        if len(p)<5:continue
        d=dt_fields(p[1],p[2])
        if d is None:continue
        try:lat=float(p[3]);lon=float(p[4])
        except:continue
        if -90<=lat<=90 and -180<=lon<=180:out.append({'utc':d,'lat':lat,'lon':lon})
    out.sort(key=lambda x:x['utc']);return out

def parse_cable(raw):
    out=[]
    for line in raw.decode('utf-8','replace').splitlines():
        p=line.split()
        if len(p)<4:continue
        d=dt_fields(p[1],p[2])
        if d is None:continue
        nums=[]
        for x in p[3:]:
            try:nums.append(float(x))
            except:pass
        if nums:out.append({'utc':d,'cable_m':nums[-1]})
    out.sort(key=lambda x:x['utc']);return out

def hav_m(a,b):
    r=6371008.8;p1,l1=map(math.radians,a);p2,l2=map(math.radians,b);dp=p2-p1;dl=l2-l1
    h=math.sin(dp/2)**2+math.cos(p1)*math.cos(p2)*math.sin(dl/2)**2
    return 2*r*math.asin(min(1.0,math.sqrt(h)))

def add_arc(nav):
    c=0.0
    for i,r in enumerate(nav):
        if i:c+=hav_m((nav[i-1]['lat'],nav[i-1]['lon']),(r['lat'],r['lon']))
        r['arc_m']=c

def interp_rows(rows,t,key):
    ts=[x['utc'].timestamp() for x in rows];q=t.timestamp();j=bisect.bisect_left(ts,q)
    if j<=0:return rows[0][key]
    if j>=len(rows):return rows[-1][key]
    a,b=rows[j-1],rows[j];f=(q-a['utc'].timestamp())/(b['utc'].timestamp()-a['utc'].timestamp())
    return a[key]+f*(b[key]-a[key])

def ship_state(nav,t):
    ts=[x['utc'].timestamp() for x in nav];q=t.timestamp();j=bisect.bisect_left(ts,q)
    if j<=0:return {'lat':nav[0]['lat'],'lon':nav[0]['lon'],'arc_m':nav[0]['arc_m']}
    if j>=len(nav):return {'lat':nav[-1]['lat'],'lon':nav[-1]['lon'],'arc_m':nav[-1]['arc_m']}
    a,b=nav[j-1],nav[j];f=(q-a['utc'].timestamp())/(b['utc'].timestamp()-a['utc'].timestamp())
    return {'lat':a['lat']+f*(b['lat']-a['lat']),'lon':a['lon']+f*(b['lon']-a['lon']),'arc_m':a['arc_m']+f*(b['arc_m']-a['arc_m'])}

def point_at_arc(nav,arc):
    arcs=[x['arc_m'] for x in nav];j=bisect.bisect_left(arcs,arc)
    if j<=0:return {'lat':nav[0]['lat'],'lon':nav[0]['lon']}
    if j>=len(nav):return {'lat':nav[-1]['lat'],'lon':nav[-1]['lon']}
    a,b=nav[j-1],nav[j]; den=b['arc_m']-a['arc_m'];f=0 if den<=0 else (arc-a['arc_m'])/den
    return {'lat':a['lat']+f*(b['lat']-a['lat']),'lon':a['lon']+f*(b['lon']-a['lon']}

def reconstruct(nav,cable,h):
    c=interp_rows(cable,h['utc'],'cable_m');p=h['pressure_m_equiv'];total=c+UMBILICAL_M
    if total<=abs(p):return {'valid':False,'cable_m':c,'pressure_m_equiv':p,'reason':'vertical_term_exceeds_total_cable'}
    horizontal=math.sqrt(total*total-p*p);ship=ship_state(nav,h['utc']);veh=point_at_arc(nav,ship['arc_m']-horizontal)
    return {'valid':True,'cable_m':c,'pressure_m_equiv':p,'total_cable_plus_umbilical_m':total,'horizontal_layback_m':horizontal,'ship':ship,'vehicle':veh}

def en_offset(origin,target):
    lat0=math.radians((origin[0]+target[0])/2); north=(target[0]-origin[0])*111195.08; east=(target[1]-origin[1])*111195.08*math.cos(lat0)
    return east,north

def map_target(h,rec):
    if not rec.get('valid') or h['heading_deg'] is None:return {'valid':False}
    e,n=en_offset((rec['vehicle']['lat'],rec['vehicle']['lon']),FROZEN);hr=math.radians(h['heading_deg'])
    along=e*math.sin(hr)+n*math.cos(hr);cross=e*math.cos(hr)-n*math.sin(hr);slant=math.sqrt(h['altitude_m']**2+cross**2);idx=int(round(slant/SAMPLE_M));side='stbd' if cross>=0 else 'port'
    return {'valid':True,'east_m':e,'north_m':n,'along_track_miss_m':along,'cross_track_m':cross,'ground_range_m':abs(cross),'slant_range_m':slant,'sample_index':idx,'side':side,'in_4000_samples':0<=idx<4000}

def local_z(vals,i,half=20):
    lo=max(0,i-half);hi=min(len(vals),i+half+1);bg=vals[lo:i]+vals[i+1:hi]
    if not bg:return None
    m=statistics.median(bg);mad=statistics.median(abs(x-m) for x in bg);sc=1.4826*mad
    return None if sc==0 else (vals[i]-m)/sc

def values_for_block(b,mapping):
    if not mapping.get('in_4000_samples'):return None
    i=mapping['sample_index'];side=mapping['side'];names=[side+'_sidescan',side+'_swath'];out={}
    for name in names:
        vals=list(struct.unpack_from('<4000h',b,OFF[name]));out[name]={'sample_index':i,'raw_signed':vals[i],'abs_amplitude':abs(vals[i]),'local_z':local_z(vals,i)}
    return out

def iso(x):return x.isoformat().replace('+00:00','Z') if x else None

def percentile(xs,p):
    if not xs:return None
    s=sorted(xs);q=(len(s)-1)*p;k=int(math.floor(q));f=q-k
    return s[k] if k==len(s)-1 else s[k]*(1-f)+s[k+1]*f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
    out={'schema':'janus.cosmos.cousteau.cd169.ground_fixed_nav_reconstruction.v1','status':'STARTED','contract':'JANUS-HANNAH-CD169-GROUND-FIXED-REPLICATION-CONTRACT-2026-08-26-v1.0','model_retuned':False,'intensity_used_for_navigation_or_ping_selection':False,'exact_prism_wireout_equivalence_claim':False}
    try:
        nr=get(NAV);cr=get(CABLE);nav=parse_nav(nr);cable=parse_cable(cr);add_arc(nav)
        out['inputs']={'nav_sha256':hashlib.sha256(nr).hexdigest(),'nav_rows':len(nav),'cable_sha256':hashlib.sha256(cr).hexdigest(),'cable_rows':len(cable)}
        headers=[];window_blocks={};rawhash=hashlib.sha256();idx=0;buf=bytearray();f=ftp();size=0
        try:
            sock=f.transfercmd('RETR '+RAW)
            try:
                while True:
                    chunk=sock.recv(1048576)
                    if not chunk:break
                    rawhash.update(chunk);size+=len(chunk);buf.extend(chunk)
                    while len(buf)>=BLOCK:
                        b=bytes(buf[:BLOCK]);del buf[:BLOCK];h=raw_header(b,idx)
                        if h['utc']:
                            headers.append(h)
                            if TARGET-WINDOW<=h['utc']<=TARGET+WINDOW:window_blocks[idx]=b
                        idx+=1
            finally:
                try:sock.close()
                except Exception:pass
        finally:
            try:f.close()
            except Exception:pass
        out['raw']={'sha256':rawhash.hexdigest(),'size_bytes':size,'records':idx,'parsed_headers':len(headers),'window_blocks':len(window_blocks)}
        hts=[h['utc'].timestamp() for h in headers]
        def nearest_header(t):
            j=bisect.bisect_left(hts,t.timestamp());cand=[]
            if j<len(headers):cand.append(headers[j])
            if j>0:cand.append(headers[j-1])
            return min(cand,key=lambda h:abs((h['utc']-t).total_seconds()))
        validations=[]
        for tim,lat,lon,logged_layback in ANCHORS:
            hh,mm,ss=map(int,tim.split(':'));t=datetime(2005,2,28,hh,mm,ss,tzinfo=timezone.utc);h=nearest_header(t);rec=reconstruct(nav,cable,h)
            residual=hav_m((rec['vehicle']['lat'],rec['vehicle']['lon']),(lat,lon)) if rec.get('valid') else None
            validations.append({'anchor_utc':iso(t),'raw_ping_utc':iso(h['utc']),'raw_ping_index':h['index'],'logged_tobi':{'lat':lat,'lon':lon,'layback_m':logged_layback},'raw_pressure_m_equiv':h['pressure_m_equiv'],'external_cable_m':rec.get('cable_m'),'reconstructed_horizontal_layback_m':rec.get('horizontal_layback_m'),'layback_delta_vs_logged_m':None if not rec.get('valid') else rec['horizontal_layback_m']-logged_layback,'reconstructed_vehicle':rec.get('vehicle'),'geodesic_residual_to_logged_tobi_m':residual})
        rs=[x['geodesic_residual_to_logged_tobi_m'] for x in validations if x['geodesic_residual_to_logged_tobi_m'] is not None]
        lds=[abs(x['layback_delta_vs_logged_m']) for x in validations if x['layback_delta_vs_logged_m'] is not None]
        out['historical_model_validation']={'role':'HISTORICAL_MODEL_REPRODUCTION_NOT_INDEPENDENT_GROUND_TRUTH','anchors':validations,'summary':{'n':len(rs),'residual_median_m':statistics.median(rs) if rs else None,'residual_p90_m':percentile(rs,0.9),'residual_max_m':max(rs) if rs else None,'abs_layback_delta_median_m':statistics.median(lds) if lds else None,'abs_layback_delta_p90_m':percentile(lds,0.9),'acceptance_threshold_m':None}}
        window=[]
        for h in headers:
            if TARGET-WINDOW<=h['utc']<=TARGET+WINDOW:
                rec=reconstruct(nav,cable,h);mp=map_target(h,rec);window.append({'index':h['index'],'utc':iso(h['utc']),'altitude_m':h['altitude_m'],'heading_deg':h['heading_deg'],'pressure_m_equiv':h['pressure_m_equiv'],'reconstruction':rec,'target_mapping':mp})
        before=[x for x in window if datetime.fromisoformat(x['utc'].replace('Z','+00:00'))<=TARGET];after=[x for x in window if datetime.fromisoformat(x['utc'].replace('Z','+00:00'))>TARGET]
        time_selected=[]
        if before:time_selected.append(max(before,key=lambda x:x['utc']))
        if after:time_selected.append(min(after,key=lambda x:x['utc']))
        geom_sorted=sorted([x for x in window if x['target_mapping'].get('valid')],key=lambda x:(abs(x['target_mapping']['along_track_miss_m']),abs(datetime.fromisoformat(x['utc'].replace('Z','+00:00')).timestamp()-TARGET.timestamp()),x['index']))[:2]
        selected_ids=sorted(set([x['index'] for x in time_selected+geom_sorted]));selected=[]
        for x in window:
            if x['index'] in selected_ids:
                y=json.loads(json.dumps(x));b=window_blocks.get(x['index']);y['measurements_after_geometry_freeze']=values_for_block(b,x['target_mapping']) if b else None;selected.append(y)
        target_h=nearest_header(TARGET);target_rec=reconstruct(nav,cable,target_h);target_map=map_target(target_h,target_rec)
        out['target_ping']={'index':target_h['index'],'utc':iso(target_h['utc']),'delta_to_frozen_target_s':(target_h['utc']-TARGET).total_seconds(),'reconstruction':target_rec,'mapping':target_map}
        out['same_line_multiping']={'window_utc':[iso(TARGET-WINDOW),iso(TARGET+WINDOW)],'time_only_selected_indices':[x['index'] for x in time_selected],'geometry_only_two_smallest_along_miss_indices':[x['index'] for x in geom_sorted],'selected_with_measurements':selected,'authority':'INTERNAL_SAME_PASS_REPLICATION_ONLY','independent_crossing':False}
        out['navigation_uncertainty_statement']='Residuals to cruise-log TOBI positions quantify reproduction of a historical layback-model product only; they do not establish absolute per-ping ground accuracy.'
        out['status']='WIREOUT_TRACK_RECONSTRUCTION_READY__G7A_GEOMETRY_MAPPED__G7B_NOT_TESTED'
    except Exception as e:
        out['status']='GROUND_FIXED_NAV_RECONSTRUCTION_FAILED';out['error_type']=type(e).__name__;out['error']=str(e)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    s=(out.get('historical_model_validation') or {}).get('summary',{});print(json.dumps({'status':out['status'],'validation_summary':s,'target_ping':out.get('target_ping'),'same_line_indices':{k:v for k,v in (out.get('same_line_multiping') or {}).items() if k.endswith('indices')}},indent=2));return 0 if out['status'].startswith('WIREOUT_TRACK_RECONSTRUCTION_READY') else 2
if __name__=='__main__':raise SystemExit(main())
