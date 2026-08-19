from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, ICRS, SkyCoord, get_body
from astropy.time import Time


def _time_grid(start: str, stop: str, step_seconds: int) -> Time:
    t0 = Time(start, scale='utc'); t1 = Time(stop, scale='utc')
    span_s = float((t1 - t0).to_value(u.s))
    n = int(math.floor(span_s / step_seconds)) + 1
    return t0 + np.arange(n, dtype=float) * step_seconds * u.s


def _unit_from_altaz_arrays(alt_deg, az_deg):
    alt = np.deg2rad(np.asarray(alt_deg, dtype=float)); az = np.deg2rad(np.asarray(az_deg, dtype=float))
    return np.column_stack([np.cos(alt)*np.sin(az), np.cos(alt)*np.cos(az), np.sin(alt)])


def _unit_from_altaz(alt_deg, az_deg):
    return _unit_from_altaz_arrays([alt_deg],[az_deg])[0]


def _altaz_from_unit(vec):
    v=np.asarray(vec,dtype=float); v=v/np.linalg.norm(v)
    alt=math.degrees(math.asin(float(np.clip(v[2],-1,1))))
    az=math.degrees(math.atan2(float(v[0]),float(v[1])))%360.0
    return alt,az


def _angle_deg_rows(a,b):
    dot=np.sum(a*b,axis=1)/(np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1))
    return np.rad2deg(np.arccos(np.clip(dot,-1,1)))


def _radec_strings(coord):
    c=coord.icrs
    return {'ra_deg':float(c.ra.deg),'dec_deg':float(c.dec.deg),
            'ra_hms':c.ra.to_string(unit=u.hourangle,sep=':',precision=3,pad=True),
            'dec_dms':c.dec.to_string(unit=u.deg,sep=':',precision=3,alwayssign=True,pad=True)}


def _evaluate(times,cfg,location,love):
    frame=AltAz(obstime=times,location=location,pressure=0*u.hPa)
    love_altaz=love.transform_to(frame)
    moon=get_body('moon',times,location=location)
    moon_altaz=moon.transform_to(frame)
    love_alt=np.asarray(love_altaz.alt.deg,float); love_az=np.asarray(love_altaz.az.deg,float)
    moon_alt=np.asarray(moon_altaz.alt.deg,float); moon_az=np.asarray(moon_altaz.az.deg,float)
    love_vec=_unit_from_altaz_arrays(love_alt,love_az); moon_vec=_unit_from_altaz_arrays(moon_alt,moon_az)
    normal_alt=float(cfg['khufu_faces']['face_normal_altitude_deg'])
    out=[]
    for face,face_az in cfg['khufu_faces']['outward_normal_azimuth_deg'].items():
        n=_unit_from_altaz(normal_alt,float(face_az))
        love_dot=love_vec@n
        edem_vec=2.0*love_dot[:,None]*n[None,:]-love_vec
        edem_vec/=np.linalg.norm(edem_vec,axis=1)[:,None]
        edem_alt=np.rad2deg(np.arcsin(np.clip(edem_vec[:,2],-1,1)))
        moon_dot=moon_vec@n; edem_dot=edem_vec@n
        sep=_angle_deg_rows(edem_vec,moon_vec)
        valid=(love_alt>0)&(moon_alt>0)&(edem_alt>0)&(love_dot>0)&(moon_dot>0)&(edem_dot>0)
        if not np.any(valid):
            out.append({'face':face,'status':'NO_VALID_SAMPLES'}); continue
        ids=np.flatnonzero(valid); j=int(ids[np.argmin(sep[ids])])
        ealt,eaz=_altaz_from_unit(edem_vec[j])
        local=AltAz(obstime=times[j],location=location,pressure=0*u.hPa)
        eicrs=SkyCoord(az=eaz*u.deg,alt=ealt*u.deg,frame=local).transform_to(ICRS())
        out.append({'face':face,'status':'OK','time_utc':times[j].utc.isot,
                    'angular_error_to_moon_deg':float(sep[j]),
                    'love_alt_deg':float(love_alt[j]),'love_az_deg':float(love_az[j]),
                    'moon_alt_deg':float(moon_alt[j]),'moon_az_deg':float(moon_az[j]),
                    'edem_alt_deg':float(ealt),'edem_az_deg':float(eaz),'edem_icrs':_radec_strings(eicrs)})
    return out


def _best(rows):
    ok=[r for r in rows if r.get('status')=='OK']
    if not ok: raise RuntimeError('No physically valid samples')
    return min(ok,key=lambda r:r['angular_error_to_moon_deg'])


def run(prereg_path:Path,output_path:Path):
    cfg=json.loads(prereg_path.read_text(encoding='utf-8'))
    o=cfg['observer']
    loc=EarthLocation.from_geodetic(lon=float(o['longitude_deg_east'])*u.deg,lat=float(o['latitude_deg'])*u.deg,height=float(o['height_m'])*u.m)
    t=cfg['love_target']
    love=SkyCoord(ra=float(t['ra_deg_icrs'])*u.deg,dec=float(t['dec_deg_icrs'])*u.deg,distance=float(t['distance_pc'])*u.pc,frame='icrs')
    tw=cfg['time_window']
    coarse=_evaluate(_time_grid(tw['start_utc'],tw['stop_utc'],int(tw['coarse_step_seconds'])),cfg,loc,love)
    cb=_best(coarse)
    center=Time(cb['time_utc'],scale='utc'); half=int(tw['refine_half_window_seconds']); step=int(tw['refine_step_seconds'])
    rb=_best(_evaluate(center+np.arange(-half,half+step,step,dtype=float)*u.s,cfg,loc,love))
    err=float(rb['angular_error_to_moon_deg']); gates=cfg['reverse_spear_geometry']
    if err<=float(gates['exact_center_gate_deg']): gate='EXACT_CENTER_GATE_PASS'
    elif err<=float(gates['strong_gate_deg']): gate='STRONG_GATE_PASS'
    elif err<=float(gates['primary_gate_deg']): gate='PRIMARY_GATE_PASS'
    else: gate='GEOMETRY_GATE_FAIL'
    base=cfg.get('comparison_baseline') or {}
    giza=None
    if 'giza_edem_ra_deg' in base:
        c1=SkyCoord(ra=rb['edem_icrs']['ra_deg']*u.deg,dec=rb['edem_icrs']['dec_deg']*u.deg)
        c0=SkyCoord(ra=float(base['giza_edem_ra_deg'])*u.deg,dec=float(base['giza_edem_dec_deg'])*u.deg)
        giza={'angular_separation_deg':float(c1.separation(c0).deg),'delta_ra_deg':float(rb['edem_icrs']['ra_deg']-float(base['giza_edem_ra_deg'])),'delta_dec_deg':float(rb['edem_icrs']['dec_deg']-float(base['giza_edem_dec_deg']))}
    result={'schema':'janus.cosmos.edem.observer_reverse_spear.result.v1','experiment_id':cfg['experiment_id'],
            'generated_utc':datetime.now(timezone.utc).isoformat(),'observer':o,
            'frame_policy':{'all_reflection_geometry_observer':o['name'],'telescope_locations_used_for_geometry':False,'love_projected_topocentrically_at_finite_distance':True,'moon_projected_topocentrically_from_observer':True,'reported_celestial_coordinate_frame':'ICRS'},
            'coarse_best':cb,'refined_best':rb,'edem_geometry_direction_candidate':rb['edem_icrs'],'gate_status':gate,
            'comparison_to_giza':giza,'edem_identity_confirmed':False,'claim_ceiling':cfg['claim_ceiling']}
    output_path.parent.mkdir(parents=True,exist_ok=True); output_path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'gate_status':gate,'face':rb['face'],'time_utc':rb['time_utc'],'error_deg':err,'edem':rb['edem_icrs'],'comparison_to_giza':giza},indent=2))
    return result


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--prereg',required=True); ap.add_argument('--output',required=True); a=ap.parse_args(); run(Path(a.prereg),Path(a.output))

if __name__=='__main__': main()
