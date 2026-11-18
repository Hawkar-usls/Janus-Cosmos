#!/usr/bin/env python3
"""R37F2A: archived natural OH sky control; conditional LSRK transform and repeat access audit."""
import collections,hashlib,io,json,math,os,pathlib,re,urllib.parse,zipfile,warnings
import numpy as np
import requests,bs4,h5py,hdf5plugin
import astropy,astropy_iers_data
from astropy import units as u
from astropy.coordinates import EarthLocation,SkyCoord,SpectralCoord,solar_system_ephemeris
from astropy.time import Time
from astropy.utils import iers
from astropy.constants import c as light
import r4_chunk_replay as base
import r37f2_w3oh as parent
import r37f1de_archive as fil
import r9_expanded_search as eng
P=pathlib.Path("research/seti/r37f2a_natural_lsrk_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r37f2a");O.mkdir(parents=True,exist_ok=True)
def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,x):(O/n).write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def previous():
    p=S["parent"];r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['artifact']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=(15,90));r.raise_for_status()
    b=r.content;assert len(b)<16*base.MIB and sha(b)==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        def load(n):
            a=[i for i in z.infolist() if pathlib.PurePosixPath(i.filename).name==n];assert len(a)==1 and a[0].file_size<8*base.MIB
            return json.loads(z.read(a[0]))
        m=load("frozen_manifest.json");d=load("result.json")
    assert m["frozen_before_power_read"] and d["candidate_dates"]==1 and len(m["pairs"])==1
    return m["pairs"][0],d["results"][0]
def transform_setup():
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error"
    table=iers.IERS_B.open(iers.IERS_B_FILE)
    iers.earth_orientation_table.set(table);solar_system_ephemeris.set("builtin")
    s=S["source"];g=s["location"]
    site=EarthLocation.from_geodetic(g["lon_deg"]*u.deg,g["lat_deg"]*u.deg,g["height_m"]*u.m)
    target=SkyCoord(ra=s["ra_deg"]*u.deg,dec=s["dec_deg"]*u.deg,distance=2*u.kpc,radial_velocity=0*u.km/u.s,frame="icrs")
    return table,site,target,{"astropy":astropy.__version__,"iers_data":astropy_iers_data.__version__,"iers_b_sha256":sha(pathlib.Path(iers.IERS_B_FILE).read_bytes()),"ephemeris":"builtin","coordinate_reference":g["reference"]}
def factors(mjd,table,site,target):
    assert float(table["MJD"][0].value)<=mjd<=float(table["MJD"][-1].value)
    t=Time(mjd,format="mjd",scale="utc");obs=site.get_itrs(obstime=t)
    nu=1667.359*u.MHz
    x=SpectralCoord(nu,observer=obs,target=target)
    direct=x.with_observer_stationary_relative_to("lsrk")
    icrs=x.with_observer_stationary_relative_to("icrs")
    chained=icrs.with_observer_stationary_relative_to("lsrk")
    back=direct.with_observer_stationary_relative_to(obs)
    diff=abs((direct-chained).to_value(u.Hz));rt=abs((back-x).to_value(u.Hz))
    fac=float(direct.to_value(u.MHz)/nu.to_value(u.MHz))
    assert diff<=S["transform"]["tolerance_hz"] and rt<=S["transform"]["tolerance_hz"],"frame consistency check"
    assert 0.9997<fac<1.0003
    # Independent approximate classical correction: barycentric + projection of LSRK solar apex.
    apex=SkyCoord(ra=270*u.deg,dec=30*u.deg,frame="fk4",equinox="B1900").icrs
    bary=target.radial_velocity_correction("barycentric",obstime=t,location=site).to_value(u.km/u.s)
    lsr_approx=float(bary+20*math.cos(apex.separation(target).rad))
    dv=float(light.to_value(u.km/u.s)*(1-fac))
    assert abs(dv-lsr_approx)<=S["transform"]["approx_velocity_tolerance_km_s"],"LSRK approximate sign/scale check"
    return {"mjd":mjd,"utc":t.isot,"topo_to_lsrk_frequency_factor":fac,"radio_velocity_correction_at_rest_km_s":dv,"barycentric_plus_solar_apex_approx_km_s":lsr_approx,"approx_difference_km_s":dv-lsr_approx,"direct_vs_chained_hz":float(diff),"roundtrip_hz":float(rt)}
def access():
    vals,meta=parent.catalog();days=S["source"]["days"];summary=[]
    for day in days:
        report={"day":day,"probes":[]};q=[v for v in vals if v["day"]==day and v["ext"]==".fil"]
        report["catalogue_fil_count"]=len(q);report["roles"]=dict(collections.Counter(v["role"] for v in q))
        # A different official listing may advertise migrated URLs. Exact basename only.
        mirrors={}
        for alias in ("DIAG_W3OH","DIAG_W3OH_OFF"):
            try:
                b,url=base.fetch(S["archive"]["official_https_index"],{"target":alias,"project":"GBT","file_type":"filterbank","perPage":100,"search":"Search","order":"asc","sort":"utc"},limit=8*base.MIB)
                soup=bs4.BeautifulSoup(b,"html.parser")
                for a in soup.select("a[href]"):
                    v=urllib.parse.urljoin(url,a["href"])
                    if base.archive_url(v) and v.lower().endswith(".fil"):mirrors.setdefault(base.basename(v),set()).add(v)
            except Exception as e:report.setdefault("index_errors",[]).append(str(e)[:180])
        def rank(v):
            try:band=abs(float(v["metadata"]["center_freq"])-1666.0)
            except (ValueError,TypeError,KeyError):band=1e9
            return (band,".gpuspec.0000." not in v["url"],v["mjd"],v["url"])
        selected=[]
        for role in ("ON","OFF"):
            selected.extend(sorted([v for v in q if v["role"]==role],key=rank)[:3])
        selected=selected[:S["archive"]["max_probes_per_day"]]
        for v in selected:
            original=v["url"];p=urllib.parse.urlparse(original)
            routes=sorted(mirrors.get(base.basename(original),set()))
            secure=urllib.parse.urlunparse(p._replace(scheme="https"))
            if secure not in routes:routes.append(secure)
            if p.scheme=="http":routes.append(original)
            r={"catalogue_target":v["target"],"role":v["role"],"mjd":v["mjd"],"attempts":[],"status":"HOLD"}
            for url in routes:
                a={"url":url}
                try:
                    g=fil.get_fil(url,eng)
                    assert str(g["source"]).strip().upper()==v["target"]
                    assert abs(g["mjd"]-v["mjd"])<=3/86400
                    a.update(status="HEADER_READ",channel_hz=abs(g["df"])*1e6,integration_s=g["ts"],rows=g["shape"][0],covers_oh=all(g["bounds"][0]<=x<=g["bounds"][1] for x in S["source"]["rest_mhz"]),transport="HTTPS" if urllib.parse.urlparse(url).scheme=="https" else "HTTP_UNAUTHENTICATED")
                    r["status"]="HEADER_ONLY" if a["transport"]=="HTTPS" else "UNAUTHENTICATED_HEADER_ONLY"
                    r["attempts"].append(a);break
                except Exception as e:a.update(status="HOLD",error=type(e).__name__+": "+str(e)[:220]);r["attempts"].append(a)
            report["probes"].append(r)
        report["https_headers"]=sum(v["status"]=="HEADER_ONLY" for v in report["probes"])
        report["independent_natural_pair_qualified"]=False
        summary.append(report)
    return {"catalogue":meta,"dates":summary}
def replay(q,old,table,site,target):
    hs=q["headers"];windows=[];kshift=round((hs[0]["f0"]-hs[1]["f0"])/hs[0]["df"])
    for rest in parent.C["line_rest_mhz"]:
        lo,hi=parent.window(rest);a,b=parent.bounds(hs[0],lo,hi);a=max(a,-kshift);b=min(b,hs[1]["shape"][-1]-kshift);b=a+((b-a)//64)*64
        assert b-a>=64*8;windows.append((rest,a,b))
    flags=[];arrays=[]
    for k,url in enumerate((q["on"],q["off"])):
        h=hs[k]
        with base.RangeFile(url,S["transform"]["max_bytes_per_file"]) as remote:
            with h5py.File(remote,"r",rdcc_nbytes=128*base.MIB) as f:
                d=f["data"];a=d.attrs;assert list(d.shape)==h["shape"] and abs(float(a["fch1"])-h["f0"])<1e-10 and abs(float(a["foff"])-h["df"])<1e-12 and abs(float(a["tstart"])-h["mjd"])<1e-10
                flag={key:str(a[key]) for key in a.keys() if any(x in key.lower() for x in ("bary","pulsar","frame","specsys","velref","doppler"))};flags.append(flag)
                shift=k*kshift
                w=[(lo+shift,hi+shift,0) for _,lo,hi in windows]
                plan=base.preflight(d,w,d.shape[0]);assert not plan or plan["compressed_bytes"]+remote.total+8*base.MIB<=S["transform"]["max_bytes_per_file"]
                aa=[]
                for _,lo,hi in windows:
                    x=np.asarray(d[:,0,lo+shift:hi+shift],dtype=np.float64);assert np.isfinite(x).all() and (x>=0).all()
                    med=np.median(x,axis=1);assert (med>0).all();aa.append((x/med[:,None]).reshape(x.shape[0],-1,64).mean(2))
                arrays.append(aa)
    # Conservative frame authority: explicit topocentric keyword only; zero flags alone do not establish it.
    explicit=all(any(k.lower() in ("specsys","freq_frame","frequency_frame") and v.upper().strip("B'\" ") in ("TOPOCENT","TOPOCENTRIC","TOPO") for k,v in f.items()) and not any(k.lower() in ("barycentric","pulsarcentric") and v not in ("0","False","0.0") for k,v in f.items()) for f in flags)
    mid=hs[0]["mjd"]+hs[0]["shape"][0]*hs[0]["ts"]/2/86400
    fc=factors(mid,table,site,target);out=[]
    for i,(rest,lo,hi) in enumerate(windows):
        a,b=arrays[0][i],arrays[1][i];diff=a.mean(0)-b.mean(0);center=float(np.median(diff));sig=1.4826*float(np.median(abs(diff-center)));assert sig>0
        ix=int(np.argmax(diff));z=float((diff[ix]-center)/sig);assert abs(z-old["lines"][i]["descriptive_peak_z"])<1e-6
        f0,df=hs[0]["f0"],hs[0]["df"];edges=sorted([f0+(lo+64*ix-.5)*df,f0+(lo+64*(ix+1)-.5)*df])
        fac=fc["topo_to_lsrk_frequency_factor"]
        vel=sorted([light.to_value(u.km/u.s)*(1-f*fac/rest) for f in edges])
        out.append({"rest_transition_mhz":rest,"replayed_descriptive_z":z,"block_center_observed_mhz":float(np.mean(edges)),"conditional_lsrk_radio_velocity_km_s":float(np.mean(vel)),"conditional_block_velocity_interval_km_s":vel,"within_broad_literature_velocity_band":bool(-60<=np.mean(vel)<=-30)})
    return {"frame_header_flags":flags,"explicit_topocentric_frame":explicit,"transform":fc,"lines":out,"frame_status":"EXPLICIT_TOPOCENTRIC_HEADER" if explicit else "CONDITIONAL_ONLY_FRAME_NOT_EXPLICIT","independent_natural_origin_confirmed":False}
def main():
    q,old=previous();table,site,target,provenance=transform_setup()
    d={"status":"R37F2A_AUDIT_ONLY","provenance":provenance,"production_adoption":False,"independent_natural_origin_confirmed":False}
    try:d["parent_replay"]=replay(q,old,table,site,target)
    except Exception as e:d["parent_replay"]={"status":"HOLD","error":type(e).__name__+": "+str(e)[:300]}
    try:d["repeat_access"]=access()
    except Exception as e:d["repeat_access"]={"status":"HOLD","error":type(e).__name__+": "+str(e)[:300]}
    save("result.json",d);print("R37F2A_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
