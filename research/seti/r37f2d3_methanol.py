#!/usr/bin/env python3
"""Historical natural methanol astronomy: checksum-gated spectral reduction."""
import hashlib,io,json,math,os,pathlib,sys,time,zipfile,urllib.parse
import numpy as np
import requests
import astropy
import astropy.units as u
from astropy.coordinates import EarthLocation,SkyCoord,SpectralCoord,solar_system_ephemeris
from astropy.time import Time
from astropy.utils import iers
import astropy_iers_data
from r37f1de_archive import fil_header,fil_grid
from r37f2d2_integrity import locator
P=pathlib.Path("research/seti/r37f2d3_methanol_protocol.json");S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f2d3");O.mkdir(parents=True,exist_ok=True)
M=S["method"];REST=S["rest_mhz"];C=299792.458
def save(name,d):(O/name).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def refs():
    p=S["parent"];r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['artifact']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90)
    r.raise_for_status();b=r.content;assert len(b)<16*1024**2 and hashlib.sha256(b).hexdigest()==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        n=[v for v in z.infolist() if pathlib.PurePosixPath(v.filename).name=="result.json"]
        assert len(n)==1 and n[0].file_size<8*1024**2
        return json.loads(z.read(n[0]))["files"]
def choose(src,rr):
    a=[v for v in rr if v["source"]==src["name"] and v["day"]==src["day"] and v["status"]=="HTTPS_PUBLISHED_MD5_JOIN" and v["reference"]["md5"]==src["md5"]]
    assert len(a)==1
    v=a[0];ref=v["reference"]
    assert locator(v["url"]) and locator(v["url"])==locator(ref["url"])
    assert urllib.parse.urlparse(ref["portal_page"]).scheme=="https" and urllib.parse.urlparse(ref["portal_page"]).hostname in ("breakthroughinitiatives.org","www.breakthroughinitiatives.org")
    assert 0<ref["size"]<=M["full_bytes_cap_each"]
    return v
def download(v,src,path):
    ref=v["reference"];md5=hashlib.md5();sha=hashlib.sha256();n=0;t=time.monotonic()
    with requests.get(v["url"],stream=True,timeout=(15,90),allow_redirects=False,headers={"Accept-Encoding":"identity"}) as r:
        r.raise_for_status();assert r.status_code==200 and r.headers.get("Content-Encoding","identity")=="identity"
        if r.headers.get("Content-Length","").isdigit():assert int(r.headers["Content-Length"])==ref["size"]
        with path.open("wb") as f:
            for chunk in r.iter_content(4*1024**2):
                n+=len(chunk);assert n<=ref["size"] and n<=M["full_bytes_cap_each"]
                md5.update(chunk);sha.update(chunk);f.write(chunk)
    assert n==ref["size"] and md5.hexdigest()==ref["md5"],"archive byte count or published MD5 mismatch"
    if src.get("expected_sha256"):assert sha.hexdigest()==src["expected_sha256"],"previous SHA256 mismatch"
    return {"bytes":n,"md5":md5.hexdigest(),"sha256":sha.hexdigest(),"seconds":time.monotonic()-t}
def pointing(d):
    def angle(q,ra):
        q=float(q);sign=1 if q>=0 else -1;q=abs(q)
        aa=int(q//10000);bb=int(q//100)%100;cc=q%100
        assert 0<=bb<60 and 0<=cc<60 and (aa<24 if ra else aa<=90)
        return sign*(aa+bb/60+cc/3600)*(15 if ra else 1)
    return SkyCoord(ra=angle(d["src_raj"],True)*u.deg,dec=angle(d["src_dej"],False)*u.deg,frame="icrs")
def inspect(path,v,src):
    with path.open("rb") as f:h,n=fil_header(f.read(65536))
    g=fil_grid(h,n,path.stat().st_size)
    assert g["source"].strip().upper()==v["target"] and abs(g["mjd"]-v["mjd"])*86400<=3
    assert h["nbits"]==32 and g["shape"][1]==1 and g["shape"][0]>=M["minimum_rows"]
    assert abs(abs(g["df"])*1e6-M["native_df_hz"])<1e-7 and abs(g["ts"]-M["native_ts_s"])<1e-7
    assert h.get("barycentric") in (0,None) and h.get("pulsarcentric") in (0,None)
    p=pointing(h);target=SkyCoord(ra=src["ra_deg"]*u.deg,dec=src["dec_deg"]*u.deg,frame="icrs")
    sep=float(p.separation(target).deg);assert sep<=.05,"header pointing inconsistent with published source"
    loc=EarthLocation.from_geodetic(S["site"]["longitude_deg"]*u.deg,S["site"]["latitude_deg"]*u.deg,S["site"]["height_m"]*u.m)
    # Calculated at start, middle and end; conditional on topocentric input.
    tar=SkyCoord(ra=target.ra,dec=target.dec,distance=2*u.kpc,pm_ra_cosdec=0*u.mas/u.yr,pm_dec=0*u.mas/u.yr,radial_velocity=0*u.km/u.s,frame="icrs")
    factors=[]
    for x in (0,.5,1):
        t=Time(g["mjd"]+x*g["ts"]*g["shape"][0]/86400,format="mjd",scale="utc")
        obs=loc.get_itrs(obstime=t)
        a=SpectralCoord(REST*u.MHz,observer=obs,target=tar,doppler_convention="radio",doppler_rest=REST*u.MHz)
        b=a.with_observer_stationary_relative_to("lsrk");back=b.with_observer_stationary_relative_to(obs)
        k=float(b.to_value(u.MHz))/REST
        assert .999<k<1.001 and abs(float(back.to_value(u.MHz))-REST)<1e-5
        factors.append(k)
    k=factors[1];vlo,vhi=src["lsr_window"];margin=M["baseline_margin_km_s"]
    fmin,fmax=sorted([REST*(1-(vlo-margin)/C)/k,REST*(1-(vhi+margin)/C)/k])
    assert g["bounds"][0]<=fmin<fmax<=g["bounds"][1],"science window outside native band"
    a,b=sorted([(fmin-h["fch1"])/g["df"],(fmax-h["fch1"])/g["df"]]);lo=math.ceil(a);hi=math.floor(b)+1
    B=M["block_bins"];hi=lo+(hi-lo)//B*B
    assert 0<=lo<hi<=g["shape"][-1] and hi-lo>=B*100
    mm=np.memmap(path,dtype="<f4",mode="r",offset=n,shape=tuple(g["shape"]))
    x=np.array(mm[:,0,lo:hi],dtype=np.float64);del mm
    assert np.isfinite(x).all() and (x>=0).all()
    vc=C*(1-(h["fch1"]+np.arange(lo,hi,dtype=float)*g["df"])*k/REST)
    off=(vc<vlo-margin/2)|(vc>vhi+margin/2)
    scale=np.median(x[:,off],axis=1);assert (scale>0).all()
    y=(x/scale[:,None]).reshape(len(x),-1,B).mean(2)
    vv=vc.reshape(-1,B).mean(1);order=np.argsort(vv);vv=vv[order];y=y[:,order]
    signal=(vv>=vlo)&(vv<=vhi);control=(vv<vlo-margin/2)|(vv>vhi+margin/2)
    assert signal.sum()>=10 and control.sum()>=100
    def stats(a):
        z=a.mean(0);med=float(np.median(z[control]));sig=1.4826*float(np.median(abs(z[control]-med)));assert sig>0
        return z,med,sig
    z,med,sig=stats(y);j=int(np.flatnonzero(signal)[np.argmax(z[signal])])
    dz=float((z[j]-med)/sig);ex=float(z[j]-med)
    halves=[]
    for part in np.array_split(y,2,axis=0):
        zz,mm,ss=stats(part);halves.append({"rows":len(part),"descriptive_z_at_mean_peak":float((zz[j]-mm)/ss),"relative_excess":float(zz[j]-mm)})
    adj=bool((j>0 and signal[j-1] and z[j-1]>med+3*sig) or (j+1<len(z) and signal[j+1] and z[j+1]>med+3*sig))
    support=bool(dz>=M["threshold_descriptive_z"] and ex>=M["minimum_relative_excess"] and adj and all(a["descriptive_z_at_mean_peak"]>=M["threshold_descriptive_z"] and a["relative_excess"]>=M["minimum_relative_excess"] for a in halves))
    # Retain only coarse astronomical spectrum, no fine-grid fingerprints.
    fname=src["name"].replace("+","_").replace(".","_")+"_science_profile.npz"
    np.savez_compressed(O/fname,velocity_km_s=vv,normalized_mean=z)
    return {"status":"CONDITIONAL_NATURAL_LINE_SUPPORT" if support else "HOLD_LINE_CRITERIA","native_df_hz":abs(g["df"])*1e6,"native_ts_s":g["ts"],"rows":g["shape"][0],"duration_s":g["ts"]*g["shape"][0],"pointing_separation_deg":sep,"barycentric":h.get("barycentric"),"pulsarcentric":h.get("pulsarcentric"),"frame_explicitly_confirmed":False,"lsrk_factor_start_mid_end":factors,"lsr_window":src["lsr_window"],"peak_conditional_lsrk_km_s":float(vv[j]),"descriptive_z":dz,"relative_excess":ex,"adjacent_support":adj,"halves":halves,"science_bin_hz":B*abs(g["df"])*1e6,"science_bin_km_s":float(np.median(np.diff(vv))),"spatial_off_verified":False,"profile":fname}
def main():
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
    rr=refs();selected=[choose(s,rr) for s in S["sources"]];assert len({v["day"] for v in selected})==len(selected)==2
    save("frozen_selection.json",{"files":selected,"frozen_before_power_read":True})
    rows=[]
    for src,v in zip(S["sources"],selected):
        r={"source":src["name"],"day":src["day"],"status":"HOLD"};p=O/"temporary.fil"
        try:
            r["integrity"]=download(v,src,p)
            r.update(inspect(p,v,src))
        except Exception as e:
            import traceback
            traceback.print_exc();r["error"]=type(e).__name__+": "+str(e)[:400]
        finally:p.unlink(missing_ok=True)
        rows.append(r);print("R37F2D3_SOURCE",json.dumps(r),flush=True)
    n=sum(r["status"]=="CONDITIONAL_NATURAL_LINE_SUPPORT" for r in rows)
    d={"status":"CONDITIONAL_NATURAL_METHANOL_DEVELOPMENT_CONTROLS" if n==len(rows) else "NATURAL_METHANOL_VALIDATION_HOLD","sources":rows,"full_integrity_matches":sum("integrity" in r for r in rows),"line_support_count":n,"distinct_dates_supported":len({r["day"] for r in rows if r["status"]=="CONDITIONAL_NATURAL_LINE_SUPPORT"}),"astropy":astropy.__version__,"iers_data":astropy_iers_data.__version__,"spatial_off_verified":False,"exact_bl_frame_verified":False,"fine_grid_protection_verified":False,"production_adoption":False}
    save("result.json",d);print("R37F2D3_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
