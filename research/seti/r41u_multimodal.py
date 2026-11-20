"""R41U: nonradio transfer pilots; no origin classifier."""
import sys, json, pathlib, time, hashlib, io, csv, math
import numpy as np
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astroquery.skyview import SkyView
from scipy.ndimage import gaussian_filter
O=pathlib.Path("runtime/r41u");O.mkdir(parents=True,exist_ok=True)
P=pathlib.Path("research/seti/r41u_multimodal_protocol.json");S=json.loads(P.read_text())
def sha(x):return hashlib.sha256(x).hexdigest()
def write(name,d):
    (O/name).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def read(name):return json.loads((O/name).read_text())
def get(url,params):
    r=requests.get(url,params=params,timeout=75);r.raise_for_status();return r
def num(x):
    try:
        v=float(x);return v if np.isfinite(v) else None
    except (ValueError,TypeError):return None
def psf_fit(data,xx,yy,fwhm):
    # Fixed-center vs alternative source centers, with the same Gaussian/background family.
    good=np.isfinite(data)&(xx**2+yy**2<=18**2)
    coords=np.argwhere(good);y=data[good];x=xx[good];z=yy[good]
    if len(y)<40:return None
    sig=fwhm/2.354820045
    grid=[(0.,0.)]+[(float(i),float(j)) for i in range(-6,7) for j in range(-6,7) if i*i+j*j<=36 and (i or j)]
    D=np.array([np.exp(-((x-i)**2+(z-j)**2)/(2*sig**2)) for i,j in grid]).T
    parity=(coords[:,0]+coords[:,1])%2;out=[]
    for fold in (0,1):
        a=parity==fold;t=~a
        if min(a.sum(),t.sum())<12:return None
        d=D[a];yc=y[a]-y[a].mean();dc=d-d.mean(0)
        amp=np.maximum(0,(dc.T@yc)/np.maximum((dc**2).sum(0),1e-30));bg=y[a].mean()-d.mean(0)*amp
        err=((d*amp+bg-y[a,None])**2).mean(0);idx=int(err.argmin())
        ev=((D[t]*amp+bg-y[t,None])**2).mean(0)
        out.append({"offset_arcsec":list(grid[idx]),"offset_norm_arcsec":float(np.hypot(*grid[idx])),"fixed_loss":float(ev[0]),"rival_loss":float(ev[idx]),"heldout_gain":float(ev[0]-ev[idx]),"fixed_nonnegative_amplitude":float(amp[0]),"rival_nonnegative_amplitude":float(amp[idx])})
    return {"folds":out,"both_gain_positive":all(x["heldout_gain"]>0 for x in out),"fold_offset_agreement_arcsec":float(np.linalg.norm(np.array(out[0]["offset_arcsec"])-out[1]["offset_arcsec"])),"scope":"correlated checkerboard pixels; not independent exposures or calibrated significance"}
def controls(data,xx,yy):
    def ap(x,y):
        rr=np.hypot(xx-x,yy-y);a=(rr<=5)&np.isfinite(data);b=(rr>=10)&(rr<=15)&np.isfinite(data)
        return float((data[a]-np.median(data[b])).sum()) if a.sum()>=3 and b.sum()>=12 else None
    c=[ap(45*np.cos(t),45*np.sin(t)) for t in np.arange(24)*2*np.pi/24];v=np.array([x for x in c if x is not None])
    target=ap(0,0);med=float(np.median(v)) if len(v) else None;sig=float(1.4826*np.median(abs(v-med))) if len(v) else None
    return {"target":target,"controls":c,"n":len(v),"robust_scale":sig,"descriptive_control_z":float((target-med)/sig) if target is not None and sig and len(v)>=12 else None}
def selftest():
    rng=np.random.default_rng(412701);yy,xx=np.mgrid[-30:31,-30:31].astype(float);sig=6.5/2.354820045
    a=rng.normal(0,.03,xx.shape)+4*np.exp(-((xx-4)**2+(yy+2)**2)/(2*sig**2))
    b=rng.normal(0,.03,xx.shape)+4*np.exp(-(xx**2+yy**2)/(2*sig**2))
    q=psf_fit(a,xx,yy,6.5);z=psf_fit(b,xx,yy,6.5)
    assert q["both_gain_positive"] and all(x["offset_arcsec"]==[4.,-2.] for x in q["folds"])
    assert all(x["offset_arcsec"]==[0.,0.] for x in z["folds"])
    write("selftest.json",{"shifted":q,"centered":z})
def gaia():
    ids=S["infrared"]["gaia_source_ids"];query="SELECT source_id,ra,dec,pmra,pmdec,ref_epoch,phot_g_mean_mag,bp_rp,parallax,ruwe FROM gaiadr3.gaia_source WHERE source_id IN ("+",".join(ids.values())+")"
    r=get("https://gea.esac.esa.int/tap-server/tap/sync",{"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":query})
    (O/"gaia.csv").write_bytes(r.content);rows={v["source_id"]:v for v in csv.DictReader(io.StringIO(r.text))}
    assert set(rows)==set(ids.values()),"not all Gaia sources"
    write("gaia.json",{"sha256":sha(r.content),"url":r.url,"rows":rows})
def infrared(label):
    ids=S["infrared"]["gaia_source_ids"];out={"label":label,"gaia_id":ids[label],"bands":{}}
    try:
        d=read("gaia.json")["rows"][ids[label]]
        ra,dec,pmra,pmdec,epoch=[num(d[x]) for x in ("ra","dec","pmra","pmdec","ref_epoch")]
        assert all(x is not None for x in (ra,dec,pmra,pmdec,epoch))
        c=SkyCoord(ra*u.deg,dec*u.deg);dt=2010.5-epoch
        c=c.spherical_offsets_by(pmra*dt*u.mas,pmdec*dt*u.mas)
        out.update(gaia=d,approx_wise_epoch=2010.5,approx_position_deg=[float(c.ra.deg),float(c.dec.deg)])
        SkyView.TIMEOUT=75
        for band,fwhm in (("WISE 3.4",6.1),("WISE 12",6.5)):
            v={"status":"HOLD"}
            try:
                ims=SkyView.get_images(position=c,survey=[band],radius=1.3*u.arcmin,pixels=160,cache=False,show_progress=False)
                if not ims:raise RuntimeError("no image")
                h=ims[0][0];data=np.asarray(h.data,float).squeeze();w=WCS(h.header).celestial;assert data.ndim==2
                y,x=np.indices(data.shape);world=w.pixel_to_world(x,y);dx,dy=c.spherical_offsets_to(world);xx=dx.arcsec;yy=dy.arcsec
                q=psf_fit(data,xx,yy,fwhm);ctl=controls(data,xx,yy)
                path=O/(label+("-W1.fits" if band=="WISE 3.4" else "-W3.fits"));h.writeto(path,overwrite=True)
                v={"status":"OK" if q is not None else "HOLD_INSUFFICIENT_PIXELS","fits_sha256":sha(path.read_bytes()),"shape":list(data.shape),"pixel_scale_arcsec":(proj_plane_pixel_scales(w)*3600).tolist(),"unit":h.header.get("BUNIT"),"psf_fit":q,"controls":ctl}
            except Exception as e:v={"status":"HOLD_SERVICE_OR_DATA","error":str(e)[:500]}
            out["bands"][band]=v
    except Exception as e:out["error"]=str(e)[:500]
    out["status"]="MEASURED_EXPLORATORY" if len(out["bands"])==2 and all(v["status"]=="OK" for v in out["bands"].values()) else "HOLD"
    write("ir-"+label+".json",out);print("R41U_IR",label,out["status"],flush=True)
def orbit():
    rows=[]
    for target in S["orbit"]["targets"]+S["orbit"]["natural_controls"]:
        v={"target":target,"role":"interstellar" if target in S["orbit"]["targets"] else "natural_control"}
        try:
            r=get("https://ssd-api.jpl.nasa.gov/sbdb.api",{"sstr":target,"full-prec":"true","cov":"mat","phys-par":"true","anc-data":"true"})
            d=r.json();assert "orbit" in d and "object" in d and "NASA/JPL" in d["signature"]["source"]
            o=d["orbit"];e={p["name"]:p for p in o["elements"]};model={p["name"]:p for p in o.get("model_pars",[]) or []}
            path=O/("jpl-"+str(len(rows)).zfill(2)+".json");path.write_bytes(r.content)
            v.update(status="MEASURED_CATALOG_ONLY",url=r.url,sha256=sha(r.content),designation=d["object"].get("fullname"),orbit_id=o.get("orbit_id"),epoch=o.get("epoch"),data_arc=o.get("data_arc"),n_obs_used=o.get("n_obs_used"),elements={k:e.get(k) for k in ("e","q","i","a")},model={k:model.get(k) for k in ("A1","A2","A3")},has_covariance=o.get("covariance") is not None,physical_parameters=d.get("phys_par",[]),interpretation="fitted catalog parameters, not independent observed residuals; missing is not zero")
        except Exception as e:v.update(status="HOLD",error=str(e)[:500])
        rows.append(v);print("R41U_ORBIT",target,v["status"],flush=True);time.sleep(.4)
    write("orbit_catalog.json",rows)
def toy():
    # Toy residual analogue of Rivals+Runs: natural smooth vs hinge vs periodic.
    rng=np.random.default_rng(412702);t=np.linspace(0,1,48);base=np.c_[np.ones(48),t,t*t]
    arms={"smooth":[base],"hinge":[np.c_[base,np.maximum(t-h,0)] for h in (.2,.3,.4)],"periodic":[np.c_[base,np.sin(2*np.pi*f*t),np.cos(2*np.pi*f*t)] for f in (1.,2.,3.)]}
    cases=[]
    for i in range(64):
        noise=rng.normal(0,.08,48);bg=.2*t+.5*t*t
        for name,y in (("noise",noise),("smooth",bg+noise),("impulse",bg+3*np.maximum(t-.3,0)+noise),("periodic",bg+.35*np.sin(4*np.pi*t)+noise)):
            report={}
            for split,train in (("future",np.arange(48)<24),("interleaved",np.arange(48)%2==0)):
                scores={}
                for arm,bank in arms.items():
                    fits=[]
                    for j,X in enumerate(bank):
                        c=np.linalg.lstsq(X[train],y[train],rcond=None)[0]
                        pred=X@c;fits.append((float(np.mean((pred[train]-y[train])**2)),j,float(np.mean((pred[~train]-y[~train])**2))))
                    fit=min(fits);scores[arm]={"training_loss":fit[0],"index":fit[1],"heldout_loss":fit[2]}
                report[split]=scores
            cases.append({"pair":i,"simulated_class":name,"splits":report})
    summary={cl:{split:{arm:float(np.median([x["splits"][split][arm]["heldout_loss"] for x in cases if x["simulated_class"]==cl])) for arm in arms} for split in ("future","interleaved")} for cl in ("noise","smooth","impulse","periodic")}
    write("toy_residual.json",{"scope":"synthetic scalar residual only; no orbital dynamics or real probe labels","cases":cases,"median_heldout_mse":summary})
def prepare():
    selftest()
    try:gaia()
    except Exception as e:write("gaia_hold.json",{"status":"HOLD","error":str(e)[:500]})
    toy()
def summary():
    ir=[read("ir-"+l+".json") for l in S["infrared"]["gaia_source_ids"]];orb=read("orbit_catalog.json");toy0=read("toy_residual.json")
    result={"status":"EXPLORATORY_ADAPTERS_ONLY","infrared_qualified":sum(x["status"]=="MEASURED_EXPLORATORY" for x in ir),"infrared_total":len(ir),"orbital_catalog_qualified":sum(x["status"]=="MEASURED_CATALOG_ONLY" for x in orb),"orbital_total":len(orb),"infrared":ir,"orbit_catalog":orb,"toy_median_heldout_mse":toy0["median_heldout_mse"],"radio_run":S["radio"]["run"],"radio_weights_reused_nonradio":False,"real_probe_detected":False,"real_astrometry_residual_tested":False,"scientific_adoption":False,"production_adoption":False,"protocol_sha256":sha(P.read_bytes()),"artifact_hashes":{p.name:sha(p.read_bytes()) for p in sorted(O.iterdir()) if p.is_file()}}
    write("result.json",result);print("R41U_SUMMARY",json.dumps({k:v for k,v in result.items() if k not in ("infrared","orbit_catalog","artifact_hashes")}),flush=True)
if __name__=="__main__":
    cmd=sys.argv[1]
    if cmd=="prepare":prepare()
    elif cmd=="ir":infrared(sys.argv[2])
    elif cmd=="orbit":orbit()
    elif cmd=="summary":summary()
    else:raise ValueError(cmd)
