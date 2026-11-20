"""Bounded official native AllWISE retrieval; prior timeout retained."""
import sys,os,io,zipfile,pathlib,json,subprocess,time,urllib.parse
import requests,numpy as np
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord
import astropy.units as u
import r41u_multimodal as old
O=pathlib.Path("runtime/r41u1");O.mkdir(parents=True,exist_ok=True)
P=pathlib.Path("research/seti/r41u1_wise_bounded_protocol.json");L=json.loads(P.read_text())
def save(k,d):(O/k).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def sha(x):return old.sha(x)
def get(url,params=None,headers=None,cap=32*1024**2):
    with requests.get(url,params=params,headers=headers,timeout=(10,25),stream=True) as r:
        r.raise_for_status();data=bytearray()
        for block in r.iter_content(65536):
            data.extend(block)
            if len(data)>cap:raise RuntimeError("response exceeds limit")
        return bytes(data)
def sia_table(raw):
    tab=parse_single_table(io.BytesIO(raw)).to_table(use_names_over_ids=True)
    required={"sia_bp_id","sia_url","coadd_id","sia_ra","sia_dec","sia_scale","sia_fmt"}
    assert required.issubset(tab.colnames),f"missing SIA named columns: {required-set(tab.colnames)}"
    return tab
def prepare():
    p=L["source"];h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"};root="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions"
    md=json.loads(get(f"{root}/artifacts/{p['artifact']}",headers=h));assert md["workflow_run"]["id"]==p["run"] and md["digest"]=="sha256:"+p["archive_sha256"]
    raw=get(f"{root}/artifacts/{p['artifact']}/zip",headers=h,cap=64*1024**2);assert sha(raw)==p["archive_sha256"]
    z=zipfile.ZipFile(io.BytesIO(raw));names=z.namelist();assert len(names)==len(set(names)) and sum(i.file_size for i in z.infolist())<128*1024**2
    for k in ["gaia.json","gaia.csv","orbit_catalog.json","selftest.json","toy_residual.json"]+[f"jpl-{i:02d}.json" for i in range(11)]:
        matches=[s for s in names if pathlib.PurePosixPath(s).name==k];assert len(matches)==1,k
        (O/k).write_bytes(z.read(matches[0]))
    g=json.loads((O/"gaia.json").read_text());assert sha((O/"gaia.csv").read_bytes())==g["sha256"] and set(g["rows"])==set(old.S["infrared"]["gaia_source_ids"].values())
    oc=json.loads((O/"orbit_catalog.json").read_text());assert len(oc)==11
    for i,v in enumerate(oc):assert v["status"]=="MEASURED_CATALOG_ONLY" and sha((O/f"jpl-{i:02d}.json").read_bytes())==v["sha256"]
    prev=L["amendment_r2"];meta=json.loads(get(f"{root}/artifacts/{prev['artifact']}",headers=h));assert meta["workflow_run"]["id"]==prev["previous_run"] and meta["digest"]=="sha256:"+prev["archive_sha256"]
    prior=get(f"{root}/artifacts/{prev['artifact']}/zip",headers=h,cap=64*1024**2);assert sha(prior)==prev["archive_sha256"]
    zz=zipfile.ZipFile(io.BytesIO(prior));rr=json.loads(zz.read("result.json"));assert rr["qualified"]==0 and rr["gaia_preserved"]==7 and rr["orbit_preserved"]==11
    for label in old.S["infrared"]["gaia_source_ids"]:
        name=f"{label}-sia.xml";raw=zz.read(name);assert sha(raw)==rr["hashes"][name];tab=sia_table(raw);assert {"W1","W3"}.issubset({str(x).strip() for x in tab["sia_bp_id"]});(O/f"{label}-prior-sia.xml").write_bytes(raw)
    save("prepare.json",{"gaia":7,"orbit":11,"source":p,"prior_sia_replayed":7,"protocol_sha256":sha(P.read_bytes())})
    fix=L["amendment_r3"];md=json.loads(get(f"{root}/artifacts/{fix['artifact']}",headers=h));assert md["workflow_run"]["id"]==fix["previous_run"] and md["digest"]=="sha256:"+fix["archive_sha256"]
    raw=get(f"{root}/artifacts/{fix['artifact']}/zip",headers=h,cap=64*1024**2);assert sha(raw)==fix["archive_sha256"]
    zz=zipfile.ZipFile(io.BytesIO(raw));assert len(zz.namelist())==len(set(zz.namelist())) and sum(i.file_size for i in zz.infolist())<128*1024**2
    result=zz.read("result.json");assert sha(result)==fix["result_sha256"];rr=json.loads(result)
    for name,digest in rr["hashes"].items():assert sha(zz.read(name))==digest
    restored=[]
    for label in old.S["infrared"]["gaia_source_ids"]:
        name=f"{label}-sia.xml";raw=zz.read(name);assert sha(raw)==rr["hashes"][name];sia_table(raw);(O/name).write_bytes(raw)
        for band in ("W1","W3"):
            name=f"{label}-{band}-native-ibe.fits"
            if name in rr["hashes"]:
                (O/name).write_bytes(zz.read(name));restored.append(name)
    assert len(restored)==6
    save("r3-restore.json",{"source":fix,"verified_hashes":len(rr["hashes"]),"restored_fits":restored})
def retry_503(url,params):
    for attempt in (1,2):
        try:return get(url,params=params,cap=2*1024**2),attempt
        except requests.HTTPError as e:
            if attempt==2 or e.response is None or e.response.status_code!=503:raise
            time.sleep(2)
def measure(label):
    sid=old.S["infrared"]["gaia_source_ids"][label];d=json.loads((O/"gaia.json").read_text())["rows"][sid]
    ra,de,pmra,pmdec,ep=[old.num(d[k]) for k in ("ra","dec","pmra","pmdec","ref_epoch")]
    assert all(x is not None for x in (ra,de,pmra,pmdec,ep))
    c=SkyCoord(ra*u.deg,de*u.deg).spherical_offsets_by(pmra*(2010.5-ep)*u.mas,pmdec*(2010.5-ep)*u.mas)
    out={"label":label,"gaia_id":sid,"position_deg":[float(c.ra.deg),float(c.dec.deg)],"epoch_approx":2010.5,"bands":{}}
    try:
        raw=(O/f"{label}-sia.xml").read_bytes()
        table=sia_table(raw)
        def val(row,k):
            x=row[k];return (x.decode() if isinstance(x,bytes) else str(x)).strip()
        grouped={}
        for row in table:
            band=val(row,"sia_bp_id");url=val(row,"sia_url");id=val(row,"coadd_id")
            if band in ("W1","W3") and urllib.parse.urlparse(url).hostname=="irsa.ipac.caltech.edu" and url.endswith("-int-3.fits"):
                grouped.setdefault(id,{})[band]=(url,float(row["sia_ra"]),float(row["sia_dec"]),np.abs(np.asarray(row["sia_scale"],float))*3600)
        choices=[(float(c.separation(SkyCoord(v["W1"][1]*u.deg,v["W1"][2]*u.deg)).deg),k,v) for k,v in grouped.items() if set(v)=={"W1","W3"}]
        assert choices,"no common W1/W3 tile";_,cid,bands=min(choices,key=lambda x:(x[0],x[1]));out["coadd_id"]=cid
        for k,fwhm in (("W1",6.1),("W3",6.5)):
            try:
                url=bands[k][0];params={"center":f"{c.ra.deg:.10f},{c.dec.deg:.10f}deg","size":"200arcsec","gzip":"false"};path=O/f"{label}-{k}-native-ibe.fits";cached=path.exists()
                if cached:raw=path.read_bytes();attempts=0
                else:raw,attempts=retry_503(url,params);path.write_bytes(raw)
                with fits.open(io.BytesIO(raw),memmap=False) as hd:
                    h=hd[0];data=np.asarray(h.data,float).squeeze();w=WCS(h.header).celestial;assert data.ndim==2
                    assert min(data.shape)>=115 and max(data.shape)<=200 and np.allclose(proj_plane_pixel_scales(w)*3600,bands[k][3],rtol=1e-3,atol=1e-5),"native scale/size mismatch"
                    cut=Cutout2D(data,c,(160*u.arcsec,160*u.arcsec),wcs=w,mode="strict")
                    y,x=np.indices(cut.data.shape);world=cut.wcs.pixel_to_world(x,y)
                    rx,ry=cut.wcs.world_to_pixel(world);assert np.max(abs(rx-x))<1e-5 and np.max(abs(ry-y))<1e-5,"WCS roundtrip mismatch"
                    same=world.transform_to(c.frame);dx,dy=c.spherical_offsets_to(same);xx=dx.arcsec;yy=dy.arcsec
                    assert np.isfinite(xx).all() and np.isfinite(yy).all() and np.min(np.hypot(xx,yy))<2.0,"target not within native pixel grid"
                    head=cut.wcs.to_header();head["BUNIT"]=str(h.header.get("BUNIT","UNKNOWN"));fits.PrimaryHDU(cut.data,header=head).writeto(O/f"{label}-{k}-cutout.fits",overwrite=True)
                    q=old.psf_fit(cut.data,xx,yy,fwhm);ctl=old.controls(cut.data,xx,yy)
                    out["bands"][k]={"status":"OK" if q is not None else "HOLD_INSUFFICIENT_PIXELS","url":url,"query":params,"source_kind":"native IBE server cutout, not full atlas","archive_reused":cached,"http_attempts":attempts,"wcs_native_frame":world.frame.name,"offset_frame":same.frame.name,"source_shape":list(data.shape),"source_sha256":sha(raw),"cutout_sha256":sha((O/f"{label}-{k}-cutout.fits").read_bytes()),"shape":list(cut.shape),"unit":h.header.get("BUNIT"),"pixel_scale_arcsec":(proj_plane_pixel_scales(cut.wcs)*3600).tolist(),"psf_fit":q,"controls":ctl}
            except Exception as e:out["bands"][k]={"status":"HOLD","error":str(e)[:400]}
            save("ir-"+label+".json",out)
    except Exception as e:out["error"]=str(e)[:400]
    out["status"]="MEASURED_EXPLORATORY" if set(out["bands"])=={"W1","W3"} and all(v["status"]=="OK" for v in out["bands"].values()) else "HOLD"
    save("ir-"+label+".json",out)
def bounded(label):
    t=time.monotonic();p=O/f"ir-{label}.json"
    try:
        subprocess.run([sys.executable,__file__,"measure",label],timeout=150,check=True)
    except (subprocess.TimeoutExpired,subprocess.CalledProcessError) as e:
        v=json.loads(p.read_text()) if p.exists() else {"label":label,"bands":{}}
        v.update(status="HOLD_CHILD_TIMEOUT_OR_ERROR",error=type(e).__name__);save(p.name,v)
    v=json.loads(p.read_text());v["wall_s"]=time.monotonic()-t;save(p.name,v);print("R41U1_IR",label,v["status"],flush=True)
def summary():
    rows=[json.loads((O/f"ir-{k}.json").read_text()) for k in old.S["infrared"]["gaia_source_ids"]]
    d={"status":"PASS_NATIVE_ATLAS_INTAKE__REVIEW_ONLY" if all(v["status"]=="MEASURED_EXPLORATORY" for v in rows) else "PARTIAL_OR_HOLD_NATIVE_ATLAS","qualified":sum(v["status"]=="MEASURED_EXPLORATORY" for v in rows),"total":7,"bands_ok":sum(b["status"]=="OK" for v in rows for b in v["bands"].values()),"rows":rows,"gaia_preserved":7,"orbit_preserved":11,"hashes":{p.name:sha(p.read_bytes()) for p in O.iterdir() if p.is_file()},"external_validation":False,"production_adoption":False}
    save("result.json",d);print("R41U1_SUMMARY",json.dumps({k:v for k,v in d.items() if k not in ("rows","hashes")}),flush=True)
if __name__=="__main__":
    cmd=sys.argv[1]
    if cmd=="prepare":prepare()
    elif cmd=="measure":measure(sys.argv[2])
    elif cmd=="bounded":bounded(sys.argv[2])
    elif cmd=="summary":summary()
