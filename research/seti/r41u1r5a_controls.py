"""R5A: preselected matched star controls on native atlas; not independent exposures."""
import csv,io,json,os,pathlib,sys,subprocess,zipfile
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import r41u1_wise_bounded as prior
import r41u_multimodal as old
O=pathlib.Path("runtime/r41u1r5a");O.mkdir(parents=True,exist_ok=True)
P=pathlib.Path("research/seti/r41u1r5a_protocol.json");L=json.loads(P.read_text())
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def read(n):return json.loads((O/n).read_text())
def pos(r):
    ra,de,pmra,pmde,ep=[old.num(r[k]) for k in ("ra","dec","pmra","pmdec","ref_epoch")]
    if None in (ra,de,pmra,pmde,ep):return None
    return SkyCoord(ra*u.deg,de*u.deg).spherical_offsets_by(pmra*(2010.5-ep)*u.mas,pmde*(2010.5-ep)*u.mas)
def prepare():
    p=L["source"];h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    root="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions"
    md=json.loads(prior.get(f"{root}/artifacts/{p['artifact']}",headers=h))
    assert md["workflow_run"]["id"]==p["run"] and md["digest"]=="sha256:"+p["archive_sha256"]
    raw=prior.get(f"{root}/artifacts/{p['artifact']}/zip",headers=h,cap=64*1024**2)
    assert old.sha(raw)==p["archive_sha256"]
    z=zipfile.ZipFile(io.BytesIO(raw));assert len(z.namelist())==len(set(z.namelist())) and sum(i.file_size for i in z.infolist())<128*1024**2
    res=z.read("result.json");assert old.sha(res)==p["result_sha256"]
    r=json.loads(res);assert r["qualified"]==7 and r["bands_ok"]==14
    for n,d in r["hashes"].items():assert old.sha(z.read(n))==d
    for n in ["gaia.json"]+[f"ir-{k}.json" for k in old.S["infrared"]["gaia_source_ids"]]:
        (O/n).write_bytes(z.read(n))
    save("prepare.json",{"source":p,"verified_files":len(r["hashes"]),"protocol_hash":old.sha(P.read_bytes())})
def select(label):
    out={"label":label,"controls":[]}
    try:
        ids=old.S["infrared"]["gaia_source_ids"]; t=read("gaia.json")["rows"][ids[label]]
        c=SkyCoord(float(t["ra"])*u.deg,float(t["dec"])*u.deg)
        g,col=[old.num(t[k]) for k in ("phot_g_mean_mag","bp_rp")]
        assert g is not None and col is not None
        q=("SELECT TOP 50000 source_id,ra,dec,pmra,pmdec,ref_epoch,phot_g_mean_mag,bp_rp,parallax,ruwe FROM gaiadr3.gaia_source "
           f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{c.ra.deg:.10f},{c.dec.deg:.10f},0.12)) ORDER BY source_id")
        raw=prior.get("https://gea.esac.esa.int/tap-server/tap/sync",params={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":q},cap=16*1024**2)
        (O/f"{label}-cone.csv").write_bytes(raw)
        rows=list(csv.DictReader(io.StringIO(raw.decode())));assert 0<len(rows)<50000 and "source_id" in rows[0]
        assert ids[label] in {x["source_id"] for x in rows}
        tp=pos(t);assert tp is not None
        eligible=[]
        for r in rows:
            if r["source_id"] in ids.values():continue
            cp=pos(r);a,b=[old.num(r[k]) for k in ("phot_g_mean_mag","bp_rp")]
            if cp is None or a is None or b is None or cp.separation(tp).arcsec<90:continue
            if abs(a-g)>1.5 or abs(b-col)>0.5:continue
            eligible.append((abs(a-g)/1.5+abs(b-col)/0.5,r["source_id"],r,cp))
        chosen=[]
        for score,_,r,cp in sorted(eligible,key=lambda x:(x[0],x[1])):
            if any(cp.separation(p).arcsec<45 for p in chosen):continue
            chosen.append(cp);out["controls"].append({"row":r,"match_score":float(score),"pos":[float(cp.ra.deg),float(cp.dec.deg)]})
            if len(chosen)==3:break
        out.update(status="SELECTED" if chosen else "HOLD_NO_MATCH",target=t,target_position=[float(tp.ra.deg),float(tp.dec.deg)],cone_count=len(rows),eligible=len(eligible),catalog_sha256=old.sha(raw),query=q)
    except Exception as e:out.update(status="HOLD_SELECTION",error=str(e)[:500])
    save(f"select-{label}.json",out)
def freeze():
    rows=[read(f"select-{l}.json") for l in old.S["infrared"]["gaia_source_ids"]]
    save("selection_lock.json",{"protocol_hash":old.sha(P.read_bytes()),"selections":rows,"hashes":{r["label"]:old.sha((O/f"select-{r['label']}.json").read_bytes()) for r in rows}})
def near(c,rows,sid):
    # Gaia neighbors at catalogue epoch (incomplete for optically dark sources).
    pairs=[]
    for r in rows:
        if r["source_id"]==sid:continue
        cp=SkyCoord(float(r["ra"])*u.deg,float(r["dec"])*u.deg)
        sep=float(c.separation(cp).arcsec)
        if sep<=18:pairs.append({"source_id":r["source_id"],"sep_arcsec_2016":sep,"g":old.num(r["phot_g_mean_mag"])})
    return sorted(pairs,key=lambda x:x["sep_arcsec_2016"])
def measure(label):
    lock=read("selection_lock.json");sel=read(f"select-{label}.json")
    assert old.sha((O/f"select-{label}.json").read_bytes())==lock["hashes"][label]
    out={"label":label,"controls":[],"status":"HOLD"}
    if sel["status"]!="SELECTED":save(f"measure-{label}.json",out);return
    rows=list(csv.DictReader(io.StringIO((O/f"{label}-cone.csv").read_text())))
    assert old.sha((O/f"{label}-cone.csv").read_bytes())==sel["catalog_sha256"]
    target=read(f"ir-{label}.json")
    t=sel["target"];out["target_neighbors"]=near(SkyCoord(float(t["ra"])*u.deg,float(t["dec"])*u.deg),rows,t["source_id"])
    for j,k in enumerate(sel["controls"]):
        r=k["row"];c=SkyCoord(k["pos"][0]*u.deg,k["pos"][1]*u.deg)
        v={"source_id":r["source_id"],"match_score":k["match_score"],"neighbors":near(SkyCoord(float(r["ra"])*u.deg,float(r["dec"])*u.deg),rows,r["source_id"]),"bands":{}}
        for band,fwhm in (("W1",6.1),("W3",6.5)):
            try:
                url=target["bands"][band]["url"]
                from urllib.parse import urlparse
                assert urlparse(url).hostname=="irsa.ipac.caltech.edu" and url.endswith("-int-3.fits")
                params={"center":f"{c.ra.deg:.10f},{c.dec.deg:.10f}deg","size":"200arcsec","gzip":"false"}
                raw=prior.get(url,params=params,cap=2*1024**2);name=f"{label}-{j}-{band}.fits";(O/name).write_bytes(raw)
                with fits.open(io.BytesIO(raw),memmap=False) as hd:
                    data=np.asarray(hd[0].data,float).squeeze();w=WCS(hd[0].header).celestial
                    assert data.ndim==2 and min(data.shape)>=115 and max(data.shape)<=200
                    assert np.allclose(proj_plane_pixel_scales(w)*3600,target["bands"][band]["pixel_scale_arcsec"],rtol=1e-3,atol=1e-5)
                    cut=Cutout2D(data,c,(160*u.arcsec,160*u.arcsec),wcs=w,mode="strict")
                    y,x=np.indices(cut.shape);world=cut.wcs.pixel_to_world(x,y);rx,ry=cut.wcs.world_to_pixel(world)
                    assert np.max(abs(rx-x))<1e-5 and np.max(abs(ry-y))<1e-5
                    xx,yy=prior.offset_icrs(c,world);assert np.min(np.hypot(xx,yy))<2
                    q=old.psf_fit(cut.data,xx,yy,fwhm);ctl=old.controls(cut.data,xx,yy)
                    assert q is not None
                    v["bands"][band]={"status":"OK","sha256":old.sha(raw),"url":url,"params":params,"psf_fit":q,"controls":ctl,"boundary":any(f["offset_norm_arcsec"]>=6 for f in q["folds"])}
            except Exception as e:v["bands"][band]={"status":"HOLD","error":str(e)[:400]}
        out["controls"].append(v);save(f"measure-{label}.json",out)
    out["status"]="MEASURED" if any(all(v["bands"].get(b,{}).get("status")=="OK" for b in ("W1","W3")) for v in out["controls"]) else "HOLD"
    save(f"measure-{label}.json",out)
def bounded(cmd,label):
    try:subprocess.run([sys.executable,__file__,cmd,label],check=True,timeout=75 if cmd=="select" else 190)
    except (subprocess.TimeoutExpired,subprocess.CalledProcessError) as e:
        n=f"{cmd}-{label}.json";r=read(n) if (O/n).exists() else {"label":label,"controls":[]};r.update(status="HOLD_TIMEOUT_OR_CHILD",error=type(e).__name__);save(n,r)
def summary():
    rows=[read(f"measure-{l}.json") for l in old.S["infrared"]["gaia_source_ids"]];lock=read("selection_lock.json")
    pairs=[v for r in rows for v in r["controls"] if all(v["bands"].get(b,{}).get("status")=="OK" for b in ("W1","W3"))]
    def d(v):return {b:{"fold_offsets":[f["offset_norm_arcsec"] for f in v["bands"][b]["psf_fit"]["folds"]],"both_gain":v["bands"][b]["psf_fit"]["both_gain_positive"],"flux_proxy":v["bands"][b]["controls"]["target"]} for b in ("W1","W3")}
    out={"status":"PASS_MATCHED_CONTROL_INTAKE__HOLD_INDEPENDENT_EXPOSURES" if len(pairs)>0 else "HOLD_NO_MEASURED_CONTROLS","selected":sum(len(r["controls"]) for r in lock["selections"]),"measured":len(pairs),"groups_measured":sum(r["status"]=="MEASURED" for r in rows),"target_diagnostics":[d(read(f"ir-{l}.json")) for l in old.S["infrared"]["gaia_source_ids"]],"control_diagnostics":[d(v) for v in pairs],"rows":rows,"scope":"Gaussian fits on coadds, not independent exposures. No calibrated significance, no origin claim or production adoption.","hashes":{p.name:old.sha(p.read_bytes()) for p in O.iterdir() if p.is_file()}}
    save("result.json",out);print("R5A",out["status"],out["selected"],out["measured"],out["groups_measured"],flush=True)
if __name__=="__main__":
    cmd=sys.argv[1]
    if cmd=="prepare":prepare()
    elif cmd=="select":select(sys.argv[2])
    elif cmd=="freeze":freeze()
    elif cmd=="measure":measure(sys.argv[2])
    elif cmd in ("bounded-select","bounded-measure"):bounded(cmd.split("-")[1],sys.argv[2])
    elif cmd=="summary":summary()
    else:raise ValueError(cmd)
