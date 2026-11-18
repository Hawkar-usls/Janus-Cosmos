#!/usr/bin/env python3
"""Natural W3OH astronomy reference-frame audit, no emitter classification."""
import collections,hashlib,io,json,math,os,pathlib,sys,urllib.parse,zipfile,re
import requests
import numpy as np
P=pathlib.Path("research/seti/r37f2a_reference_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r37f2a");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def official(u):
    p=urllib.parse.urlparse(u)
    return p.scheme in ("http","https") and (p.hostname=="bldata.berkeley.edu" or bool(re.fullmatch(r"blpd\d+\.ssl\.berkeley\.edu",p.hostname or ""))) and not p.username and not p.password
def parent():
    p=S["parent"];r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['artifact']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90)
    r.raise_for_status();b=r.content;assert len(b)<16*1024**2 and sha(b)==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        a=[v for v in z.infolist() if pathlib.PurePosixPath(v.filename).name=="frozen_manifest.json"]
        assert len(a)==1 and a[0].file_size<8*1024**2
        d=json.loads(z.read(a[0]))
    assert d["frozen_before_power_read"] and len(d["pairs"])==1
    return d["pairs"][0]
def access():
    from r37f1de_archive import fil_header,fil_grid
    r=requests.get(S["catalogue"],params={"target":"W3OH","telescope":"GBT","file-type":"HDF5"},timeout=60)
    r.raise_for_status();b=r.content;assert len(b)<=S["metadata_limit"]
    j=json.loads(b);rr=j if isinstance(j,list) else next(j[k] for k in ("data","files","results") if isinstance(j.get(k),list))
    out=[];files={}
    for v in rr:
        if not isinstance(v,dict):continue
        u=v.get("url");t=str(v.get("target","")).strip().upper()
        if not isinstance(u,str) or not official(u) or t not in ("DIAG_W3OH","DIAG_W3OH_OFF","W3OH","W3OH_OFF"):continue
        try:m=float(v["mjd"]);assert int(m) in S["days"]
        except (ValueError,TypeError,AssertionError,KeyError):continue
        ext=pathlib.PurePosixPath(urllib.parse.urlparse(u).path).suffix.lower()
        if ext not in (".fil",".h5",".raw"):continue
        files[u]={"url":u,"mjd":m,"target":t,"role":"OFF" if t.endswith("_OFF") else "ON","ext":ext}
    def probe(v):
        u=v["url"];schemes=[u];p=urllib.parse.urlparse(u)
        if p.scheme=="http":schemes.append(urllib.parse.urlunparse(p._replace(scheme="https")))
        attempts=[]
        for url in schemes:
            a={"scheme":urllib.parse.urlparse(url).scheme}
            try:
                cap=S["header_probe_bytes"]
                with requests.get(url,headers={"Range":f"bytes=0-{cap-1}","Accept-Encoding":"identity"},stream=True,timeout=(10,25),allow_redirects=False) as r:
                    a["http_status"]=r.status_code
                    if r.status_code!=206:raise IOError("exact range unavailable")
                    m=re.fullmatch(r"bytes 0-(\d+)/(\d+)",r.headers.get("Content-Range",""));assert m and int(m[1])<cap
                    buf=bytearray()
                    for c in r.iter_content(8192):buf.extend(c);assert len(buf)<=cap
                    assert len(buf)==int(m[1])+1
                    h,n=fil_header(bytes(buf));g=fil_grid(h,n,int(m[2]))
                    assert g["source"].strip().upper()==v["target"] and abs(g["mjd"]-v["mjd"])<=3/86400
                    assert math.isfinite(g["df"]) and math.isfinite(g["ts"])
                    a.update(status="HEADER_PARSED",header_sha256=sha(bytes(buf[:n])),channel_hz=abs(g["df"])*1e6,integration_s=g["ts"],rows=g["shape"][0],covers_oh=all(g["bounds"][0]<=f<=g["bounds"][1] for f in S["source"]["rest_mhz"]),authenticated_transport=a["scheme"]=="https",frame_flags={k:h.get(k) for k in ("barycentric","pulsarcentric")})
            except Exception as e:a.update(status="HOLD",error=type(e).__name__+": "+str(e)[:200])
            attempts.append(a)
            if a["status"]=="HEADER_PARSED" and a["authenticated_transport"]:break
        return {**v,"attempts":attempts}
    for day in S["days"]:
        fs=[v for v in files.values() if int(v["mjd"])==day]
        sel=[]
        # Explicit cap: first and last exact ON/OFF observations, prefer fine product.
        for role in ("ON","OFF"):
            q=sorted([v for v in fs if v["role"]==role and v["ext"]==".fil"],key=lambda v:(".gpuspec.0000." not in v["url"],v["mjd"],v["url"]))
            seen=set()
            for v in q:
                if v["mjd"] in seen:continue
                sel.append(v);seen.add(v["mjd"])
                if len(seen)>=S["header_probe_max_per_day"]//2:break
        results=[probe(v) for v in sel]
        out.append({"day":day,"catalogue_unique":len(fs),"formats":dict(collections.Counter(v["ext"] for v in fs)),"roles":dict(collections.Counter(v["role"] for v in fs)),"probes":results,"authenticated_headers":sum(any(a.get("authenticated_transport") and a["status"]=="HEADER_PARSED" for a in v["attempts"]) for v in results),"provisional_http_headers":sum(any(a.get("authenticated_transport") is False and a["status"]=="HEADER_PARSED" for a in v["attempts"]) for v in results)})
    d={"status":"ARCHIVE_ACCESS_AUDIT_ONLY","catalogue_sha256":sha(b),"days":out,"independent_power_validated":False}
    save("access.json",d);print("R37F2A_ACCESS_SUMMARY",json.dumps(d),flush=True)
def frame():
    import astropy,astropy_iers_data
    import astropy.units as u
    from astropy.coordinates import EarthLocation,SkyCoord,SpectralCoord,solar_system_ephemeris
    from astropy.time import Time
    from astropy.utils import iers
    import h5py,hdf5plugin
    import r4_chunk_replay as base
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error"
    solar_system_ephemeris.set("builtin")
    q=parent();hs=q["headers"];cfg=S["frame"];src=S["source"];m=S["measure"]
    loc=EarthLocation.from_geodetic(cfg["longitude_deg"]*u.deg,cfg["latitude_deg"]*u.deg,cfg["approx_height_m"]*u.m)
    target=SkyCoord(ra=src["ra_deg"]*u.deg,dec=src["dec_deg"]*u.deg,distance=2*u.kpc,pm_ra_cosdec=0*u.mas/u.yr,pm_dec=0*u.mas/u.yr,radial_velocity=0*u.km/u.s,frame="icrs")
    def factor(mjd):
        t=Time(mjd,format="mjd",scale="utc");obs=loc.get_itrs(obstime=t)
        a=SpectralCoord(1*u.GHz,observer=obs,target=target,doppler_convention="radio",doppler_rest=1*u.GHz)
        b=a.with_observer_stationary_relative_to("lsrk")
        back=b.with_observer_stationary_relative_to(obs)
        k=float(b.to_value(u.GHz));err=abs(float(back.to_value(u.GHz))-1)
        assert 0.999<k<1.001 and err<1e-8
        return {"factor":k,"roundtrip_relative_error":err}
    factors=[]
    for h in hs:
        start=h["mjd"];dur=h["ts"]*h["shape"][0]/86400
        f=[factor(start+x*dur) for x in (0,.5,1)]
        factors.append(f)
    # Correct each pointing's time separately, interpolate only normalized OFF spectrum
    # onto ON velocity centers; do not claim raw native-grid equality after correction.
    windows=[];rd=[]
    for k,(url,h) in enumerate(zip((q["on"],q["off"]),hs)):
        assert base.archive_url(url)
        with base.RangeFile(url,m["read_budget_bytes"]) as remote:
            with h5py.File(remote,"r",rdcc_nbytes=128*base.MIB) as f:
                d=f["data"];a=d.attrs;assert list(d.shape)==h["shape"]
                assert abs(float(a["fch1"])-h["f0"])<1e-10 and abs(float(a["foff"])-h["df"])<1e-12 and abs(float(a["tstart"])-h["mjd"])<1e-9
                flags={x:int(a[x]) if x in a else None for x in ("barycentric","pulsarcentric")}
                assert flags["barycentric"] in (None,0) and flags["pulsarcentric"] in (None,0)
                rows=[];kk=factors[k][1]["factor"]
                for rest,(vlo,vhi) in zip(src["rest_mhz"],src["published_lsr_ranges_km_s"]):
                    # Broad fixed surrounding baseline of +/-20 km/s beyond published endpoints.
                    va,vb=vlo-20,vhi+20
                    fm=sorted([rest*(1-va/299792.458)/kk,rest*(1-vb/299792.458)/kk])
                    lo,hi=sorted([(x-h["f0"])/h["df"] for x in fm])
                    lo=max(0,math.ceil(lo));hi=min(d.shape[-1],math.floor(hi)+1)
                    hi=lo+(hi-lo)//m["block_bins"]*m["block_bins"];assert hi-lo>m["block_bins"]*10
                    x=np.asarray(d[:,0,lo:hi],dtype=np.float64);assert np.isfinite(x).all() and (x>=0).all()
                    med=np.median(x,axis=1);assert (med>0).all()
                    y=(x/med[:,None]).reshape(len(x),-1,m["block_bins"]).mean(2)
                    ch=np.arange(lo,hi,dtype=float).reshape(-1,m["block_bins"]).mean(1)
                    v=299792.458*(1-(h["f0"]+ch*h["df"])*kk/rest)
                    order=np.argsort(v);rows.append((v[order],y[:,order]))
                windows.append(rows);rd.append({"flags":flags,"bytes_read":remote.total})
    lines=[]
    for i,(rest,(lo,hi)) in enumerate(zip(src["rest_mhz"],src["published_lsr_ranges_km_s"])):
        v,on=windows[0][i];w,off=windows[1][i]
        keep=(v>=w.min())&(v<=w.max());v=v[keep];on=on[:,keep]
        control=np.interp(v,w,off.mean(0));diff=on.mean(0)-control
        med=float(np.median(diff));sigma=1.4826*float(np.median(abs(diff-med)));assert sigma>0
        tol=src["tolerance_km_s"];ok=(v>=lo-tol)&(v<=hi+tol);assert ok.any()
        j=int(np.where(ok)[0][np.argmax(diff[ok])]);z=float((diff[j]-med)/sigma)
        pos=int(sum(on[:,j]>control[j]));adj=bool((j>0 and diff[j-1]>med+3*sigma) or (j+1<len(v) and diff[j+1]>med+3*sigma))
        lines.append({"rest_mhz":rest,"published_lsr_interval":[lo,hi],"peak_lsrk_km_s":float(v[j]),"descriptive_z":z,"relative_excess":float(diff[j]),"positive_rows":pos,"rows":len(on),"adjacent":adj,"criterion_met":bool(z>=m["threshold_z"] and diff[j]>=m["min_excess"] and pos==len(on) and adj)})
    iers_table=iers.earth_orientation_table.get()
    # Record deterministic installed data-file hashes, not an internet ephemeris guess.
    data_files=[]
    for name in ("IERS_A_FILE","IERS_B_FILE","IERS_LEAP_SECOND_FILE"):
        path=getattr(astropy_iers_data,name,None)
        if path and pathlib.Path(path).is_file():data_files.append({"name":name,"sha256":sha(pathlib.Path(path).read_bytes())})
    d={"status":"CONDITIONAL_LSRK_LINE_AUDIT","astropy":astropy.__version__,"iers_data_version":astropy_iers_data.__version__,"iers_table_class":type(iers_table).__name__,"iers_file_hashes":data_files,"ephemeris":"builtin","day":q["day"],"factors_start_mid_end":factors,"reads":rd,"lines":lines,"both_criteria_met":all(v["criterion_met"] for v in lines),"archive_frame_independently_verified":False,"independent_epoch_confirmed":False,"production_adoption":False}
    save("frame.json",d);print("R37F2A_FRAME_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":
    {"access":access,"frame":frame}[sys.argv[1]]()
