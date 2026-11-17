#!/usr/bin/env python3
"""R10 frozen candidate follow-up. No candidate promotion or broadcast."""
import hashlib,json,pathlib,math
import numpy as np
import h5py,hdf5plugin
from astropy import units as u
from astropy.constants import c
from astropy.coordinates import SkyCoord,EarthLocation
from astropy.time import Time
from astropy.utils import iers
import r4_chunk_replay as base
import r9_expanded_search as r9

S=json.loads(pathlib.Path("research/seti/r10_frozen_candidate.json").read_text())
A=S["anchor"]; F=A["frequency_mhz"]; D=A["drift_hz_s"]; T=A["mjd"]
OUT=pathlib.Path("runtime/r10");OUT.mkdir(parents=True,exist_ok=True)
iers.conf.auto_download=False
loc=EarthLocation.from_geodetic(S["barycentric"]["location_lon_deg"]*u.deg,S["barycentric"]["location_lat_deg"]*u.deg,S["barycentric"]["height_m"]*u.m)
coord=SkyCoord(ra=A["ra_hour"]*u.hourangle,dec=A["dec_deg"]*u.deg,frame="icrs")
def factor(mjd):
    return 1+float((coord.radial_velocity_correction(obstime=Time(mjd,format="mjd",scale="utc"),location=loc)/c).decompose().value)
FB=F/factor(T)
def models(h,same=False):
    t=h["mjd"];span=(h["shape"][0]-1)*h["ts"];dt=(t-T)*86400
    p=FB*factor(t); drift=(FB*factor(t+span/86400)-p)*1e6/span
    m=[("fixed_original_drift",F,D),("fixed_zero_drift",F,0.)]
    if same and abs(dt)<3600:m.append(("continuous",F+D*dt/1e6,D))
    m.append(("constant_barycentric",p,drift))
    return m
def score(x,tile,f,dr):
    a=r9.matched(r9.normalized(x),tile,f,dr)
    return {"predicted_mhz":f,"predicted_drift_hz_s":dr,"best":a,"strong":r9.present(a)}
def load_prior():
    root=pathlib.Path("prior-r9"); data=json.loads((root/"r9-result.json").read_text())
    assert data["targets_completed"]==4
    cand=[a for a in data["candidates"] if a["target"]==S["target"] and abs(a["frequency_mhz"]-F)<1e-12]
    assert len(cand)==1 and cand[0]["status"].startswith("PRELIMINARY")
    rows=[]
    for q in data["results"]:
        for o in q["observations"]:
            tile=next(t for t in o["tiles"] if t["center_mhz"]==3000)
            x=np.load(root/tile["artifact"],allow_pickle=False)
            assert list(x.shape)==tile["shape"] and hashlib.sha256(x.tobytes()).hexdigest()==tile["slice_sha256"]
            assert np.isfinite(x).all()
            m=models(o["header"],q["target"]==S["target"])
            if q["target"]!=S["target"]:m=m[:2]
            row={"target":q["target"],"role":o["role"],"name":o["name"],"mjd":o["header"]["mjd"],"tests":{k:score(x,tile,f,d) for k,f,d in m}}
            rows.append(row);print("R10_PRIOR",json.dumps(row),flush=True)
    return rows,cand[0]
def freeze():
    fs,q=r9.query(S["target"])
    fs.sort(key=lambda f:(f["day"],f["seconds"],f["url"]))
    orig=set(S["source_cadence"]);selected=[]
    # First subsequent R and S after source ON2; same day, within 30 minutes.
    for p in ("R","S"):
        pool=[f for f in fs if f["day"]==int(T) and f["position"]==p and 13929<f["seconds"]<=13929+1800 and f["name"] not in orig]
        if pool:selected.append({**pool[0],"role":"EXTRA_SAME_DAY_"+p})
    cad=[]
    for i in range(len(fs)-2):
        a,b,d=fs[i:i+3]
        if [a["position"],b["position"],d["position"]]!=["S","R","S"] or len({a["day"],b["day"],d["day"]})!=1 or a["day"]==int(T):continue
        if not 0<b["seconds"]-a["seconds"]<=600 or not 0<d["seconds"]-b["seconds"]<=600:continue
        cad.append((a,b,d))
    for side in (-1,1):
        pool=sorted((x for x in cad if (x[0]["day"]-int(T))*side>0),key=lambda x:(abs(x[0]["day"]-int(T)),x[0]["seconds"],x[0]["url"]))
        if pool:
            for f,role in zip(pool[0],("ON1","OFF","ON2")):selected.append({**f,"role":("EARLIER_" if side<0 else "LATER_")+role})
    result={"query":q,"selected":selected,"metadata_only":True}
    (OUT/"manifest.json").write_text(json.dumps(result,indent=2)+"\n")
    print("R10_MANIFEST",json.dumps({"query":q,"selected":[{"name":f["name"],"role":f["role"],"day":f["day"],"seconds":f["seconds"]} for f in selected]}),flush=True)
    return selected
def inspect(f):
    row={"name":f["name"],"role":f["role"],"url":f["url"],"status":"HOLD","tests":{}}
    try:
        h=r9.hdr(f);row["header"]=h
        if not r9.valid(h,f,S["target"]):raise ValueError("ineligible header")
        m=models(h,f["day"]==int(T))
        with base.RangeFile(f["url"],S["limits"]["budget_mib_per_file"]*base.MIB) as remote:
            with h5py.File(remote,"r",rdcc_nbytes=128*base.MIB) as z:
                d=z["data"];a=d.attrs
                if float(a["tstart"])!=h["mjd"] or float(a["foff"])!=h["df"] or list(d.shape)!=h["shape"]:raise ValueError("header drift")
                for k,freq,dr in m:
                    w=base.window(d,h["f0"],h["df"],freq,S["window_half_hz"])
                    if w is None:row["tests"][k]={"status":"HOLD_OUTSIDE_COVERAGE"};continue
                    lo,hi,_=w
                    x=np.asarray(d[:,0,lo:hi],dtype=np.float32)
                    if not np.isfinite(x).all() or x.shape[0]!=h["shape"][0]:raise ValueError("bad spectral slice")
                    tile={"channel0_mhz":h["f0"]+lo*h["df"],"df_hz":h["df"]*1e6,"ts":h["ts"]}
                    name=hashlib.sha256((f["url"]+"|"+str(freq)).encode()).hexdigest()[:24]+".npy"
                    np.save(OUT/name,x,allow_pickle=False)
                    row["tests"][k]={**score(x,tile,freq,dr),"artifact":name,"sha256":hashlib.sha256(x.tobytes()).hexdigest(),"tile":tile}
                row["bytes_read"]=remote.total;row["status"]="SCREEN_COMPLETE"
    except Exception as e:row["error"]=str(e)
    print("R10_EXTRA",json.dumps(row),flush=True)
    return row
def main():
    selected=freeze();prior,candidate=load_prior();extra=[]
    for f in selected:
        extra.append(inspect(f))
        (OUT/"partial.json").write_text(json.dumps(extra,indent=2)+"\n")
    cross=[r for r in prior if r["target"]!=S["target"] and any(t.get("strong") for t in r["tests"].values())]
    off=[r for r in extra if (r["role"].endswith("OFF") or r["role"].endswith("_R")) and any(t.get("strong") for t in r["tests"].values())]
    summary={"schema":"janus.seti.r10.followup.v1","prior_scans":len(prior),"extra_selected":len(selected),"extra_completed":sum(r["status"]=="SCREEN_COMPLETE" for r in extra),"cross_target_strong_tracks":len(cross),"extra_off_strong_tracks":len(off),"barycentric_reference_mhz":FB,"status":"RFI_CONCERN_REVIEW" if cross or off else "HOLD_INDEPENDENT_CONFIRMATION","confirmed_technosignatures":0,"independently_confirmed":0,"limitations":["R9 OFF already had score 12.83, persistence .25; absence is not established.","Fixed corridors and existing selection are descriptive; no absolute flux or blind false-alarm calibration.","Constant barycentric frequency is only one restrictive emission model; no correction for unknown source acceleration.","Archive S/R positions are not new independent observations; related frequency tracks do not establish common emitter."]}
    (OUT/"result.json").write_text(json.dumps({**summary,"prior":prior,"extra":extra,"r9_candidate":candidate},indent=2)+"\n")
    print("R10_SUMMARY",json.dumps(summary),flush=True)
if __name__=="__main__":main()
