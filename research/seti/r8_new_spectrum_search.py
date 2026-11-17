#!/usr/bin/env python3
"""R8 prospectively specified, metadata-selected search. No catalog hit input.
Only fixed narrow tiles; descriptive scores and bounded follow-up.
"""
import hashlib,math,json,pathlib,re,sys,collections
from urllib.parse import urlparse,unquote
import h5py,hdf5plugin,numpy as np
import r4_chunk_replay as base

SPEC=json.loads(pathlib.Path("research/seti/r8_preregistered_search.json").read_text())
OUT=pathlib.Path("runtime/r8");OUT.mkdir(parents=True,exist_ok=True)
API="https://seti.berkeley.edu/opendata/api/query-files"

def allowed(u):
    p=urlparse(u)
    return p.scheme in ("http","https") and (p.hostname=="bldata.berkeley.edu" or bool(re.fullmatch(r"blpd\d+\.ssl\.berkeley\.edu",p.hostname or "")))

def query(target):
    q={"target":target,"telescope":"Parkes","file-type":"HDF5"}
    raw,url=base.fetch(API,q,limit=16*base.MIB)
    d=json.loads(raw); rows=d if isinstance(d,list) else next((d[k] for k in ("data","files","results") if isinstance(d.get(k),list)),[])
    (OUT/f"index-{target}.json").write_bytes(raw)
    fs={}
    for row in rows:
        if not isinstance(row,dict):continue
        u=row.get("url","")
        if not isinstance(u,str) or not allowed(u):continue
        m=re.fullmatch(r"Parkes_(\d{5})_(\d{5})_"+target+r"_([SR])_fine\.h5",unquote(base.basename(u)),re.I)
        if not m:continue
        f={"url":u,"day":int(m[1]),"seconds":int(m[2]),"position":m[3].upper(),"name":unquote(base.basename(u)),"metadata":row}
        fs[u]=f
    return list(fs.values()),{"url":url,"rows":len(rows),"fine":len(fs),"sha256":hashlib.sha256(raw).hexdigest()}

def hdr(f):
    with base.RangeFile(f["url"],8*base.MIB) as remote:
        with h5py.File(remote,"r") as h:
            d=h["data"];a=d.attrs
            f0=float(a["fch1"]);df=float(a["foff"]);ts=float(a["tsamp"]);mjd=float(a["tstart"]);name=a.get("source_name","")
            if isinstance(name,bytes):name=name.decode(errors="replace")
            bounds=sorted([f0,f0+(d.shape[-1]-1)*df])
            return {"source":str(name),"shape":list(d.shape),"f0":f0,"df":df,"ts":ts,"mjd":mjd,"bounds":bounds,"chunks":list(d.chunks) if d.chunks else None,"src_raj":str(a.get("src_raj","")),"src_dej":str(a.get("src_dej",""))}
def valid(h,f,target):
    return (h["source"].upper()==target+"_"+f["position"] and len(h["shape"])==3 and 12<=h["shape"][0]<=64 and 0<abs(h["df"])*1e6<=10 and 0<h["ts"]<=60 and
            abs((h["mjd"]-f["day"])*86400-f["seconds"])<=3 and all(h["bounds"][0]<x-.021 and h["bounds"][1]>x+.021 for x in SPEC["window_centers_mhz"]))

def freeze():
    chosen=[];attempts=[]
    for target in SPEC["targets_in_fixed_order"]:
        try:
            fs,q=query(target);attempt={"target":target,**q,"cadence_attempts":[]};attempts.append(attempt)
            fs.sort(key=lambda f:(f["day"],f["seconds"],f["url"]))
            cad=[]
            for i in range(len(fs)-2):
                a,b,c=fs[i:i+3]
                if [a["position"],b["position"],c["position"]]!=["S","R","S"]:continue
                if len({a["day"],b["day"],c["day"]})!=1:continue
                if not 0<b["seconds"]-a["seconds"]<=600 or not 0<c["seconds"]-b["seconds"]<=600:continue
                cad.append((a,b,c))
            cad.sort(key=lambda c:(c[0]["day"]<58575,c[0]["day"],c[0]["seconds"],c[0]["url"]))
            # At most five metadata-only attempts per target, no power data read.
            for a,b,c in cad[:5]:
                try:
                    hs=[hdr(f) for f in (a,b,c)]
                    if not all(valid(h,f,target) for h,f in zip(hs,(a,b,c))):raise ValueError("ineligible source header")
                    if len({(h["shape"][0],round(h["ts"],7),round(h["df"],12)) for h in hs})!=1:raise ValueError("resolution mismatch")
                    if hs[0]["src_dej"]==hs[1]["src_dej"] and hs[0]["src_raj"]==hs[1]["src_raj"]:raise ValueError("no header pointing difference")
                    chosen.append({"target":target,"cohort":"post_2019_march" if a["day"]>=58575 else "archival_reanalysis","files":[{**f,"header":h,"role":role} for f,h,role in zip((a,b,c),hs,("ON1","OFF","ON2"))]})
                    attempt["cadence_attempts"].append({"day":a["day"],"seconds":a["seconds"],"status":"SELECTED"});break
                except Exception as e:attempt["cadence_attempts"].append({"day":a["day"],"seconds":a["seconds"],"error":str(e)})
            if len(chosen)==SPEC["limits"]["targets"]:break
        except Exception as e:attempts.append({"target":target,"error":str(e)})
    result={"schema":"janus.seti.r8.frozen_archive_manifest.v1","spec_sha256":hashlib.sha256(pathlib.Path("research/seti/r8_preregistered_search.json").read_bytes()).hexdigest(),"selected":chosen,"attempts":attempts,"no_intensity_read_during_selection":True}
    (OUT/"manifest.json").write_text(json.dumps(result,indent=2)+"\n")
    print("R8_MANIFEST",json.dumps({"selected":[{"target":c["target"],"cohort":c["cohort"],"files":[{"role":f["role"],"name":f["name"],"header":f["header"]} for f in c["files"]]} for c in chosen],"attempts":attempts}),flush=True)
    if not chosen:print("HOLD_NO_ELIGIBLE_CADENCE")

def get_tiles(f):
    out=[];h0=f["header"]
    with base.RangeFile(f["url"],SPEC["limits"]["spectrum_budget_per_file_mib"]*base.MIB) as remote:
        with h5py.File(remote,"r",rdcc_nbytes=128*base.MIB,rdcc_nslots=10007) as h:
            d=h["data"];a=d.attrs;df=float(a["foff"]);f0=float(a["fch1"]);ts=float(a["tsamp"])
            if df!=h0["df"] or f0!=h0["f0"] or ts!=h0["ts"] or float(a["tstart"])!=h0["mjd"] or list(d.shape)!=h0["shape"]:raise ValueError("frozen header mismatch")
            w=[base.window(d,f0,df,x,SPEC["window_half_hz"]) for x in SPEC["window_centers_mhz"]]
            if any(v is None for v in w):raise ValueError("missing frequency tile")
            plan=base.preflight(d,w,d.shape[0])
            if plan and plan["compressed_bytes"]+remote.total+8*base.MIB>remote.budget:raise ValueError("full scan exceeds budget")
            for fc,(lo,hi,center) in zip(SPEC["window_centers_mhz"],w):
                x=np.asarray(d[:,0,lo:hi],dtype=np.float32)
                if x.shape[0]!=d.shape[0] or not np.isfinite(x).all():raise ValueError("bad tile")
                name=hashlib.sha256((f["url"]+"|"+str(fc)).encode()).hexdigest()[:24]+".npy"
                np.save(OUT/name,x,allow_pickle=False)
                out.append({"center_mhz":fc,"channel0_mhz":f0+lo*df,"df_hz":df*1e6,"ts":ts,"shape":list(x.shape),"artifact":name,"slice_sha256":hashlib.sha256(x.tobytes()).hexdigest()})
            print("R8_READ",json.dumps({"name":f["name"],"role":f["role"],"bytes":remote.total,"tiles":out}),flush=True)
    return out

def normalized(x):
    m=np.median(x,axis=1);s=1.4826*np.median(np.abs(x-m[:,None]),axis=1)
    if (s<=0).any():raise ValueError("zero scale rows")
    return (x-m[:,None])/s[:,None]

def scores(y,shift):
    n,m=y.shape
    ds=np.rint(np.arange(n)*shift/(n-1)).astype(int)
    lo=max(0,-int(ds.min()));hi=min(m,m-int(ds.max()))
    if hi<=lo:return None
    v=np.zeros(hi-lo,dtype=np.float32)
    for i,k in enumerate(ds):v+=y[i,lo+k:hi+k]
    return v/np.sqrt(n),lo,ds

def top(y,df,ts,nkeep,threshold):
    n,m=y.shape;lim=int(math.floor(4*(n-1)*ts/abs(df)));pool=[]
    # All integer endpoint shifts; keep local maxima per drift with guards.
    for shift in range(-lim,lim+1):
        r=scores(y,shift)
        if r is None:continue
        v,lo,ds=r
        if len(v)<5:continue
        idx=np.flatnonzero((v[1:-1]>=v[:-2])&(v[1:-1]>=v[2:])&(v[1:-1]>=threshold))+1
        if len(idx):
            idx=idx[np.argsort(v[idx])[-max(8,2*nkeep):]]
            for k in idx:
                c=lo+int(k);vals=y[np.arange(n),c+ds]
                pers=float(np.mean(vals>3))
                if pers>=SPEC["anchor_min_row_excess_fraction"]:
                    pool.append({"bin":c,"shift":shift,"score":float(v[k]),"persistence":pers,"drift_hz_s":shift*df/((n-1)*ts)})
    pool.sort(key=lambda a:(-a["score"],a["bin"],a["shift"]))
    keep=[]
    for a in pool:
        # Collapse nearby endpoint trajectories, including parallel drift templates.
        if any(abs(a["bin"]-b["bin"])*abs(df)<28 and abs((a["bin"]+a["shift"])-(b["bin"]+b["shift"]))*abs(df)<28 for b in keep):continue
        keep.append(a)
        if len(keep)==nkeep:break
    return keep,{"drift_trials":2*lim+1,"channels":m,"upper_bound_templates":(2*lim+1)*m,"eligible_local_peak_pool":len(pool)}

def matched(y,tile,freq,drift,prior_amp=None):
    df=tile["df_hz"];ts=tile["ts"];n,m=y.shape
    c=(freq-tile["channel0_mhz"])*1e6/df
    s0=int(round(drift*(n-1)*ts/df))
    choices=[]
    for s in range(s0-1,s0+2):
        r=scores(y,s)
        if r is None:continue
        v,lo,ds=r
        for k in range(round(c)-10,round(c)+11):
            if k<lo or k>=lo+len(v):continue
            vals=y[np.arange(n),k+ds]
            choices.append({"score":float(v[k-lo]),"persistence":float(np.mean(vals>3)),"frequency_mhz":tile["channel0_mhz"]+k*df/1e6,"drift_hz_s":s*df/((n-1)*ts)})
    if not choices:return None
    return max(choices,key=lambda a:a["score"])

def present(a):return a is not None and a["score"]>=SPEC["off_and_repeat_threshold_score"] and a["persistence"]>=SPEC["anchor_min_row_excess_fraction"]

def selftest():
    rng=np.random.default_rng(801);df=-2.8;ts=18.;n=16;m=1024
    x=rng.normal(0,1,(n,m));shift=124;c=330
    for i,k in enumerate(np.rint(np.arange(n)*shift/(n-1)).astype(int)):x[i,c+k]+=14
    a,_=top(normalized(x),df,ts,5,12)
    assert any(abs(b["bin"]-c)<=2 and abs(b["shift"]-shift)<=2 for b in a),a
    noise=0
    for _ in range(8):
        a,_=top(normalized(rng.normal(0,1,(n,m))),df,ts,5,12);noise+=len(a)
    assert noise==0,noise
    print("R8_SELFTEST",json.dumps({"synthetic_injection_recovered":1,"injected":1,"gaussian_noise_tiles":8,"noise_candidates":noise,"note":"Synthetic smoke, not empirical global false-alarm calibration."}),flush=True)

def replay():
    selftest()
    manifest=json.loads((OUT/"manifest.json").read_text())
    assert manifest["spec_sha256"]==hashlib.sha256(pathlib.Path("research/seti/r8_preregistered_search.json").read_bytes()).hexdigest()
    results=[];all_candidates=[];off_count=0
    for q in manifest["selected"]:
        result={"target":q["target"],"cohort":q["cohort"],"status":"HOLD","observations":[],"tiles":[]};results.append(result)
        try:
            for f in q["files"]:result["observations"].append({**f,"tiles":get_tiles(f)})
            for i,fc in enumerate(SPEC["window_centers_mhz"]):
                obs=result["observations"];tiles=[o["tiles"][i] for o in obs]
                arrays=[]
                for tile in tiles:
                    x=np.load(OUT/tile["artifact"],allow_pickle=False)
                    assert hashlib.sha256(x.tobytes()).hexdigest()==tile["slice_sha256"]
                    arrays.append(normalized(x))
                y,off,rep=arrays;a,b,c=tiles
                hits,counts=top(y,a["df_hz"],a["ts"],SPEC["top_candidates_per_tile"],SPEC["anchor_min_score"])
                oh,oc=top(off,b["df_hz"],b["ts"],SPEC["top_candidates_per_tile"],SPEC["anchor_min_score"])
                off_count+=len(oh)
                result["tiles"].append({"center_mhz":fc,"on_search":counts,"off_search":oc,"off_candidates":oh,"on_candidate_count":len(hits)})
                for h in hits:
                    freq=a["channel0_mhz"]+h["bin"]*a["df_hz"]/1e6;dr=h["drift_hz_s"]
                    candidate={"target":q["target"],"center_mhz":fc,"frequency_mhz":freq,**h,"comparisons":[],"status":"HOLD_NO_DIRECTIONAL_REPEAT"}
                    for o,t,z in zip(obs[1:],tiles[1:],arrays[1:]):
                        dt=(o["header"]["mjd"]-obs[0]["header"]["mjd"])*86400
                        pred=freq+dr*dt/1e6
                        cv=matched(z,t,pred,dr);fv=matched(z,t,freq,dr)
                        # Normalized injections near, but not at, predicted frequency.
                        n,m=y.shape;ds=np.rint(np.arange(n)*h["shift"]/(n-1)).astype(int)
                        amp=float(np.median(y[np.arange(n),h["bin"]+ds]))
                        inject=[]
                        for offset in (-70,70):
                            ftest=pred+offset*abs(t["df_hz"])/1e6
                            k=int(round((ftest-t["channel0_mhz"])*1e6/t["df_hz"]))
                            s=int(round(dr*(len(z)-1)*t["ts"]/t["df_hz"]))
                            dsi=np.rint(np.arange(len(z))*s/(len(z)-1)).astype(int)
                            if k+dsi.min()<0 or k+dsi.max()>=z.shape[1]:continue
                            zz=z.copy();zz[np.arange(len(z)),k+dsi]+=amp
                            inject.append(present(matched(zz,t,ftest,dr)))
                        candidate["comparisons"].append({"role":o["role"],"dt_seconds":dt,"continuous":cv,"fixed":fv,"injections_passed":sum(inject),"injections":len(inject)})
                    offc,repc=candidate["comparisons"]
                    if present(offc["continuous"]) or present(offc["fixed"]):candidate["status"]="CONTROL_TRACK_PRESENT_RFI_CONCERN"
                    elif offc["continuous"] is not None and offc["fixed"] is not None and offc["injections"]==2 and offc["injections_passed"]==2 and present(repc["continuous"]):
                        candidate["status"]="PRELIMINARY_ON_OFF_ON_HOLD_CALIBRATION_AND_INDEPENDENT_REPEAT"
                    all_candidates.append(candidate)
                    print("R8_CANDIDATE",json.dumps(candidate),flush=True)
                print("R8_TILE",json.dumps({"target":q["target"],**result["tiles"][-1]}),flush=True)
            result["status"]="FIXED_TILE_SEARCH_COMPLETE"
        except Exception as e:result["error"]=str(e);print("R8_HOLD",q["target"],str(e),flush=True)
        (OUT/"r8-result.json").write_text(json.dumps({"results":results,"candidates":all_candidates},indent=2)+"\n")
    counts=dict(collections.Counter(a["status"] for a in all_candidates))
    final={"schema":"janus.seti.r8.fixed_tile_search.v1","frozen_manifest_sha256":hashlib.sha256((OUT/"manifest.json").read_bytes()).hexdigest(),"targets_selected":len(manifest["selected"]),"targets_completed":sum(q["status"]=="FIXED_TILE_SEARCH_COMPLETE" for q in results),"on_tiles_searched":sum(len(q["tiles"]) for q in results),"off_tiles_searched":sum(len(q["tiles"]) for q in results),"on_candidates":len(all_candidates),"off_candidates":off_count,"candidate_status_counts":counts,"independently_confirmed":0,"confirmed_technosignatures":0,"limitations":["Small fixed frequency tiles not full survey; top-five cap, descriptive search maxima.","Metadata-selected archived observations may have been analyzed before.","No absolute flux calibration or verified global false alarm rate.","No independent cross-epoch/barycentric repeat. Public HTTP archive partial reads, not full file digest checks."]}
    (OUT/"r8-result.json").write_text(json.dumps({**final,"results":results,"candidates":all_candidates},indent=2)+"\n")
    print("R8_SUMMARY",json.dumps(final),flush=True)

if __name__=="__main__":
    {"freeze":freeze,"scan":replay}[sys.argv[1]]()
