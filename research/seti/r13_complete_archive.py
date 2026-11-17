#!/usr/bin/env python3
"""R13 full archive census; fixed follow-up models, 1 MHz morphology and independent metadata audit."""
import hashlib,json,pathlib,sys,collections,re,math
import numpy as np
import h5py,hdf5plugin
import r12_candidate_followup as prev
import r9_expanded_search as eng
import r4_chunk_replay as base
P=pathlib.Path("research/seti/r13_archive_protocol.json")
S=json.loads(P.read_text());ROOT=pathlib.Path("runtime/r13")
ROOT.mkdir(parents=True,exist_ok=True);eng.OUT=ROOT
F=prev.F;REV=S["reverse_control_mhz"];TARGET=S["target"]
def save(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+"\n")
def emit(k,v):print(k,json.dumps(v,ensure_ascii=False),flush=True)
def freeze():
    fs,q=eng.query(TARGET);fs.sort(key=lambda f:(f["day"],f["seconds"],f["url"]))
    assert len({f["name"] for f in fs})==len(fs),"duplicate archive filename"
    assert 0<len(fs)<=S["limits"]["max_files"],"archive count exceeds bounded census"
    external=[]
    for target in (TARGET,TARGET+"_S",TARGET+"_R"):
        for tel in ("ALL","GBT","MeerKAT","ATA"):
            params={"target":target,"file-type":"HDF5"}
            if tel!="ALL":params["telescope"]=tel
            row={"target":target,"query_telescope":tel}
            try:
                raw,url=base.fetch(eng.API,params,limit=16*base.MIB)
                d=json.loads(raw);rs=d if isinstance(d,list) else next((d[k] for k in ("data","files","results") if isinstance(d.get(k),list)),[])
                (ROOT/f"external-{target}-{tel}.json").write_bytes(raw)
                exact=[r for r in rs if isinstance(r,dict) and str(r.get("target","")).upper() in (TARGET,TARGET+"_S",TARGET+"_R")]
                counts=collections.Counter(str(r.get("telescope","unknown")) for r in exact)
                row.update(url=url,sha256=hashlib.sha256(raw).hexdigest(),rows=len(rs),exact_rows=len(exact),telescope_counts=dict(counts),non_parkes_products=[{k:r.get(k) for k in ("target","telescope","url","center_freq","file_type","mjd","size")} for r in exact if str(r.get("telescope","")).lower()!="parkes"][:100])
            except Exception as e:row["error"]=str(e)
            external.append(row)
    m={"schema":"janus.seti.r13.manifest.v1","protocol_sha256":hashlib.sha256(P.read_bytes()).hexdigest(),"query":q,"files":fs,"days":sorted({f["day"] for f in fs}),"position_counts":dict(collections.Counter(f["position"] for f in fs)),"external_queries":external,"selection":"complete metadata census, no intensity selection"}
    save(ROOT/"manifest.json",m)
    emit("R13_MANIFEST",{"files":len(fs),"days":m["days"],"position_counts":m["position_counts"],"query":q,"external_queries":[{k:v for k,v in x.items() if k!="non_parkes_products"} for x in external]})
def tile(h,lo):
    return {"channel0_mhz":h["f0"]+lo*h["df"],"df_hz":h["df"]*1e6,"ts":h["ts"]}
def valid(h,f):
    return eng.valid(h,f,TARGET) and h["bounds"][0]<F-.5 and h["bounds"][1]>F+.5
def injections(x,t,f,dr):
    n,m=x.shape;df=t["df_hz"];ts=t["ts"];scale=1.4826*np.median(np.abs(x-np.median(x,axis=1)[:,None]),axis=1)
    s=int(round(dr*(n-1)*ts/df));ds=np.rint(np.arange(n)*s/(n-1)).astype(int)
    results=[]
    for amp in (4,8):
        for shift in (-70,70):
            ft=f+shift*abs(df)/1e6;k=round((ft-t["channel0_mhz"])*1e6/df)
            if k+min(ds)<0 or k+max(ds)>=m:continue
            z=x.copy();z[np.arange(n),k+ds]+=amp*scale
            a=eng.matched(eng.normalized(z),t,ft,dr)
            results.append({"amplitude":amp,"offset_bins":shift,"recovered":eng.present(a)})
    return results
def morphology(x,t,out,key):
    y=eng.normalized(x);n,m=y.shape
    med=np.median(y,axis=0); persistence=np.mean(y>3,axis=0)
    ix=np.flatnonzero((med[1:-1]>=med[:-2])&(med[1:-1]>=med[2:])&(persistence[1:-1]>=.5)&(med[1:-1]>=3))+1
    # Greedy separated peaks, report exact pre-dedup number and retained count.
    keep=[]
    for i in sorted(ix,key=lambda k:(-float(med[k]),int(k))):
        if all(abs(int(i)-j)*abs(t["df_hz"])>=28 for j in keep):keep.append(int(i))
        if len(keep)==100:break
    freq=lambda i:t["channel0_mhz"]+i*t["df_hz"]/1e6
    peaks=[{"frequency_mhz":freq(i),"offset_from_candidate_hz":(freq(i)-F)*1e6,"median_row_excess":float(med[i]),"persistence":float(persistence[i])} for i in keep]
    path=out/f"{key}-wide.npz"
    np.savez_compressed(path,power=x.astype(np.float32),median_row_excess=med.astype(np.float32),persistence=persistence.astype(np.float32),channel0_mhz=t["channel0_mhz"],df_hz=t["df_hz"],ts=t["ts"])
    return {"shape":list(x.shape),"bandwidth_hz":m*abs(t["df_hz"]),"stationary_local_peaks_before_dedup":len(ix),"retained_peaks":len(keep),"cap_reached":len(keep)==100,"peaks":peaks,"artifact":path.name,"sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
def inspect(f,out):
    r={"name":f["name"],"url":f["url"],"day":f["day"],"position":f["position"],"status":"HOLD","tests":{}}
    key=hashlib.sha256(f["url"].encode()).hexdigest()[:24]
    try:
        h=eng.hdr(f);r["header"]=h
        if not valid(h,f):raise ValueError("ineligible header or coverage")
        models=prev.models(h,f["day"]==int(prev.T));wide=f["day"]==int(prev.T)
        budget=S["limits"]["wide_file_mib"] if wide else S["limits"]["per_file_mib"]
        with base.RangeFile(f["url"],budget*base.MIB) as remote:
            with h5py.File(remote,"r",rdcc_nbytes=128*base.MIB,rdcc_nslots=10007) as z:
                d=z["data"];a=d.attrs
                if float(a["tstart"])!=h["mjd"] or float(a["foff"])!=h["df"] or float(a["fch1"])!=h["f0"] or list(d.shape)!=h["shape"]:raise ValueError("header mismatch")
                windows=[]
                for k,fr,dr in models:
                    w=base.window(d,h["f0"],h["df"],fr,S["narrow_half_hz"])
                    if w is None:raise ValueError("predicted window missing")
                    windows.append((k,fr,dr,w))
                ww=base.window(d,h["f0"],h["df"],F,S["wide_half_hz"]) if wide else None
                plan=base.preflight(d,[v[3] for v in windows]+([ww] if ww else []),d.shape[0])
                if plan and plan["compressed_bytes"]+remote.total+8*base.MIB>remote.budget:raise ValueError("preflight budget")
                arrays={}
                for k,fr,dr,(lo,hi,_) in windows:
                    if (lo,hi) not in arrays:
                        x=np.asarray(d[:,0,lo:hi],dtype=np.float32)
                        if not np.isfinite(x).all():raise ValueError("nonfinite window")
                        name=f"{key}-{lo}.npy";np.save(out/name,x,allow_pickle=False)
                        arrays[lo,hi]=(x,name)
                    x,name=arrays[lo,hi];t=tile(h,lo)
                    result=prev.score(x,t,fr,dr)
                    if k in ("fixed_original_drift","constant_barycentric"):result["injections"]=injections(x,t,fr,dr)
                    r["tests"][k]={**result,"tile":t,"artifact":name,"slice_sha256":hashlib.sha256(x.tobytes()).hexdigest()}
                if ww:
                    lo,hi,_=ww;x=np.asarray(d[:,0,lo:hi],dtype=np.float32)
                    if not np.isfinite(x).all():raise ValueError("nonfinite wide window")
                    r["morphology"]=morphology(x,tile(h,lo),out,key)
                r["bytes_read"]=remote.total;r["status"]="COMPLETE"
    except Exception as e:r["error"]=str(e)
    emit("R13_SCAN",{**r,"morphology":{k:v for k,v in r["morphology"].items() if k!="peaks"} if "morphology" in r else None})
    return r
def shard(i):
    m=json.loads((ROOT/"manifest.json").read_text());assert m["protocol_sha256"]==hashlib.sha256(P.read_bytes()).hexdigest()
    out=ROOT/f"shard-{i}";out.mkdir(exist_ok=True)
    fs=m["files"][i::S["limits"]["shards"]];rows=[]
    for f in fs:
        rows.append(inspect(f,out));save(out/"result.json",{"shard":i,"manifest_sha256":hashlib.sha256((ROOT/"manifest.json").read_bytes()).hexdigest(),"rows":rows})
    emit("R13_SHARD",{"shard":i,"planned":len(fs),"completed":sum(r["status"]=="COMPLETE" for r in rows)})
def aggregate():
    m=json.loads((ROOT/"manifest.json").read_text());digest=hashlib.sha256((ROOT/"manifest.json").read_bytes()).hexdigest()
    rows=[];missing=[]
    for i in range(S["limits"]["shards"]):
        p=ROOT/f"shard-{i}"/"result.json"
        if not p.exists():missing.append(i);continue
        x=json.loads(p.read_text());assert x["manifest_sha256"]==digest and x["shard"]==i
        assert [r["name"] for r in x["rows"]]==[f["name"] for f in m["files"][i::8][:len(x["rows"])]]
        rows.extend(x["rows"])
    done=[r for r in rows if r["status"]=="COMPLETE"]
    def strong(r,ks):return any(r["tests"].get(k,{}).get("strong",False) for k in ks)
    fixed=("fixed_original_drift","fixed_zero_drift","continuous")
    same=[r for r in done if r["day"]==int(prev.T)];other=[r for r in done if r["day"]!=int(prev.T)]
    tallies={}
    for name,group in (("same_day",same),("other_days",other)):
        for pos in ("S","R"):
            g=[r for r in group if r["position"]==pos]
            tallies[name+"_"+pos]={"scans":len(g),"fixed_strong":sum(strong(r,fixed) for r in g),"barycentric_strong":sum(strong(r,("constant_barycentric",)) for r in g),"reverse_line_strong":sum(strong(r,("reverse_control_fixed",)) for r in g)}
    inj={}
    for k in ("fixed_original_drift","constant_barycentric"):
        inj[k]={}
        for amp in (4,8):
            vals=[q for r in done for q in r["tests"].get(k,{}).get("injections",[]) if q["amplitude"]==amp]
            inj[k][str(amp)]={"passed":sum(v["recovered"] for v in vals),"trials":len(vals)}
    wide=[r for r in done if "morphology" in r]
    pairs=[]
    for a in wide:
        for b in wide:
            if a["position"]!="S" or b["position"]!="R":continue
            ap=a["morphology"]["peaks"];bp=b["morphology"]["peaks"]
            pairs.append({"on":a["name"],"off":b["name"],"on_peaks":len(ap),"off_peaks":len(bp),"on_peaks_matched_within_3_channels":sum(any(abs(p["frequency_mhz"]-q["frequency_mhz"])*1e6<=3*abs(a["header"]["df"]*1e6) for q in bp) for p in ap)})
    off_hit=any(r["position"]=="R" and strong(r,fixed) for r in done)
    s={"schema":"janus.seti.r13.summary.v1","manifest_sha256":digest,"planned":len(m["files"]),"reported":len(rows),"completed":len(done),"days_completed":len({r["day"] for r in done}),"days_planned":len(m["days"]),"missing_shards":missing,"tallies":tallies,"wide_scans":len(wide),"morphology_pairs":pairs,"targeted_relative_injection":inj,"status":"HOLD_INCOMPLETE_ARCHIVE" if len(done)!=len(m["files"]) else ("RFI_CONCERN_REVIEW" if off_hit else "HOLD_INDEPENDENT_CONFIRMATION"),"independently_confirmed":0,"confirmed_technosignatures":0,"limitations":["Local fixed model matches and stationary peak similarity do not establish common emitter. Barycentric hypothesis assumes constant emitted frequency.","Matched raw-power injections scaled by robust row noise are not blind full-window calibration, nor Jy/EIRP or beam sensitivity.","Public archive coverage only, not every physical observation; no independent new telescope time.","External metadata query is not actual independent coverage verification; no observations commissioned."]}
    save(ROOT/"summary.json",{**s,"rows":rows,"external_queries":m["external_queries"]});emit("R13_SUMMARY",s)
if __name__=="__main__":
    if sys.argv[1]=="freeze":freeze()
    elif sys.argv[1]=="shard":shard(int(sys.argv[2]))
    elif sys.argv[1]=="aggregate":aggregate()
