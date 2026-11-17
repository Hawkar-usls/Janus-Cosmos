#!/usr/bin/env python3
"""Bounded SHA-verified replay, not a fresh independent observation."""
import hashlib,json,pathlib,zipfile,os,sys,time,struct,io,math
import numpy as np
import requests
import r9_expanded_search as eng
import r4_chunk_replay as base
import r15c_r16c_spectral as sig
import r11_control_calibration as r11
P=json.loads(pathlib.Path("research/seti/r15d_r16d_r17d_protocol.json").read_text())
lane=sys.argv[1];OUT=pathlib.Path("runtime/"+lane);OUT.mkdir(parents=True,exist_ok=True)
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def emit(k,x):print(k,json.dumps(x),flush=True)
def archive(key):
    a=P["artifacts"][key];path=OUT/(key+".zip");h=hashlib.sha256();n=0
    with requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},stream=True,timeout=(15,120)) as r:
        r.raise_for_status()
        with path.open("wb") as f:
            for b in r.iter_content(1024*1024):
                n+=len(b);assert n<=500*1024**2;h.update(b);f.write(b)
    assert h.hexdigest()==a["digest"]
    z=zipfile.ZipFile(path);assert sum(i.file_size for i in z.infolist())<=1200*1024**2
    names={}
    for i in z.infolist():
        if i.is_dir():continue
        b=pathlib.PurePosixPath(i.filename).name
        assert b not in names;names[b]=i.filename
    return z,names,path
def get(z,n,k):return json.loads(z.read(n[k]))
def tile(z,n,t,k=None):
    x=np.load(io.BytesIO(z.read(n[k or t["artifact"]])),allow_pickle=False)
    assert list(x.shape)==t["shape"] and hashlib.sha256(x.tobytes()).hexdigest()==t["slice_sha256"]
    return x
def compare_reader(f,x,t):
    from blimpy.io import sigproc
    from blimpy import Waterfall
    h=f["header"];p=OUT/"independent_anchor_crop.fil";lo=round((t["channel0_mhz"]-h["fch1"])/h["foff"])
    with base.RangeFile(f["url"],8*base.MIB) as r:
        r.seek(0);raw=r.read(h["offset"]);p.write_bytes(raw)
        external=sigproc.read_header(str(p))
        for k in ("nchans","nifs","nbits","fch1","foff","tsamp","tstart"):assert float(external[k])==float(h[k]),k
        assert sigproc.len_header(str(p))==h["offset"]
        sigproc.fix_header(str(p),"nchans",x.shape[1])
        sigproc.fix_header(str(p),"fch1",t["channel0_mhz"])
        with p.open("ab") as w:
            for j in range(x.shape[0]):
                r.seek(h["offset"]+4*(j*h["nchans"]+lo));b=r.read(x.shape[1]*4)
                assert len(b)==x.shape[1]*4
                assert np.array_equal(np.asarray(struct.unpack("<"+str(x.shape[1])+"f",b),np.float32),x[j])
                w.write(b)
    wf=Waterfall(str(p));y=np.asarray(wf.data)[:,0,:]
    assert y.shape==x.shape
    equal=bool(np.array_equal(y,x));reverse=bool(np.array_equal(y[:,::-1],x))
    assert equal or reverse
    return {"status":"PASS","header_fields":7,"payload_samples":int(x.size),"blimpy_equal":equal,"blimpy_reversed":reverse}
def search(x,t,sham=None):
    s=P["r15d"];y=eng.normalized(x)
    if sham is not None:y=y[np.random.default_rng(sham).permutation(len(y))]
    nr,nc=y.shape;df=t["df_hz"];ts=t["ts"];best=[]
    for width in sorted(set(min(nr,w) for w in s["windows"])):
        ans=None
        for start in sorted(set(list(range(0,nr-width+1,max(1,width//2)))+[nr-width])):
            tt=np.arange(width)*ts
            for d in s["drift_grid_hz_s"]:
                for a in s["drift_change_grid_hz_s"]:
                    # d is initial rate of segment; a is change over full scan.
                    shifts=np.rint((d*tt+0.5*a*tt**2/max((nr-1)*ts,1))/df).astype(int)
                    lo=max(0,-int(shifts.min()));hi=min(nc,nc-int(shifts.max()))
                    if hi<=lo:continue
                    vals=sum(y[start+i,lo+k:hi+k] for i,k in enumerate(shifts))
                    med=np.median(vals);sd=1.4826*np.median(np.abs(vals-med))
                    if sd<=0:continue
                    scores=(vals-med)/sd
                    freq=t["channel0_mhz"]+np.arange(lo,hi)*df/1e6
                    ix=np.flatnonzero(np.abs(freq-s["center_mhz"])*1e6<=s["half_hz"])
                    if not len(ix):continue
                    j=int(ix[np.argmax(scores[ix])]);v=float(scores[j])
                    record={"score":v,"frequency_mhz":float(freq[j]),"drift_hz_s":d,"drift_change_hz_s":a,"start":start,"rows":width}
                    if ans is None or v>ans["score"]:ans=record
        best.append(ans)
    return best
def r15():
    z,n,p=archive("blc");m,mn,mp=archive("blc_manifest")
    manifest=get(m,mn,"manifest.json");result=get(z,n,"result.json");rows=[];ind={"status":"HOLD"}
    for r in result["rows"]:
        if "tile" not in r:continue
        x=tile(z,n,r["tile"],r["artifact"]);best=search(x,r["tile"])
        shams=[search(x,r["tile"],12000+i) for i in range(P["r15d"]["sham_repeats"])]
        item={"name":r["name"],"role":r["role"],"best":best,"sham":shams}
        if r["name"]==P["r15d"]["anchor"]:
            f=next(f for f in manifest["selected"] if f["name"]==r["name"])
            try:ind=compare_reader(f,x,r["tile"])
            except Exception as e:ind={"status":"HOLD","error":str(e)}
        rows.append(item);emit("R15D_SCAN",item)
    a=next(r for r in rows if r["name"]==P["r15d"]["anchor"])
    d={"rows":rows,"anchor":a,"independent_decoder":ind,"scope":"Known-signal recovery pilot; multiple grid/window trials and sham; score not published SNR; no origin classification","production_adoption":False}
    save("result.json",d);emit("R15D_SUMMARY",{"anchor":a,"independent_decoder":ind,"read":len(rows)})
    z.close();m.close();p.unlink();mp.unlink()
def r16():
    z,n,p=archive("atlas");d=get(z,n,"result.json");rows=[]
    for w in d["windows"]:
        if w["status"]!="COMPLETE" or w["center_mhz"] not in P["r16d"]["centers"]:continue
        ts=[r["tile"] for r in w["rows"]];ys=[eng.normalized(tile(z,n,t)) for t in ts]
        for key,indices in (("forward",[0,2]),("reverse",[1,3])):
            for a in w[key]:
                if a["status"]!=P["r16d"]["status"]:continue
                # Existing candidates only; compare equivalent fixed-frequency and continuous predictions on all six scans.
                # Candidate list does not carry explicit anchor index; match recorded bin/shift to stored original hits.
                possible=[i for i in indices if any(h["bin"]==a["bin"] and h["shift"]==a["shift"] and h["score"]==a["score"] for h in w["rows"][i]["hits"])]
                if len(possible)!=1:
                    rows.append({"direction":key,"candidate":a,"status":"HOLD_ANCHOR_AMBIGUOUS"});continue
                i=possible[0];f=a["frequency_mhz"];dr=a["drift_hz_s"];matches=[]
                for k,(t,y) in enumerate(zip(ts,ys)):
                    if k==i:continue
                    # R16C stores dt for adjacent comparisons only. Fixed frequency is valid for all scans;
                    # use recorded dt solely where it exists, never fabricate absolute time.
                    c=a["comparisons"][k-i-1] if k-i in (1,2) else None
                    pred=f+dr*c["dt_seconds"]/1e6 if c else None
                    fixed=r11.track_match(y,t,f,dr)
                    continuous=r11.track_match(y,t,pred,dr) if pred is not None else None
                    rng=np.random.default_rng(16000+10*i+k)
                    sham=[r11.track_match(y[rng.permutation(len(y))],t,pred if pred is not None else f,dr) for _ in range(P["r16d"]["sham_repeats"])]
                    matches.append({"role":w["rows"][k]["role"],"fixed":fixed,"continuous":continuous,"sham":sham})
                item={"center_mhz":w["center_mhz"],"direction":key,"anchor":w["rows"][i]["role"],"candidate":a,"matches":matches,"status":"HOLD_SAME_ARCHIVE_ONLY"}
                rows.append(item);emit("R16D_CANDIDATE",item)
    save("result.json",{"candidates":rows,"n":len(rows),"confirmed_technosignatures":0,"scope":"Fixed candidate replay, same archive only; no independent observation","production_adoption":False})
    emit("R16D_SUMMARY",{"n":len(rows),"confirmed_technosignatures":0})
    z.close();p.unlink()
def r17():
    z,n,p=archive("hybrid");d=get(z,n,"result.json");cal=[];test=[]
    for k in sorted(n):
        if k.startswith("target-") and k.endswith(".json"):
            a=get(z,n,k);cal.extend(r for r in a if r["split"]=="calibration");test.extend(r for r in a if r["split"]=="evaluation")
    assert len(cal)==512 and len(test)==1024
    ct=d["lock"]["thresholds"]["classical"];ht=d["lock"]["thresholds"]["hybrid"]
    def met(a,band):
        # Preserve every classical alert; consult AI only on below-threshold ambiguous cases.
        route=[ct-band<=r["classical"]<ct for r in a]
        pred=[r["classical"]>=ct or (q and r["hybrid"]>=ht) for q,r in zip(route,a)]
        pos=sum(r["label"]==1 for r in a);neg=len(a)-pos
        tp=sum(v and r["label"]==1 for v,r in zip(pred,a));fp=sum(v and r["label"]==0 for v,r in zip(pred,a))
        return {"n":len(a),"detected":tp,"control_alerts":fp,"recall":tp/pos if pos else None,"ai_calls":sum(route),"estimated_seconds":sum(r["classical_seconds"]+(r["hybrid_seconds"] if q else 0) for q,r in zip(route,a))}
    options=[{"band":b,**met(cal,b)} for b in [0]+P["r17d"]["bands"]]
    chosen=min(options,key=lambda r:(-r["recall"],r["control_alerts"],r["estimated_seconds"]))
    lock={"band":chosen["band"],"classical_threshold":ct,"hybrid_threshold":ht,"options":options,"scope":"Calibration-only selection maximizing recall first; no classical alert veto"}
    save("lock.json",lock)
    out={"lock":lock,"test_classical":met(test,0),"test_selective":met(test,chosen["band"]),"hard_classical":met([r for r in test if r["hard"]],0),"hard_selective":met([r for r in test if r["hard"]],chosen["band"]),"test_full_hybrid":d["prospective"]["hybrid"],"scope":"Previously examined R17C data, estimated additive timing from stored per-method timings; development only; fresh run and actual routing timing required","production_adoption":False}
    save("result.json",out);emit("R17D_SUMMARY",out);z.close();p.unlink()
{"r15d":r15,"r16d":r16,"r17d":r17}[lane]()
