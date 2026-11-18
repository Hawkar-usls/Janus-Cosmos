#!/usr/bin/env python3
"""R39 bounded, uninjected real archive search; no automatic origin judgement."""
import hashlib, io, json, os, pathlib, sys, time, zipfile
import numpy as np
import requests
import r28_hybrid as h

P=pathlib.Path("research/seti/r39_discovery_protocol.json")
S=json.loads(P.read_text())
L=S["detector"]["lock"]
TARGET=sys.argv[1]
OUT=pathlib.Path("runtime/r39")/TARGET
OUT.mkdir(parents=True,exist_ok=True)
h.OUT=OUT; h.r.OUT=OUT; h.eng.OUT=OUT
assert h.eng.SPEC["window_centers_mhz"]==S["windows_mhz"]
assert h.eng.SPEC["window_half_hz"]==S["half_hz"]
assert h.eng.SPEC["limits"]["spectrum_budget_per_file_mib"]==S["max_mib_per_scan"]
def sha(b):return hashlib.sha256(b).hexdigest()
def save(name,d):
    (OUT/name).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def emit(name,d):print(name,json.dumps(d,allow_nan=False),flush=True)
def manifest():
    a=S["source_manifest"]
    r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['artifact']}/zip",
        headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90)
    r.raise_for_status();b=r.content;assert len(b)<=32*1024**2 and sha(b)==a["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        names=[n for n in z.infolist() if pathlib.PurePosixPath(n.filename).name=="manifest.json"]
        assert len(names)==1 and names[0].file_size<16*1024**2
        d=json.loads(z.read(names[0]))
    assert d["new_obs_power_read"] is False
    selected=[q for q in d["selected"]["Parkes"] if q["target"]==TARGET]
    assert len(selected)==1
    q=selected[0]
    assert [TARGET,q["day"]] in S["targets"] and TARGET not in S["excluded"]
    assert [v["role"] for v in q["files"]]==S["roles"]
    assert len({v["name"] for v in q["files"]})==4
    assert len({v["url"] for v in q["files"]})==4
    save("manifest.json",{"protocol_sha256":sha(P.read_bytes()),"source":a,"group":q,"lock":L})
    emit("R39_LOCK",{"target":TARGET,"day":q["day"],"manifest_sha256":sha((OUT/"manifest.json").read_bytes())})
    return q
def score(a,X,b,R,role,window,mjd):
    X,R,offset=h.r.common(a,b,X,R)
    assert X.shape[0]==R.shape[0]==16 and X.shape[1]>=128
    w=S["patch"];stride=S["stride"];starts=np.array(sorted(set(list(range(0,X.shape[1]-w+1,stride))+[X.shape[1]-w])),dtype=int)
    scores={k:[] for k in S_MODES};paths={k:[] for k in S_MODES}
    for lo in range(0,len(starts),S["chunk"]):
        st=starts[lo:lo+S["chunk"]]
        xs=np.stack([X[:,i:i+w].astype(float) for i in st])
        rs=np.stack([R[:,i:i+w].astype(float) for i in st])
        q,_,_,_=h.prep([(x,) for x in xs],rs)
        for k in S_MODES:
            z,p=q[k][:2] if k in h.BASE else h.fuse(q,L,k)[:2]
            assert np.isfinite(z).all() and p.shape==(len(st),16)
            scores[k].append(z);paths[k].append(p+st[:,None])
    z={k:np.concatenate(v) for k,v in scores.items()}
    p={k:np.concatenate(v) for k,v in paths.items()}
    df=float(a["df_hz"]);f0=float(a["channel0_mhz"])*1e6+offset*df
    out={"starts":starts,"df_hz":df,"frequency0_hz":f0,"ts":float(a["ts"]),"mjd":mjd}
    for k in S_MODES:out["score_"+k]=z[k];out["path_"+k]=p[k]
    fname=f"{role}_{window}_scores.npz"
    np.savez_compressed(OUT/fname,**out)
    rec={"role":role,"window_mhz":window,"patches":len(starts),"native_channels":X.shape[1],"frequency0_hz":f0,"df_hz":df,"ts":float(a["ts"]),"mjd":mjd,"scores_file":fname,"sha256":sha((OUT/fname).read_bytes()),"alerts":{k:int((z[k]>=L["thresholds"][k]).sum()) for k in S_MODES}}
    emit("R39_SCAN",rec)
    return rec,z,p
def candidates(part,mode):
    rec,z,p=part
    ix=np.flatnonzero(z[mode]>=L["thresholds"][mode])
    freqs=rec["frequency0_hz"]+p[mode][ix]*rec["df_hz"]
    return ix,freqs,z[mode][ix]
def correspondence(anchor,other,mode):
    ar,af=anchor
    r=other[0];ix,ff,zz=candidates(other,mode)
    if not len(ix):return None
    elapsed=(r["mjd"]-ar["mjd"])*86400
    drift=(af[-1]-af[0])/(15*ar["ts"])
    expected0=af[0]+elapsed*drift;expected1=expected0+15*r["ts"]*drift
    delta=np.maximum(abs(ff[:,0]-expected0),abs(ff[:,-1]-expected1))
    valid=np.flatnonzero(delta<=4*abs(ar["df_hz"]))
    if not len(valid):return None
    j=valid[np.argmax(zz[valid])]
    return {"mode":mode,"patch_index":int(ix[j]),"score":float(zz[j]),"endpoint_error_hz":float(delta[j])}
S_MODES=h.S["modes"]
def main():
    q=manifest();mh=sha((OUT/"manifest.json").read_bytes());start=time.perf_counter()
    h.a.init_tree()
    tiles={}
    for f in q["files"]:
        assert sha((OUT/"manifest.json").read_bytes())==mh
        tiles[f["role"]]=h.r.tiles(f)
    records=[];shortlist=[]
    paired={"ON1":"OFF1","ON2":"OFF2","OFF1":"OFF2","OFF2":"OFF1"}
    for j,center in enumerate(S["windows_mhz"]):
        parts={}
        for f in q["files"]:
            role=f["role"];a,X=tiles[role][j];b,R=tiles[paired[role]][j]
            assert a["center_mhz"]==b["center_mhz"]==center
            parts[role]=score(a,X,b,R,role,center,float(f["header"]["mjd"]))
            records.append(parts[role][0])
        ar,z,p=parts["ON1"];ix,ff,zz=candidates(parts["ON1"],"registered_hybrid")
        order=np.argsort(-zz,kind="stable")
        used=[]
        for k in order:
            i=int(ix[k]);af=ff[k]
            if any(max(abs(af[0]-v[0]),abs(af[-1]-v[1]))<=4*abs(ar["df_hz"]) for v in used):continue
            used.append((float(af[0]),float(af[-1])))
            repeat=correspondence((ar,af),parts["ON2"],"registered_hybrid")
            controls={r:{m:correspondence((ar,af),parts[r],m) for m in ("raw_linear","raw_odonto")} for r in ("OFF1","OFF2")}
            comparisons={}
            for m in ("registered_linear","registered_odonto","registered_taylor"):
                # compare same ON1 scan, at the same physical trajectory, including overlapping patches
                comparisons[m]=correspondence((ar,af),parts["ON1"],m)
            item={"target":TARGET,"window_mhz":center,"score":float(zz[k]),"patch_index":i,"frequency_hz":float(af[0]),"end_frequency_hz":float(af[-1]),"path_hz":af.tolist(),"linearized_drift_hz_s":float((af[-1]-af[0])/(15*ar["ts"])),"on2":repeat,"off":controls,"component_matches":comparisons,"status":"REVIEW_ONLY","origin_verified":False}
            shortlist.append(item)
            if len(used)>=S["shortlist_per_window"]:break
        save("partial.json",{"records":records,"review":shortlist})
    assert sha((OUT/"manifest.json").read_bytes())==mh
    repeated=[v for v in shortlist if v["on2"] is not None and not any(m is not None for rr in v["off"].values() for m in rr.values())]
    exclusive=[v for v in repeated if not any(x is not None for x in v["component_matches"].values())]
    result={"status":"ARCHIVE_SEARCH_COMPLETE_REVIEW_ONLY","target":TARGET,"day":q["day"],"records":records,"patches_total":sum(v["patches"] for v in records),"shortlisted":len(shortlist),"repeat_without_matched_off":len(repeated),"repeat_without_component_match":len(exclusive),"review":shortlist,"elapsed_seconds":time.perf_counter()-start,"original_tile_hashes":[{"role":r,"tile":t} for r,vv in tiles.items() for t,x in vv],"no_injection":True,"global_false_alarm_calibrated":False,"confirmed_new_signals":0,"production_adoption":False,"scope":S["limits"]}
    save("result.json",result);emit("R39_SUMMARY",{k:v for k,v in result.items() if k not in ("records","review","original_tile_hashes")})
if __name__=="__main__":
    try:main()
    except Exception as e:
        save("hold.json",{"status":"HOLD_INCOMPLETE_SEARCH","target":TARGET,"error":type(e).__name__+": "+str(e)[:500]})
        raise
