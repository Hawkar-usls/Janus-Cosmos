#!/usr/bin/env python3
"""Measured-power resolution transfer check: derived observations, no new labels."""
import hashlib,io,json,os,pathlib,zipfile
import numpy as np
import requests
import r30c_hannah_memory as mem
P=pathlib.Path("research/seti/r37a_hannah_resolution_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r37a");O.mkdir(parents=True,exist_ok=True)
def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def retrieve(cfg):
    u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{cfg['artifact']}/zip"
    r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90);r.raise_for_status()
    raw=r.content;assert len(raw)<100*1024**2 and sha(raw)==cfg["digest"]
    z=zipfile.ZipFile(io.BytesIO(raw));assert sum(i.file_size for i in z.infolist())<256*1024**2
    names={pathlib.PurePosixPath(i.filename).name:i.filename for i in z.infolist() if not i.is_dir()}
    assert len(names)==len([i for i in z.infolist() if not i.is_dir()])
    b=pathlib.Path(cfg["sidecar"]).read_bytes();assert sha(b)==cfg["sidecar_sha256"]
    return z,names,json.loads(b)
def records():
    for k,cfg in S["source"].items():
        z,n,d=retrieve(cfg);r=json.loads(z.read(n["raw_receipts.json"]))
        if k=="stationary":
            rr=r;ps=[(e["family_group"],e["raw_evidence"]["passport"],"EXTERNAL_CENTRAL_REPLAY") for e in d["entries"]]
        else:
            rr=[m for q in r for m in q["measurements"] if "tile" in m]
            ps=[(e["published"]["id"],p,"UNKNOWN_REVIEW_ONLY") for e in d["entries"] for p in [e["passport"]]+[v["passport"] for v in e["off_matches"]]]
        for group,pr,label in ps:
            mem.validate(pr);m=pr["provenance"];h=m["source_tile_sha256"]
            matches=[v for v in rr if v["tile"]["slice_sha256"]==h and v["scan"]==m["source_scan"]];assert len(matches)==1
            v=matches[0];nm=v["artifact"] if k=="stationary" else v["raw_artifact"]
            x=np.load(io.BytesIO(z.read(n[nm])),allow_pickle=False);assert list(x.shape)==v["tile"]["shape"] and sha(x.tobytes())==h
            start=m["start_row"];p=np.asarray(pr["identity"]["track"],int);x=x[start:start+len(p)]
            if k=="drifting":assert sha(x.tobytes())==m["row_slice_sha256"]
            assert mem.make_record(x,p,m,verify_full=(k=="stationary"))==pr
            yield {"id":pr["id"],"group":group,"scope":label,"hash":h,"x":x[:16].astype(float),"p":p[:16],"df":m["df_hz"],"ts":m["ts"]}
        z.close()
def weights(src,dst):
    q=np.maximum(0,np.minimum(dst[1:,None],src[None,1:])-np.maximum(dst[:-1,None],src[None,:-1]))
    cov=q.sum(1);expected=np.diff(dst)
    assert np.allclose(cov,expected,rtol=1e-10,atol=1e-9),"incomplete physical support"
    return q/expected[:,None]
def edges(src,step,phase):
    lo=src[0]+phase*step;count=int(np.floor((src[-1]-lo)/step+1e-10))
    assert count>=1;return lo+step*np.arange(count+1)
def rebin(x,t,f,newt,newf):
    wt=weights(t,newt);wf=weights(f,newf)
    y=wt@x@wf.T
    # Compare integrated power over exactly the represented rectangular support.
    old=np.sum(x*(wt.T@np.diff(newt))[:,None]*(wf.T@np.diff(newf))[None,:])
    new=np.sum(y*np.diff(newt)[:,None]*np.diff(newf)[None,:])
    err=abs(old-new)/max(1.,abs(old));assert err<=1e-10
    return y,err
def descriptor(x,t,f,trajectory_t,trajectory_f,offset=0.):
    y=mem.norm(x.astype(float));time_centers=(t[:-1]+t[1:])/2
    track=np.interp(time_centers,trajectory_t,trajectory_f)+offset
    target=np.asarray(S["physical_grid"]["frequency_edges_hz"],float)
    a=np.stack([weights(f,target+fc)@row for fc,row in zip(track,y)])
    cells=weights(t,np.asarray(S["physical_grid"]["time_edges_seconds"],float))@a
    # Retain unscaled contrast and shape separately. Never compare absolute SNR across grids.
    contrast=cells-(cells[:,[0,-1]].mean(1)[:,None]);v=np.maximum(contrast,0).ravel()
    length=np.linalg.norm(v)
    if not np.isfinite(length) or length<=1e-12:raise ValueError("no supported shape")
    return v/length,float(length),cells.tolist()
def sim(a,b):return float(np.clip(a@b,-1.,1.))
def main():
    u=np.arange(10,dtype=float);x=np.full((9,9),4.);t=np.array([.5,2.5,4.5,6.5,8.5]);y,e=rebin(x,u,u,t,t)
    assert np.array_equal(y,np.full((4,4),4.)) and e==0
    rows=list(records());assert len(rows)==10 and len({r["hash"] for r in rows})==10
    out=[];descs=[]
    for r in rows:
        x=r["x"];p=r["p"];df=r["df"];ts=r["ts"];assert x.shape[0]==16
        if df<0:x=x[:,::-1];p=x.shape[1]-1-p
        df=abs(df);te=np.arange(17)*ts;fe=(np.arange(x.shape[1]+1)-.5)*df
        tt=(te[1:]+te[:-1])/2;ff=p*df
        a,amp,_=descriptor(x,te,fe,tt,ff);descs.append((r,a))
        for test in S["physical_grid"]["transforms"]:
            t=edges(te,ts*test["time_factor"],test["phase"]);f=edges(fe,df*test["frequency_factor"],test["phase"])
            y,e=rebin(x,te,fe,t,f)
            row={"id":r["id"],"group":r["group"],"scope":r["scope"],"transform":test,"source_hash":r["hash"],"power_relative_error":e}
            try:
                b,ap,c=descriptor(y,t,f,tt,ff)
                row.update({"status":"COMPARED","cosine":sim(a,b),"raw_contrast_norm":amp,"degraded_contrast_norm":ap,"descriptor":b.tolist()})
            except ValueError as err:row.update({"status":"ABSTAIN","reason":str(err)})
            out.append(row)
    assert len(out)==40
    pairs=[]
    for i,(r,a) in enumerate(descs):
        for q,b in descs[i+1:]:
            pairs.append({"a":r["id"],"b":q["id"],"same_reported_group":r["group"]==q["group"],"scope":"DISTINCT_RAW_SAME_CAMPAIGN_NOT_INDEPENDENT","cosine":sim(a,b)})
    valid=[v for v in out if v["status"]=="COMPARED"];ok=sum(v["cosine"]>=S["gate"]["similarity_min"] for v in valid)
    # Retrieval task, not physical-origin classification: candidates remain UNKNOWN.
    same=[v["cosine"] for v in pairs if v["same_reported_group"]];other=[v["cosine"] for v in pairs if not v["same_reported_group"]]
    d={"status":"PASS_RESOLUTION_TRANSFER_DEVELOPMENT" if ok>=36 and len(valid)==40 else "HOLD_RESOLUTION_TRANSFER","measured_passports":len(rows),"distinct_measured_tiles":len(rows),"derived_pairs":len(out),"compared":len(valid),"above_prespecified_threshold":ok,"min_cosine":min([v["cosine"] for v in valid],default=None),"median_cosine":float(np.median([v["cosine"] for v in valid])) if valid else None,"distinct_scan_same_group_median":float(np.median(same)),"different_group_median":float(np.median(other)),"different_group_above_threshold":sum(v>=.9 for v in other),"different_group_pairs":len(other),"power_conservation":"PASS","verified_promoted":0,"scope":"Same measured raw tiles with controlled coarsening; given trajectory, not recovered detector path; one campaign, no independent labels or protected real control. No eligibility for radio origin or automatic suppression.","production_adoption":False}
    save("transfer_pairs.json",out);save("distinct_pairs.json",pairs);save("result.json",d);print("R37A_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
