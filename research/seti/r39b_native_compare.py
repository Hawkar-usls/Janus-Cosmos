#!/usr/bin/env python3
"""R39B: native comparator replay, same archived inputs, historical null lock."""
import hashlib,io,json,os,pathlib,sys,time,zipfile
import numpy as np
import requests
import r28_hybrid as h
import r35g_native_compare as native
MODE,TARGET=sys.argv[1:3]
assert MODE in ("bliss","fastdd","turbo")
P=pathlib.Path("research/seti/r39b_protocol.json")
S=json.loads(P.read_text());C=S["native_compare"];L=S["detector"]["lock"]
OUT=pathlib.Path("runtime/r39b_native")/MODE/TARGET
OUT.mkdir(parents=True,exist_ok=True)
WORK=OUT/"work";WORK.mkdir(exist_ok=True);native.OUT=WORK
# Restore original detector constants after importing legacy adapters.
h.S=json.loads(pathlib.Path("research/seti/r28_hybrid_protocol.json").read_text())
h.r.S=h.S;h.od.S=h.S
def sha(x):return hashlib.sha256(x).hexdigest()
def save(n,v):(OUT/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def emit(k,v):print(k,json.dumps(v,allow_nan=False),flush=True)
def artifact(a,cap=128):
    u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{a['id']}/zip"
    r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=120)
    r.raise_for_status();b=r.content
    assert len(b)<cap*1024**2 and sha(b)==a["sha256"]
    z=zipfile.ZipFile(io.BytesIO(b));info=[i for i in z.infolist() if not i.is_dir()]
    assert sum(i.file_size for i in info)<512*1024**2
    n={pathlib.PurePosixPath(i.filename).name:i.filename for i in info};assert len(n)==len(info)
    return z,n
def r28(V):
    X=V.astype(np.float64)
    a=h.a.straight(X);b=h.od.sliced(X,0,None)[:2]
    q={"registered_linear":(*a,0.),"registered_odonto":(*b,0.)}
    z,p,_,_=h.fuse(q,L,"registered_hybrid")
    return z,p
def cut(a):
    v=np.sort(np.asarray(a,dtype=float))[::-1]
    assert len(v)==256 and np.isfinite(v).all()
    t=float(np.nextafter(v[int(len(v)*C["rate"])],np.inf))
    assert int((v>=t).sum())<=12
    return t,int((v>=t).sum())
def calibration(fn):
    z,n=artifact(C["calibration"])
    rec=json.loads(z.read(n["receipt.json"]));ms=json.loads(z.read(n["manifest.json"]))
    assert rec["status"]=="PASS_COMMON_REPLAY"
    assert not {q["target"] for q in ms["selected"]}&{q[0] for q in S["targets"]}
    assert not {int(q["day"]) for q in ms["selected"]}&{q[1] for q in S["targets"]}
    a=[];b=[]
    for r in rec["records"]:
        if not r["name"].startswith("cal_"):continue
        raw=z.read(n[r["name"]]);assert sha(raw)==r["sha256"]
        with np.load(io.BytesIO(raw),allow_pickle=False) as q:
            V=q["V"];assert V.dtype==np.float32 and sha(V.tobytes())==r["array_sha256"]
            m=q["y"]==0
            ns,_,_=fn(V,"cal_"+pathlib.Path(r["name"]).stem)
            rs,_=r28(V);a.extend(ns[m]);b.extend(rs[m])
    z.close()
    nt,nc=cut(a);rt,rc=cut(b)
    d={"status":"FROZEN_HISTORICAL_COMMON_NULL","native_threshold":nt,"r28_threshold":rt,"native_cal_alerts":nc,"r28_cal_alerts":rc,"controls":len(a),"calibration":C["calibration"],"protocol_sha256":sha(P.read_bytes())}
    save("lock.json",d);emit("R39B_NATIVE_LOCK",d);return d
def main():
    begin=time.perf_counter();fn,setup=native.init()
    lock=calibration(fn);lh=sha((OUT/"lock.json").read_bytes())
    art=S["parent_search_artifacts"][TARGET];z,n=artifact(art)
    m=json.loads(z.read(n["manifest.json"]));r=json.loads(z.read(n["result.json"]))
    assert r["status"]=="ARCHIVE_SEARCH_COMPLETE_REVIEW_ONLY" and r["target"]==TARGET and r["no_injection"]
    assert r["patches_total"]==8944 and len(r["records"])==16
    assert [TARGET,r["day"]] in S["targets"] and m["lock"]==L
    tiles={}
    for v in r["original_tile_hashes"]:
        t=v["tile"];raw=z.read(n[t["artifact"]]);x=np.load(io.BytesIO(raw),allow_pickle=False)
        assert list(x.shape)==t["shape"] and sha(x.tobytes())==t["slice_sha256"]
        tiles[v["role"],t["center_mhz"]]=(t,x)
    assert len(tiles)==16
    paired={"ON1":"OFF1","ON2":"OFF2","OFF1":"OFF2","OFF2":"OFF1"}
    records=[];candidate_review=[]
    for v in r["records"]:
        role=v["role"];w=v["window_mhz"];tag=f"{role}_{w}"
        assert sha((OUT/"lock.json").read_bytes())==lh
        a,X=tiles[role,w];b,R=tiles[paired[role],w]
        X,R,offset=h.r.common(a,b,X,R)
        raw=z.read(n[v["scores_file"]]);assert sha(raw)==v["sha256"]
        with np.load(io.BytesIO(raw),allow_pickle=False) as q:
            starts=q["starts"].copy();orig=q["score_registered_hybrid"].copy();origp=q["path_registered_hybrid"].copy()
            assert float(q["df_hz"])==v["df_hz"] and float(q["frequency0_hz"])==v["frequency0_hz"]
        expected=np.array(sorted(set(list(range(0,X.shape[1]-128+1,64))+[X.shape[1]-128])))
        assert np.array_equal(starts,expected) and len(starts)==v["patches"]
        assert abs(a["channel0_mhz"]*1e6+offset*a["df_hz"]-v["frequency0_hz"])<1e-5
        ss=[];pp=[];rr=[];rp=[];hh=[];tim=[];err=0.
        for lo in range(0,len(starts),C["chunk"]):
            st=starts[lo:lo+C["chunk"]]
            xs=np.stack([X[:,i:i+128].astype(float) for i in st])
            rs=np.stack([R[:,i:i+128].astype(float) for i in st])
            reg=h.r.register(h.r.norm(xs),h.r.norm(rs))[0]
            V=np.ascontiguousarray(reg,dtype=np.float32);hh.append(sha(V.tobytes()))
            t=time.perf_counter();rz,pr=r28(V);rt=time.perf_counter()-t
            t=time.perf_counter();nz,pn,nt=fn(V,tag+"_"+str(lo));tt=time.perf_counter()-t
            assert len(nz)==len(rz)==len(st) and pn.shape==pr.shape==(len(st),16)
            assert np.isfinite(nz).all() and np.isfinite(rz).all()
            err=max(err,float(np.max(abs(rz-orig[lo:lo+len(st)]))))
            ss.append(nz);pp.append(pn);rr.append(rz);rp.append(pr);tim.append({"native_adapter_wall":tt,"r28_wall":rt,"native":nt})
        ns=np.concatenate(ss);np0=np.concatenate(pp);rs=np.concatenate(rr);rp0=np.concatenate(rp)
        na=ns>=lock["native_threshold"];ra=rs>=lock["r28_threshold"]
        valid=((np0>=0)&(np0<128)).all(1);ng=na&valid
        np1=np0+starts[:,None];rp1=rp0+starts[:,None]
        # Do not count edge-truncated native paths as evidence; preserve their actual alerts.
        delta=np.abs(np1-rp1)
        agree=(delta<=C["tolerance_channels"]).sum(1)>=12
        row={"role":role,"window_mhz":w,"patches":len(starts),"native_alerts":int(na.sum()),"native_complete_alerts":int(ng.sum()),"native_edge_alerts":int((na&~valid).sum()),"r28_alerts":int(ra.sum()),"native_only_patch":int((na&~ra).sum()),"r28_only_patch":int((ra&~na).sum()),"both_alert_same_patch_path":int((ng&ra&agree).sum()),"r28_original_max_score_difference":err,"common_float32_chunk_sha256":hh,"timing":tim}
        np.savez_compressed(OUT/(tag+"_scores.npz"),starts=starts,native_score=ns,native_path=np1,native_alert=na,native_complete_path=valid,r28_score=rs,r28_path=rp1,r28_alert=ra,original_r39_score=orig,original_r39_path=origp)
        if role=="ON1":
            for c in r["review"]:
                if c["window_mhz"]!=w:continue
                i=c["patch_index"];path=(np.asarray(c["path_hz"])-v["frequency0_hz"])/v["df_hz"]
                assert np.allclose(path,origp[i],atol=1e-5)
                match=lambda p:np.flatnonzero((abs(p-path)<=C["tolerance_channels"]).sum(1)>=12)
                ni=match(np1);ri=match(rp1)
                candidate_review.append({"window_mhz":w,"patch_index":i,"frequency_hz":c["frequency_hz"],"original_score":c["score"],"native_matching_alerts":ni[ng[ni]].tolist(),"r28_matching_alerts":ri[ra[ri]].tolist(),"original_on2":c["on2"] is not None,"original_off":any(x is not None for a in c["off"].values() for x in a.values())})
        records.append(row);save("partial.json",{"records":records,"review":candidate_review})
        emit("R39B_NATIVE_SCAN",{k:x for k,x in row.items() if k not in ("timing","common_float32_chunk_sha256")})
    assert len(candidate_review)==r["shortlisted"] and sha((OUT/"lock.json").read_bytes())==lh
    def total(role):
        rr=[x for x in records if x["role"].startswith(role)]
        return {k:sum(x[k] for x in rr) for k in ("patches","native_alerts","native_complete_alerts","native_edge_alerts","r28_alerts","native_only_patch","r28_only_patch","both_alert_same_patch_path")}
    out={"status":"NATIVE_ARCHIVE_REPLAY_REVIEW_ONLY","target":TARGET,"mode":MODE,"source":art,"lock_sha256":lh,"on":total("ON"),"off":total("OFF"),"shortlist":len(candidate_review),"shortlist_native_match":sum(bool(v["native_matching_alerts"]) for v in candidate_review),"shortlist_r28_match":sum(bool(v["r28_matching_alerts"]) for v in candidate_review),"elapsed_seconds":time.perf_counter()-begin,"setup_seconds":setup,"records":records,"review":candidate_review,"production_adoption":False,"scope":C["scope"]}
    save("result.json",out);emit("R39B_NATIVE_SUMMARY",{k:v for k,v in out.items() if k not in ("records","review")})
if __name__=="__main__":
    try:main()
    except Exception as e:
        save("hold.json",{"status":"HOLD_INCOMPLETE_NATIVE_REPLAY","error":type(e).__name__+": "+str(e)[:500]})
        raise
