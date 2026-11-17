#!/usr/bin/env python3
"""Known Voyager 1 carrier calibration with JANUS native drift bank."""
import hashlib,json,pathlib,time
from urllib.parse import urlparse
import numpy as np
import requests,h5py,hdf5plugin
import r9_expanded_search as eng
P=json.loads(pathlib.Path("research/seti/r26f_voyager_protocol.json").read_text())
OUT=pathlib.Path("runtime/r26f");OUT.mkdir(parents=True,exist_ok=True)
def download():
    dest=OUT/"voyager.h5";h=hashlib.sha256();size=0
    start=time.perf_counter()
    try:
        response=requests.get(P["url"],stream=True,timeout=(20,90),allow_redirects=False)
    except requests.exceptions.SSLError:
        print("R26F_TRANSPORT official public HTTP fallback; SHA256 verification required",flush=True)
        response=requests.get(P["public_http_fallback"],stream=True,timeout=(20,90),allow_redirects=False)
    with response:
        response.raise_for_status()
        assert response.status_code==200 and urlparse(response.url).hostname=="blpd0.ssl.berkeley.edu"
        assert response.url in (P["url"],P["public_http_fallback"])
        with dest.open("wb") as f:
            for b in response.iter_content(1024**2):
                size+=len(b);assert size<=P["bytes"];h.update(b);f.write(b)
    assert size==P["bytes"] and h.hexdigest()==P["sha256"],"source size/hash mismatch"
    return dest,time.perf_counter()-start
def scalar(a):
    return a.item() if isinstance(a,np.generic) else a
def main():
    path,load=download();start=time.perf_counter()
    with h5py.File(path,"r") as f:
        d=f["data"];a=d.attrs;f0=float(a["fch1"]);df=float(a["foff"])*1e6;ts=float(a["tsamp"])
        x=np.asarray(d[:,0,:],dtype=np.float32)
        source=scalar(a["source_name"]);source=source.decode() if isinstance(source,bytes) else str(source)
        assert source.lower()=="voyager1" and x.shape==(16,1048576) and np.isfinite(x).all()
        assert abs(df+2.7939677238464355)<1e-8 and abs(ts-18.253611008)<1e-6
    norm=eng.normalized(x);del x
    assert np.isfinite(norm).all()
    w=norm.shape[1];n=norm.shape[0];step=abs(df)/((n-1)*ts)
    lim=int(np.floor(P["search"]["max_abs_drift"]/step))
    block=P["search"]["chunk_bins"];guard=P["search"]["guard_bins"];assert guard>lim
    candidates=[];anchor=None;target=P["expected"]
    for lo in range(0,w,block):
        hi=min(w,lo+block);a=max(0,lo-guard);b=min(w,hi+guard)
        y=norm[:,a:b];pool=[]
        for shift in range(-lim,lim+1):
            z,offset,ds=eng.scores(y,shift)
            aa=max(lo,a+offset);bb=min(hi,a+offset+len(z))
            if aa>=bb:continue
            zz=z[aa-a-offset:bb-a-offset];dr=shift*df/((n-1)*ts)
            def rec(j):
                k=aa+int(j);freq=f0+k*df/1e6
                q=k+ds;assert np.all((q>=0)&(q<w))
                vals=norm[np.arange(n),q]
                score=float(zz[j])
                assert abs(float(vals.sum()/np.sqrt(n))-score)<=max(1e-4,abs(score)*1e-6)
                return {"frequency_mhz":freq,"drift_hz_s":dr,"score":score,"start_bin":k,"end_bin":int(q[-1]),"rows":n}
            k=min(4,len(zz))
            for j in np.argpartition(zz,-k)[-k:]:
                pool.append(rec(j))
            if abs(dr-target["drift_hz_s"])<=target["drift_tolerance"]:
                center=(target["frequency_mhz"]-f0)*1e6/df;rr=target["frequency_tolerance_hz"]/abs(df)
                l=max(0,int(np.ceil(center-rr))-aa);h=min(len(zz),int(np.floor(center+rr))+1-aa)
                if l<h:
                    q=rec(l+int(np.argmax(zz[l:h])))
                    if anchor is None or q["score"]>anchor["score"]:anchor=q
        keep=[]
        for q in sorted(pool,key=lambda z:-z["score"]):
            if any(abs(q["start_bin"]-u["start_bin"])<8 for u in keep):continue
            keep.append(q)
            if len(keep)==P["search"]["top_per_chunk"]:break
        candidates+=keep
    best=[]
    for q in sorted(candidates,key=lambda z:-z["score"]):
        if any(abs(q["start_bin"]-u["start_bin"])<8 for u in best):continue
        best.append(q)
        if len(best)==16:break
    recovered=anchor is not None and anchor["score"]>=P["search"]["minimum_robust_score"]
    d={"status":"KNOWN_CARRIER_RECOVERED" if recovered else "KNOWN_CARRIER_NOT_RECOVERED","file_sha256":P["sha256"],"source_bytes":P["bytes"],"header":{"source":source,"fch1_mhz":f0,"df_hz":df,"ts_s":ts,"shape":[n,w]},"anchor":anchor,"top16":best,"download_seconds":load,"search_seconds":time.perf_counter()-start,"drift_trials":2*lim+1,"scope":P["scope"],"production_adoption":False}
    (OUT/"result.json").write_text(json.dumps(d,indent=2)+"\n");print("R26F_SUMMARY",json.dumps(d),flush=True)
    assert recovered
if __name__=="__main__":main()
