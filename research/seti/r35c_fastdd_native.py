#!/usr/bin/env python3
"""R35C: author's unmodified fastDD CPU functions, pinned source and reconstructed path witness."""
import hashlib,json,pathlib,subprocess,ctypes,time,resource,os,io,zipfile
import numpy as np
import requests
OUT=pathlib.Path("runtime/r35c");OUT.mkdir(parents=True,exist_ok=True)
SOURCE="2ba67ea833b330c2f507a790c0f45e864f76126e"
FILES={"fastDD_cpu.h":"e3caa8d3a02dca37f74b66904c05b349d4c63667","fastDD_cpu_fns.cu":"77b45fcfe54eb5297530b729530cc18132f635d9"}
ART={"id":14405617460,"sha256":"fa60d60f8c961d521fe9bfb0f915447b327a05f03c5756d5b743a965ce46a9b0"}
def save(p,v):(OUT/p).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def build():
    src=OUT/"src";src.mkdir(exist_ok=True)
    for name,digest in FILES.items():
        r=requests.get(f"https://raw.githubusercontent.com/khouston22/seticore2/{SOURCE}/{name}",timeout=30);r.raise_for_status();b=r.content
        assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()==digest
        (src/name).write_bytes(b)
    # Only dependency declarations are shimmed; neither native CPU source nor native header is edited.
    (src/"cuda_util.h").write_text("#pragma once\n#define __global__\nusing uint=unsigned int;\n")
    (src/"bridge.cpp").write_text(r"""
#include <vector>
#include <cstring>
#include <cmath>
#include "fastDD_cpu_fns.cu"
extern "C" int native(float* x,int n,int w,float* dst,int cap,double* rates) {
  DD_metadata dd;
  gen_fastDD_metadata(&dd,0.,1.,w,1.,n,-32./15.,32./15.,1.,n,8);
  if(dd.Nr>cap || dd.Nt!=n) return -1;
  std::vector<float> a(dd.Nr_ext*w,0.f),b(dd.Nr_ext*w,0.f);
  float* work[2]={a.data(),b.data()};
  float* z=fastDD_cpu(x,work,&dd);
  std::memcpy(dst,z,size_t(dd.Nr)*w*sizeof(float));
  for(int q=0;q<dd.Nr;q++) rates[q]=(dd.m_min+q)*dd.d_dfdt;
  return dd.Nr;
}
""")
    subprocess.run(["g++","-std=c++17","-O3","-fPIC","-shared","-o",str(src/"libfastdd.so"),str(src/"bridge.cpp")],check=True)
    lib=ctypes.CDLL(str((src/"libfastdd.so").resolve()));fn=lib.native
    fp=np.ctypeslib.ndpointer(np.float32,flags="C_CONTIGUOUS");dp=np.ctypeslib.ndpointer(np.float64,flags="C_CONTIGUOUS")
    fn.argtypes=[fp,ctypes.c_int,ctypes.c_int,fp,ctypes.c_int,dp];fn.restype=ctypes.c_int
    return fn
def plane(fn,x):
    x=np.ascontiguousarray(x,np.float32);n,w=x.shape;z=np.empty((128,w),np.float32);rates=np.empty(128,np.float64)
    k=fn(x,n,w,z,128,rates);assert 0<k<=128
    return z[:k],rates[:k]
def recover_paths(fn):
    w=128;n=16;probe=64;paths=[];rates0=None
    for i in range(n):
        x=np.zeros((n,w),np.float32);x[i]=np.arange(w,dtype=np.float32)
        a,rates=plane(fn,x)
        if rates0 is None:rates0=rates.copy()
        else:assert np.array_equal(rates,rates0)
        q=a[:,probe];assert np.all(q==np.rint(q));paths.append(q.astype(int)-probe)
    D=np.stack(paths,axis=1)
    use=(np.abs(rates0)<=32/15+1e-9)&(np.abs(D[:,-1]-D[:,0])<=32)
    assert np.any(use&(rates0==0)) and np.any(use&(rates0>0)) and np.any(use&(rates0<0))
    return D,rates0,use
def score(fn,V,D,R,use):
    b=[];p=[];compute=bridge=0.
    for x in V:
        start=time.perf_counter();a,r=plane(fn,x);compute+=time.perf_counter()-start
        start=time.perf_counter();assert np.array_equal(r,R)
        best=-np.inf;bp=None;w=x.shape[1]
        for i,d in enumerate(D):
            if not use[i]:continue
            lo=max(0,-int(d.min()));hi=min(w,w-int(d.max()))
            if hi<=lo:continue
            z=a[i,lo:hi]/4
            # Native plane compared with independently reconstructed complete sums; do not use truncated edge paths.
            sums=sum(x[j,lo+k:hi+k].astype(np.float64) for j,k in enumerate(d))/4
            assert np.allclose(z,sums,rtol=2e-5,atol=2e-5),"native path/plane mismatch"
            k=int(z.argmax());s=float(z[k])
            if s>best:best=s;bp=lo+k+d
        assert bp is not None;b.append(best);p.append(bp);bridge+=time.perf_counter()-start
    return np.asarray(b),np.stack(p),{"native_plus_transfer_seconds":compute,"path_validation_seconds":bridge}
def restore():
    u=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{ART['id']}/zip";r=requests.get(u,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=60);r.raise_for_status();b=r.content
    assert len(b)<64*1024**2 and hashlib.sha256(b).hexdigest()==ART["sha256"]
    z=zipfile.ZipFile(io.BytesIO(b));assert sum(i.file_size for i in z.infolist())<64*1024**2 and len(z.namelist())==6
    return z
def main():
    st=time.perf_counter();fn=build();buildsec=time.perf_counter()-st;D,R,use=recover_paths(fn)
    rng=np.random.default_rng(3512);x=rng.normal(size=(3,16,128)).astype(np.float32)
    for i,r in enumerate((-1.,0.,1.)):
        d=D[np.argmin(abs(R-r))];x[i,np.arange(16),64+d]+=40
    s,p,t=score(fn,x,D,R,use)
    assert np.all(abs(p[:,0]-64)<=1) and np.all(s>100)
    save("selftest.json",{"status":"PASS_NATIVE_FASTDD_COMPLETE_PATHS","rates":R.tolist(),"eligible_rates":R[use].tolist(),"paths":D.tolist(),"source":SOURCE,"files":FILES,"build_seconds":buildsec})
    print("R35C_SELFTEST",json.dumps({"status":"PASS","eligible_rates":int(use.sum()),"scores":s.tolist()}),flush=True)
    z=restore();cal=[];ev=[]
    for f in sorted(k for k in z.namelist() if k.startswith("cal_")):
        with np.load(io.BytesIO(z.read(f)),allow_pickle=False) as q:
            s,p,t=score(fn,q["V"],D,R,use);cal.append(s[q["y"]==0])
    assert len(cal)==2
    v=np.sort(np.concatenate(cal))[::-1];th=float(np.nextafter(v[int(len(v)*.05)],np.inf));save("lock.json",{"threshold":th,"calibration_controls":len(v)})
    for f in sorted(k for k in z.namelist() if k.startswith("eval_")):
        with np.load(io.BytesIO(z.read(f)),allow_pickle=False) as q:
            s,p,t=score(fn,q["V"],D,R,use);y=q["y"];loc=((abs(p-q["truth"])<=2)&q["active"]).sum(1)>=.75*q["active"].sum(1);pred=s>=th;hit=pred&loc&(y==1);base=q["registered_hybrid_alert"]&q["registered_hybrid_localized"]&(y==1)
            ev.append({"target":f,"localized":int(hit.sum()),"weak":int((hit&(q["amp"]<=1.25)).sum()),"very_weak":int((hit&(q["amp"]==.75)).sum()),"control_alerts":int((pred&(y==0)).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum()),"added_vs_r28":int((hit&~base).sum()),"lost_vs_r28":int((base&~hit).sum()),**t})
            save("partial.json",ev)
    z.close();assert len(ev)==4;total={k:sum(v[k] for v in ev) for k in ev[0] if k!="target"};assert total["positives"]==total["controls"]==512
    d={"status":"NATIVE_FASTDD_CPU_KERNEL_HISTORICAL_COMMON_INPUT_COMPARISON","source":SOURCE,"source_files":FILES,"input_artifact":ART,"metrics":total,"per_target":ev,"threshold":th,"peak_rss_kib_process":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"scope":"Historical R28B common registered input, native fastDD CPU core N0=8; same rate envelope, own exact native paths; no full seticore2 pipeline, GPU, origin, or independent exam claim.","production_adoption":False}
    save("result.json",d);print("R35C_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
