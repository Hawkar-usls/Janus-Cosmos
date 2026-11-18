#!/usr/bin/env python3
"""Unmodified BLISS CPU integration, HDF5 adapter, exact complete-path validation."""
import json,pathlib,time,hashlib
import numpy as np,h5py,torch,blissdedrift as bliss
OUT=pathlib.Path("runtime/r28b_bliss");OUT.mkdir(parents=True,exist_ok=True)
IN=pathlib.Path("runtime/r28b/bliss_inputs")
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def write(V,path):
    # Synthetic index coordinates for fixed drift-bank equivalence; no physical header claim.
    n,N,W=V.shape
    with h5py.File(path,"w") as f:
        f.attrs["CLASS"]="FILTERBANK";f.attrs["VERSION"]="1.0"
        d=f.create_dataset("data",data=V.transpose(1,0,2).reshape(N,1,n*W),dtype="f4")
        for k,v in {"fch1":1000.,"foff":15e-6,"tsamp":1.,"tstart":58000.,"src_raj":0.,"src_dej":1.,"source_name":"R28B_INDEX_ADAPTER","telescope_id":4,"machine_id":0,"data_type":1,"nbits":32,"nchans":n*W,"nifs":1}.items():d.attrs[k]=v
        d.attrs.create("DIMENSION_LABELS",np.asarray(["time","feed_id","frequency"],dtype=object),dtype=h5py.string_dtype())
def arr(a):return torch.utils.dlpack.from_dlpack(a.__dlpack__()).numpy().copy()
def run(V,tag):
    path=OUT/(tag+".h5");t=time.perf_counter();write(V,path);write_seconds=time.perf_counter()-t
    t=time.perf_counter();scan=bliss.scan(str(path),V.shape[2]);scan.set_device("cpu",False)
    opts=bliss.integrate_drifts_options();opts.desmear=False;opts.low_rate_Hz_per_sec=-32.;opts.high_rate_Hz_per_sec=33.;opts.resolution=1
    n,N,W=V.shape;assert scan.num_coarse_channels==n
    scores=[];paths=[];compute=0.;bridge=0.;read=0.
    for i in range(n):
        t=time.perf_counter();cc=scan.read_coarse_channel(i);x=arr(cc.data)
        assert x.shape==V[i].shape and np.array_equal(x,V[i]),"BLISS input mismatch"
        read+=time.perf_counter()-t;t=time.perf_counter()
        integrated=bliss.integrate_drifts(cc,opts);plane=integrated.integrated_drift_plane();compute+=time.perf_counter()-t
        t=time.perf_counter();a=arr(plane.integrated_drifts);rates=plane.drift_rate_info()
        assert a.shape==(65,W) and len(rates)==65
        best=-np.inf;p=None;shifts=[]
        for rate in rates:
            sh=int(round(float(rate.drift_rate_slope)*(N-1)));shifts.append(sh)
            d=np.sign(sh)*np.floor(np.abs(np.arange(N)*float(rate.drift_rate_slope))+.5).astype(int)
            lo=max(0,-int(d.min()));hi=min(W,W-int(d.max()));z=a[rate.index_in_plane,lo:hi]*np.sqrt(N)
            q=int(z.argmax());score=float(z[q]);trace=lo+q+d
            # Check all complete paths against independent sums, tolerating float32 accumulation.
            sums=sum(x[j,lo+k:hi+k].astype(np.float64) for j,k in enumerate(d))/np.sqrt(N)
            assert np.allclose(z,sums,rtol=2e-5,atol=2e-5),"native plane mismatch"
            if score>best:best=score;p=trace
        assert sorted(shifts)==list(range(-32,33)) and p is not None
        scores.append(best);paths.append(p);bridge+=time.perf_counter()-t
    del scan
    path.unlink()
    return np.asarray(scores),np.stack(paths),{"write_seconds":write_seconds,"read_and_input_check_seconds":read,"native_integration_seconds":compute,"plane_transfer_and_witness_seconds":bridge,"total_measured_seconds":write_seconds+read+compute+bridge}
def main():
    rng=np.random.default_rng(2811);V=rng.normal(size=(4,16,128)).astype(np.float32)
    for i,s in enumerate((-32,-9,0,24)):
        d=np.rint(np.arange(16)*s/15).astype(int);V[i,np.arange(16),64+d]+=15
    z,p,t=run(V,"selftest");assert np.all(abs(p[:,0]-64)<=1)
    print("R28B_BLISS_SELFTEST",json.dumps({"status":"PASS","version":bliss.__version__,"scores":z.tolist()}),flush=True)
    cal=[];ev=[]
    for f in sorted(IN.glob("cal_*.npz")):
        with np.load(f,allow_pickle=False) as q:
            s,p,t=run(q["V"],f.stem);cal.append(s[q["y"]==0])
    assert len(cal)==2;v=np.sort(np.concatenate(cal))[::-1];th=float(np.nextafter(v[int(len(v)*.05)],np.inf));save("lock.json",{"threshold":th,"controls":len(v)})
    for f in sorted(IN.glob("eval_*.npz")):
        with np.load(f,allow_pickle=False) as q:
            s,p,t=run(q["V"],f.stem);y=q["y"];amp=q["amp"];active=q["active"];loc=((abs(p-q["truth"])<=2)&active).sum(1)>=.75*active.sum(1);pred=s>=th
            d={"target":f.stem,"localized":int(np.sum(pred&loc&(y==1))),"detected":int(np.sum(pred&(y==1))),"controls":int(np.sum(pred&(y==0))),"positives":int(np.sum(y==1)),"negatives":int(np.sum(y==0)),"weak_localized":int(np.sum(pred&loc&(y==1)&(amp<=1.25))),"timing":t,"paired":{}}
            for k in ("adaptive_router","registered_hybrid","registered_linear","registered_odonto","registered_taylor"):
                bp=q[k+"_alert"];bl=q[k+"_localized"];aa=pred&loc&(y==1);bb=bp&bl&(y==1);d["paired"][k]={"bliss_added":int(np.sum(aa&~bb)),"bliss_lost":int(np.sum(bb&~aa)),"bliss_added_controls":int(np.sum(pred&~bp&(y==0))),"bliss_removed_controls":int(np.sum(bp&~pred&(y==0)))}
            ev.append(d);save("partial.json",ev)
    assert len(ev)==4
    total={k:sum(v[k] for v in ev) for k in ("localized","detected","controls","positives","negatives","weak_localized")};assert total["positives"]==total["negatives"]==512
    total["timing"]={k:sum(v["timing"][k] for v in ev) for k in ev[0]["timing"]}
    d={"status":"NATIVE_BLISS_CPU_KERNEL_COMPARISON","version":bliss.__version__,"source":"2b98afe960f13ee7e467aca499576d87ee7502f5","threshold":th,"metrics":total,"per_target":ev,"scope":"Original unmodified CPU BLISS integration only. Same normalized registered float32 input; artificial coordinate header for slope bank, 65 drift slopes, desmear disabled, complete paths and independently checked sums. Separate HDF5 and DLPack bridge time; not full BLISS pipeline or GPU comparison.","production_adoption":False};save("result.json",d);print("R28B_BLISS_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
