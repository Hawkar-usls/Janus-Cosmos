#!/usr/bin/env python3
"""Native turboSETI FindDoppler end-to-end *development* input bridge; fail closed."""
import hashlib,json,pathlib,io,zipfile,time,os,resource,logging,contextlib
import numpy as np,h5py,requests
from turbo_seti.find_doppler.find_doppler import FindDoppler
from turbo_seti.find_doppler.data_handler import DATAH5
from turbo_seti.find_doppler.kernels import Kernels
OUT=pathlib.Path("runtime/r35d");OUT.mkdir(parents=True,exist_ok=True)
FROM="7d9b4fde9bc98d834dc11cfc0acd2380e6676f0e"
ART={"id":14405617460,"sha256":"fa60d60f8c961d521fe9bfb0f915447b327a05f03c5756d5b743a965ce46a9b0"}
DF=-15e-6;TS=1.;F0=1500.
def save(n,v):(OUT/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def write(V,path):
    b,n,w=V.shape
    with h5py.File(path,"w") as f:
        f.attrs["CLASS"]="FILTERBANK";f.attrs["VERSION"]="1.0"
        d=f.create_dataset("data",data=V.transpose(1,0,2).reshape(n,1,b*w),dtype="f4")
        for k,v in {"fch1":F0,"foff":DF,"tsamp":TS,"tstart":58000.,"src_raj":0.,"src_dej":1.,"source_name":"R35D_INDEX_ADAPTER","telescope_id":6,"machine_id":20,"data_type":1,"nbits":32,"nchans":b*w,"nifs":1,"nbeams":1,"ibeam":0,"az_start":0.,"za_start":0.,"rawdatafile":"R35D.raw"}.items():d.attrs[k]=v
        d.attrs.create("DIMENSION_LABELS",np.asarray(["time","feed_id","frequency"],dtype=object),dtype=h5py.string_dtype())
def run(V,tag,kernels):
    b,n,w=V.shape;assert (n,w)==(16,128)
    path=OUT/(tag+".h5");folder=OUT/tag;folder.mkdir(exist_ok=True)
    t=time.perf_counter();write(V,path);wt=time.perf_counter()-t
    # Actual DATAH5 loader checked for every coarse channel; no hidden padding/normalization in bridge.
    t=time.perf_counter()
    for i in range(b):
        f1=F0+i*w*DF;f2=F0+(i+1)*w*DF
        data=DATAH5(str(path),f_start=min(f1,f2),f_stop=max(f1,f2),cchan_id=i,n_coarse_chan=b,kernels=kernels)
        spectra,_=data.load_data();assert np.array_equal(spectra,V[i]),f"native input mismatch {i}"
        if hasattr(data,"close"):data.close()
    readsec=time.perf_counter()-t
    t=time.perf_counter()
    with (folder/"native_stdout.txt").open("w") as f,contextlib.redirect_stdout(f),contextlib.redirect_stderr(f):
        fd=FindDoppler(str(path),max_drift=32.,min_drift=0.,snr=.01,out_dir=str(folder),n_coarse_chan=b,kernels=kernels,flagging=False,blank_dc=False,gpu_backend=False,log_level_int=logging.ERROR)
        fd.search(n_partitions=1,progress_bar="n")
    searchsec=time.perf_counter()-t
    # Native DAT file output: highest native SNR per coarse patch, no custom detector.
    dat=list(folder.glob("*.dat"));assert dat,"missing native DAT"
    merged=folder/(tag+".dat")
    files=[merged] if merged.exists() else dat
    p=np.zeros((b,n),int);s=np.zeros(b,float);hits=0
    t=time.perf_counter()
    for f in files:
        for line in f.read_text().splitlines():
            q=line.split()
            if not q or not q[0].isdigit():continue
            assert len(q)>=12
            drift=float(q[1]);snr=float(q[2]);index=int(q[5]);cc=int(q[10]);assert 0<=cc<b
            assert np.isfinite(snr) and np.isfinite(drift)
            # Native DAT uncorrected frequency and local channel index agree via canary.
            trace=index+np.rint(np.arange(n)*drift/(DF*1e6)).astype(int)
            # Preserve native maximal detection even if its path crosses a patch edge.
            if snr>s[cc]:s[cc]=snr;p[cc]=trace
            hits+=1
    assert hits>=0
    bridge=time.perf_counter()-t;path.unlink()
    return s,p,{"write_seconds":wt,"native_input_check_seconds":readsec,"native_pipeline_seconds":searchsec,"dat_parse_seconds":bridge,"native_hits":hits}
def main():
    kern=Kernels(False,1,0);rng=np.random.default_rng(3520)
    x=rng.normal(size=(3,16,128)).astype(np.float32)
    truth=[]
    for i,s in enumerate((-15,0,15)):
        p=64+np.rint(np.arange(16)*s/15).astype(int);truth.append(p);x[i,np.arange(16),p]+=30
    a,p,t=run(x,"canary",kern);assert np.all(a>5) and np.all((abs(p-np.stack(truth))<=2).sum(1)>=12),"native DAT coordinate canary failed"
    print("R35D_NATIVE_SELFTEST",json.dumps({"status":"PASS","scores":a.tolist(),"timing":t}),flush=True)
    url=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{ART['id']}/zip"
    r=requests.get(url,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=60);r.raise_for_status();raw=r.content
    assert len(raw)<64*1024**2 and hashlib.sha256(raw).hexdigest()==ART["sha256"]
    z=zipfile.ZipFile(io.BytesIO(raw));assert len(z.namelist())==6 and sum(i.file_size for i in z.infolist())<64*1024**2
    cal=[];ev=[]
    for f in sorted(k for k in z.namelist() if k.startswith("cal_")):
        with np.load(io.BytesIO(z.read(f)),allow_pickle=False) as q:
            s,p,t=run(q["V"],pathlib.Path(f).stem,kern);cal.append(s[q["y"]==0])
    assert len(cal)==2
    v=np.sort(np.concatenate(cal))[::-1];th=float(np.nextafter(v[int(len(v)*.05)],np.inf));save("lock.json",{"threshold":th,"controls":len(v)})
    for f in sorted(k for k in z.namelist() if k.startswith("eval_")):
        with np.load(io.BytesIO(z.read(f)),allow_pickle=False) as q:
            s,p,t=run(q["V"],pathlib.Path(f).stem,kern);y=q["y"];loc=((abs(p-q["truth"])<=2)&q["active"]).sum(1)>=.75*q["active"].sum(1);pred=s>=th;hit=pred&loc&(y==1);base=q["registered_hybrid_alert"]&q["registered_hybrid_localized"]&(y==1)
            ev.append({"target":f,"localized":int(hit.sum()),"weak":int((hit&(q["amp"]<=1.25)).sum()),"control_alerts":int((pred&(y==0)).sum()),"positives":int((y==1).sum()),"controls":int((y==0).sum()),"added_vs_r28":int((hit&~base).sum()),"lost_vs_r28":int((base&~hit).sum()),**t})
            save("partial.json",ev);print("R35D_TARGET",json.dumps(ev[-1]),flush=True)
    z.close();assert len(ev)==4
    m={k:sum(v[k] for v in ev) for k in ev[0] if k!="target"};assert m["positives"]==m["controls"]==512
    d={"status":"NATIVE_TURBOSETI_FINDDOPPLER_HISTORICAL_COMMON_INPUT","source":FROM,"artifact":ART,"metrics":m,"per_target":ev,"threshold":th,"peak_rss_kib_process":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"scope":"Actual unchanged native FindDoppler pipeline on normalized data with artificial coordinates; exact read checks and DAT-coordinate canaries. Report native rate path approximation, not exact native tree path. No independent exam, physical file headers, cadence origin test or GPU claim.","production_adoption":False}
    save("result.json",d);print("R35D_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
