#!/usr/bin/env python3
"""Frozen small targeted archive pilots for BLC1 and ATLAS."""
import collections,hashlib,json,math,pathlib,re,sys
from urllib.parse import urljoin,urlparse,unquote
import numpy as np
import h5py,hdf5plugin
from bs4 import BeautifulSoup
import r4_chunk_replay as base
import r9_expanded_search as eng
import r11_control_calibration as r11
S=json.loads(pathlib.Path("research/seti/r15_r16_r17_protocol.json").read_text())
def emit(k,x): print(k,json.dumps(x),flush=True)
def save(p,x): p.write_text(json.dumps(x,indent=2)+"\n")
def header(url):
    with base.RangeFile(url,8*base.MIB) as r:
        with h5py.File(r,"r") as h:
            d=h["data"];a=d.attrs
            f=float(a["fch1"]);df=float(a["foff"])
            name=a.get("source_name","")
            if isinstance(name,bytes):name=name.decode(errors="replace")
            return {"shape":list(d.shape),"source":str(name),"f0":f,"df":df,"ts":float(a["tsamp"]),"mjd":float(a["tstart"]),"bounds":sorted([f,f+(d.shape[-1]-1)*df])}
def eligible(h,fc,w):
    return len(h["shape"])==3 and 12<=h["shape"][0]<=256 and 0<abs(h["df"])*1e6<=10 and 0<h["ts"]<=60 and h["bounds"][0]<fc-w/1e6 and h["bounds"][1]>fc+w/1e6
def slice_read(u,h,fc,w,out):
    with base.RangeFile(u,S["limits"]["per_file_mib"]*base.MIB) as r:
        with h5py.File(r,"r",rdcc_nbytes=64*base.MIB) as z:
            d=z["data"];a=d.attrs
            assert list(d.shape)==h["shape"] and float(a["tstart"])==h["mjd"] and float(a["foff"])==h["df"] and float(a["fch1"])==h["f0"]
            win=base.window(d,h["f0"],h["df"],fc,w)
            if win is None:raise ValueError("window outside coverage")
            plan=base.preflight(d,[win],d.shape[0])
            if plan and plan["compressed_bytes"]+r.total+8*base.MIB>r.budget:raise ValueError("bounded spectral budget exceeded")
            lo,hi,_=win; x=np.asarray(d[:,0,lo:hi],dtype=np.float32)
            if not np.isfinite(x).all():raise ValueError("nonfinite data")
            name=hashlib.sha256((u+"|"+str(fc)).encode()).hexdigest()[:24]+".npy";np.save(out/name,x,allow_pickle=False)
            return x,{"center_mhz":fc,"channel0_mhz":h["f0"]+lo*h["df"],"df_hz":h["df"]*1e6,"ts":h["ts"],"shape":list(x.shape),"artifact":name,"slice_sha256":hashlib.sha256(x.tobytes()).hexdigest(),"bytes_read":r.total}
def blc_top(x,t,maxdrift):
    y=eng.normalized(x);n,m=y.shape;df=t["df_hz"];ts=t["ts"]
    lim=int(math.floor(maxdrift*(n-1)*ts/abs(df)));best=None
    for shift in range(-lim,lim+1):
        r=eng.scores(y,shift)
        if r is None:continue
        v,lo,ds=r
        for j in np.argsort(v)[-5:]:
            k=lo+int(j);pers=float(np.mean(y[np.arange(n),k+ds]>3))
            if pers<.5:continue
            a={"score":float(v[j]),"persistence":pers,"frequency_mhz":t["channel0_mhz"]+k*df/1e6,"drift_hz_s":shift*df/((n-1)*ts),"bin":k,"shift":shift}
            if not best or a["score"]>best["score"]:best=a
    return best
def blc_freeze(out):
    p=S["r15"];rows={};queries=[]
    for name in p["target_names"]:
        try:
            raw,url=base.fetch(eng.API,{"target":name,"telescope":"Parkes","file-type":"HDF5"},limit=16*base.MIB)
            d=json.loads(raw);rs=d if isinstance(d,list) else next((d[k] for k in ("data","files","results") if isinstance(d.get(k),list)),[])
            queries.append({"target":name,"rows":len(rs),"sha256":hashlib.sha256(raw).hexdigest()})
            for v in rs:
                if not isinstance(v,dict):continue
                u=v.get("url","")
                if not eng.allowed(u):continue
                b=unquote(base.basename(u))
                if not re.search(r"58602",b) or not re.search(r"(fine|gpuspec\.0000|rawspec\.0000)\.?(h5)$",b,re.I):continue
                label=str(v.get("target","")).upper().replace("_","").replace(" ","")
                if not label:continue
                if "PROX" in label:role="ON"
                elif "1421" in label or "1934" in label:role="CALIBRATOR"
                else:continue
                rows[u]={"url":u,"name":b,"role":role,"metadata":v}
        except Exception as e:queries.append({"target":name,"error":str(e)})
    selected=[];attempts=[]
    for role,cap in (("ON",12),("CALIBRATOR",8)):
        for f in sorted((v for v in rows.values() if v["role"]==role),key=lambda x:x["name"])[:24]:
            try:
                h=header(f["url"])
                if int(h["mjd"])!=58602 or not eligible(h,p["windows_mhz"][0],p["half_width_hz"]):raise ValueError("not original date/coverage")
                selected.append({**f,"header":h})
                if sum(v["role"]==role for v in selected)>=cap:break
            except Exception as e:attempts.append({"name":f["name"],"error":str(e)})
    m={"queries":queries,"attempts":attempts,"selected":selected};save(out/"manifest.json",m);emit("R15_MANIFEST",{"queries":queries,"attempts":attempts,"selected":[{"name":f["name"],"role":f["role"],"header":f["header"]} for f in selected]})
def atlas_freeze(out):
    p=S["r16"];raw,url=base.fetch(p["source"],limit=8*base.MIB)
    (out/"directory.html").write_bytes(raw)
    soup=BeautifulSoup(raw,"html.parser");groups=collections.defaultdict(dict)
    for a in soup.find_all("a",href=True):
        u=urljoin(url,a["href"]);b=unquote(base.basename(u))
        if urlparse(u).hostname!="bldata.berkeley.edu":continue
        m=re.search(r"^(blc\d+)_guppi_61027_\d+_DIAG_3I_ATLAS(?:_OFF)?_(\d{4})\.(.+)\.h5$",b)
        if not m or int(m[2]) not in p["scans"]:continue
        groups[(m[1],m[3])][int(m[2])]=u
    attempts=[];selected=[]
    for key,g in sorted(groups.items()):
        if set(g)!=set(p["scans"]):continue
        try:
            h=header(g[13])
            if not eligible(h,p["window_mhz"],p["half_width_hz"]):continue
            hs=[h]+[header(g[i]) for i in p["scans"][1:]]
            if not all(eligible(v,p["window_mhz"],p["half_width_hz"]) and int(v["mjd"])==p["day"] for v in hs):continue
            if len({(v["shape"][0],round(v["df"],12),round(v["ts"],7)) for v in hs})!=1:continue
            selected=[{"url":g[i],"name":unquote(base.basename(g[i])),"header":v,"role":("ON" if i%2 else "OFF")+str((i-13)//2+1),"scan":i} for i,v in zip(p["scans"],hs)]
            break
        except Exception as e:attempts.append({"backend":key,"error":str(e)})
    m={"directory_sha256":hashlib.sha256(raw).hexdigest(),"groups":len(groups),"attempts":attempts,"selected":selected};save(out/"manifest.json",m);emit("R16_MANIFEST",m)
def scan15(out):
    p=S["r15"];m=json.loads((out/"manifest.json").read_text());rows=[]
    for f in m["selected"]:
        r={"name":f["name"],"role":f["role"],"header":f["header"],"windows":[]}
        for fc in p["windows_mhz"]:
            try:
                x,t=slice_read(f["url"],f["header"],fc,p["half_width_hz"],out)
                a=blc_top(x,t,p["drift_max"]);r["windows"].append({**t,"best":a,"strong":bool(a and a["score"]>=12)})
            except Exception as e:r["windows"].append({"center_mhz":fc,"error":str(e)})
        rows.append(r);emit("R15_SCAN",r)
    summary={"selected":len(m["selected"]),"read_anchor":sum(any(v.get("center_mhz")==p["windows_mhz"][0] and "best" in v for v in r["windows"]) for r in rows),"anchor_strong":sum(any(v.get("center_mhz")==p["windows_mhz"][0] and v.get("strong") for v in r["windows"]) for r in rows),"companion_strong":sum(v.get("strong",False) for r in rows for v in r["windows"] if v["center_mhz"]!=p["windows_mhz"][0]),"status":"TARGETED_REPLAY_ONLY","confirmed_technosignatures":0}
    save(out/"result.json",{"summary":summary,"rows":rows});emit("R15_SUMMARY",summary)
def scan16(out):
    p=S["r16"];m=json.loads((out/"manifest.json").read_text());obs=[];tiles=[];ys=[];hits=[];rows=[]
    eng.SPEC=r11.S
    for f in m["selected"]:
        try:
            x,t=slice_read(f["url"],f["header"],p["window_mhz"],p["half_width_hz"],out);y=eng.normalized(x);a,c=r11.top(y,t)
            obs.append(f);tiles.append(t);ys.append(y);hits.append(a);r={"name":f["name"],"role":f["role"],"tile":t,"hits":a,"counts":c};rows.append(r);emit("R16_SCAN",r)
        except Exception as e:emit("R16_HOLD",{"name":f["name"],"error":str(e)})
    forward=[];reverse=[]
    if len(obs)==6:
        for i in (0,2):
            for h in hits[i]:forward.append(r11.candidate(h,obs[i:i+3],tiles[i:i+3],ys[i:i+3],"ON_OFF_ON","3I_ATLAS",p["window_mhz"]))
        for i in (1,3):
            for h in hits[i]:reverse.append(r11.candidate(h,obs[i:i+3],tiles[i:i+3],ys[i:i+3],"OFF_ON_OFF","3I_ATLAS",p["window_mhz"]))
    summary={"selected":len(m["selected"]),"read":len(obs),"on_alerts":sum(len(h) for f,h in zip(obs,hits) if f["role"].startswith("ON")),"off_alerts":sum(len(h) for f,h in zip(obs,hits) if f["role"].startswith("OFF")),"forward_status":dict(collections.Counter(x["status"] for x in forward)),"reverse_status":dict(collections.Counter(x["status"] for x in reverse)),"scope":"single 100 kHz L band tile, not full 1-12 GHz survey","confirmed_technosignatures":0}
    save(out/"result.json",{"summary":summary,"rows":rows,"forward":forward,"reverse":reverse});emit("R16_SUMMARY",summary)
if __name__=="__main__":
    lane=sys.argv[1];out=pathlib.Path("runtime/"+lane);out.mkdir(parents=True,exist_ok=True)
    {("r15","freeze"):blc_freeze,("r15","scan"):scan15,("r16","freeze"):atlas_freeze,("r16","scan"):scan16}[lane,sys.argv[2]](out)
