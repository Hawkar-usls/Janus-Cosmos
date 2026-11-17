#!/usr/bin/env python3
"""Bounded exact SIGPROC and fixed multi-band ATLAS replay."""
import collections,hashlib,json,math,pathlib,struct,sys,re
from urllib.parse import urljoin,urlparse,unquote
import numpy as np
from bs4 import BeautifulSoup
import r4_chunk_replay as base
import r9_expanded_search as eng
import r11_control_calibration as r11
import r15_r16_archive_pilots as old
import r15b_r16b_followup as prior
S=json.loads(pathlib.Path("research/seti/r15c_r16c_r17c_protocol.json").read_text())
def emit(k,x):print(k,json.dumps(x),flush=True)
def save(p,x):p.write_text(json.dumps(x,indent=2)+"\n")
STR={"rawdatafile","source_name"}
INT={"machine_id","telescope_id","data_type","barycentric","pulsarcentric","nbits","nifs","nchans","nbeams","ibeam","nsamples","signed"}
DBL={"az_start","za_start","src_raj","src_dej","tstart","tsamp","fch1","foff","refdm","period"}
def exact(r,n):
    b=r.read(n)
    if len(b)!=n:raise ValueError("short SIGPROC read")
    return b
def string(r):
    n=struct.unpack("<i",exact(r,4))[0]
    if not 0<n<=4096:raise ValueError("SIGPROC string length")
    return exact(r,n).decode("ascii")
def parse_header(r):
    if string(r)!="HEADER_START":raise ValueError("missing HEADER_START")
    h={}
    for _ in range(128):
        k=string(r)
        if k=="HEADER_END":break
        if k in h:raise ValueError("duplicate header field")
        if k in STR:h[k]=string(r)
        elif k in INT:h[k]=struct.unpack("<i",exact(r,4))[0]
        elif k in DBL:h[k]=struct.unpack("<d",exact(r,8))[0]
        else:raise ValueError("unsupported SIGPROC key "+k)
    else:raise ValueError("missing HEADER_END")
    offset=r.tell()
    for k in ("nchans","nifs","nbits","fch1","foff","tsamp","tstart"):assert k in h
    if h["nifs"]!=1 or h["nbits"]!=32 or not 0<h["nchans"]<=1000000000:raise ValueError("unsupported sample layout")
    if not all(math.isfinite(h[k]) for k in ("fch1","foff","tsamp","tstart")):raise ValueError("nonfinite header")
    if not 0<abs(h["foff"])*1e6<=10 or not 0<h["tsamp"]<=60:raise ValueError("not high-resolution time/frequency product")
    row=h["nchans"]*4
    count,rem=divmod(r.size-offset,row)
    if rem or not 12<=count<=512:raise ValueError("invalid payload dimensions")
    h.update(offset=offset,rows=count,size=r.size)
    return h
def fil_window(r,h,fc,half):
    i=round((fc-h["fch1"])/h["foff"]);delta=math.ceil(half/abs(h["foff"]*1e6))
    lo=i-delta;hi=i+delta+1
    if lo<0 or hi>h["nchans"]:raise ValueError("window outside band")
    x=np.zeros((h["rows"],hi-lo),np.float32)
    for j in range(h["rows"]):
        r.seek(h["offset"]+4*(j*h["nchans"]+lo))
        x[j]=np.frombuffer(exact(r,4*(hi-lo)),dtype="<f4")
    if not np.isfinite(x).all():raise ValueError("nonfinite samples")
    return x,{"channel0_mhz":h["fch1"]+lo*h["foff"],"df_hz":h["foff"]*1e6,"ts":h["tsamp"],"shape":list(x.shape),"slice_sha256":hashlib.sha256(x.tobytes()).hexdigest()}
def freeze15(out):
    p=S["r15c"];z,n=prior.archive(p["artifact"]);d=json.loads(z.read(n["manifest.json"]));z.close()
    rows=[]
    for f in d["exact_records"]:
        name=f["name"]
        if not name.startswith("blc03_") or not name.endswith(".fil") or "_ND_" in name:continue
        if re.match(r"blc03_ProxCen_S_",name,re.I):role="ON"
        elif re.match(r"blc03_1421-490_[SR]_",name,re.I):role="CONTROL"
        else:continue
        rows.append({**f,"role":role})
    rows.sort(key=lambda f:(f["name"]!=p["published_anchor"],f["name"]))
    chosen=[];holds=[]
    for f in rows:
        if sum(a["role"]==f["role"] for a in chosen)>=(9 if f["role"]=="ON" else 12):continue
        try:
            with base.RangeFile(f["url"],1024*1024) as r:
                h=parse_header(r)
                assert int(h["tstart"])==58602
                assert int(f["metadata"]["size"])==r.size
                lo,hi=sorted((h["fch1"],h["fch1"]+(h["nchans"]-1)*h["foff"]))
                assert lo<p["center_mhz"]<hi
            chosen.append({**f,"header":h})
        except Exception as e:holds.append({"name":f["name"],"error":str(e)})
    m={"indexed":len(rows),"selected":chosen,"holds":holds,"published_anchor":p["published_anchor"]}
    save(out/"manifest.json",m);emit("R15C_MANIFEST",{"indexed":len(rows),"selected":[{"name":f["name"],"role":f["role"],"header":f["header"]} for f in chosen],"holds":holds})
def scan15(out):
    p=S["r15c"];m=json.loads((out/"manifest.json").read_text());rows=[]
    for f in m["selected"]:
        row={"name":f["name"],"role":f["role"]}
        try:
            with base.RangeFile(f["url"],8*base.MIB) as r:
                h=parse_header(r);assert h==f["header"]
                x,t=fil_window(r,h,p["center_mhz"],p["half_hz"])
                a=old.blc_top(x,t,p["drift_max_hz_s"])
                path=hashlib.sha256(f["url"].encode()).hexdigest()[:24]+".npy";np.save(out/path,x,allow_pickle=False)
                row.update(tile=t,header=h,best=a,strong=bool(a and a["score"]>=12),artifact=path,bytes_read=r.total)
        except Exception as e:row["error"]=str(e)
        rows.append(row);emit("R15C_SCAN",row)
    anchor=next((r for r in rows if r["name"]==p["published_anchor"]),{})
    summary={"selected":len(m["selected"]),"read":sum("tile" in r for r in rows),"anchor_read":"tile" in anchor,"anchor_strong":anchor.get("strong",False),"anchor":anchor,"on_read":sum(r["role"]=="ON" and "tile" in r for r in rows),"control_read":sum(r["role"]=="CONTROL" and "tile" in r for r in rows),"status":"TARGETED_REPLAY_NOT_ORIGIN_CLASSIFICATION","confirmed_technosignatures":0}
    save(out/"result.json",{"summary":summary,"rows":rows});emit("R15C_SUMMARY",summary)
def freeze16(out):
    p=S["r16c"];raw,url=base.fetch("https://bldata.berkeley.edu/ATLAS/GB_ATLAS/",limit=8*base.MIB)
    (out/"directory.html").write_bytes(raw)
    groups=collections.defaultdict(dict)
    for a in BeautifulSoup(raw,"html.parser").find_all("a",href=True):
        u=urljoin(url,a["href"]);b=unquote(base.basename(u))
        if urlparse(u).hostname!="bldata.berkeley.edu":continue
        m=re.search(r"^(blc\d+)_guppi_61027_\d+_DIAG_3I_ATLAS(_OFF)?_(\d{4})\.(.+)\.h5$",b)
        if m:groups[(m[1],m[4])][int(m[3])]=(u,bool(m[2]))
    selected=[];holds=[];cache={}
    for fc in p["centers_mhz"]:
        found=False
        for key,g in sorted(groups.items()):
            for start in sorted(g):
                seq=list(range(start,start+6))
                if not all(i in g for i in seq) or [g[i][1] for i in seq]!=[False,True]*3:continue
                try:
                    uu=[g[i][0] for i in seq]
                    if uu[0] not in cache:cache[uu[0]]=old.header(uu[0])
                    if not old.eligible(cache[uu[0]],fc,p["half_hz"]):continue
                    for u in uu[1:]:
                        if u not in cache:cache[u]=old.header(u)
                    hs=[cache[u] for u in uu]
                    if not all(old.eligible(h,fc,p["half_hz"]) and int(h["mjd"])==61027 for h in hs):continue
                    if len({(h["shape"][0],round(h["df"],12),round(h["ts"],7)) for h in hs})!=1:continue
                    if not all(0<(b["mjd"]-a["mjd"])*86400<=600 for a,b in zip(hs,hs[1:])):continue
                    selected.append({"center_mhz":fc,"backend":key,"files":[{"url":u,"name":unquote(base.basename(u)),"header":h,"role":("ON" if j%2==0 else "OFF")+str(j//2+1)} for j,(u,h) in enumerate(zip(uu,hs))]})
                    found=True;break
                except Exception as e:holds.append({"center_mhz":fc,"backend":key,"error":str(e)})
            if found:break
    m={"directory_sha256":hashlib.sha256(raw).hexdigest(),"selected":selected,"unsupported":[x for x in p["centers_mhz"] if x not in [v["center_mhz"] for v in selected]],"holds":holds}
    save(out/"manifest.json",m);emit("R16C_MANIFEST",m)
def scan16(out):
    p=S["r16c"];m=json.loads((out/"manifest.json").read_text());results=[]
    for w in m["selected"]:
        fc=w["center_mhz"];rows=[];ys=[];ts=[];hits=[]
        try:
            for f in w["files"]:
                x,t=old.slice_read(f["url"],f["header"],fc,p["half_hz"],out)
                y=eng.normalized(x);a,c=r11.top(y,t);ys.append(y);ts.append(t);hits.append(a);rows.append({"role":f["role"],"name":f["name"],"tile":t,"hits":a,"counts":c})
            forward=[];reverse=[]
            for i in (0,2):
                for h in hits[i]:forward.append(r11.candidate(h,w["files"][i:i+3],ts[i:i+3],ys[i:i+3],"ON_OFF_ON","3I_ATLAS",fc))
            for i in (1,3):
                for h in hits[i]:reverse.append(r11.candidate(h,w["files"][i:i+3],ts[i:i+3],ys[i:i+3],"OFF_ON_OFF","3I_ATLAS",fc))
            r={"center_mhz":fc,"status":"COMPLETE","rows":rows,"forward":forward,"reverse":reverse,"forward_counts":dict(collections.Counter(h["status"] for h in forward)),"reverse_counts":dict(collections.Counter(h["status"] for h in reverse)),"on_alerts":sum(len(h) for h in hits[::2]),"off_alerts":sum(len(h) for h in hits[1::2])}
        except Exception as e:r={"center_mhz":fc,"status":"HOLD","rows":rows,"error":str(e)}
        results.append(r);save(out/"partial.json",results);emit("R16C_WINDOW",{k:v for k,v in r.items() if k not in ("rows","forward","reverse")})
    good=[r for r in results if r["status"]=="COMPLETE"]
    s={"planned":len(p["centers_mhz"]),"selected":len(m["selected"]),"complete":len(good),"unsupported":m["unsupported"],"scan_tiles":6*len(good),"on_alerts":sum(r["on_alerts"] for r in good),"off_alerts":sum(r["off_alerts"] for r in good),"forward":dict(collections.Counter(h["status"] for r in good for h in r["forward"])),"reverse":dict(collections.Counter(h["status"] for r in good for h in r["reverse"])),"scope":"Disjoint fixed 100 kHz windows only; same archive day; top-five cap","confirmed_technosignatures":0}
    save(out/"result.json",{"summary":s,"windows":results});emit("R16C_SUMMARY",s)
if __name__=="__main__":
    lane,phase=sys.argv[1:];out=pathlib.Path("runtime/"+lane);out.mkdir(parents=True,exist_ok=True)
    {("r15c","freeze"):freeze15,("r15c","scan"):scan15,("r16c","freeze"):freeze16,("r16c","scan"):scan16}[lane,phase](out)
