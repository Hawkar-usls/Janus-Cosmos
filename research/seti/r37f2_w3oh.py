#!/usr/bin/env python3
"""R37F2: bounded natural astronomical OH line validation only."""
import collections,hashlib,io,json,math,os,pathlib,re,urllib.parse,zipfile
import numpy as np
import h5py,hdf5plugin
import requests
import r4_chunk_replay as base
import r37b_source_header as geom
P=pathlib.Path("research/seti/r37f2_w3oh_protocol.json")
S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f2");O.mkdir(parents=True,exist_ok=True)
C=S["source"];L=S["selection"];M=S["measurement"]
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,x):(O/n).write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
def parent():
    p=S["parent"]
    r=requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['artifact']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=(15,90))
    r.raise_for_status();b=r.content;assert len(b)<16*1024**2 and sha(b)==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        n=[i for i in z.infolist() if pathlib.PurePosixPath(i.filename).name=="natural.json"]
        assert len(n)==1 and n[0].file_size<8*1024**2
        d=json.loads(z.read(n[0]))
    pairs=[v for v in d["headers"] if v["status"]=="HEADER_VERIFIED" and v["field"]=="W3OH" and v["comparable_grid"] and v["oh_coverage"]]
    assert len(pairs)==2 and {v["role"] for v in pairs}=={"ON","OFF"}
    return sorted(pairs,key=lambda v:v["role"]!="ON")
def catalog():
    b,u=base.fetch(L["query"],{"target":C["target_alias"],"telescope":L["telescope"],"file-type":"HDF5"},limit=16*base.MIB)
    assert urllib.parse.urlparse(u).hostname=="seti.berkeley.edu"
    d=json.loads(b);rr=d if isinstance(d,list) else next(d[k] for k in ("data","files","results") if isinstance(d.get(k),list))
    files={}
    for r in rr:
        if not isinstance(r,dict):continue
        u=r.get("url");t=str(r.get("target","")).strip().upper()
        if not isinstance(u,str) or t not in ("DIAG_W3OH","DIAG_W3OH_OFF","W3OH","W3OH_OFF"):continue
        p=urllib.parse.urlparse(u);ext=pathlib.PurePosixPath(p.path).suffix.lower()
        if ext not in (".h5",".fil",".raw"):continue
        try:
            mjd=float(r["mjd"]);assert 50000<=mjd<=65000
        except (ValueError,TypeError,AssertionError,KeyError):continue
        f={"url":u,"target":t,"role":"OFF" if t.endswith("_OFF") else "ON","mjd":mjd,"day":int(mjd),"ext":ext,"metadata":r}
        if u in files and files[u]!=f:raise ValueError("conflicting catalogue rows")
        files[u]=f
    vals=list(files.values())
    summary={"sha256":sha(b),"rows":len(rr),"unique":len(vals),"dates_by_format":{x:sorted({v["day"] for v in vals if v["ext"]==x}) for x in (".h5",".fil",".raw")},"counts_by_format":dict(collections.Counter(v["ext"] for v in vals))}
    return vals,summary
def hdr(f):
    assert base.archive_url(f["url"])
    import r9_expanded_search as eng
    h=eng.hdr(f)
    assert str(h["source"]).strip().upper()==f["target"]
    assert abs(h["mjd"]-f["mjd"])<=3/86400
    assert len(h["shape"])==3 and 2<=h["shape"][0]<=L["max_rows"] and h["shape"][1]>=1
    assert 0<abs(h["df"])*1e6<=10 and 0<h["ts"]<=60
    p=geom.pointing(h,f["metadata"])
    return h,p
def qualify(a,b):
    assert a["role"]=="ON" and b["role"]=="OFF" and a["day"]==b["day"]
    gap=(b["mjd"]-a["mjd"])*86400;assert 0<gap<=L["max_gap_s"]
    h,p=hdr(a);g,q=hdr(b);ref=(C["coordinate_ra_deg"],C["coordinate_dec_deg"])
    sep_on=geom.sep(p,ref);sep_off=geom.sep(q,ref);sep_pair=geom.sep(p,q)
    assert sep_on<=C["max_on_pointing_separation_deg"],"ON pointing not close to literature coordinate"
    assert sep_off>=C["min_off_pointing_separation_deg"] and sep_pair>=C["min_off_pointing_separation_deg"],"OFF pointing insufficient"
    assert abs(h["df"]-g["df"])<1e-12 and abs(h["ts"]-g["ts"])<1e-7
    delta=(h["f0"]-g["f0"])/h["df"];assert abs(delta-round(delta))<.001,"grid alignment"
    for f in C["line_rest_mhz"]:
        lo,hi=window(f)
        assert all(v["bounds"][0]<=lo<hi<=v["bounds"][1] for v in (h,g))
    return {"day":a["day"],"on":a["url"],"off":b["url"],"headers":[h,g],"pointing_separations_deg":{"on_to_literature":sep_on,"off_to_literature":sep_off,"on_to_off":sep_pair},"gap_s":gap,"status":"HEADER_QUALIFIED"}
def window(rest):
    c=299792.458;v0,v1=C["topocentric_velocity_search_km_s"]
    return rest*(1-v1/c),rest*(1-v0/c)
def select(vals,pa):
    by={v["url"]:v for v in vals}
    chosen=[];attempts=[]
    def check(a,b,primary=False):
        try:
            q=qualify(a,b);q["primary"]=primary;chosen.append(q);attempts.append({"day":a["day"],"status":"HEADER_QUALIFIED"})
        except Exception as e:
            attempts.append({"day":a["day"],"on":a["url"],"off":b["url"],"error":type(e).__name__+": "+str(e)[:240]})
    try:check(by[pa[0]["url"]],by[pa[1]["url"]],True)
    except KeyError:attempts.append({"status":"HOLD_PARENT_CATALOGUE_JOIN"})
    fs=[v for v in vals if v["ext"]==".h5" and base.archive_url(v["url"]) and v["url"].endswith(".gpuspec.0000.h5")]
    for day in sorted({v["day"] for v in fs}):
        if day==int(pa[0]["mjd"]):continue
        if len(chosen)>=L["max_pairs"]:break
        on=sorted([v for v in fs if v["day"]==day and v["role"]=="ON"],key=lambda v:(v["mjd"],v["url"]))
        off=sorted([v for v in fs if v["day"]==day and v["role"]=="OFF"],key=lambda v:(v["mjd"],v["url"]))
        for a in on:
            if len(attempts)>=L["max_header_candidates"]:break
            bs=[b for b in off if 0<(b["mjd"]-a["mjd"])*86400<=L["max_gap_s"]]
            if bs:
                check(a,bs[0])
                if any(v["day"]==day for v in chosen):break
    return chosen,attempts
def bounds(h,lo,hi):
    f0=h["f0"];df=h["df"];a,b=sorted(((lo-f0)/df,(hi-f0)/df))
    return max(0,math.ceil(a-1e-7)),min(h["shape"][-1],math.floor(b+1e-7)+1)
def read_pair(q):
    hs=q["headers"];windows=[]
    for rest in C["line_rest_mhz"]:
        lo,hi=window(rest);a,b=bounds(hs[0],lo,hi)
        # Use exact overlapping grid indices, not interpolation.
        shift=round((hs[0]["f0"]-hs[1]["f0"])/hs[0]["df"])
        a=max(a,-shift);b=min(b,hs[1]["shape"][-1]-shift)
        b=a+((b-a)//64)*64;assert b-a>=64*8
        windows.append({"rest":rest,"indices":[[a,b],[a+shift,b+shift]]})
    arrays=[];receipts=[]
    for k,url in enumerate((q["on"],q["off"])):
        h=hs[k]
        with base.RangeFile(url,L["max_bytes_per_file"]) as remote:
            with h5py.File(remote,"r",rdcc_nbytes=128*base.MIB) as f:
                d=f["data"];assert list(d.shape)==h["shape"]
                assert abs(float(d.attrs["fch1"])-h["f0"])<1e-10 and abs(float(d.attrs["foff"])-h["df"])<1e-12
                ww=[(*v["indices"][k],0) for v in windows]
                plan=base.preflight(d,ww,d.shape[0])
                if plan and plan["compressed_bytes"]+remote.total+8*base.MIB>L["max_bytes_per_file"]:raise ValueError("compressed plan exceeds per-file budget")
                aa=[]
                for v in windows:
                    lo,hi=v["indices"][k];x=np.asarray(d[:,0,lo:hi],dtype=np.float64)
                    assert np.isfinite(x).all() and (x>=0).all()
                    med=np.median(x,axis=1);assert (med>0).all()
                    y=(x/med[:,None]).reshape(x.shape[0],-1,64).mean(2)
                    aa.append(y)
                arrays.append(aa);receipts.append({"bytes_read":remote.total,"file_bytes":remote.size,"chunk_plan":plan,"rows":int(d.shape[0])})
    report=[]
    for ix,v in enumerate(windows):
        on,off=arrays[0][ix],arrays[1][ix]
        diff=on.mean(0)-off.mean(0);center=float(np.median(diff));sigma=1.4826*float(np.median(abs(diff-center)))
        assert sigma>0 and math.isfinite(sigma)
        i=int(np.argmax(diff));z=float((diff[i]-center)/sigma)
        excess=float(diff[i]);positive=int(np.sum(on[:,i]>off[:,i].mean()))
        adjacent=bool(i>0 and diff[i-1]>center+3*sigma or i+1<len(diff) and diff[i+1]>center+3*sigma)
        ok=z>=M["minimum_descriptive_z"] and excess>=M["minimum_on_relative_excess"] and positive==len(on) and adjacent
        report.append({"rest_transition_mhz":v["rest"],"descriptive_peak_z":z,"relative_on_minus_off":excess,"on_rows_positive":positive,"on_rows":len(on),"adjacent_support":adjacent,"blocks":len(diff),"criterion_met":bool(ok)})
    return {"day":q["day"],"status":"NATURAL_LINE_CANDIDATE" if all(v["criterion_met"] for v in report) else "HOLD_LINE_CRITERIA","lines":report,"io":receipts,"scope":"Descriptive ON/OFF OH window excess; exact velocity frame not established."}
def main():
    pa=parent();vals,summary=catalog();pairs,attempts=select(vals,pa)
    save("frozen_manifest.json",{"catalogue":summary,"pairs":pairs,"attempts":attempts,"frozen_before_power_read":True})
    print("R37F2_METADATA",json.dumps({"catalogue":summary,"qualified_pairs":len(pairs),"attempts":attempts}),flush=True)
    results=[]
    for q in pairs:
        try:results.append(read_pair(q))
        except Exception as e:results.append({"day":q["day"],"status":"HOLD_POWER","error":type(e).__name__+": "+str(e)[:300]})
    n=sum(v["status"]=="NATURAL_LINE_CANDIDATE" for v in results)
    d={"status":"NATURAL_CONTROL_DEVELOPMENT_CANDIDATE" if n else "HOLD_NATURAL_LINE_VALIDATION","catalogue":summary,"qualified_pairs":len(pairs),"independent_qualified_dates":len({v["day"] for v in pairs}),"candidate_dates":n,"results":results,"actual_natural_origin_confirmed":False,"independent_scientific_confirmation":False,"exact_lsr_correction_verified":False,"production_adoption":False}
    save("result.json",d);print("R37F2_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
