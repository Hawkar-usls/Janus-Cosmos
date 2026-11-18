#!/usr/bin/env python3
"""R37F3A: read-only astronomical profile preservation, never a signal exclusion policy."""
import ast,hashlib,io,json,os,pathlib,zipfile
import numpy as np
import requests
P=pathlib.Path("research/seti/r37f3a_natural_preservation_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r37f3a");O.mkdir(parents=True,exist_ok=True)
L=S["limits"]
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def pinned_rebin():
    path=pathlib.Path(S["implementation"]["path"]);b=path.read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
    assert blob=="e9010d9905a6100fecc719281d30581b74826d68","R37A source changed"
    funcs=[x for x in ast.parse(b).body if isinstance(x,ast.FunctionDef) and x.name in ("weights","rebin")]
    assert len(funcs)==2
    # Use only exact, pinned mathematical functions, no module initialization or memory writes.
    ns={"np":np};exec(compile(ast.Module(body=funcs,type_ignores=[]),str(path),"exec"),ns)
    return ns["rebin"],blob
def fetch():
    p=S["parent"];url=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['artifact']}/zip"
    with requests.get(url,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90,stream=True) as r:
        r.raise_for_status();b=bytearray()
        for c in r.iter_content(65536):b.extend(c);assert len(b)<16*1024**2
    assert sha(b)==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        names={pathlib.PurePosixPath(x.filename).name:x for x in z.infolist() if not x.is_dir()}
        assert len(names)==len([x for x in z.infolist() if not x.is_dir()]) and sum(x.file_size for x in names.values())<32*1024**2
        d=json.loads(z.read(names["result.json"]));assert d["line_support_count"]==2 and d["full_integrity_matches"]==2
        for c in S["controls"]:
            r=[x for x in d["sources"] if x["source"]==c["source"] and x["day"]==c["day"]]
            assert len(r)==1 and r[0]["profile"]==c["profile"] and r[0]["status"]=="CONDITIONAL_NATURAL_LINE_SUPPORT"
            b=z.read(names[c["profile"]])
            with np.load(io.BytesIO(b),allow_pickle=False) as a:
                v=np.asarray(a["velocity_km_s"],float);y=np.asarray(a["normalized_mean"],float)
            assert v.ndim==y.ndim==1 and v.shape==y.shape and 100<len(v)<10000
            assert np.isfinite(v).all() and np.isfinite(y).all() and (y>=0).all()
            step=float(np.median(np.diff(v)));assert step>0 and np.allclose(np.diff(v),step,rtol=1e-6,atol=1e-10)
            e=np.r_[v[0]-step/2,(v[:-1]+v[1:])/2,v[-1]+step/2]
            yield c,r[0],v,y,e,step,sha(b)
def validate(e,y):
    assert e.ndim==y.ndim==1 and len(e)==len(y)+1
    assert np.isfinite(e).all() and np.isfinite(y).all() and (y>=0).all() and (np.diff(e)>0).all()
def metrics(e,y,bg,line):
    validate(e,y);lo,hi=line;assert e[0]<=lo<hi<=e[-1]
    a=np.maximum(e[:-1],lo);b=np.minimum(e[1:],hi);w=np.maximum(0,b-a)
    assert abs(w.sum()-(hi-lo))<1e-8
    x=y-bg;positive=np.maximum(0,x);mass=float(positive@w)
    area=float(x@w);assert mass>0 and area>0
    return {"area":area,"centroid":float(positive@(w*(a+b)/2)/mass),"peak":float(x[w>0].max()),"positive_mass":mass}
def compare(e,y,f,z,bg,line):
    try:
        validate(e,y);validate(f,z)
        assert abs(e[0]-f[0])<1e-9 and abs(e[-1]-f[-1])<1e-9,"incomplete physical support"
        a=metrics(e,y,bg,line);b=metrics(f,z,bg,line)
        lo,hi=line;g=np.unique(np.r_[lo,e[(e>lo)&(e<hi)],f[(f>lo)&(f<hi)],hi]);mid=(g[:-1]+g[1:])/2
        iy=np.searchsorted(e,mid,side="right")-1;iz=np.searchsorted(f,mid,side="right")-1
        w=np.diff(g);rel=float(np.abs(y[iy]-z[iz])@w/max(float(np.abs(y[iy]-bg)@w),1e-15))
        m={"area_rel":abs(b["area"]-a["area"])/a["area"],"centroid_km_s":abs(b["centroid"]-a["centroid"]),"peak_rel":abs(b["peak"]-a["peak"])/a["peak"],"profile_l1_rel":rel}
        failed=[k for k,v in m.items() if not np.isfinite(v) or v>L[k]]
        return {"status":"PRESERVED" if not failed else "REVIEW","metrics":m,"failed":failed}
    except (AssertionError,ValueError,FloatingPointError) as err:return {"status":"REVIEW","error":str(err) or type(err).__name__}
def main():
    rebin,blob=pinned_rebin();rows=[];canaries=[];immutable=[]
    for c,r,v,y,e,step,receipt in fetch():
        src=c["source"];lo,hi=c["line"];bg=float(np.median(y[(v<lo-10)|(v>hi+10)]))
        assert np.isfinite(bg) and bg>0
        original=sha(y.tobytes());y.setflags(write=False)
        ref=metrics(e,y,bg,c["line"])
        for t in S["transforms"]:
            delta=step*t["factor"];start=e[0]+delta*t["phase"]
            inner=np.arange(start,e[-1],delta);inner=inner[(inner>e[0]+1e-9)&(inner<e[-1]-1e-9)]
            f=np.r_[e[0],inner,e[-1]]
            x,err=rebin(y[None,:],np.array([0.,1.]),e,np.array([0.,1.]),f)
            out=compare(e,y,f,x[0],bg,c["line"])
            assert err<=L["conservation_rel"]
            rows.append({"source":src,"day":c["day"],"transform":t,"input_profile_sha256":receipt,"reference":ref,"whole_support_conservation_rel":float(err),**out})
        mask=(v>=lo)&(v<=hi);peak=v[np.flatnonzero(mask)[np.argmax(y[mask])]]
        variants={"unchanged":(e,y),"attenuation_50pct":(e,bg+(y-bg)*.5),"peak_notch":(e,np.where(abs(v-peak)<=.2,bg,y)),"shift_0.5km_s":(e,np.interp(v-.5,v,y,left=bg,right=bg)),"nonfinite":(e,np.where(abs(v-peak)<=.02,np.nan,y)),"incomplete_support":(e[1:],y[1:])}
        for name in S["canaries"]:
            f,z=variants[name];res=compare(e,y,f,z,bg,c["line"]);expected="PRESERVED" if name=="unchanged" else "REVIEW"
            canaries.append({"source":src,"name":name,"expected":expected,"correct":res["status"]==expected,**res})
        assert sha(y.tobytes())==original;immutable.append({"source":src,"original_unchanged":True,"input_profile_sha256":receipt})
    assert len(rows)==16 and len(canaries)==12
    passed=sum(x["status"]=="PRESERVED" for x in rows);canary=sum(x["correct"] for x in canaries)
    d={"status":"PASS_ISOLATED_NATURAL_PRESERVATION" if passed==16 and canary==12 else "HOLD_TRANSFORM_OR_CANARY_PRESERVATION","transforms":rows,"canaries":canaries,"preserved":passed,"transform_count":16,"canaries_correct":canary,"canary_count":12,"immutable":immutable,"r37a_git_blob":blob,"real_input_profiles":2,"input_science_bin_hz":715.2557373046875,"native_fine_grid_tested":False,"end_to_end_hannah_tested":False,"spatial_off_verified":False,"exact_bl_frame_verified":False,"production_adoption":False}
    save("result.json",d);print("R37F3A_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
