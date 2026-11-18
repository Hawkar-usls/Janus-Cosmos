#!/usr/bin/env python3
"""Isolated astronomy archive authority/replay; display derivatives never authoritative."""
import dataclasses,hashlib,json,pathlib,stat
import numpy as np
P=pathlib.Path("research/seti/r37f3b1_native_authority_protocol.json");S=json.loads(P.read_text())
b=(P.parent/"r37f3b_native_preservation.py").read_bytes()
assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==S["pinned_native_blob"]
import r37f3b_native_preservation as prev
science=prev.science;guard=prev.guard
from astropy.utils import iers
from astropy.coordinates import solar_system_ephemeris
O=pathlib.Path("runtime/r37f3b1");O.mkdir(parents=True,exist_ok=True)
@dataclasses.dataclass(frozen=True)
class Ref:
    role:str
    source:str
    day:int
    path:pathlib.Path
    sha256:str
    size:int
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(4*1024**2),b""):h.update(b)
    return h.hexdigest()
def authorize(ref,src,selection):
    if ref.role!="ORIGINAL_ARCHIVE":raise ValueError("display-only or unknown role")
    if (ref.source,ref.day)!=(src["name"],src["day"]):raise ValueError("observation mismatch")
    if ref.sha256!=src["expected_sha256"] or ref.size!=selection["reference"]["size"]:raise ValueError("reference mismatch")
    if ref.path.is_symlink() or not ref.path.is_file():raise ValueError("archive path invalid")
    if ref.path.stat().st_size!=ref.size or digest(ref.path)!=ref.sha256:raise ValueError("bytes do not match reference")
    if ref.path.stat().st_mode & (stat.S_IWUSR|stat.S_IWGRP|stat.S_IWOTH):raise ValueError("local original not read-only")
def measure(ref,src,selected,frozen):
    authorize(ref,src,selected)
    r=science.inspect(ref.path,selected,src)
    assert r["status"]=="CONDITIONAL_NATURAL_LINE_SUPPORT"
    with np.load(science.O/r["profile"],allow_pickle=False) as a:
        coarse_exact=bool(np.array_equal(a["velocity_km_s"],frozen[0]) and np.array_equal(a["normalized_mean"],frozen[1]))
    vc,x,e,step=prev.native(ref.path,selected,src,r)
    metrics=[]
    lo,hi=src["lsr_window"]
    for group,part in zip(prev.S["groups"],[x,x[:3],x[3:]]):
        y=part.mean(0);bg=float(np.median(y[(vc<lo-10)|(vc>hi+10)]));metrics.append({"group":group,**guard.metrics(e,y,bg,[lo,hi])})
    native_hash=hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()
    assert digest(ref.path)==ref.sha256
    return {"coarse_reference_exact":coarse_exact,"native_summary":metrics,"native_replay_hash":native_hash,"rows":len(x),"native_channels":len(vc)},(vc,x,e,step)
def rejected(ref,src,selected):
    try:authorize(ref,src,selected)
    except ValueError as exc:return {"rejected":True,"reason":str(exc)}
    return {"rejected":False,"reason":"unexpected authorization"}
def main():
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
    refs=science.refs();frozen={c["source"]:(v,y) for c,r,v,y,e,step,receipt in guard.fetch()}
    rows=[]
    for s in science.S["sources"]:
        src={**s,"expected_sha256":prev.S["expected_sha256"][s["name"]]};selected=science.choose(src,refs)
        p=O/"original.fil";fake=O/"forged_view.bin";row={"source":src["name"],"day":src["day"],"status":"HOLD"}
        try:
            row["integrity"]=science.download(selected,src,p);p.chmod(0o444)
            ref=Ref("ORIGINAL_ARCHIVE",src["name"],src["day"],p,src["expected_sha256"],selected["reference"]["size"])
            first,(vc,x,e,step)=measure(ref,src,selected,frozen[src["name"]])
            views=[];checks=[]
            y=x.mean(0);lo,hi=src["lsr_window"];bg=float(np.median(y[(vc<lo-10)|(vc>hi+10)]))
            for t in prev.S["transforms"]:
                f=prev.grid(e,step,t);z,err=prev.fast_rebin(y,e,f)
                verdict=guard.compare(e,y,f,z,bg,[lo,hi])
                # Derivative remains ephemeral, no sample arrays published.
                ch=rejected(dataclasses.replace(ref,role="DISPLAY_ONLY"),src,selected)
                checks.append({"name":"display_role",**ch})
                views.append({"transform":t,"science_authority":False,"fidelity":verdict["status"],"conservation_rel":err})
            del x,y,vc,e
            fake.write_bytes(b"derived display placeholder");fake.chmod(0o444)
            for name,bad in [("forged_original_role",dataclasses.replace(ref,path=fake)),("wrong_observation",dataclasses.replace(ref,day=ref.day+1)),("wrong_digest",dataclasses.replace(ref,sha256="0"*64))]:
                checks.append({"name":name,**rejected(bad,src,selected)})
            second,_=measure(ref,src,selected,frozen[src["name"]])
            # Native replay hashes compared internally, not exported.
            same=first==second;orig=digest(p)==ref.sha256
            row.update(first_reference_exact=first["coarse_reference_exact"],second_reference_exact=second["coarse_reference_exact"],native_replay_exact=same,original_hash_unchanged=orig,views=views,provenance_checks=checks,native_channels=first["native_channels"],rows=first["rows"])
            if same and orig and first["coarse_reference_exact"] and second["coarse_reference_exact"] and len(views)==8 and all(c["rejected"] for c in checks):
                row["status"]="PASS_ORIGINAL_AUTHORITY"
        except Exception as exc:
            import traceback
            traceback.print_exc();row["error"]=type(exc).__name__+": "+str(exc)[:300]
        finally:p.unlink(missing_ok=True);fake.unlink(missing_ok=True)
        rows.append(row);print("R37F3B1_SOURCE",json.dumps(row),flush=True)
    d={"status":"PASS_ISOLATED_ORIGINAL_AUTHORITY" if len(rows)==2 and all(r["status"]=="PASS_ORIGINAL_AUTHORITY" for r in rows) else "HOLD_ORIGINAL_AUTHORITY","sources":rows,"original_authority_pass":sum(r["status"]=="PASS_ORIGINAL_AUTHORITY" for r in rows),"views":sum(len(r.get("views",[])) for r in rows),"negative_checks":sum(len(r.get("provenance_checks",[])) for r in rows),"negative_rejected":sum(c["rejected"] for r in rows for c in r.get("provenance_checks",[])),"parent_native_coarsening_pass":4,"parent_native_coarsening_total":48,"persistent_local_archive":False,"end_to_end_hannah_tested":False,"production_adoption":False}
    (O/"result.json").write_text(json.dumps(d,indent=2,allow_nan=False)+"\n");print("R37F3B1_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
