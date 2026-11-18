#!/usr/bin/env python3
"""Paired read-only astronomy science-report regression, no candidate processing."""
import hashlib,io,json,os,pathlib,time,zipfile
import requests
P=pathlib.Path("research/seti/r37f3cb1_paired_science_protocol.json");S=json.loads(P.read_text())
O=pathlib.Path("runtime/r37f3cb1");O.mkdir(parents=True,exist_ok=True)
def blob(p):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def check_pins():
    for n,h in S["pinned_blobs"].items():assert blob(P.parent/n)==h,n
check_pins()
import r37f3ca_science_provenance as ca
import r37f3cb_science_sidecar as cb
auth=ca.auth;science=auth.science
from astropy.utils import iers
from astropy.coordinates import solar_system_ephemeris
def parent():
    p=S["parent"]
    with requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['artifact']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90,stream=True) as r:
        r.raise_for_status();b=bytearray()
        for x in r.iter_content(65536):b.extend(x);assert len(b)<16*1024**2
    assert hashlib.sha256(b).hexdigest()==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        ns=[n for n in z.infolist() if pathlib.PurePosixPath(n.filename).name=="result.json"]
        assert len(ns)==1 and ns[0].file_size<8*1024**2
        d=json.loads(z.read(ns[0]))
    assert d["status"]=="PASS_ISOLATED_SCIENCE_REPORT_BOUNDARY" and d["g4_unchanged"] and d["report_references"]==2
def armRun(arm,ledger,src,sel,ref,frozen,raw,expected):
    start=time.perf_counter()
    hook=cb.ObservationAccumulator(ledger,{"id":cb.S["parent"]["artifact"],"digest":cb.S["parent"]["sha256"]})
    side=None
    if arm=="B":
        side=cb.ScienceReportSidecar(hook)
        assert side.attach(raw,expected)=="ATTACHED"
    m,arrays=auth.measure(ref,src,sel,frozen);del arrays
    env=json.loads(raw)
    # Compare historical report's measurements only; same-runtime equivalence
    # additionally checks the internal native replay digest, never exported.
    fresh=ca.report(src,sel,ref,m,env["report"]["method"])
    historic=ca.verify(env,fresh)
    # Baseline and interface arms both check frozen original report against fresh measurement.
    assert historic and m["coarse_reference_exact"]
    refs=0 if side is None else len(side.export()["references"])
    assert refs==(1 if arm=="B" else 0)
    final,delta,replay=hook.seal()
    unchanged=ca.canonical(final)==ca.canonical(ledger) and delta["added"]==0 and replay["added"]==0
    assert unchanged and not hook.records and not hook.observations and not hook.seen
    return m,{"arm":arm,"seconds":time.perf_counter()-start,"frozen_report_exact":historic,"g4_unchanged":unchanged,"reference_count":refs,"candidate_records_added":0}
def main():
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
    parent();cb.pinned();raws=cb.fetch();rref=science.refs()
    frozen={c["source"]:(v,y) for c,r,v,y,e,step,receipt in auth.guard.fetch()}
    ledger=json.loads((P.parent/"memory/hannah_library_v0004.json").read_bytes())
    before=ca.canonical(ledger);rows=[]
    for s in science.S["sources"]:
        src={**s,"expected_sha256":auth.prev.S["expected_sha256"][s["name"]]}
        sel=science.choose(src,rref);p=O/"original.fil";row={"source":src["name"],"day":src["day"],"status":"HOLD"}
        try:
            expected=next(x for x in cb.S["reports"] if x["source"]==src["name"] and x["day"]==src["day"])
            raw=raws[expected["file"]];assert cb.sha(raw)==expected["sha256"]
            row["integrity"]=science.download(sel,src,p);p.chmod(0o444)
            ref=auth.Ref("ORIGINAL_ARCHIVE",src["name"],src["day"],p,src["expected_sha256"],sel["reference"]["size"])
            repeats=[];runs=[]
            for arm in S["order"][src["name"]]:
                m,receipt=armRun(arm,ledger,src,sel,ref,frozen[src["name"]],raw,expected)
                repeats.append(ca.canonical(m));runs.append(receipt)
            exact=all(x==repeats[0] for x in repeats)
            unchanged=auth.digest(p)==ref.sha256 and ca.canonical(ledger)==before
            row.update(runs=runs,measurement_exact=exact,original_unchanged=unchanged,matched_measurements=len(repeats) if exact else 0)
            assert exact and unchanged and len(runs)==4
            row["status"]="PASS_PAIRED_NATURAL_SCIENCE"
        except Exception as e:
            import traceback
            traceback.print_exc();row["error"]=type(e).__name__+": "+str(e)[:400]
        finally:p.unlink(missing_ok=True)
        rows.append(row);print("R37F3CB1_SOURCE",json.dumps(row),flush=True)
    check_pins();cb.pinned()
    n=sum(x["status"]=="PASS_PAIRED_NATURAL_SCIENCE" for x in rows)
    times={a:[r["seconds"] for x in rows for r in x.get("runs",[]) if r["arm"]==a] for a in ("A","B")}
    d={"status":"PASS_PAIRED_NATURAL_SCIENCE_REPORTING" if n==2 else "HOLD_PAIRED_NATURAL_SCIENCE_REPORTING","sources":rows,"sources_pass":n,"measurements_exact":sum(x.get("matched_measurements",0) for x in rows),"measurements_planned":8,"wall_seconds":times,"g4_entries":len(ledger["entries"]),"g4_unchanged":ca.canonical(ledger)==before,"production_adoption":False,"end_to_end_hannah_tested":False,"detector_comparison":False,"persistent_local_archive":False,"independent_signature":False,"externally_replicated":False}
    (O/"result.json").write_bytes(ca.canonical(d)+b"\n");print("R37F3CB1_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
