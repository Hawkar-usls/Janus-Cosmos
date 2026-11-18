#!/usr/bin/env python3
"""Isolated historical astronomical archive-to-report integrity audit."""
import copy,hashlib,io,json,os,pathlib,platform,zipfile
import numpy as np
import requests
import astropy
import astropy_iers_data
from astropy.utils import iers
from astropy.coordinates import solar_system_ephemeris
P=pathlib.Path("research/seti/r37f3ca_science_provenance_protocol.json")
S=json.loads(P.read_text())
for name,expected in S["pinned_blobs"].items():
    b=(P.parent/name).read_bytes()
    assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==expected
import r37f3b1_native_authority as auth
science=auth.science
O=pathlib.Path("runtime/r37f3ca")
O.mkdir(parents=True,exist_ok=True)
def canonical(d):
    return json.dumps(d,sort_keys=True,separators=(",",":"),ensure_ascii=False,allow_nan=False).encode()
def sha(b):return hashlib.sha256(b).hexdigest()
def seal(body):return {"report":body,"report_sha256":sha(canonical(body))}
def save(path,d):path.write_bytes(canonical(d)+b"\n")
def parent():
    p=S["parent"]
    with requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['artifact']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90,stream=True) as r:
        r.raise_for_status();b=bytearray()
        for x in r.iter_content(65536):
            b.extend(x);assert len(b)<16*1024**2
    assert sha(b)==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        names=[n for n in z.infolist() if pathlib.PurePosixPath(n.filename).name=="result.json"]
        assert len(names)==1 and names[0].file_size<8*1024**2
        d=json.loads(z.read(names[0]))
    assert d["status"]=="PASS_ISOLATED_ORIGINAL_AUTHORITY" and d["original_authority_pass"]==2 and d["negative_rejected"]==22
    return d
def report(src,sel,ref,measurement,context):
    return {"schema":"janus.cosmos.natural_science_report.v1","source":src["name"],"day":src["day"],
        "original":{"role":ref.role,"archive_url":sel["url"],"bytes":ref.size,"sha256":ref.sha256,
        "published_md5":sel["reference"]["md5"],"catalogue_https":sel["reference"]["portal_page"]},
        "method":context,"measurement":{"coarse_reference_exact":measurement["coarse_reference_exact"],
        "native_channels":measurement["native_channels"],"rows":measurement["rows"],
        "native_summary":measurement["native_summary"],"native_df_hz":science.M["native_df_hz"]},
        "limitations":{"conditional_lsrk":True,"absolute_flux_calibrated":False,"spatial_off_verified":False,
        "closed_independent_test":False,"persistent_local_archive":False,"end_to_end_hannah_tested":False,
        "production_adoption":False}}
def verify(envelope,expected):
    try:
        assert set(envelope)=={"report","report_sha256"}
        assert envelope["report_sha256"]==sha(canonical(envelope["report"]))
        # Expected fields come from separately repeated authorized original measurement,
        # not from any values supplied by the report under examination.
        assert canonical(envelope["report"])==canonical(expected)
        return True
    except (AssertionError,TypeError,ValueError,KeyError):return False
def mutations(body):
    cases=[
        ("source",["source"],"OTHER"),
        ("date",["day"],body["day"]+1),
        ("role",["original","role"],"DISPLAY_ONLY"),
        ("archive_url",["original","archive_url"],"https://example.invalid/derived"),
        ("bytes",["original","bytes"],body["original"]["bytes"]-1),
        ("sha256",["original","sha256"],"0"*64),
        ("md5",["original","published_md5"],"0"*32),
        ("parent",["method","parent","sha256"],"0"*64),
        ("code",["method","authority_blob"],"0"*40),
        ("resolution",["measurement","native_df_hz"],715.2557373046875),
        ("measurement",["measurement","native_summary",0,"area"],0.0),
        ("false_claim",["limitations","end_to_end_hannah_tested"],True)]
    for name,path,value in cases:
        bad=copy.deepcopy(body);q=bad
        for k in path[:-1]:q=q[k]
        q[path[-1]]=value
        yield name,seal(bad)
def main():
    iers.conf.auto_download=False;iers.conf.iers_degraded_accuracy="error";solar_system_ephemeris.set("builtin")
    p=parent();refs=science.refs()
    frozen={c["source"]:(v,y) for c,r,v,y,e,step,receipt in auth.guard.fetch()}
    context={"parent":S["parent"],"authority_blob":S["pinned_blobs"]["r37f3b1_native_authority.py"],
        "pinned_measurement_blobs":auth.prev.S["pinned_blobs"],"protocol_sha256":sha(P.read_bytes()),
        "python":platform.python_version(),"numpy":np.__version__,"astropy":astropy.__version__,"iers_data":astropy_iers_data.__version__}
    rows=[]
    for source in science.S["sources"]:
        src={**source,"expected_sha256":auth.prev.S["expected_sha256"][source["name"]]}
        sel=science.choose(src,refs);path=O/"original.fil"
        row={"source":src["name"],"day":src["day"],"status":"HOLD"}
        try:
            pr=[r for r in p["sources"] if r["source"]==src["name"] and r["day"]==src["day"]]
            assert len(pr)==1 and pr[0]["status"]=="PASS_ORIGINAL_AUTHORITY"
            row["integrity"]=science.download(sel,src,path);path.chmod(0o444)
            ref=auth.Ref("ORIGINAL_ARCHIVE",src["name"],src["day"],path,src["expected_sha256"],sel["reference"]["size"])
            first,_=auth.measure(ref,src,sel,frozen[src["name"]]);del _
            body=report(src,sel,ref,first,context);del first
            env=seal(body);file=O/(("w3oh" if src["name"]=="W3OH" else "g42")+"_science_report.json")
            save(file,env)
            # Full independent replay within the same runtime, not external replication.
            second,_=auth.measure(ref,src,sel,frozen[src["name"]]);del _
            expected=report(src,sel,ref,second,context)
            parsed=json.loads(file.read_bytes())
            replay=verify(parsed,expected)
            checks=[{"name":name,"rejected":not verify(bad,expected)} for name,bad in mutations(body)]
            assert len(checks)==12 and replay and all(c["rejected"] for c in checks)
            assert auth.digest(path)==ref.sha256 and second["coarse_reference_exact"]
            row.update(status="PASS_ARCHIVE_TO_REPORT",report_file=file.name,report_sha256=env["report_sha256"],report_file_sha256=auth.digest(file),replay_exact=replay,checks=checks,original_unchanged=True)
        except Exception as e:
            import traceback
            traceback.print_exc();row["error"]=type(e).__name__+": "+str(e)[:400]
        finally:path.unlink(missing_ok=True)
        rows.append(row)
        print("R37F3CA_SOURCE",json.dumps(row),flush=True)
    n=sum(r["status"]=="PASS_ARCHIVE_TO_REPORT" for r in rows)
    k=sum(c["rejected"] for r in rows for c in r.get("checks",[]))
    d={"status":"PASS_ISOLATED_ARCHIVE_TO_REPORT" if n==2 and k==24 else "HOLD_ARCHIVE_TO_REPORT",
        "sources":rows,"reports_verified":n,"tamper_rejected":k,"tamper_count":sum(len(r.get("checks",[])) for r in rows),
        "parent_verified":True,"persistent_local_archive":False,"independent_signature":False,
        "end_to_end_hannah_tested":False,"spatial_off_verified":False,"production_adoption":False}
    save(O/"result.json",d);print("R37F3CA_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
