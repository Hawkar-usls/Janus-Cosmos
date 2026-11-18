#!/usr/bin/env python3
"""Read-only astronomical report sidecar at HANNAH's research boundary."""
import copy,hashlib,io,json,os,pathlib,zipfile
import requests
P=pathlib.Path("research/seti/r37f3cb_science_sidecar_protocol.json")
S=json.loads(P.read_text());O=pathlib.Path("runtime/r37f3cb");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def canonical(x):return json.dumps(x,sort_keys=True,separators=(",",":"),ensure_ascii=False,allow_nan=False).encode()
def blob(p):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def pinned():
    for name,h in S["pinned_blobs"].items():assert blob(P.parent/name)==h,name
pinned()
from hannah_observation_hook import ObservationAccumulator
def fetch():
    p=S["parent"];url=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{p['artifact']}/zip"
    with requests.get(url,headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},timeout=90,stream=True) as r:
        r.raise_for_status();b=bytearray()
        for x in r.iter_content(65536):b.extend(x);assert len(b)<16*1024**2
    assert sha(b)==p["sha256"]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        names={pathlib.PurePosixPath(x.filename).name:x for x in z.infolist() if not x.is_dir()}
        assert len(names)==len([x for x in z.infolist() if not x.is_dir()]) and sum(x.file_size for x in names.values())<16*1024**2
        d=json.loads(z.read(names["result.json"]));assert d["status"]=="PASS_ISOLATED_ARCHIVE_TO_REPORT" and d["reports_verified"]==2 and d["tamper_rejected"]==24
        return {r["file"]:z.read(names[r["file"]]) for r in S["reports"]}
class ScienceReportSidecar:
    """The report namespace cannot be written to candidate memory or used for decisions."""
    def __init__(self,observation_boundary):
        self.boundary=observation_boundary
        self.references={}
    def attach(self,raw,expected):
        if expected not in S["reports"]:raise ValueError("untrusted authority")
        if not isinstance(raw,bytes) or len(raw)>1024**2 or sha(raw)!=expected["sha256"]:raise ValueError("frozen report bytes mismatch")
        try:
            env=json.loads(raw);b=env["report"]
            if set(env)!={"report","report_sha256"} or env["report_sha256"]!=sha(canonical(b)):raise ValueError("report digest")
            if canonical(env)+b"\n"!=raw:raise ValueError("noncanonical")
            if b["schema"]!="janus.cosmos.natural_science_report.v1" or (b["source"],b["day"])!=(expected["source"],expected["day"]):raise ValueError("identity")
            if b["original"]["role"]!="ORIGINAL_ARCHIVE" or b["measurement"]["rows"]!=6 or not b["measurement"]["coarse_reference_exact"]:raise ValueError("source contract")
            if b["limitations"]["production_adoption"] or b["limitations"]["end_to_end_hannah_tested"]:raise ValueError("unsupported authority")
        except (KeyError,TypeError,AssertionError,UnicodeError) as e:raise ValueError("invalid report") from e
        key=(expected["source"],expected["day"])
        record={"source":key[0],"day":key[1],"report_file_sha256":expected["sha256"],"report_sha256":env["report_sha256"],"artifact":S["parent"]["artifact"],"artifact_sha256":S["parent"]["sha256"],"role":"SCIENCE_REPORT_REFERENCE_ONLY"}
        if key in self.references:
            if self.references[key]!=record:raise ValueError("conflicting reference")
            return "DUPLICATE"
        self.references[key]=record
        return "ATTACHED"
    def export(self):
        return {"schema":"janus.cosmos.science_references.v1","references":list(self.references.values()),"candidate_memory_adoption":False,"production_adoption":False}
def mutate(raw,path,value):
    env=json.loads(raw);q=env["report"]
    for k in path[:-1]:q=q[k]
    q[path[-1]]=value
    env["report_sha256"]=sha(canonical(env["report"]))
    return canonical(env)+b"\n"
def main():
    reports=fetch()
    parent=json.loads((P.parent/"memory/hannah_library_v0004.json").read_bytes())
    before=canonical(parent)
    hook=ObservationAccumulator(parent,{"id":S["parent"]["artifact"],"digest":S["parent"]["sha256"]})
    first,delta,rep=hook.seal()
    assert canonical(first)==before and delta["added"]==0 and rep["added"]==0
    side=ScienceReportSidecar(hook);rows=[]
    for expected in S["reports"]:
        raw=reports[expected["file"]]
        assert side.attach(raw,expected)=="ATTACHED"
        assert side.attach(raw,expected)=="DUPLICATE"
        cases=[("malformed",b"{",expected),("wrong_authority",raw,{"source":"OTHER"}),("wrong_source",mutate(raw,["source"],"OTHER"),expected),("wrong_date",mutate(raw,["day"],expected["day"]+1),expected),("derived",mutate(raw,["original","role"],"DISPLAY_ONLY"),expected),("changed_measurement",mutate(raw,["measurement","native_df_hz"],715.2557373046875),expected),("false_claim",mutate(raw,["limitations","end_to_end_hannah_tested"],True),expected),("extra_bytes",raw+b" ",expected)]
        checks=[]
        for name,b,e in cases:
            try:side.attach(b,e);ok=False
            except ValueError:ok=True
            checks.append({"name":name,"rejected":ok})
        rows.append({"source":expected["source"],"day":expected["day"],"attached":True,"duplicate_idempotent":True,"checks":checks})
    second,d2,r2=hook.seal()
    unchanged=canonical(second)==before and canonical(parent)==before and d2["added"]==0 and r2["added"]==0
    assert not hook.records and not hook.observations and not hook.seen
    pinned()
    k=sum(x["rejected"] for a in rows for x in a["checks"])
    out=side.export();assert len(out["references"])==2
    (O/"science_references.json").write_bytes(canonical(out)+b"\n")
    d={"status":"PASS_ISOLATED_SCIENCE_REPORT_BOUNDARY" if unchanged and k==16 else "HOLD_SCIENCE_REPORT_BOUNDARY","sources":rows,"report_references":2,"duplicates_idempotent":2,"tamper_rejected":k,"tamper_count":16,"g4_entries":len(parent["entries"]),"g4_unchanged":unchanged,"candidate_records_added":0,"actual_hook_empty_seal_tested":True,"new_candidate_ingest_tested":False,"original_measurement_replayed_in_this_run":False,"end_to_end_hannah_tested":False,"persistent_local_archive":False,"independent_signature":False,"production_adoption":False}
    (O/"result.json").write_bytes(canonical(d)+b"\n")
    print("R37F3CB_SUMMARY",json.dumps(d),flush=True)
if __name__=="__main__":main()
