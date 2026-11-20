"""R42A4B: authenticate and reuse two frozen pilot groups, replay nine remaining groups."""
import hashlib,io,json,os,pathlib,zipfile
import requests
import r42a4a_noise_coverage as source
P=pathlib.Path("research/seti/r42a4b_full_noise_coverage_protocol.json")
Q=json.loads(P.read_text())
O=pathlib.Path("runtime/r42a4b")
O.mkdir(parents=True,exist_ok=True)
source.O=O
source.P=P
source.S=Q
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def authenticated_pilot():
    q=Q["pilot"];h={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
    u="https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/"+str(q["artifact"])
    def get(url):
        r=requests.get(url,headers=h,timeout=120);r.raise_for_status();return r
    m=get(u).json()
    assert m["workflow_run"]["id"]==q["run"] and m["digest"]=="sha256:"+q["digest"] and not m["expired"]
    b=get(u+"/zip").content
    assert len(b)<128*1024**2 and sha(b)==q["digest"]
    z=zipfile.ZipFile(io.BytesIO(b))
    files=[i for i in z.infolist() if not i.is_dir()]
    assert sum(i.file_size for i in files)<256*1024**2
    names={pathlib.PurePosixPath(i.filename).name:i for i in files}
    assert len(names)==len(files)
    return lambda name:z.read(names[name])
def prepare():
    source.prepare()
    r=authenticated_pilot()
    summary=json.loads(r("result.json"))
    assert summary["status"]=="PASS_ALL_CANDIDATE_NOISE_COVERAGE_PILOT__HOLD_GLOBAL_SCORE_AND_PROSPECTIVE_DATA"
    assert summary["pilot_groups"]==Q["pilot"]["groups"]
    assert summary["rows"]==Q["pilot"]["result"]["rows"] and summary["score_supported"]==Q["pilot"]["result"]["scores"]
    assert summary["original_qualified"]==Q["pilot"]["result"]["original_qualified"] and summary["rejected_scored"]==Q["pilot"]["result"]["rejected_scored"]
    for g in Q["pilot"]["groups"]:
        num=f"{g:02d}"
        d=json.loads(r("group-"+num+".json"))
        assert d["group"]==g and d["rows"]==288
        recs=json.loads(r("rows-"+num+".json"))
        assert len(recs)==288 and all(x["group"]==g for x in recs)
        for record in d["records"]:
            name=f'{num}-{record["kind"]}-noise-all.npz'
            raw=r(name)
            assert sha(raw)==record["file_sha256"]
            (O/name).write_bytes(raw)
        save("group-"+num+".json",d|{"pilot_reused":True})
        save("rows-"+num+".json",recs)
    save("pilot-receipt.json",{"status":"AUTHENTICATED_PILOT_REUSE","source_run":Q["pilot"]["run"],"source_artifact":Q["pilot"]["artifact"],"pilot_digest":Q["pilot"]["digest"],"pilot_groups":Q["pilot"]["groups"],"script_sha256":sha(pathlib.Path(__file__).read_bytes()),"protocol_sha256":sha(P.read_bytes())})
    print("R42A4B_PREPARED",json.dumps({"groups_reused":Q["pilot"]["groups"],"rows_reused":576}),flush=True)
def group(i):
    assert (i+1) not in Q["pilot"]["groups"] and 0<=i<11
    assert (O/"pilot-receipt.json").is_file()
    source.group(i)
def summary():
    assert (O/"pilot-receipt.json").is_file()
    source.summary()
    r=json.loads((O/"result.json").read_text())
    assert len(r["per_group"])==11 and r["rows"]==Q["expected"]["total_rows"] and r["pilot_groups"]==Q["pilot_groups"]
    assert sum(g["fit"] for g in r["per_group"])==6
    assert all(g["rows"]==288 for g in r["per_group"])
    reasons={}
    for g in r["per_group"]:
        for record in g["records"]:
            for reason,count in record["unsupported_reasons"].items():
                reasons[reason]=reasons.get(reason,0)+count
            raw=(O/f'{g["group"]:02d}-{record["kind"]}-noise-all.npz').read_bytes()
            assert sha(raw)==record["file_sha256"]
    assert reasons.get("SCORED")==r["score_supported"]
    assert sum(reasons.values())==r["rows"]
    rejected_total=sum(x["rejected_total"] for g in r["per_group"] for x in g["records"])
    assert r["rejected_scored"]<=rejected_total
    stratified={}
    for kind,subset in (("fit",[g for g in r["per_group"] if g["fit"]]),("audit",[g for g in r["per_group"] if not g["fit"]])):
        stratified[kind]={"groups":len(subset),"rows":sum(g["rows"] for g in subset),"supported":sum(g["score_supported"] for g in subset),"rejected_scored":sum(g["rejected_scored"] for g in subset),"original_qualified":sum(g["original_qualified"] for g in subset)}
    assert stratified["fit"]["rows"]==1728 and stratified["audit"]["rows"]==1440
    r.update(status="PASS_FULL_HISTORICAL_NOISE_FEATURE_COVERAGE__HOLD_GLOBAL_SCORE_AND_EXTERNAL_DATA",pilot_reused=True,reason_counts=reasons,rejected_total=rejected_total,rejected_coverage_fraction=r["rejected_scored"]/rejected_total,strata=stratified,continuous_combined_janus_noise=False,full_pipeline_roc=False,production_adoption=False,external_validation=False)
    save("result.json",r)
    print("R42A4B_SUMMARY",json.dumps({k:v for k,v in r.items() if k!="per_group"}),flush=True)
