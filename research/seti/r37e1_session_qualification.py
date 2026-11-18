#!/usr/bin/env python3
"""R37E1: qualify entire dated capture groups using author SigMF metadata only."""
import collections,datetime,hashlib,json,math,pathlib,re,zipfile,urllib.parse
import requests
P=pathlib.Path("research/seti/r37e1_session_qualification_protocol.json");S=json.loads(P.read_text());SRC=S["source"]
O=pathlib.Path("runtime/r37e1");O.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,v):(O/n).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def get(url,path,expected):
    h=hashlib.md5();s=hashlib.sha256();n=0
    with requests.get(url,stream=True,timeout=(20,120)) as r:
        r.raise_for_status();assert urllib.parse.urlparse(r.url).hostname=="zenodo.org"
        with path.open("wb") as f:
            for b in r.iter_content(1048576):
                n+=len(b);assert n<=expected["bytes"];h.update(b);s.update(b);f.write(b)
    assert n==expected["bytes"] and h.hexdigest()==expected["md5"]
    return {"bytes":n,"md5":h.hexdigest(),"sha256":s.hexdigest()}
def audit(cfg):
    name=cfg["name"];env=cfg["environment"];path=O/name
    url=f"https://zenodo.org/records/{SRC['record']}/files/{name}?download=1"
    try:
        receipt=get(url,path,cfg)
        rows=[];issues=[];other=collections.Counter()
        with zipfile.ZipFile(path) as z:
            f={}
            for i in z.infolist():
                if i.is_dir():continue
                p=pathlib.PurePosixPath(i.filename)
                assert not p.is_absolute() and ".." not in p.parts and "\\" not in i.filename and i.filename not in f
                assert not i.flag_bits&1
                f[i.filename]=i
            metas={n[:-len(".sigmf-meta")]:n for n in f if n.endswith(".sigmf-meta")}
            data={n[:-len(".sigmf-data")]:n for n in f if n.endswith(".sigmf-data")}
            assert metas and set(metas)==set(data)
            for stem,n in sorted(metas.items()):
                i=f[n];j=f[data[stem]]
                assert 0<i.file_size<=65536
                raw=z.read(i);m=json.loads(raw);g=m["global"];caps=m["captures"];assert len(caps)>=1
                k=re.fullmatch(r"(\d{3})_(\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2})_(\d+)_"+re.escape(env),pathlib.PurePosixPath(stem).name)
                if not k:issues.append({"name":n,"reason":"unexpected capture filename"});continue
                date=datetime.datetime.strptime(k[2],"%Y_%m_%d_%H_%M_%S")
                freqs=[v.get("core:frequency") for v in caps];rate=g.get("core:sample_rate");typ=g.get("core:datatype")
                times=[v.get("core:datetime") for v in caps]
                valid=all(isinstance(v,(int,float)) and math.isfinite(v) and v>0 for v in [rate,*freqs])
                # Only declared SigMF datatype widths, never decode the sample array.
                width={"cf32_le":8,"cf32_be":8,"cf64_le":16,"cf64_be":16,"ci16_le":4,"ci16_be":4,"ci8":2,"cu8":2}.get(typ)
                integral=bool(width and j.file_size%width==0 and j.file_size>=width)
                capture_time_agreement=[]
                for v in times:
                    if v is None:capture_time_agreement.append("MISSING");continue
                    try:
                        dt=datetime.datetime.fromisoformat(v.replace("Z","+00:00"))
                        delta=(dt.replace(tzinfo=None)-date).total_seconds()
                        capture_time_agreement.append({"wall_clock_equal":abs(delta)<=1,"offset_seconds_uninterpreted":delta})
                    except (ValueError,AttributeError):capture_time_agreement.append("MALFORMED")
                # Record metadata key names only for extra identifiers; no waveform fields.
                labels={key:g[key] for key in ("core:description","core:hw","core:recorder","core:version") if key in g and isinstance(g[key],str) and len(g[key])<=1024}
                rows.append({"environment":env,"group_token":k[1],"session":date.isoformat(),"day":date.date().isoformat(),"replicate":int(k[3]),"meta":n,"data":data[stem],"meta_sha256":sha(raw),"data_crc":f"{j.CRC:08x}","data_bytes":j.file_size,"datatype":typ,"sample_rate":rate,"receiver_center_frequencies":freqs,"capture_times":times,"time_comparison":capture_time_agreement,"valid_physical_metadata":valid,"integral_sample_count":integral,"global_keys":sorted(g),"capture_keys":sorted({key for v in caps for key in v}),"annotation_keys":sorted({key for v in m.get("annotations",[]) for key in v}),"descriptive_metadata":labels})
            for n in f:
                if n not in metas.values() and n not in data.values():other[pathlib.PurePosixPath(n).suffix]+=1
        save(env+"_inventory.json",rows)
        return {"status":"ARCHIVE_INVENTORY_VERIFIED","receipt":receipt,"rows":rows,"issues":issues,"other_representations":dict(other)}
    finally:path.unlink(missing_ok=True)
def main():
    results=[]
    for f in SRC["files"]:
        try:
            r=audit(f);r["environment"]=f["environment"];results.append(r)
        except Exception as e:results.append({"environment":f["environment"],"status":"HOLD_ARCHIVE","error":type(e).__name__+": "+str(e)[:400]})
    rows=[v for d in results for v in d.get("rows",[])]
    by=collections.defaultdict(list);session=collections.defaultdict(set)
    for r in rows:
        by[(r["environment"],r["day"])].append(r);session[(r["environment"],r["session"])].add(r["group_token"])
    days=[{"environment":k[0],"day":k[1],"recordings":len(v),"group_tokens":sorted({r["group_token"] for r in v}),"sessions":len({r["session"] for r in v}),"all_five_tokens":len({r["group_token"] for r in v})==5} for k,v in sorted(by.items())]
    freqs=sorted({float(v) for r in rows if r["valid_physical_metadata"] for v in r["receiver_center_frequencies"]})
    rates=sorted({float(r["sample_rate"]) for r in rows if r["valid_physical_metadata"]})
    crc=collections.defaultdict(list)
    for r in rows:crc[(r["data_crc"],r["data_bytes"])].append((r["environment"],r["session"],r["group_token"],r["replicate"]))
    duplicates=[{"crc":k[0],"bytes":k[1],"members":v} for k,v in crc.items() if len(v)>1]
    chamber=[d for d in days if d["environment"]=="anechoic_chamber" and d["all_five_tokens"]]
    candidate_split=None
    if len(chamber)>=3:
        candidate_split={"fit":chamber[0],"calibration":chamber[1],"evaluation":chamber[-1],"separate_environment":[d for d in days if d["environment"]=="office_room"],"authority":"CANDIDATE_ONLY_NO_WAVEFORM_ACCESS"}
    result={"status":"HOLD_AUTHOR_LABEL_AND_TRANSMIT_CHANNEL_REVIEW","archives":[{k:v for k,v in r.items() if k!="rows"} for r in results],"recordings":len(rows),"days":days,"sessions":len(session),"sessions_with_all_five_filename_tokens":sum(len(v)==5 for v in session.values()),"receiver_center_frequencies_hz":freqs,"sample_rates_hz":rates,"valid_metadata":sum(r["valid_physical_metadata"] for r in rows),"integral_sample_count":sum(r["integral_sample_count"] for r in rows),"crc_size_duplicate_candidates":duplicates,"candidate_split":candidate_split,"capture_metadata_keys":sorted({k for r in rows for k in r["capture_keys"]}),"global_metadata_keys":sorted({k for r in rows for k in r["global_keys"]}),"descriptive_metadata_examples":list({json.dumps(r["descriptive_metadata"],sort_keys=True) for r in rows})[:10],"raw_waveform_samples_read":0,"classifier_trained":False,"exam_started":False,"production_adoption":False,"limitations":"Filename tokens are not independently verified device labels; common receiver tuning is not proof of identical emitted center frequency. CRC/size duplicate checks are not independent sample-hash verification. No classification."}
    save("result.json",result);print("R37E1_SUMMARY",json.dumps(result),flush=True)
if __name__=="__main__":main()
