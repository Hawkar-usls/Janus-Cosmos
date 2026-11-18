#!/usr/bin/env python3
"""R37E0: authenticate laboratory archive and qualify independent capture units; no classifier fit."""
import hashlib,json,pathlib,re,zipfile,urllib.parse
import numpy as np
import requests
P=pathlib.Path("research/seti/r37e_identity_qualification_protocol.json")
S=json.loads(P.read_text());O=pathlib.Path("runtime/r37e");O.mkdir(parents=True,exist_ok=True)
def save(n,d):(O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
def canaries():
    x=np.array([1+2j,3-1j,-2+4j,1-3j,2+1j,-4-1j,2-2j,3+2j])
    y=1j*x
    def stft_power(v):
        return np.stack([abs(np.fft.fft(v[i:i+4]))**2 for i in range(5)])
    e=float(np.max(abs(stft_power(x)-stft_power(y))));assert e<=1e-12 and not np.array_equal(x,y)
    # Different relative phases, identical integrated spectrum; not just global phase.
    u=np.array([1,2,3,4,5,6,7,8],complex);v=u*np.array([1,1j,-1,-1j,1,-1,1j,-1j])
    a=np.fft.ifft(u);b=np.fft.ifft(v)
    err=float(np.max(abs(abs(np.fft.fft(a))**2-abs(np.fft.fft(b))**2)))
    assert err<1e-10 and not np.allclose(a,b) and not np.allclose(a/a[0],b/b[0])
    return {"global_phase_stft_power_max_error":e,"distinct_relative_phase_spectrum_max_error":err,"status":"PASS_NONINJECTIVE_POWER_OBSERVATION","scope":"Mathematical toy examples only; not observations of distinct actual transmitters."}
def url_ok(u):
    x=urllib.parse.urlparse(u);return x.scheme=="https" and x.hostname=="zenodo.org"
def metadata(record):
    u=f"https://zenodo.org/api/records/{record}"
    with requests.get(u,timeout=(15,60),stream=True) as r:
        r.raise_for_status();assert url_ok(r.url)
        b=bytearray()
        for x in r.iter_content(65536):
            b.extend(x);assert len(b)<=4*1024**2
    d=json.loads(b);assert int(d["id"])==record
    (O/f"metadata-{record}.json").write_bytes(b)
    return {"record":record,"doi":d["doi"],"title":d["metadata"]["title"],"license":d["metadata"].get("license"),"files":[{"key":v["key"],"size":v["size"],"checksum":v["checksum"]} for v in d["files"]],"metadata_sha256":hashlib.sha256(b).hexdigest()}
def inventory(cfg,md):
    assert md["doi"]==cfg["doi"]
    f=next(x for x in md["files"] if x["key"]==cfg["filename"]);assert f["checksum"]=="md5:"+cfg["md5"] and 0<f["size"]<=cfg["max_bytes"]
    u=f"https://zenodo.org/records/{cfg['record']}/files/{cfg['filename']}?download=1"
    path=O/"source.zip";h=hashlib.md5();s=hashlib.sha256();size=0
    try:
        with requests.get(u,timeout=(20,120),stream=True) as r:
            r.raise_for_status();assert url_ok(r.url)
            with path.open("wb") as out:
                for b in r.iter_content(1048576):
                    size+=len(b);assert size<=f["size"] and size<=cfg["max_bytes"]
                    h.update(b);s.update(b);out.write(b)
        assert size==f["size"] and h.hexdigest()==cfg["md5"]
        with zipfile.ZipFile(path) as z:
            ins=z.infolist();assert len(ins)<=50000
            rows=[];names=set()
            for i in ins:
                n=i.filename;p=pathlib.PurePosixPath(n)
                assert not p.is_absolute() and ".." not in p.parts and "\\" not in n and n not in names
                names.add(n)
                if i.is_dir():continue
                rows.append({"name":n,"bytes":i.file_size,"crc32":f"{i.CRC:08x}","encrypted":bool(i.flag_bits&1)})
            counts={}
            for i in rows:
                suffix=pathlib.PurePosixPath(i["name"]).suffix.lower();counts[suffix]=counts.get(suffix,0)+1
            # Device/capture hints are not accepted as labels until author schema review.
            hints=sorted(set(m.group(0) for i in rows for m in re.finditer(r"(?i)(?:device|node|trial|session|capture)[_-]?\d+",i["name"])))
        save("archive_inventory.json",rows)
        return {"status":"CHECKSUM_VERIFIED_INVENTORY_ONLY","bytes":size,"sha256":s.hexdigest(),"published_md5":h.hexdigest(),"members":len(rows),"extensions":counts,"device_session_hints":hints[:200],"example_names":[v["name"] for v in rows[:20]],"raw_signal_samples_read":0,"extraction":False}
    finally:path.unlink(missing_ok=True)
def main():
    c=canaries();save("observability.json",c);print("R37E_OBSERVABILITY",json.dumps(c),flush=True)
    rows=[]
    for cfg in S["sources"]:
        if "record" not in cfg:continue
        try:
            m=metadata(cfg["record"]);assert m["doi"]==cfg["doi"];rows.append({"name":cfg["name"],"status":"METADATA_VERIFIED","metadata":m})
        except Exception as e:rows.append({"name":cfg["name"],"status":"HOLD_METADATA","error":type(e).__name__+": "+str(e)[:400]})
    save("source_metadata.json",rows);print("R37E_METADATA",json.dumps(rows),flush=True)
    cfg=S["sources"][0];m=next(x for x in rows if x["name"]=="INRIA")
    try:
        assert m["status"]=="METADATA_VERIFIED"
        inv=inventory(cfg,m["metadata"])
    except Exception as e:inv={"status":"HOLD_ARCHIVE","error":type(e).__name__+": "+str(e)[:400]}
    result={"status":"HOLD_CAPTURE_UNIT_QUALIFICATION","observability":c,"metadata":rows,"archive":inv,"exam_started":False,"classifier_trained":False,"verified_promoted":0,"production_adoption":False}
    save("result.json",result);print("R37E_SUMMARY",json.dumps(result),flush=True)
if __name__=="__main__":main()
