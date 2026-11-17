#!/usr/bin/env python3
"""R20C: sealed cross-run tile cache, verified on independent range replay."""
import hashlib,io,json,os,pathlib,re,time,statistics,zipfile,sys
import numpy as np
import requests
import r9_expanded_search as eng
import juxtapose_undina_r18bc as old
import r17c_frozen_holdout as hold
P=pathlib.Path("research/seti/r20c_cache_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r20c");OUT.mkdir(parents=True,exist_ok=True)
CACHE=OUT/"cache";CACHE.mkdir(exist_ok=True)
COLD=OUT/"cold";COLD.mkdir(exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
def emit(k,d):print(k,json.dumps(d),flush=True)
def download(spec):
    t=time.perf_counter();b=bytearray();h=hashlib.sha256()
    with requests.get(f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{spec['id']}/zip",headers={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"},stream=True,timeout=(15,120)) as r:
        r.raise_for_status()
        for c in r.iter_content(1024*1024):
            b.extend(c);h.update(c)
            if len(b)>S["limits"]["zip_bytes"]:raise ValueError("ZIP_BUDGET")
    if h.hexdigest()!=spec["digest"]:raise ValueError("ZIP_DIGEST")
    if "size" in spec and len(b)!=spec["size"]:raise ValueError("ZIP_SIZE")
    return b,{"seconds":time.perf_counter()-t,"bytes":len(b)}
def safe_members(z):
    result={}
    for m in z.infolist():
        if m.is_dir():continue
        p=pathlib.PurePosixPath(m.filename)
        if p.is_absolute() or ".." in p.parts or "\\" in m.filename:raise ValueError("BAD_MEMBER_PATH")
        if p.name in result:raise ValueError("DUPLICATE_BASENAME")
        result[p.name]=m
    return result
def tile_bytes(t,p):
    if not re.fullmatch(r"[0-9a-f]{24}\.npy",t["artifact"]):raise ValueError("BAD_TILE_NAME")
    if not p.is_file() or p.stat().st_size>S["limits"]["member_bytes"]:raise ValueError("MISSING_OR_OVERSIZED_TILE")
    x=np.load(p,allow_pickle=False)
    if x.dtype!=np.float32 or list(x.shape)!=t["shape"] or not np.isfinite(x).all():raise ValueError("TILE_SCHEMA")
    if sha(x.tobytes())!=t["slice_sha256"]:raise ValueError("TILE_DIGEST")
    return x
def identity(f,t):
    if hashlib.sha256((f["url"]+"|"+str(t["center_mhz"])).encode()).hexdigest()[:24]+".npy"!=t["artifact"]:raise ValueError("SOURCE_IDENTITY")
    h=f["header"]
    if t["df_hz"]!=h["df"]*1e6 or t["ts"]!=h["ts"] or t["shape"][0]!=h["shape"][0]:raise ValueError("STALE_HEADER")
    q=(t["channel0_mhz"]-h["f0"])*1e6/t["df_hz"]
    if abs(q-round(q))>1e-4 or round(q)<0 or round(q)+t["shape"][1]>h["shape"][-1]:raise ValueError("CHANNEL_ALIGNMENT")
    if abs(t["center_mhz"]-t["channel0_mhz"])*1e6>60000:raise ValueError("WINDOW_ALIGNMENT")
def restore():
    b,mt=download(S["source_manifest"])
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        ms=safe_members(z);m=ms["manifest.json"]
        assert m.file_size<1024*1024
        manifest=json.loads(z.read(m))
    sources={f["name"]:f for q in manifest["selected"] for f in q["files"] if f["role"].startswith("OFF")}
    by={r["name"]:r for r in S["files"]}
    assert len(sources)==len(by)==12 and set(sources)==set(by)
    alltiles={}
    for name,r in by.items():
        f=sources[name]
        assert r["role"]==f["role"] and len(r["tiles"])==4
        for t in r["tiles"]:
            identity(f,t)
            assert t["artifact"] not in alltiles
            alltiles[t["artifact"]]=t
    b,at=download(S["source_artifact"])
    start=time.perf_counter();total=0
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        ms=safe_members(z)
        assert {k for k in ms if k.endswith(".npy")}==set(alltiles)
        for name,t in alltiles.items():
            m=ms[name];total+=m.file_size
            if m.file_size>S["limits"]["member_bytes"] or total>S["limits"]["total_uncompressed_tiles"]:raise ValueError("EXTRACT_BUDGET")
            payload=z.read(m);p=CACHE/name
            tmp=p.with_suffix(".tmp");tmp.write_bytes(payload)
            tile_bytes(t,tmp)
            tmp.replace(p)
    return sources,by,manifest,{"manifest_download":mt,"artifact_download":at,"verified_extract_seconds":time.perf_counter()-start,"tile_bytes":total,"tiles":len(alltiles)}
def cached_get(f,sources,by):
    expected=sources.get(f["name"])
    if expected is None or expected["url"]!=f["url"] or expected["header"]!=f["header"] or expected["role"]!=f["role"]:raise ValueError("CACHE_SOURCE_MISMATCH")
    rec=by[f["name"]]
    for t in rec["tiles"]:
        identity(f,t);tile_bytes(t,CACHE/t["artifact"])
    return rec["tiles"]
def fail(fn,reason):
    try:fn()
    except ValueError:return reason
    raise AssertionError("negative control accepted "+reason)
def selftests(sources,by):
    f=next(iter(sources.values()));t=by[f["name"]]["tiles"][0];src=CACHE/t["artifact"]
    p=OUT/"corruption-test.npy";b=bytearray(src.read_bytes());b[-1]^=1;p.write_bytes(b)
    checks=[fail(lambda:tile_bytes(t,p),"CORRUPT_TILE_REJECTED")]
    p.unlink();checks.append(fail(lambda:tile_bytes(t,p),"MISSING_TILE_REJECTED"))
    wrong={**f,"header":{**f["header"],"ts":f["header"]["ts"]+1}}
    checks.append(fail(lambda:cached_get(wrong,sources,by),"STALE_HEADER_REJECTED"))
    wrong={**f,"url":f["url"]+"?changed=1"}
    checks.append(fail(lambda:cached_get(wrong,sources,by),"STALE_SOURCE_REJECTED"))
    checks.append(fail(lambda:tile_bytes({**t,"shape":[1,1]},src),"SHAPE_MISMATCH_REJECTED"))
    return checks
def read_all(sources,by):
    a=[]
    for f in sources.values():
        ts=cached_get(f,sources,by)
        a.extend(t["slice_sha256"] for t in ts)
    return a
def main():
    sources,by,manifest,restore_cost=restore();checks=selftests(sources,by)
    expected=sorted(t["slice_sha256"] for r in by.values() for t in r["tiles"])
    warm=[]
    for _ in range(S["repeats"]):
        t=time.perf_counter();a=read_all(sources,by);warm.append(time.perf_counter()-t)
        assert sorted(a)==expected
    original=eng.get_tiles;eng.OUT=COLD;old.OUT=COLD
    old.S={**old.S,"seed":"R20C-replay-seed"}
    cold_rows=[];actual_cold=[]
    for name in S["cold_replay_files"]:
        f=sources[name];t=time.perf_counter();tiles=original(f);sec=time.perf_counter()-t
        expected_tiles=by[name]["tiles"]
        assert tiles==expected_tiles
        for x in tiles:assert np.array_equal(tile_bytes(x,COLD/x["artifact"]),tile_bytes(x,CACHE/x["artifact"]))
        cold_rows.append({"name":name,"seconds":sec,"verified_tiles":len(tiles)})
        actual_cold.append(f)
    # Recreate the exact same patch stream twice, from independently fetched tiles and sealed cache.
    q={"files":actual_cold}
    eng.get_tiles=lambda f:by[f["name"]]["tiles"]
    old.OUT=COLD
    t=time.perf_counter();a=old.gather(q);cold_patch_seconds=time.perf_counter()-t
    eng.get_tiles=lambda f:cached_get(f,sources,by)
    old.OUT=CACHE
    t=time.perf_counter();b=old.gather(q);cache_patch_seconds=time.perf_counter()-t
    assert len(a)==len(b)==256
    for x,y in zip(a,b):
        assert x[1:]==y[1:] and np.array_equal(x[0],y[0])
    cold_seconds=sum(r["seconds"] for r in cold_rows)
    # Fair same-two-files validated-read comparator, not whole-cache/local versus one cold file.
    pair=[]
    for _ in range(S["repeats"]):
        t=time.perf_counter()
        for f in actual_cold:cached_get(f,sources,by)
        pair.append(time.perf_counter()-t)
    d={"status":"VERIFIED_REPLAY","source_tiles":48,"cold_replayed_tiles":8,"identical_patches":256,"mismatches":0,"negative_controls":checks,"restore":restore_cost,"warm_all_seconds":warm,"warm_all_median_seconds":statistics.median(warm),"warm_pair_seconds":pair,"warm_pair_median_seconds":statistics.median(pair),"cold_pair_seconds":cold_seconds,"cold_pair_details":cold_rows,"cold_patch_build_seconds":cold_patch_seconds,"cache_patch_build_seconds":cache_patch_seconds,"historical_source_range_bytes":sum(r["bytes"] for r in S["files"]),"scope":"Same R19E data, fresh remote replay of 2 predetermined OFF scans; all 48 pins checked. Cache restore cost included separately; warm ratio does not mean new archive acquisition accelerated. Cache scoped to pinned manifest, not arbitrary latest inputs. No production adoption.","production_adoption":False}
    save("result.json",d);emit("R20C_SUMMARY",d)
if __name__=="__main__":main()
