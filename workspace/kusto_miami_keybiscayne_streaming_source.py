#!/usr/bin/env python3
from __future__ import annotations
import json, struct
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out";OUT.mkdir(parents=True,exist_ok=True)
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-MIAMI-KEYBISCAYNE-STREAMING-SOURCE-PREREG-2026-09-23-v1.0.json").read_text())
BASE=f"https://www.sciencebase.gov/catalog/file/get/{PRE['dataset']['sciencebase_item_id']}"
UA={"User-Agent":"JANUS-KUSTO-Miami-KeyBiscayne-streaming-source/1.0"}
TAIL=int(PRE["range_safety"]["tail_probe_bytes"])
MAX_META=int(PRE["range_safety"]["maximum_nonpayload_metadata_bytes"])

def head(name):
    r=requests.head(BASE,params={"name":name},headers=UA,timeout=90,allow_redirects=True)
    return {"status":r.status_code,"resolved_url":r.url,"content_length":r.headers.get("content-length"),"accept_ranges":r.headers.get("accept-ranges"),"content_type":r.headers.get("content-type")}

def range_get(name,spec):
    h=dict(UA);h["Range"]=spec
    r=requests.get(BASE,params={"name":name},headers=h,timeout=120,allow_redirects=True,stream=True)
    meta={"status":r.status_code,"resolved_url":r.url,"content_range":r.headers.get("content-range"),"content_length":r.headers.get("content-length"),"accept_ranges":r.headers.get("accept-ranges")}
    if r.status_code!=206:
        r.close()
        return meta,None
    cl=r.headers.get("content-length")
    if cl and int(cl)>MAX_META:
        r.close();meta["rejected"]="metadata response exceeds frozen cap";return meta,None
    b=r.content
    if len(b)>MAX_META:
        raise RuntimeError("metadata response exceeded frozen cap after read")
    return meta,b

def parse_eocd(tail):
    sig=b"PK\x05\x06";i=tail.rfind(sig)
    if i<0:return None
    if i+22>len(tail):return None
    v=struct.unpack_from("<4s4H2LH",tail,i)
    return {
      "disk_number":v[1],"central_directory_start_disk":v[2],
      "records_this_disk":v[3],"records_total":v[4],
      "central_directory_size":v[5],"central_directory_offset":v[6],
      "comment_length":v[7],"tail_eocd_offset":i
    }

def parse_cd(data):
    out=[];p=0
    while p+46<=len(data):
        if data[p:p+4]!=b"PK\x01\x02":break
        v=struct.unpack_from("<4s6H3I5H2I",data,p)
        flags=v[3];method=v[4];crc=v[7];csize=v[8];usize=v[9];n=v[10];e=v[11];cm=v[12];loff=v[16]
        q=p+46
        nameb=data[q:q+n]
        enc="utf-8" if flags & 0x800 else "cp437"
        name=nameb.decode(enc,errors="replace")
        out.append({"name":name,"compressed_size":csize,"uncompressed_size":usize,"compression_method":method,"crc32":f"{crc:08x}","local_header_offset":loff})
        p=q+n+e+cm
    return out,p

results=[]
for a in PRE["dataset"]["archives"]:
    name=a["name"]
    h=head(name)
    tail_meta,tail=range_get(name,f"bytes=-{TAIL}")
    rec={**a,"head":h,"tail_range":tail_meta,"numeric_payload_read":False}
    if tail is None:
        rec["status"]="RANGE_UNSUPPORTED_OR_UNAVAILABLE"
        results.append(rec);continue
    e=parse_eocd(tail)
    rec["eocd"]=e
    if e is None:
        rec["status"]="EOCD_NOT_FOUND_IN_FROZEN_TAIL"
        results.append(rec);continue
    if e["central_directory_size"]>MAX_META:
        rec["status"]="CENTRAL_DIRECTORY_EXCEEDS_FROZEN_METADATA_CAP"
        results.append(rec);continue
    start=e["central_directory_offset"];end=start+e["central_directory_size"]-1
    cd_meta,cd=range_get(name,f"bytes={start}-{end}")
    rec["central_directory_range"]=cd_meta
    if cd is None:
        rec["status"]="CENTRAL_DIRECTORY_RANGE_UNAVAILABLE"
        results.append(rec);continue
    members,consumed=parse_cd(cd)
    rec["member_count_parsed"]=len(members)
    rec["central_directory_bytes_consumed"]=consumed
    rec["members"]=members
    rec["status"]="STREAMABLE_ZIP_INDEX_BOUND" if members else "CENTRAL_DIRECTORY_PARSE_EMPTY"
    results.append(rec)
    print(name,rec["status"],"members",len(members),flush=True)

out={
 "artifact_id":"JANUS-KUSTO-MIAMI-KEYBISCAYNE-STREAMING-SOURCE-RUN-2026-09-23-v1.0",
 "prereg":PRE["artifact_id"],
 "depth_values_read":False,
 "numeric_member_payload_read":False,
 "archives":results,
 "streamable_archive_count":sum(1 for r in results if r["status"]=="STREAMABLE_ZIP_INDEX_BOUND"),
 "next_gate":"FREEZE_MEMBER_SELECTION_AND_NATIVE_RESOLUTION_FROM_ARCHIVE_METADATA_BEFORE_NUMERIC_READ",
 "claim_ceiling":"REMOTE_ARCHIVE_STRUCTURE_AND_ACCESS_CAPABILITY_ONLY"
}
p=OUT/"JANUS-KUSTO-MIAMI-KEYBISCAYNE-STREAMING-SOURCE-RUN-2026-09-23-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "streamable_archive_count":out["streamable_archive_count"],
 "archives":[{"id":r["id"],"status":r["status"],"members":r.get("member_count_parsed"),"head":r["head"],"tail_range":r["tail_range"]} for r in results],
 "depth_values_read":False
},indent=2))
