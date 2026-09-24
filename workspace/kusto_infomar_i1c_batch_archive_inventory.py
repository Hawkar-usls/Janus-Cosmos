#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, struct
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
REC=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1B-PREOUTCOME-BATCH-METADATA-AND-TRANSPORT-RECEIPT-2026-09-24-v1.0.json").read_text())
PRE=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I1C-BATCH-ARCHIVE-MEMBER-INVENTORY-PREREG-2026-09-24-v1.0.json").read_text())
S=requests.Session(); S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I1C/1.0","Accept-Encoding":"identity"})
EOCD=bytes.fromhex("504b0506"); CDSIG=bytes.fromhex("504b0102")

def ranged(url,start,end):
    r=S.get(url,headers={"Range":f"bytes={start}-{end}"},timeout=180,allow_redirects=True)
    r.raise_for_status()
    if r.status_code!=206: raise RuntimeError(f"range status {r.status_code}: {url}")
    return r.content

def size_probe(url):
    r=S.get(url,headers={"Range":"bytes=0-0"},timeout=120,allow_redirects=True,stream=True)
    if r.status_code!=206: raise RuntimeError(f"size probe status {r.status_code}: {url}")
    cr=r.headers.get("Content-Range") or r.headers.get("content-range") or ""
    if "/" not in cr: raise RuntimeError(f"missing Content-Range: {url}")
    return int(cr.rsplit("/",1)[1])

def inventory(url,expected_size):
    size=size_probe(url)
    if size!=int(expected_size): raise RuntimeError(f"archive size drift {size} vs {expected_size}: {url}")
    n=min(size,65557); start=size-n
    tail=ranged(url,start,size-1)
    p=tail.rfind(EOCD)
    if p<0 or p+22>len(tail): raise RuntimeError(f"EOCD not found: {url}")
    disk,cd_disk,disk_n,total_n,cd_size,cd_offset,comment_len=struct.unpack_from("<HHHHIIH",tail,p+4)
    if disk!=0 or cd_disk!=0 or disk_n!=total_n: raise RuntimeError("multidisk ZIP unsupported")
    if total_n==0xffff or cd_size==0xffffffff or cd_offset==0xffffffff: raise RuntimeError("ZIP64 requires explicit gate")
    cd=ranged(url,cd_offset,cd_offset+cd_size-1)
    members=[]; off=0
    while off<len(cd):
        if cd[off:off+4]!=CDSIG: raise RuntimeError(f"bad CD signature at {off}")
        fields=struct.unpack_from("<6H3I5H2I",cd,off+4)
        _vm,_vn,flags,method,_mt,_md,crc,csize,usize,nlen,xlen,clen,_disk,_ia,_ea,lhoff=fields
        a=off+46; b=a+nlen; c=b+xlen; d=c+clen
        raw=cd[a:b]; enc="utf-8" if flags&0x800 else "cp437"; name=raw.decode(enc,errors="replace")
        members.append({"filename":name,"flags":flags,"compression_method":method,"crc32":f"{crc:08x}",
                        "compressed_size":csize,"uncompressed_size":usize,"local_header_offset":lhoff})
        off=d
    if len(members)!=total_n: raise RuntimeError(f"entry count mismatch {len(members)} vs {total_n}")
    exact_tif=[m for m in members if m["filename"].lower().endswith(".tif")]
    return {"archive_size_bytes":size,"member_count":len(members),"members":members,
            "exact_tif_member_count":len(exact_tif),"exact_tif_members":exact_tif,
            "tail_sha256":hashlib.sha256(tail).hexdigest(),"central_directory_sha256":hashlib.sha256(cd).hexdigest(),
            "local_file_payload_bytes_consumed":0,"raster_pixels_read":False}

out_members=[]
for item in REC["prospective"]:
    sid=item["SURVEY_ID"]
    if item["batch_state"]!="TRANSPORT_USABLE":
        out_members.append({"rank":item["rank"],"SURVEY_ID":sid,"batch_state":item["batch_state"],
                            "inventory_state":"PRESERVED_TRANSPORT_BLOCK__NOT_REPLACED"})
        continue
    b=inventory(item["bathymetry_url"],item["bathymetry_archive_size_bytes"])
    k=inventory(item["backscatter_url"],item["backscatter_archive_size_bytes"])
    pair_ok=b["exact_tif_member_count"]==1 and k["exact_tif_member_count"]==1
    out_members.append({"rank":item["rank"],"SURVEY_ID":sid,"batch_state":item["batch_state"],
                        "inventory_state":"PASS_UNIQUE_EXACT_TIF_PAIR" if pair_ok else "MEMBER_SELECTION_BLOCKED",
                        "bathymetry":b,"backscatter":k})

out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I1C-BATCH-ARCHIVE-MEMBER-INVENTORY-RECEIPT-2026-09-24-v1.0",
 "prereg":PRE["artifact_id"],"batch":out_members,
 "usable_unique_tif_pair_count":sum(x.get("inventory_state")=="PASS_UNIQUE_EXACT_TIF_PAIR" for x in out_members),
 "transport_blocked_count":sum("TRANSPORT_BLOCK" in x.get("inventory_state","") for x in out_members),
 "local_file_payload_bytes_consumed":0,"raster_pixels_read":False,
 "next_gate":"I1D_BATCH_PAIRED_RASTER_MEMBER_AND_HEADER_GEOMETRY",
 "claim_ceiling":PRE["claim_ceiling"]
}
p=OUT/"JANUS-KUSTO-INFOMAR-I1C-BATCH-ARCHIVE-MEMBER-INVENTORY-RECEIPT-2026-09-24-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps({
 "artifact_id":out["artifact_id"],
 "summary":[{"rank":x["rank"],"SURVEY_ID":x["SURVEY_ID"],"inventory_state":x["inventory_state"],
             "bathy_tif":(x.get("bathymetry") or {}).get("exact_tif_members"),
             "back_tif":(x.get("backscatter") or {}).get("exact_tif_members")} for x in out_members],
 "raster_pixels_read":False
},indent=2,ensure_ascii=False))
