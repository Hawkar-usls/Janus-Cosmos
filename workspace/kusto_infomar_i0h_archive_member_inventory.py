#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, struct
from pathlib import Path
from urllib.parse import urlparse
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"workspace"/"kusto_global_groundtruth_out"; OUT.mkdir(parents=True,exist_ok=True)
PR=json.loads((ROOT/"data/cousteau/JANUS-KUSTO-INFOMAR-I0H-ARCHIVE-MEMBER-INVENTORY-PREREG-2026-09-24-v1.0.json").read_text())
S=requests.Session()
S.headers.update({"User-Agent":"JANUS-KUSTO-INFOMAR-I0H/1.0","Accept-Encoding":"identity"})

EOCD_SIG=b"PK\\x05\\x06"
CD_SIG=b"PK\\x01\\x02"

def ranged(url,start,end):
    r=S.get(url,headers={"Range":f"bytes={start}-{end}"},timeout=180,allow_redirects=True)
    r.raise_for_status()
    if r.status_code != 206:
        raise RuntimeError(f"range request did not return 206: {r.status_code} {url}")
    return r.content, dict(r.headers), r.url

def archive_size(url):
    r=S.get(url,headers={"Range":"bytes=0-0"},timeout=120,allow_redirects=True,stream=True)
    if r.status_code != 206:
        raise RuntimeError(f"size probe did not return 206: {r.status_code} {url}")
    cr=r.headers.get("Content-Range") or r.headers.get("content-range") or ""
    if "/" not in cr:
        raise RuntimeError(f"missing Content-Range total: {cr}")
    return int(cr.rsplit("/",1)[1]), dict(r.headers), r.url

def parse_archive(label,spec):
    url=spec["url"]
    size,h0,resolved=archive_size(url)
    expected=int(spec["expected_size_bytes"])
    if size != expected:
        raise RuntimeError(f"{label} size drift: got {size}, frozen {expected}")
    tail_n=min(size, 65557)
    tail_start=size-tail_n
    tail,ht,_=ranged(url,tail_start,size-1)
    p=tail.rfind(EOCD_SIG)
    if p < 0 or p+22 > len(tail):
        raise RuntimeError(f"{label}: EOCD not found")
    disk_no,cd_disk,disk_entries,total_entries,cd_size,cd_offset,comment_len=struct.unpack_from("<HHHHIIH",tail,p+4)
    if disk_no!=0 or cd_disk!=0 or disk_entries!=total_entries:
        raise RuntimeError(f"{label}: multi-disk ZIP unsupported")
    if total_entries==0xffff or cd_size==0xffffffff or cd_offset==0xffffffff:
        raise RuntimeError(f"{label}: ZIP64 requires explicit next gate")
    if p+22+comment_len > len(tail):
        raise RuntimeError(f"{label}: incomplete EOCD comment")
    cd, hc, _ = ranged(url, cd_offset, cd_offset+cd_size-1)
    members=[]
    off=0
    while off < len(cd):
        if cd[off:off+4] != CD_SIG:
            raise RuntimeError(f"{label}: bad central-directory signature at {off}")
        if off+46 > len(cd):
            raise RuntimeError(f"{label}: truncated central-directory header")
        fields=struct.unpack_from("<6H3I5H2I",cd,off+4)
        _ver_made,_ver_need,flags,method,_mtime,_mdate,crc,csize,usize,nlen,xlen,clen,_disk,_iattr,_eattr,lhoff=fields
        a=off+46; b=a+nlen; c=b+xlen; d=c+clen
        if d>len(cd):
            raise RuntimeError(f"{label}: truncated central-directory variable fields")
        rawname=cd[a:b]
        enc="utf-8" if (flags & 0x800) else "cp437"
        name=rawname.decode(enc,errors="replace")
        suffixes=[s.lower() for s in Path(name).suffixes]
        members.append({
            "filename":name,
            "compression_method":method,
            "flags":flags,
            "crc32":f"{crc:08x}",
            "compressed_size":csize,
            "uncompressed_size":usize,
            "local_header_offset":lhoff,
            "filename_suffix":"".join(suffixes[-2:]) if len(suffixes)>=2 and suffixes[-1] in (".xml",".ovr") else (suffixes[-1] if suffixes else "")
        })
        off=d
    if len(members)!=total_entries:
        raise RuntimeError(f"{label}: entry count mismatch parsed={len(members)} eocd={total_entries}")
    consumed=len(tail)+len(cd)+1
    return {
      "requested_url":url,
      "resolved_url":resolved,
      "archive_size_bytes":size,
      "archive_size_matches_frozen":True,
      "eocd":{"absolute_offset":tail_start+p,"entry_count":total_entries,"central_directory_size":cd_size,"central_directory_offset":cd_offset,"comment_length":comment_len},
      "members":members,
      "member_count":len(members),
      "range_bytes_consumed":consumed,
      "range_fraction_of_archive":consumed/size,
      "tail_sha256":hashlib.sha256(tail).hexdigest(),
      "central_directory_sha256":hashlib.sha256(cd).hexdigest(),
      "local_file_payload_bytes_consumed":0,
      "raster_pixels_read":False
    }

results={k:parse_archive(k,v) for k,v in PR["archives"].items()}
passed=all(x["member_count"]>0 and x["local_file_payload_bytes_consumed"]==0 for x in results.values())
out={
 "artifact_id":"JANUS-KUSTO-INFOMAR-I0H-ARCHIVE-MEMBER-INVENTORY-RECEIPT-2026-09-24-v1.0",
 "prereg":PR["artifact_id"],
 "survey_id":PR["survey_id"],
 "result":"PASS" if passed else "FAIL",
 "archives":results,
 "archive_full_bodies_consumed":False,
 "local_file_payload_bytes_consumed":0,
 "raster_pixels_read":False,
 "next_gate":PR["next_gate_on_pass"] if passed else PR["next_gate_on_fail"],
 "claim_ceiling":PR["claim_ceiling"]
}
p=OUT/"JANUS-KUSTO-INFOMAR-I0H-ARCHIVE-MEMBER-INVENTORY-RECEIPT-2026-09-24-v1.0.json"
p.write_text(json.dumps(out,indent=2,ensure_ascii=False))
print(json.dumps(out,indent=2,ensure_ascii=False))
