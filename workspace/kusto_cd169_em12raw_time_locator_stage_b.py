#!/usr/bin/env python3
import ftplib, hashlib, json, struct, zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

HOST="livftp.noc.ac.uk"
ZIP_PATH="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11280/EM12/EM12raw.zip"
EXPECTED_ZIP_SHA="060b65865312330e866f842959a89674b744ff44569dcb1fa13daa7b5b0cefb0"
ENTRIES=["0030_260205_082242_raw.all","0031_260205_191055_raw.all","0032_270205_144041_raw.all"]
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
LOCAL=Path("/tmp/CD169_EM12raw.zip")
TARGET_START=datetime(2005,2,26,21,16,0,tzinfo=timezone.utc).timestamp()
TARGET_END=datetime(2005,2,26,21,21,0,tzinfo=timezone.utc).timestamp()

h=hashlib.sha256()
f=ftplib.FTP(timeout=180);f.connect(HOST,21);f.login("anonymous","janus-kusto@example.invalid");f.voidcmd("TYPE I")
with LOCAL.open("wb") as out:
    def cb(b): out.write(b);h.update(b)
    f.retrbinary("RETR "+ZIP_PATH,cb,blocksize=1024*1024)
try:f.quit()
except:f.close()
sha=h.hexdigest()
if sha!=EXPECTED_ZIP_SHA:raise RuntimeError(f"ZIP SHA mismatch {sha}")

def date_epoch(date,time_ms):
    y=date//10000;m=(date//100)%100;d=date%100
    hh=time_ms//3600000; rem=time_ms%3600000
    mm=rem//60000; rem%=60000
    ss=rem//1000; ms=rem%1000
    dt=datetime(y,m,d,hh,mm,ss,ms*1000,tzinfo=timezone.utc)
    return dt.timestamp(),dt.isoformat().replace("+00:00","Z")

def plausible(date,time_ms,model):
    y=date//10000;m=(date//100)%100;d=date%100
    return 2000<=y<=2010 and 1<=m<=12 and 1<=d<=31 and 0<=time_ms<86400000 and 1<=model<=10000

def parse_entry(data,name):
    pos=0; records=[]; errors=[]; endian_counts=Counter()
    n=len(data)
    while pos+16<=n:
        if data[pos+4]!=2:
            # bounded bytewise resync looking for [len][STX]
            found=None
            for q in range(pos+1,min(n-16,pos+4096)):
                if data[q+4]==2:
                    for endian in ("<",">"):
                        L=struct.unpack_from(endian+"I",data,q)[0]
                        if 12<=L<=2_000_000 and q+4+L<=n:
                            model=struct.unpack_from(endian+"H",data,q+6)[0]
                            date=struct.unpack_from(endian+"I",data,q+8)[0]
                            tm=struct.unpack_from(endian+"I",data,q+12)[0]
                            if plausible(date,tm,model):
                                found=(q,endian,L);break
                    if found:break
            if not found:
                errors.append({"offset":pos,"reason":"no_resync_within_4096"});break
            errors.append({"offset":pos,"reason":"resync","next_offset":found[0]})
            pos=found[0]
        choices=[]
        for endian in ("<",">"):
            L=struct.unpack_from(endian+"I",data,pos)[0]
            if not (12<=L<=2_000_000 and pos+4+L<=n):continue
            if data[pos+4]!=2:continue
            typ=data[pos+5]
            model=struct.unpack_from(endian+"H",data,pos+6)[0]
            date=struct.unpack_from(endian+"I",data,pos+8)[0]
            tm=struct.unpack_from(endian+"I",data,pos+12)[0]
            if plausible(date,tm,model):
                choices.append((endian,L,typ,model,date,tm))
        if not choices:
            errors.append({"offset":pos,"reason":"no_plausible_header"});pos+=1;continue
        endian,L,typ,model,date,tm=choices[0]
        epoch,iso=date_epoch(date,tm)
        records.append({"offset":pos,"length":L,"type_byte":typ,"type_chr":chr(typ) if 32<=typ<=126 else None,
                        "model":model,"date":date,"time_ms":tm,"epoch":epoch,"utc":iso})
        endian_counts[endian]+=1
        pos+=4+L
    ts=[r["epoch"] for r in records]
    window=[r for r in records if TARGET_START<=r["epoch"]<=TARGET_END]
    # unique timestamp gaps, independent of datagram multiplicity at same ping/time
    uts=sorted(set(ts))
    gaps=[]
    for a,b in zip(uts,uts[1:]):
        if b-a>=60:
            gaps.append({"from_utc":datetime.fromtimestamp(a,tz=timezone.utc).isoformat().replace("+00:00","Z"),
                         "to_utc":datetime.fromtimestamp(b,tz=timezone.utc).isoformat().replace("+00:00","Z"),
                         "gap_s":b-a})
    type_counts=Counter((r["type_chr"] or f"0x{r['type_byte']:02x}") for r in records)
    return {
      "filename":name,"bytes":n,
      "datagram_count":len(records),"endian_counts":dict(endian_counts),
      "first_utc":records[0]["utc"] if records else None,
      "last_utc":records[-1]["utc"] if records else None,
      "min_utc":datetime.fromtimestamp(min(ts),tz=timezone.utc).isoformat().replace("+00:00","Z") if ts else None,
      "max_utc":datetime.fromtimestamp(max(ts),tz=timezone.utc).isoformat().replace("+00:00","Z") if ts else None,
      "type_counts":dict(type_counts),
      "target_window_datagram_count":len(window),
      "target_window_first":window[0] if window else None,
      "target_window_last":window[-1] if window else None,
      "gaps_ge_60s":gaps,
      "parse_errors":errors[:100],
      "parse_error_count":len(errors),
      "depth_fields_read":False
    }

with zipfile.ZipFile(LOCAL,"r") as z:
    results=[]
    for name in ENTRIES:
        b=z.read(name)
        results.append(parse_entry(b,name))

r31=next(x for x in results if x["filename"].startswith("0031_"))
coverage=r31["target_window_datagram_count"]>0
out={
 "artifact_id":"JANUS-KUSTO-CD169-EM12RAW-ZIP-TIME-LOCATOR-STAGE-B-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CD169-EM12RAW-ZIP-TIME-LOCATOR-PREREG-2026-09-23-v1.0.json",
 "zip":{"path":ZIP_PATH,"sha256":sha},
 "frozen_window_utc":["2005-02-26T21:16:00Z","2005-02-26T21:21:00Z"],
 "entries_scanned":results,
 "raw_target_window_present_in_0031":coverage,
 "locator_verdict":"PASS_RAW_FILE_0031_CONTAINS_FROZEN_INTERVAL" if coverage else "FAIL_RAW_0031_DOES_NOT_CONTAIN_FROZEN_INTERVAL",
 "depth_or_beam_values_read":False,
 "next":"If PASS, freeze exact raw datagram types/timestamps and preregister beam-level raw forensics for the two dominant cells. If FAIL, preserve missing/corrupt archival interval.",
 "claim_ceiling":"RAW_DATAGRAM_TIME_AND_TYPE_PROVENANCE_ONLY"
}
p=OUT/"JANUS-KUSTO-CD169-EM12RAW-ZIP-TIME-LOCATOR-STAGE-B-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
