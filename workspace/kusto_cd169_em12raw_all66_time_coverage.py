#!/usr/bin/env python3
import ftplib, hashlib, json, zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

HOST="livftp.noc.ac.uk"
ZIP_PATH="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11280/EM12/EM12raw.zip"
ZIP_SHA="060b65865312330e866f842959a89674b744ff44569dcb1fa13daa7b5b0cefb0"
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
LOCAL=Path("/tmp/CD169_EM12raw.zip")
T0=datetime(2005,2,26,21,16,0,tzinfo=timezone.utc).timestamp()
T1=datetime(2005,2,26,21,21,0,tzinfo=timezone.utc).timestamp()

TYPES={
  0x0285:("START",424,"parameter"),
  0x0286:("STOP",424,"parameter"),
  0x0287:("PARAMETER",424,"parameter"),
  0x0293:("POS",93,"parameter"),
  0x029A:("SVP",419,"survey"),
  0x0294:("EM12DS_BATH",926,"survey"),
  0x0295:("EM12DP_BATH",926,"survey"),
  0x0296:("EM12S_BATH",926,"survey"),
  0x02C8:("EM12DP_SS",554,None),
  0x02C9:("EM12DS_SS",554,None),
  0x02CA:("EM12S_SS",554,None),
  0x02CB:("EM12DP_SSP",1468,None),
  0x02CC:("EM12DS_SSP",1468,None),
  0x02CD:("EM12S_SSP",1468,None)
}
BATH_TYPES={"EM12DS_BATH","EM12DP_BATH","EM12S_BATH"}

h=hashlib.sha256()
f=ftplib.FTP(timeout=180); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I")
with LOCAL.open("wb") as out:
    def cb(b): out.write(b); h.update(b)
    f.retrbinary("RETR "+ZIP_PATH,cb,blocksize=1024*1024)
try:f.quit()
except:f.close()
sha=h.hexdigest()
if sha!=ZIP_SHA: raise RuntimeError(f"ZIP SHA mismatch {sha}")

def digits(b,a,z):
    s=b[a:z]
    if len(s)!=(z-a) or not all(48<=x<=57 for x in s): return None
    return int(s.decode("ascii"))

def parse_time(payload,mode):
    try:
        dd=digits(payload,0,2); mm=digits(payload,2,4); yy=digits(payload,4,6)
        if None in (dd,mm,yy): return None
        if mode=="parameter":
            hh=digits(payload,7,9); mi=digits(payload,9,11); ss=digits(payload,11,13); cs=digits(payload,13,15)
        else:
            hh=digits(payload,6,8); mi=digits(payload,8,10); ss=digits(payload,10,12); cs=digits(payload,12,14)
        if None in (hh,mi,ss,cs): return None
        year=2000+yy if yy<80 else 1900+yy
        dt=datetime(year,mm,dd,hh,mi,ss,cs*10000,tzinfo=timezone.utc)
        return dt.timestamp(),dt.isoformat().replace("+00:00","Z")
    except Exception:
        return None

def scan_entry(data,name):
    records=[]
    # The old format is a label+fixed-payload stream. Search labels and validate timestamps.
    for i in range(len(data)-18):
        if data[i] != 0x02: continue
        typ=(data[i]<<8)|data[i+1]
        spec=TYPES.get(typ)
        if spec is None: continue
        tname,size,mode=spec
        if mode is None or i+2+size>len(data): continue
        parsed=parse_time(data[i+2:i+2+32],mode)
        if parsed is None: continue
        epoch,utc=parsed
        if not (datetime(2005,2,1,tzinfo=timezone.utc).timestamp() <= epoch <= datetime(2005,4,1,tzinfo=timezone.utc).timestamp()):
            continue
        records.append({"offset":i,"type":tname,"epoch":epoch,"utc":utc})
    # exact duplicate candidate labels are not expected, but de-duplicate identical offset/type/time defensively.
    uniq={(r["offset"],r["type"],r["epoch"]):r for r in records}
    records=sorted(uniq.values(),key=lambda r:(r["epoch"],r["offset"],r["type"]))
    bath=[r for r in records if r["type"] in BATH_TYPES]
    win=[r for r in records if T0<=r["epoch"]<=T1]
    bwin=[r for r in bath if T0<=r["epoch"]<=T1]
    types=Counter(r["type"] for r in records)
    return {
      "filename":name,
      "bytes":len(data),
      "timestamped_record_count":len(records),
      "record_type_counts":dict(types),
      "bath_record_count":len(bath),
      "first_any_utc":records[0]["utc"] if records else None,
      "last_any_utc":records[-1]["utc"] if records else None,
      "first_bath_utc":bath[0]["utc"] if bath else None,
      "last_bath_utc":bath[-1]["utc"] if bath else None,
      "records_in_frozen_interval":len(win),
      "bath_records_in_frozen_interval":len(bwin),
      "frozen_interval_records":[{"type":r["type"],"utc":r["utc"],"offset":r["offset"]} for r in win],
      "_records":records
    }

with zipfile.ZipFile(LOCAL,"r") as z:
    infos=z.infolist()
    scans=[]
    all_records=[]
    for info in infos:
        s=scan_entry(z.read(info.filename),info.filename)
        scans.append(s)
        for r in s.pop("_records"):
            all_records.append({"filename":info.filename,**r})

all_records.sort(key=lambda r:(r["epoch"],r["filename"],r["offset"]))
inside=[r for r in all_records if T0<=r["epoch"]<=T1]
bath_inside=[r for r in inside if r["type"] in BATH_TYPES]
before=[r for r in all_records if r["epoch"]<T0]
after=[r for r in all_records if r["epoch"]>T1]
prev=before[-1] if before else None
nxt=after[0] if after else None
gap=None
if prev and nxt:
    gap={
      "last_record_before_interval":{"filename":prev["filename"],"type":prev["type"],"utc":prev["utc"],"offset":prev["offset"]},
      "first_record_after_interval":{"filename":nxt["filename"],"type":nxt["type"],"utc":nxt["utc"],"offset":nxt["offset"]},
      "record_to_record_gap_s":nxt["epoch"]-prev["epoch"],
      "missing_from_last_record_to_interval_start_s":T0-prev["epoch"],
      "missing_from_interval_end_to_first_record_s":nxt["epoch"]-T1
    }

intersecting=sorted({r["filename"] for r in inside})
bath_intersecting=sorted({r["filename"] for r in bath_inside})
verdict="RAW_INTERVAL_PRESENT" if inside else "RAW_INTERVAL_MISSING_FROM_ACCESSIBLE_ARCHIVE"
bath_verdict="RAW_BATH_INTERVAL_PRESENT" if bath_inside else "RAW_BATH_INTERVAL_MISSING_FROM_ACCESSIBLE_ARCHIVE"

out={
 "artifact_id":"JANUS-KUSTO-CD169-EM12RAW-ALL66-TIME-COVERAGE-AUDIT-RUN-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CD169-EM12RAW-ALL66-TIME-COVERAGE-AUDIT-PREREG-2026-09-23-v1.0.json",
 "source":{"ftp_path":ZIP_PATH,"sha256":sha,"entry_count":len(scans)},
 "parser_authority":{"format":"MBF_EMOLDRAW 51","source":"MB-System mbsys_simrad.h + mbr_emoldraw.c"},
 "frozen_interval_utc":["2005-02-26T21:16:00Z","2005-02-26T21:21:00Z"],
 "entries":scans,
 "all_entries_intersecting_frozen_interval":intersecting,
 "all_bath_entries_intersecting_frozen_interval":bath_intersecting,
 "frozen_interval_record_count_all66":len(inside),
 "frozen_interval_bath_record_count_all66":len(bath_inside),
 "bracketing_raw_archive_gap":gap,
 "verdict":verdict,
 "bath_verdict":bath_verdict,
 "depth_or_beam_values_read":False,
 "modern_all_parser_run_35792379449":"VOID_PARSER_MISMATCH",
 "claim_ceiling":"COMPLETE_ACCESSIBLE_RAW_ARCHIVE_TIME_COVERAGE_AUDIT_ONLY"
}
p=OUT/"JANUS-KUSTO-CD169-EM12RAW-ALL66-TIME-COVERAGE-AUDIT-RUN-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "source":out["source"],
 "frozen_interval_utc":out["frozen_interval_utc"],
 "all_entries_intersecting_frozen_interval":intersecting,
 "all_bath_entries_intersecting_frozen_interval":bath_intersecting,
 "frozen_interval_record_count_all66":len(inside),
 "frozen_interval_bath_record_count_all66":len(bath_inside),
 "bracketing_raw_archive_gap":gap,
 "verdict":verdict,
 "bath_verdict":bath_verdict,
 "entry_envelopes":[{k:s[k] for k in ["filename","timestamped_record_count","bath_record_count","first_any_utc","last_any_utc","first_bath_utc","last_bath_utc"]} for s in scans]
},indent=2))
