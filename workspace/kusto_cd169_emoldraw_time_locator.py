#!/usr/bin/env python3
import ftplib, hashlib, json, zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

HOST="livftp.noc.ac.uk"
ZIP_PATH="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11280/EM12/EM12raw.zip"
ZIP_SHA="060b65865312330e866f842959a89674b744ff44569dcb1fa13daa7b5b0cefb0"
ENTRIES=["0030_260205_082242_raw.all","0031_260205_191055_raw.all","0032_270205_144041_raw.all"]
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
LOCAL=Path("/tmp/CD169_EM12raw.zip")
T0=datetime(2005,2,26,21,16,0,tzinfo=timezone.utc).timestamp()
T1=datetime(2005,2,26,21,21,0,tzinfo=timezone.utc).timestamp()

# MB-System mbsys_simrad.h / mbr_emoldraw.c source-bound constants.
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

h=hashlib.sha256()
f=ftplib.FTP(timeout=180);f.connect(HOST,21);f.login("anonymous","janus-kusto@example.invalid");f.voidcmd("TYPE I")
with LOCAL.open("wb") as out:
    def cb(b):out.write(b);h.update(b)
    f.retrbinary("RETR "+ZIP_PATH,cb,blocksize=1024*1024)
try:f.quit()
except:f.close()
if h.hexdigest()!=ZIP_SHA:raise RuntimeError("ZIP SHA mismatch")

def n2(b,a,z):
    s=b[a:z]
    if len(s)!=(z-a) or not all(48<=x<=57 for x in s):return None
    return int(s.decode("ascii"))

def parse_time(payload,mode):
    try:
        dd=n2(payload,0,2);mm=n2(payload,2,4);yy=n2(payload,4,6)
        if None in (dd,mm,yy):return None
        if mode=="parameter":
            hh=n2(payload,7,9);mi=n2(payload,9,11);ss=n2(payload,11,13);cs=n2(payload,13,15)
        else:
            hh=n2(payload,6,8);mi=n2(payload,8,10);ss=n2(payload,10,12);cs=n2(payload,12,14)
        if None in (hh,mi,ss,cs):return None
        year=2000+yy if yy<80 else 1900+yy
        dt=datetime(year,mm,dd,hh,mi,ss,cs*10000,tzinfo=timezone.utc)
        return dt.timestamp(),dt.isoformat().replace("+00:00","Z")
    except:return None

def scan(data,name):
    records=[]
    # Scan valid big-endian two-byte labels; time validation suppresses payload false positives.
    for i in range(len(data)-16):
        if data[i]!=0x02:continue
        typ=(data[i]<<8)|data[i+1]
        if typ not in TYPES:continue
        tname,size,mode=TYPES[typ]
        if mode is None:continue
        if i+2+size>len(data):continue
        parsed=parse_time(data[i+2:i+2+32],mode)
        if parsed is None:continue
        epoch,iso=parsed
        if not (datetime(2005,2,1,tzinfo=timezone.utc).timestamp() <= epoch <= datetime(2005,4,1,tzinfo=timezone.utc).timestamp()):
            continue
        records.append({"offset":i,"type_hex":f"0x{typ:04x}","type":tname,"record_size_after_label":size,"epoch":epoch,"utc":iso})
    records.sort(key=lambda r:(r["epoch"],r["offset"]))
    bath=[r for r in records if "BATH" in r["type"]]
    window=[r for r in records if T0<=r["epoch"]<=T1]
    bwin=[r for r in bath if T0<=r["epoch"]<=T1]
    tc=Counter(r["type"] for r in records)
    return {
      "filename":name,"bytes":len(data),"valid_timestamped_record_count":len(records),
      "record_type_counts":dict(tc),
      "first_utc":records[0]["utc"] if records else None,
      "last_utc":records[-1]["utc"] if records else None,
      "first_bath_utc":bath[0]["utc"] if bath else None,
      "last_bath_utc":bath[-1]["utc"] if bath else None,
      "target_window_any_record_count":len(window),
      "target_window_bath_record_count":len(bwin),
      "target_window_any_first":window[0] if window else None,
      "target_window_any_last":window[-1] if window else None,
      "target_window_bath_first":bwin[0] if bwin else None,
      "target_window_bath_last":bwin[-1] if bwin else None,
      "depth_fields_read":False,
      "bath_records_in_window_metadata_only":[{"offset":r["offset"],"type":r["type"],"utc":r["utc"]} for r in bwin]
    }

with zipfile.ZipFile(LOCAL,"r") as z:
    results=[scan(z.read(n),n) for n in ENTRIES]

r31=next(r for r in results if r["filename"].startswith("0031_"))
passed=r31["target_window_bath_record_count"]>0
out={
 "artifact_id":"JANUS-KUSTO-CD169-EM12RAW-EMOLDRAW-TIME-LOCATOR-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CD169-EM12RAW-ZIP-TIME-LOCATOR-PREREG-2026-09-23-v1.0.json",
 "supersedes_stage_b_parser_interpretation":{
   "run_id":35792379449,
   "reason":"Technical parser mismatch: attempted modern Kongsberg ALL layout on old Simrad MBF_EMOLDRAW bytes.",
   "scientific_verdict_from_superseded_run":"VOID"
 },
 "parser_authority":{
   "source":"MB-System src/mbio/mbsys_simrad.h + src/mbio/mbr_emoldraw.c",
   "format":"MBF_EMOLDRAW 51",
   "em12s_bath_label":"0x0296",
   "em12_bath_record_size_bytes":926,
   "timestamp":"ASCII DDMMYYHHMMSScc at start of bath payload after two-byte label"
 },
 "zip":{"path":ZIP_PATH,"sha256":ZIP_SHA},
 "frozen_window_utc":["2005-02-26T21:16:00Z","2005-02-26T21:21:00Z"],
 "entries":results,
 "raw_bath_records_present_in_0031_frozen_window":passed,
 "verdict":"PASS_RAW_EM12_BATH_RECORDS_PRESENT_IN_FROZEN_WINDOW" if passed else "FAIL_NO_RAW_EM12_BATH_RECORDS_IN_FROZEN_WINDOW",
 "depth_or_beam_values_read":False,
 "next":"If PASS: freeze exact bath ping offsets/timestamps and preregister raw beam-level cell forensics. If FAIL: preserve missing raw interval.",
 "claim_ceiling":"SOURCE_EXACT_RAW_EM12_TIMESTAMP_AND_RECORD_TYPE_PROVENANCE_ONLY"
}
p=OUT/"JANUS-KUSTO-CD169-EM12RAW-EMOLDRAW-TIME-LOCATOR-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
