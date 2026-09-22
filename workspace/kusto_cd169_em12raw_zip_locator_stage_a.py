#!/usr/bin/env python3
import ftplib, hashlib, json, os, re, zipfile
from pathlib import Path

HOST="livftp.noc.ac.uk"
ZIP_PATH="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11280/EM12/EM12raw.zip"
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
LOCAL=Path("/tmp/CD169_EM12raw.zip")
WINDOW=("2005-02-26T21:16:00Z","2005-02-26T21:21:00Z")

# download exact bytes and hash as they arrive
h=hashlib.sha256(); total=0
f=ftplib.FTP(timeout=180); f.connect(HOST,21); f.login("anonymous","janus-kusto@example.invalid"); f.voidcmd("TYPE I")
with LOCAL.open("wb") as out:
    def cb(chunk):
        global total
        out.write(chunk); h.update(chunk)
    f.retrbinary("RETR "+ZIP_PATH, cb, blocksize=1024*1024)
try:f.quit()
except:f.close()
total=LOCAL.stat().st_size
sha=h.hexdigest()

with zipfile.ZipFile(LOCAL,"r") as z:
    infos=z.infolist()
    entries=[]
    for i,info in enumerate(infos):
        entries.append({
          "index":i,"filename":info.filename,
          "compressed_size":info.compress_size,
          "uncompressed_size":info.file_size,
          "crc32":f"{info.CRC:08x}",
          "zip_datetime":"%04d-%02d-%02dT%02d:%02d:%02d" % info.date_time,
          "compression_type":info.compress_type
        })
    bad=z.testzip()

# filename-only date/time extraction; deliberately no entry contents read for selection.
patterns=[
 re.compile(r'(?P<y>2005)[^0-9]?(?P<m>0?2)[^0-9]?(?P<d>26)[^0-9]?(?P<h>21)[^0-9]?(?P<mi>1[6-9]|20|21)',re.I),
 re.compile(r'(?P<y>05)(?P<m>02)(?P<d>26)[^0-9]?(?P<h>21)(?P<mi>[0-5][0-9])',re.I),
 re.compile(r'(?P<doy>0?57)[^0-9]?(?P<h>21)[^0-9]?(?P<mi>[0-5][0-9])',re.I)
]
filename_time_candidates=[]
for e in entries:
    hits=[]
    for p in patterns:
        m=p.search(e["filename"])
        if m:hits.append(m.groupdict())
    if hits:
        filename_time_candidates.append({"entry":e,"filename_time_matches":hits})

# Also report likely numbering structure and gaps purely from names.
numbers=[]
for e in entries:
    stem=Path(e["filename"]).name
    ms=re.findall(r'(?<!\d)(\d{1,4})(?!\d)',stem)
    numbers.append({"filename":e["filename"],"numeric_tokens":[int(x) for x in ms]})

receipt={
 "artifact_id":"JANUS-KUSTO-CD169-EM12RAW-ZIP-TIME-LOCATOR-STAGE-A-2026-09-23-v1.0",
 "prereg":"data/cousteau/JANUS-KUSTO-CD169-EM12RAW-ZIP-TIME-LOCATOR-PREREG-2026-09-23-v1.0.json",
 "source":{"ftp_path":ZIP_PATH,"bytes":total,"sha256":sha},
 "zip":{"entry_count":len(entries),"testzip_bad_entry":bad,"integrity_pass":bad is None},
 "entries":entries,
 "filename_time_candidates_for_frozen_window":filename_time_candidates,
 "filename_numeric_tokens":numbers,
 "depth_or_bathymetry_values_inspected":False,
 "raw_entry_contents_inspected_for_selection":False,
 "frozen_locator_window_utc":list(WINDOW),
 "stage_A_status":"PASS_ZIP_DIRECTORY_FROZEN" if bad is None else "FAIL_ZIP_INTEGRITY",
 "next":"If filenames encode target time, freeze exact entry. Otherwise preregister header/time-only scan of entry records without reading depth.",
 "claim_ceiling":"RAW_ARCHIVE_DIRECTORY_AND_FILENAME_TIME_PROVENANCE_ONLY"
}
p=OUT/"JANUS-KUSTO-CD169-EM12RAW-ZIP-TIME-LOCATOR-STAGE-A-2026-09-23-v1.0.json"
p.write_text(json.dumps(receipt,indent=2))
print(json.dumps({
 "source":receipt["source"],"zip":receipt["zip"],
 "entries":entries,
 "filename_time_candidates_for_frozen_window":filename_time_candidates,
 "stage_A_status":receipt["stage_A_status"]
},indent=2))
