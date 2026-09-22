#!/usr/bin/env python3
import ftplib, hashlib, json
from pathlib import Path
HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281"
OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
FILES=[
 ROOT+"/README.txt",
 ROOT+"/Nav/cd169leg1_1min.listit",
 ROOT+"/Nav/cd169leg1_5min.listit"
]
def fetch(path):
    f=ftplib.FTP(timeout=120);f.connect(HOST,21);f.login("anonymous","janus-kusto@example.invalid");f.voidcmd("TYPE I")
    chunks=[]
    try:f.retrbinary("RETR "+path,chunks.append)
    finally:
        try:f.quit()
        except:f.close()
    b=b"".join(chunks)
    return {"path":path,"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest(),"preview":b.decode("utf-8","replace")[:30000]}
rows=[fetch(p) for p in FILES]
out={"artifact_id":"JANUS-KUSTO-CD169-README-NAV-FORMAT-PROBE-2026-09-22-v1.0","files":rows,"claim_ceiling":"FORMAT_AND_METADATA_PROBE_ONLY"}
p=OUT/"JANUS-KUSTO-CD169-README-NAV-FORMAT-PROBE-2026-09-22-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
