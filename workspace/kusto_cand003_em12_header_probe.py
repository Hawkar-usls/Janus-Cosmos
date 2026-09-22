#!/usr/bin/env python3
import ftplib, hashlib, json
from pathlib import Path

HOST="livftp.noc.ac.uk"
ROOT="/bodc/bodc/data/BODCREQ-9406/CD169_TOBI/sd11281/EM12"
FILES=["B1-81-1_Acceptl15-22.xyz.ascii","B1-81-1_Acceptl28-33.xyz.ascii"]
OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)

f=ftplib.FTP(timeout=120);f.connect(HOST,21);f.login("anonymous","janus-kusto@example.invalid");f.voidcmd("TYPE I")
rows={}
try:
    for name in FILES:
        data=bytearray();f.retrbinary("RETR "+ROOT+"/"+name,data.extend)
        b=bytes(data);txt=b.decode("utf-8","replace")
        lines=txt.splitlines()
        rows[name]={
          "bytes":len(b),
          "sha256":hashlib.sha256(b).hexdigest(),
          "first_40_lines":lines[:40],
          "header_lines":[x for x in lines[:100] if (x.strip().startswith("***") or ":" in x or "accepted" in x.lower() or "cell size" in x.lower())]
        }
finally:
    try:f.quit()
    except:f.close()
out={"artifact_id":"JANUS-KUSTO-CAND003-CONFLICTING-EM12-PRODUCT-HEADER-PROBE-2026-09-22-v1.0","files":rows,"claim_ceiling":"SOURCE_HEADER_COMPARISON_ONLY"}
p=OUT/"JANUS-KUSTO-CAND003-CONFLICTING-EM12-PRODUCT-HEADER-PROBE-2026-09-22-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
