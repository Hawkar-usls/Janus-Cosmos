#!/usr/bin/env python3
import hashlib,json,requests
from pathlib import Path
from kusto_mbldeoih_reader import summarize

OUT=Path("workspace/kusto_cross_survey_out");OUT.mkdir(parents=True,exist_ok=True)
FILES={
 "PF0501":"https://data.ngdc.noaa.gov/platforms/ocean/ships/pathfinder/PF0501/multibeam/data/version1/MB/generated/Atlantic_line_269_raw.all.mb51.fbt",
 "RB1202":"https://data.ngdc.noaa.gov/platforms/ocean/ships/ronald_h._brown/RB1202/multibeam/data/version1/MB/generated/0040_20120707_085106_RHB.all.mb58.fbt"
}
res={}
for sid,u in FILES.items():
    r=requests.get(u,timeout=180,headers={"User-Agent":"JANUS-KUSTO-MBLDEOIH-selftest/1.0"});r.raise_for_status()
    b=r.content;s=summarize(b)
    s.update({"url":u,"sha256":hashlib.sha256(b).hexdigest()})
    # source-format implementation checks only, not morphology thresholds
    nav=s["nav_bounds"];depth=s["good_depth_range_m"]
    plausible=bool(s["eof_exact"] and s["data_records"]>0 and nav and
                   -180<=nav[0]<=180 and -90<=nav[1]<=90 and -180<=nav[2]<=180 and -90<=nav[3]<=90 and
                   depth and -12000<depth[0]<12000 and -12000<depth[1]<12000)
    s["implementation_selftest_pass"]=plausible
    res[sid]=s
out={
 "artifact_id":"JANUS-KUSTO-MBLDEOIH-V1-V5-PARSER-SELFTEST-2026-09-23-v1.0",
 "authority":"MB-System src/mbio/mbr_mbldeoih.c",
 "scientific_scoring_performed":False,
 "results":res,
 "all_pass":all(x["implementation_selftest_pass"] for x in res.values()),
 "claim_ceiling":"PARSER_IMPLEMENTATION_SELFTEST_ONLY__NO_UNSEEN_SCIENTIFIC_RESULT"
}
p=OUT/"JANUS-KUSTO-MBLDEOIH-V1-V5-PARSER-SELFTEST-2026-09-23-v1.0.json";p.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
if not out["all_pass"]:raise SystemExit(2)
