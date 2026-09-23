#!/usr/bin/env python3
import io, json, re, tarfile, requests, time, xml.etree.ElementTree as ET
from urllib.parse import urljoin
from pathlib import Path

OUT=Path("workspace/kusto_cross_survey_out"); OUT.mkdir(parents=True,exist_ok=True)
UA={"User-Agent":"JANUS-KUSTO-unseen-ex2105-file-metadata/1.0"}
REPORT="https://www.ngdc.noaa.gov/ships/okeanos_explorer/EX2105_mb.html"

def get(u,stream=False):
    last=None
    for i in range(6):
        try:
            r=requests.get(u,headers=UA,timeout=90,stream=stream)
            if r.status_code==429:
                time.sleep(min(20,2**i)); continue
            r.raise_for_status(); return r
        except Exception as e:
            last=e
            if i==5: raise
            time.sleep(min(20,2**i))
    raise last

html=get(REPORT).text
m=re.search(r'href=["\']([^"\']*v1_file_metadata\.tar\.gz)["\']',html,re.I)
if not m: raise RuntimeError("v1_file_metadata.tar.gz link not found")
url=urljoin(REPORT,m.group(1))
blob=get(url).content

members=[]
records=[]
with tarfile.open(fileobj=io.BytesIO(blob),mode="r:gz") as tf:
    for ti in tf.getmembers():
        if not ti.isfile(): continue
        b=tf.extractfile(ti).read()
        text=b.decode("utf-8","replace")
        members.append({"name":ti.name,"bytes":len(b)})
        rec={"member":ti.name,"text_preview":text[:6000]}
        # Generic extraction from XML/text metadata, no depth reading.
        vals={}
        patterns={
          "west":[r'westBoundLongitude[^>]*>([-+0-9.eE]+)<',r'west[^:=]*[:=]\s*([-+0-9.eE]+)'],
          "east":[r'eastBoundLongitude[^>]*>([-+0-9.eE]+)<',r'east[^:=]*[:=]\s*([-+0-9.eE]+)'],
          "south":[r'southBoundLatitude[^>]*>([-+0-9.eE]+)<',r'south[^:=]*[:=]\s*([-+0-9.eE]+)'],
          "north":[r'northBoundLatitude[^>]*>([-+0-9.eE]+)<',r'north[^:=]*[:=]\s*([-+0-9.eE]+)'],
        }
        for k,ps in patterns.items():
            for p in ps:
                mm=re.search(p,text,re.I|re.S)
                if mm:
                    try: vals[k]=float(mm.group(1)); break
                    except: pass
        # Raw/source filename candidates.
        fnames=sorted(set(re.findall(r'([0-9]{4}_[0-9]{8}_[0-9]{6}_EX2105_MB\.kmall(?:\.gz)?)',text,re.I)))
        if not fnames:
            fnames=sorted(set(re.findall(r'([A-Za-z0-9_.-]+\.kmall(?:\.gz)?)',text,re.I)))
        rec["raw_filenames"]=fnames
        rec["bounds"]=vals if len(vals)==4 else None
        # time strings only
        times=sorted(set(re.findall(r'20(?:21)[-T:0-9Z+.]+',text)))
        rec["time_tokens"]=times[:20]
        records.append(rec)

out={
 "artifact_id":"JANUS-KUSTO-UNSEEN-EX2105-FILE-METADATA-PROBE-2026-09-23-v1.0",
 "parent_prereg":"data/cousteau/JANUS-KUSTO-UNSEEN-MULTIBEAM-PF0501-X-EX2105-SOURCE-INVENTORY-PREREG-2026-09-23-v1.0.json",
 "metadata_url":url,
 "archive_bytes":len(blob),
 "depth_values_read":False,
 "member_count":len(members),
 "members":members,
 "records":records,
 "records_with_bounds":sum(r["bounds"] is not None for r in records),
 "records_with_raw_filenames":sum(bool(r["raw_filenames"]) for r in records),
 "claim_ceiling":"UNSEEN_EX2105_FILE_LEVEL_METADATA_ONLY"
}
p=OUT/"JANUS-KUSTO-UNSEEN-EX2105-FILE-METADATA-PROBE-2026-09-23-v1.0.json"
p.write_text(json.dumps(out,indent=2))
print(json.dumps({
 "metadata_url":url,
 "archive_bytes":len(blob),
 "member_count":len(members),
 "records_with_bounds":out["records_with_bounds"],
 "records_with_raw_filenames":out["records_with_raw_filenames"],
 "sample_records":[{"member":r["member"],"raw_filenames":r["raw_filenames"],"bounds":r["bounds"],"time_tokens":r["time_tokens"][:4]} for r in records[:20]]
},indent=2))
