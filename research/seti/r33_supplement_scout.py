#!/usr/bin/env python3
"""Bounded official supplementary-source retrieval; no inferred membership labels."""
import hashlib,json,pathlib,re,requests
from pypdf import PdfReader
P=pathlib.Path("research/seti/r33_verified_families_protocol.json");S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/r33");OUT.mkdir(parents=True,exist_ok=True)
def main():
    result={"url":S["publication"]["supplement"],"status":"BLOCKED","admission":"NO_AUTOMATIC_LABELS"}
    try:
        f=OUT/"official_supplement.pdf";limit=S["supplement_scouter"]["limit_bytes"];n=0;sha=hashlib.sha256()
        with requests.get(result["url"],stream=True,timeout=(15,60),headers={"User-Agent":"JANUS-COSMOS/R33 academic source verification"}) as r:
            r.raise_for_status()
            assert r.url.startswith("https://static-content.springer.com/")
            with f.open("wb") as w:
                for b in r.iter_content(1024**2):
                    n+=len(b);assert n<=limit;sha.update(b);w.write(b)
        assert f.read_bytes()[:4]==b"%PDF"
        reader=PdfReader(f);pages=[]
        for i,p in enumerate(reader.pages):
            t=p.extract_text() or ""
            hits=re.findall(r"(?i)\b(?:supplementary\s+)?(?:figure|fig\.?|table)\s*(?:[Ss])?\s*(\d+)",t)
            floats=re.findall(r"(?<!\d)(?:[789]\d{2}|[123]\d{3})\.\d{3,9}(?!\d)",t)
            if hits or floats:pages.append({"page":i+1,"caption_numbers":sorted(set(hits)),"numeric_frequency_like_tokens":sorted(set(floats))[:100],"text_chars":len(t),"text_sha256":hashlib.sha256(t.encode()).hexdigest()})
        result.update(status="SOURCE_RETRIEVED_METADATA_ONLY",bytes=n,sha256=sha.hexdigest(),pages=len(reader.pages),candidate_pages=pages)
        # The supplement remains an external source, not a source of new labels absent exact per-signal joins.
    except Exception as e:result["reason"]=str(e)[:400]
    (OUT/"supplement_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("R33_SUPPLEMENT",json.dumps({k:v for k,v in result.items() if k!="candidate_pages"}),flush=True)
if __name__=="__main__":main()
