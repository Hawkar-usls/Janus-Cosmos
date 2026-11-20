#!/usr/bin/env python3
"""Read-only Gaia DR3 and AllWISE recheck of the source-verified JANUS point 3773."""
import hashlib,io,json,os,pathlib,zipfile
import requests
import collapse_3773 as original

RUN=44836632590
ARTIFACT=14679400570
DIGEST="0806e5559b068ff4eca699bd42872c9986d6b780f9aed400cbad00d86d74f34b"
FINGERPRINT="1c42dc53a9591f477cfba9678121c7ab331cc2afa2286b7d5187ef376fced080"
O=pathlib.Path("runtime/numeric_sky/singularity")
O.mkdir(parents=True,exist_ok=True)
H={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
def sha(b):return hashlib.sha256(b).hexdigest()
def get(url):
    r=requests.get(url,headers=H,timeout=120);r.raise_for_status();return r
def main():
    uri=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{ARTIFACT}"
    source=get(uri).json()
    assert source["workflow_run"]["id"]==RUN and source["digest"]=="sha256:"+DIGEST and not source["expired"]
    raw=get(uri+"/zip").content;assert sha(raw)==DIGEST
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        names=z.namelist();assert len(names)==len(set(names))
        c=json.loads(z.read("3773_coordinate.json"))
        s=json.loads(z.read("3773_catalog_checks.json"))
        assert c["fingerprint_sha256"]==FINGERPRINT==s["source_fingerprint"]
        assert c["sources"]==3772 and c["direction"]["ordinal"]==3773
        assert s["counts"]["singularity"]["near_30arcsec"]==0
        assert s["counts"]["singularity"]["near_60arcsec"]==1
        assert s["counts"]["random_control"]["available"]==256
        assert s["counts"]["random_control"]["near_60arcsec"]==75
    coord=c["direction"]
    results=original.independent_multiband(coord["ra_deg"],coord["dec_deg"])
    receipt={"schema":"janus.cosmos.zerkalo_janusa.3773.multiband.v1","date":"2026-11-21",
        "source_run":RUN,"source_artifact":ARTIFACT,"source_zip_sha256":DIGEST,
        "source_fingerprint":FINGERPRINT,"coordinate":coord,
        "simbad_preserved":{"radius_arcsec":60,"target_found":True,"nearest_objects":
          s["observations"][0]["catalog_entries_top40"][:3],
          "random_control_matches":s["counts"]["random_control"]},
        "independent_catalog":results,"conclusion":"Catalog crossmatches are not evidence of a new sky object or technological source.",
        "production_change":False,"discovery_claim":False}
    (O/"3773_multiband_recheck.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("JANUS_3773_MULTIBAND",json.dumps(receipt,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
