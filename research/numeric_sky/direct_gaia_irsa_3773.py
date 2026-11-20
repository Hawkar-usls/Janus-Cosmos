#!/usr/bin/env python3
"""Direct ESA Gaia and NASA/IPAC AllWISE TAP crossmatch for JANUS singularity 3773.
VizieR unavailable: query independent archive services directly, fail closed.
"""
import csv,hashlib,io,json,math,os,pathlib,zipfile
import requests
RUN=44836632590
ARTIFACT=14679400570
DIGEST="0806e5559b068ff4eca699bd42872c9986d6b780f9aed400cbad00d86d74f34b"
FINGERPRINT="1c42dc53a9591f477cfba9678121c7ab331cc2afa2286b7d5187ef376fced080"
O=pathlib.Path("runtime/numeric_sky/singularity")
O.mkdir(parents=True,exist_ok=True)
H={"User-Agent":"JANUS-COSMOS-ZerkaloJanusa-TAP/1.0 (read-only)"}
GITHUB_H={"Authorization":"Bearer "+os.environ["GH_TOKEN"],"Accept":"application/vnd.github+json"}
def sha(b):return hashlib.sha256(b).hexdigest()
def get(url,headers=None,params=None,timeout=90):
    r=requests.get(url,headers=headers or H,params=params,timeout=timeout);r.raise_for_status();return r
def load_frozen():
    uri=f"https://api.github.com/repos/Hawkar-usls/Janus-Cosmos/actions/artifacts/{ARTIFACT}"
    m=get(uri,headers=GITHUB_H).json()
    assert m["workflow_run"]["id"]==RUN and m["digest"]=="sha256:"+DIGEST and not m["expired"]
    raw=get(uri+"/zip",headers=GITHUB_H,timeout=120).content
    assert sha(raw)==DIGEST
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        c=json.loads(z.read("3773_coordinate.json"))
        sim=json.loads(z.read("3773_catalog_checks.json"))
    assert c["fingerprint_sha256"]==sim["source_fingerprint"]==FINGERPRINT
    assert sim["counts"]["singularity"]["near_60arcsec"]==1
    assert sim["counts"]["random_control"]["available"]==256
    return c,sim
def distance(ra1,dec1,ra2,dec2):
    a,b,c,d=map(math.radians,(ra1,dec1,ra2,dec2))
    u=math.sin((b-d)/2)**2+math.cos(b)*math.cos(d)*math.sin((a-c)/2)**2
    return 2*math.degrees(math.asin(min(1,math.sqrt(max(0,u)))))*3600
def query(services,adql):
    errors=[]
    for service in services:
        for method in ("get","post"):
            try:
                args={"REQUEST":"doQuery","LANG":"ADQL","QUERY":adql,"FORMAT":"csv"}
                if method=="get":
                    r=requests.get(service,params=args,headers=H,timeout=(10,90))
                else:
                    r=requests.post(service,data=args,headers=H,timeout=(10,90))
                r.raise_for_status()
                s=r.text.strip()
                if not s or s.startswith(("<?xml","<html","<!DOCTYPE")):raise ValueError("Non-CSV/no-response body from TAP")
                rows=list(csv.DictReader([line for line in s.splitlines() if not line.startswith("#")]))
                if rows and (not "ra" in {x.lower() for x in rows[0]} or not "dec" in {x.lower() for x in rows[0]}):
                    raise ValueError("Unexpected CSV headers: "+str(list(rows[0]))[:200])
                return {"status":"OK","service":service,"method":method,"rows":rows,"errors":errors}
            except Exception as e:errors.append(service+" "+method+": "+str(e)[:190])
    return {"status":"UNAVAILABLE","service":None,"rows":[],"errors":errors}
def check(label,target,urls,adql,idfield,extra):
    result=query(urls,adql)
    if result["status"]!="OK":return {"catalog":label,**result}
    out=[]
    ra=target["ra_deg"];dec=target["dec_deg"]
    for a in result["rows"]:
        try:
            obj={k.lower():v for k,v in a.items()}
            r=float(obj["ra"]);d=float(obj["dec"])
            sep=distance(ra,dec,r,d)
            if sep>60.05:continue
            keep={k:obj.get(k.lower()) for k in extra}
            out.append({"catalog_id":obj.get(idfield.lower()),"ra_deg":r,"dec_deg":d,
               "separation_arcsec":sep,"catalog_data":keep})
        except Exception:continue
    out.sort(key=lambda o:o["separation_arcsec"])
    return {"catalog":label,"status":"OK","service":result["service"],"method":result["method"],
      "sources_within_60arcsec":len(out),
      "sources_within_2arcsec":sum(x["separation_arcsec"]<=2 for x in out),
      "sources_within_10arcsec":sum(x["separation_arcsec"]<=10 for x in out),
      "sources_within_30arcsec":sum(x["separation_arcsec"]<=30 for x in out),
      "nearest":out[:15],"truncated_possible":len(result["rows"])>=100,"source_rows":len(result["rows"])}
def main():
    source,sim=load_frozen()
    tgt=source["direction"]
    ra=tgt["ra_deg"];dec=tgt["dec_deg"];rad=60/3600
    gaia_adql=("SELECT TOP 100 source_id, ra, dec, phot_g_mean_mag, parallax, pmra, pmdec "
         "FROM gaiadr3.gaia_source "
         f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{rad}))")
    wise_adql=("SELECT TOP 100 designation, ra, dec, w1mpro, w2mpro, w3mpro, w4mpro "
         "FROM allwise_p3as_psd "
         f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{rad}))")
    gaia=check("Gaia DR3 ESA",tgt,
         ["https://gea.esac.esa.int/tap-server/tap/sync","https://gaia.ari.uni-heidelberg.de/tap/sync"],
         gaia_adql,"source_id",["phot_g_mean_mag","parallax","pmra","pmdec"])
    print("JANUS_3773_GAIA_DIRECT",json.dumps(gaia,ensure_ascii=False),flush=True)
    wise=check("AllWISE NASA/IPAC",tgt,
        ["https://irsa.ipac.caltech.edu/TAP/sync"],wise_adql,"designation",
        ["w1mpro","w2mpro","w3mpro","w4mpro"])
    print("JANUS_3773_WISE_DIRECT",json.dumps(wise,ensure_ascii=False),flush=True)
    out={"schema":"janus.cosmos.zerkalo_janusa.3773.direct_tap.v1","date":"2026-11-21",
       "source_run":RUN,"source_artifact":ARTIFACT,"source_sha256":DIGEST,
       "source_fingerprint":FINGERPRINT,"target":tgt,
       "simbad_preserved":{"known_radio_sources_within_60arcsec":sim["observations"][0]["catalog_entries_top40"],
                           "256_random_controls":sim["counts"]["random_control"]},
       "gaia":gaia,"allwise":wise,
       "scientific_conclusion":"No novel source can be claimed from catalog matching or pseudorandom hash coordinate, even if matches exist.",
       "production_change":False,"independent_new_observation":False}
    (O/"3773_direct_gaia_allwise.json").write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("JANUS_3773_DIRECT_TAP_SUMMARY",json.dumps({"coordinate":tgt,"gaia":gaia,"allwise":wise},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
