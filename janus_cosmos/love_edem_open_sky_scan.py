from __future__ import annotations

import argparse, json, math, re
from pathlib import Path
from urllib.parse import urljoin

import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

TARGETS = {
    "LOVE": (204.30267916666668, -36.78240527777778),
    "EDEM_ZAPORIZHZHIA_DIRECTION_CANDIDATE": (139.22409686590188, 30.26038779947318),
}


def sep_arcsec(ra1, dec1, ra2, dec2):
    a = SkyCoord(ra1*u.deg, dec1*u.deg)
    b = SkyCoord(float(ra2)*u.deg, float(dec2)*u.deg)
    return float(a.separation(b).arcsec)


def irsa_query(catalog, ra, dec, radius_arcsec, selcols):
    url = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
    params = {
        "catalog": catalog, "spatial": "cone", "radius": radius_arcsec,
        "radunits": "arcsec", "objstr": f"{ra} {dec}", "outfmt": 1,
        "outrows": 200, "selcols": ",".join(selcols),
    }
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    text = r.text
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("\\")]
    # IPAC tables are easiest and safest to parse with astropy.
    p = Path("/tmp/irsa.tbl"); p.write_text(text, encoding="utf-8")
    try:
        t = Table.read(p, format="ascii.ipac")
    except Exception:
        return {"status": "PARSE_ERROR", "raw_preview": text[:1000], "count": 0, "rows": []}
    rows=[]
    for row in t[:50]:
        d={}
        for c in t.colnames:
            v=row[c]
            try: v=v.item()
            except Exception: pass
            if hasattr(v, "mask") and v.mask: v=None
            d[c]=None if str(v)=="--" else v
        if "ra" in d and "dec" in d:
            try: d["separation_arcsec"] = sep_arcsec(ra, dec, d["ra"], d["dec"])
            except Exception: pass
        rows.append(d)
    rows.sort(key=lambda x: x.get("separation_arcsec", 1e99))
    return {"status":"OK", "count":len(t), "rows":rows}


def gaia_query(ra, dec, radius_arcsec):
    from astroquery.gaia import Gaia
    coord = SkyCoord(ra*u.deg, dec*u.deg)
    job = Gaia.cone_search_async(coord, radius_arcsec*u.arcsec)
    t = job.get_results()
    rows=[]
    wanted=["source_id","ra","dec","parallax","parallax_error","pmra","pmdec","phot_g_mean_mag","bp_rp","ruwe"]
    for row in t[:100]:
        d={}
        for c in wanted:
            if c not in t.colnames: continue
            v=row[c]
            try: v=v.item()
            except Exception: pass
            try:
                if math.isnan(float(v)): v=None
            except Exception: pass
            d[c]=v
        d["separation_arcsec"] = sep_arcsec(ra,dec,d["ra"],d["dec"])
        rows.append(d)
    rows.sort(key=lambda x:x["separation_arcsec"])
    return {"status":"OK", "count":len(t), "rows":rows}


def simbad_query(ra, dec, radius_arcsec):
    from astroquery.simbad import Simbad
    s = Simbad()
    try:
        s.add_votable_fields("otype", "ids", "flux(V)")
    except Exception:
        pass
    t = s.query_region(SkyCoord(ra*u.deg, dec*u.deg), radius=radius_arcsec*u.arcsec)
    if t is None:
        return {"status":"OK", "count":0, "rows":[]}
    rows=[]
    for row in t[:100]:
        d={}
        for c in t.colnames:
            v=row[c]
            if isinstance(v, bytes): v=v.decode("utf-8", "replace")
            try: v=v.item()
            except Exception: pass
            if str(v)=="--": v=None
            d[c]=v
        # modern SIMBAD returns ra/dec in degrees under ra/dec; legacy may differ.
        raval=d.get("ra"); decval=d.get("dec")
        if isinstance(raval,(int,float)) and isinstance(decval,(int,float)):
            d["separation_arcsec"] = sep_arcsec(ra,dec,raval,decval)
        rows.append(d)
    rows.sort(key=lambda x:x.get("separation_arcsec",1e99))
    return {"status":"OK", "count":len(t), "rows":rows}


def panstarrs_query(ra, dec, radius_arcsec):
    if dec < -30:
        return {"status":"OUTSIDE_PS1_3PI_DEC_COVERAGE", "count":0, "rows":[]}
    from astroquery.mast import Catalogs
    t = Catalogs.query_region(f"{ra} {dec}", radius=radius_arcsec*u.arcsec, catalog="Panstarrs", data_release="dr2")
    rows=[]
    cols=[c for c in ["objID","raMean","decMean","nDetections","gMeanPSFMag","rMeanPSFMag","iMeanPSFMag","zMeanPSFMag","yMeanPSFMag"] if c in t.colnames]
    for row in t[:100]:
        d={}
        for c in cols:
            v=row[c]
            try: v=v.item()
            except Exception: pass
            try:
                if math.isnan(float(v)): v=None
            except Exception: pass
            d[c]=v
        if d.get("raMean") is not None and d.get("decMean") is not None:
            d["separation_arcsec"] = sep_arcsec(ra,dec,d["raMean"],d["decMean"])
        rows.append(d)
    rows.sort(key=lambda x:x.get("separation_arcsec",1e99))
    return {"status":"OK", "count":len(t), "rows":rows}


def skyview_query(label, ra, dec, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    surveys=["DSS2 Red","2MASS-J","WISE 3.4"]
    result={}
    base="https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"
    for survey in surveys:
        params={"Position":f"{ra},{dec}","coordinates":"J2000","pixels":700,
                "projection":"Tan","scaling":"Sqrt","survey":survey,"Size":"0.2"}
        try:
            r=requests.get(base,params=params,timeout=120)
            r.raise_for_status()
            html=r.text
            safe=re.sub(r"[^A-Za-z0-9]+","_",survey).strip("_").lower()
            (outdir/f"{safe}.html").write_text(html,encoding="utf-8")
            hrefs=re.findall(r'href=["\']?([^"\' >]+)',html,re.I)
            srcs=re.findall(r'src=["\']?([^"\' >]+)',html,re.I)
            links=[]
            for h in hrefs+srcs:
                full=urljoin(r.url,h)
                if any(x in full.lower() for x in (".fits",".jpg",".jpeg",".png")) and full not in links:
                    links.append(full)
            downloaded=[]
            for i,link in enumerate(links[:8]):
                try:
                    rr=requests.get(link,timeout=120); rr.raise_for_status()
                    ctype=rr.headers.get("content-type","").lower()
                    ext=".bin"
                    for e in (".fits",".jpg",".jpeg",".png"):
                        if e in link.lower(): ext=e; break
                    if "jpeg" in ctype: ext=".jpg"
                    elif "png" in ctype: ext=".png"
                    elif "fits" in ctype: ext=".fits"
                    fn=f"{safe}_{i}{ext}"
                    (outdir/fn).write_bytes(rr.content)
                    downloaded.append({"file":fn,"url":link,"bytes":len(rr.content),"content_type":ctype})
                except Exception as exc:
                    downloaded.append({"url":link,"error":f"{type(exc).__name__}: {exc}"})
            result[survey]={"status":"OK","http_status":r.status_code,"links_found":links,"downloaded":downloaded}
        except Exception as exc:
            result[survey]={"status":"ERROR","error":f"{type(exc).__name__}: {exc}"}
    return result


def run(output_dir: Path, radius_arcsec: float):
    output_dir.mkdir(parents=True, exist_ok=True)
    receipt={"schema":"janus.cosmos.open_sky_scan.v1","radius_arcsec":radius_arcsec,"targets":{}}
    for label,(ra,dec) in TARGETS.items():
        tdir=output_dir/label.lower(); tdir.mkdir(exist_ok=True)
        entry={"ra_deg_icrs":ra,"dec_deg_icrs":dec}
        for name,func in [
            ("gaia_dr3", lambda:gaia_query(ra,dec,radius_arcsec)),
            ("simbad", lambda:simbad_query(ra,dec,radius_arcsec)),
            ("2mass_psc", lambda:irsa_query("fp_psc",ra,dec,radius_arcsec,["ra","dec","designation","j_m","h_m","k_m"])),
            ("allwise", lambda:irsa_query("allwise_p3as_psd",ra,dec,radius_arcsec,["ra","dec","designation","w1mpro","w2mpro","w3mpro","w4mpro"])),
            ("panstarrs_dr2", lambda:panstarrs_query(ra,dec,radius_arcsec)),
        ]:
            try: entry[name]=func()
            except Exception as exc: entry[name]={"status":"ERROR","error":f"{type(exc).__name__}: {exc}","count":0,"rows":[]}
        entry["skyview"] = skyview_query(label,ra,dec,tdir/"skyview")
        # compact summary
        entry["summary"]={k:{"status":v.get("status"),"count":v.get("count"),"nearest":(v.get("rows") or [None])[0]}
                          for k,v in entry.items() if isinstance(v,dict) and "count" in v}
        receipt["targets"][label]=entry
    p=output_dir/"open-sky-scan.json"; p.write_text(json.dumps(receipt,indent=2,default=str)+"\n",encoding="utf-8")
    print(json.dumps({k:v["summary"] for k,v in receipt["targets"].items()},indent=2,default=str))

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--output-dir",required=True); ap.add_argument("--radius-arcsec",type=float,default=120.0)
    a=ap.parse_args(); run(Path(a.output_dir),a.radius_arcsec)
