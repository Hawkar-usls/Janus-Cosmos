from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
import requests

TARGETS = {
    "TARGET_A": {
        "semantic_alias": "LOVE",
        "ra_deg": 204.30267916666668,
        "dec_deg": -36.78240527777778,
    },
    "TARGET_B": {
        "semantic_alias": "EDEM_ZAPORIZHZHIA_DIRECTION_CANDIDATE",
        "ra_deg": 139.22409686590188,
        "dec_deg": 30.26038779947318,
    },
}


def _clean(v):
    try:
        if getattr(v, "mask", False):
            return None
    except Exception:
        pass
    try:
        v = v.item()
    except Exception:
        pass
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    try:
        if math.isnan(float(v)):
            return None
    except Exception:
        pass
    if str(v) in {"--", "nan", "masked"}:
        return None
    return v


def _sep(ra1, dec1, ra2, dec2) -> float:
    a = SkyCoord(float(ra1) * u.deg, float(dec1) * u.deg, frame="icrs")
    b = SkyCoord(float(ra2) * u.deg, float(dec2) * u.deg, frame="icrs")
    return float(a.separation(b).arcsec)


def _rows(table, columns, center_ra, center_dec, ra_col="ra", dec_col="dec", limit=200):
    out = []
    for row in table[:limit]:
        d = {c: _clean(row[c]) for c in columns if c in table.colnames}
        if d.get(ra_col) is not None and d.get(dec_col) is not None:
            d["separation_from_center_arcsec"] = _sep(center_ra, center_dec, d[ra_col], d[dec_col])
        out.append(d)
    out.sort(key=lambda x: x.get("separation_from_center_arcsec", 1e99))
    return out


def gaia(ra, dec, radius_arcsec):
    from astroquery.gaia import Gaia
    coord = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    # astroquery 0.4.11 requires keyword arguments for cone_search_async.
    job = Gaia.cone_search_async(coordinate=coord, radius=radius_arcsec * u.arcsec)
    t = job.get_results()
    cols = [
        "source_id", "ra", "dec", "parallax", "parallax_error", "pmra", "pmdec",
        "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp", "ruwe",
    ]
    rows = _rows(t, cols, ra, dec)
    return {"status": "OK", "count": len(t), "rows": rows}


def irsa(catalog, ra, dec, radius_arcsec, columns):
    url = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
    params = {
        "catalog": catalog,
        "spatial": "cone",
        "radius": radius_arcsec,
        "radunits": "arcsec",
        "objstr": f"{ra} {dec}",
        "outfmt": 1,
        "outrows": 500,
        "selcols": ",".join(columns),
    }
    r = requests.get(url, params=params, timeout=90)
    r.raise_for_status()
    tmp = Path("/tmp/janus-center-irsa.tbl")
    tmp.write_text(r.text, encoding="utf-8")
    t = Table.read(tmp, format="ascii.ipac")
    rows = _rows(t, list(t.colnames), ra, dec)
    return {"status": "OK", "count": len(t), "rows": rows}


def simbad(ra, dec, radius_arcsec):
    from astroquery.simbad import Simbad
    s = Simbad()
    try:
        s.add_votable_fields("otype", "ids", "V")
    except Exception:
        try:
            s.add_votable_fields("otype", "ids")
        except Exception:
            pass
    t = s.query_region(SkyCoord(ra * u.deg, dec * u.deg, frame="icrs"), radius=radius_arcsec * u.arcsec)
    if t is None:
        return {"status": "OK", "count": 0, "rows": []}
    rows = []
    for row in t[:200]:
        d = {c: _clean(row[c]) for c in t.colnames}
        raval, decval = d.get("ra"), d.get("dec")
        if isinstance(raval, (int, float)) and isinstance(decval, (int, float)):
            d["separation_from_center_arcsec"] = _sep(ra, dec, raval, decval)
        rows.append(d)
    rows.sort(key=lambda x: x.get("separation_from_center_arcsec", 1e99))
    return {"status": "OK", "count": len(t), "rows": rows}


def sdss(ra, dec, radius_arcsec):
    from astroquery.sdss import SDSS
    coord = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    fields = ["objid", "ra", "dec", "type", "u", "g", "r", "i", "z"]
    t = SDSS.query_region(coord, radius=radius_arcsec * u.arcsec, photoobj_fields=fields, data_release=17)
    if t is None:
        return {"status": "OK", "count": 0, "rows": []}
    return {"status": "OK", "count": len(t), "rows": _rows(t, fields, ra, dec)}


def ps1_vizier(ra, dec, radius_arcsec):
    if dec < -30.0:
        return {"status": "OUTSIDE_PS1_3PI_DEC_COVERAGE", "count": 0, "rows": []}
    from astroquery.vizier import Vizier
    coord = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    v = Vizier(columns=["*", "+_r"], row_limit=200)
    tables = v.query_region(coord, radius=radius_arcsec * u.arcsec, catalog="II/349/ps1")
    if not tables:
        return {"status": "OK", "count": 0, "rows": []}
    t = tables[0]
    ra_col = "RAJ2000" if "RAJ2000" in t.colnames else "RAdeg"
    dec_col = "DEJ2000" if "DEJ2000" in t.colnames else "DEdeg"
    cols = [c for c in ["objID", ra_col, dec_col, "gmag", "rmag", "imag", "zmag", "ymag", "_r"] if c in t.colnames]
    rows = _rows(t, cols, ra, dec, ra_col=ra_col, dec_col=dec_col)
    return {"status": "OK", "count": len(t), "rows": rows}


def nearest_crossmatch(origin: dict | None, candidates: list[dict], ra_key="ra", dec_key="dec"):
    if not origin or origin.get("ra") is None or origin.get("dec") is None or not candidates:
        return None
    best = None
    for row in candidates:
        r = row.get(ra_key)
        d = row.get(dec_key)
        if r is None or d is None:
            continue
        sep = _sep(origin["ra"], origin["dec"], r, d)
        if best is None or sep < best[0]:
            best = (sep, row)
    if best is None:
        return None
    return {"separation_arcsec": float(best[0]), "row": best[1]}


def run(output: Path, radius_arcsec: float):
    result = {
        "schema": "janus.cosmos.love_edem.center_object_probe.v1",
        "radius_arcsec": radius_arcsec,
        "targets": {},
        "firewall": {
            "catalog_proximity_is_not_identity": True,
            "edem_identity_confirmed": False,
            "love_candidate_activated": False,
        },
    }
    for label, target in TARGETS.items():
        ra, dec = target["ra_deg"], target["dec_deg"]
        entry = {**target, "analysis_label": label}
        queries = {
            "gaia_dr3": lambda: gaia(ra, dec, radius_arcsec),
            "simbad": lambda: simbad(ra, dec, radius_arcsec),
            "allwise": lambda: irsa(
                "allwise_p3as_psd", ra, dec, radius_arcsec,
                ["ra", "dec", "designation", "w1mpro", "w2mpro", "w3mpro", "w4mpro",
                 "w1snr", "w2snr", "w3snr", "w4snr", "ph_qual", "cc_flags", "ext_flg", "var_flg", "nb", "na"],
            ),
            "2mass_psc": lambda: irsa(
                "fp_psc", ra, dec, radius_arcsec,
                ["ra", "dec", "designation", "j_m", "h_m", "k_m", "ph_qual", "cc_flg", "gal_contam"],
            ),
            "sdss_dr17": lambda: sdss(ra, dec, radius_arcsec),
            "panstarrs_vizier": lambda: ps1_vizier(ra, dec, radius_arcsec),
        }
        for name, fn in queries.items():
            try:
                entry[name] = fn()
            except Exception as exc:
                entry[name] = {"status": "ERROR", "error": f"{type(exc).__name__}: {exc}", "count": 0, "rows": []}

        wise0 = (entry["allwise"].get("rows") or [None])[0]
        entry["nearest_allwise_crossmatches"] = {
            "gaia_dr3": nearest_crossmatch(wise0, entry["gaia_dr3"].get("rows", [])),
            "sdss_dr17": nearest_crossmatch(wise0, entry["sdss_dr17"].get("rows", [])),
            "simbad": None,
        }
        # SIMBAD column names are service-version dependent; use numeric ra/dec when present.
        if wise0:
            candidates = []
            for row in entry["simbad"].get("rows", []):
                if isinstance(row.get("ra"), (int, float)) and isinstance(row.get("dec"), (int, float)):
                    candidates.append(row)
            entry["nearest_allwise_crossmatches"]["simbad"] = nearest_crossmatch(wise0, candidates)

        entry["summary"] = {
            name: {
                "status": value.get("status"),
                "count": value.get("count"),
                "nearest": (value.get("rows") or [None])[0],
            }
            for name, value in entry.items()
            if isinstance(value, dict) and "count" in value
        }
        result["targets"][label] = entry

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({
        label: {
            "summary": entry["summary"],
            "nearest_allwise_crossmatches": entry["nearest_allwise_crossmatches"],
        }
        for label, entry in result["targets"].items()
    }, indent=2, default=str))
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default="results/love_edem_center_object_probe/result.json")
    ap.add_argument("--radius-arcsec", type=float, default=30.0)
    args = ap.parse_args()
    run(Path(args.output), args.radius_arcsec)


if __name__ == "__main__":
    main()
