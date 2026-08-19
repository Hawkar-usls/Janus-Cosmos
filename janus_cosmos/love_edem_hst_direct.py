from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS


def _v(row, key, default=""):
    try:
        v = row[key]
    except Exception:
        return default
    if v is None:
        return default
    return v


def _score_product(row, preferred):
    subgroup = str(_v(row, "productSubGroupDescription", "")).upper()
    fn = str(_v(row, "productFilename", "")).lower()
    science = 1 if str(_v(row, "productType", "")).upper() == "SCIENCE" else 0
    fits_ok = 1 if fn.endswith((".fits", ".fits.gz")) else 0
    try:
        pref = len(preferred) - preferred.index(subgroup)
    except ValueError:
        pref = 0
    try:
        calib = int(_v(row, "calib_level", 0) or 0)
    except Exception:
        calib = 0
    return (science, fits_ok, pref, calib, fn)


def _find_containing_hdu(path: Path, coord: SkyCoord):
    with fits.open(path, memmap=False) as hdul:
        for idx, hdu in enumerate(hdul):
            data = getattr(hdu, "data", None)
            if data is None or np.ndim(data) < 2:
                continue
            arr = np.asarray(data)
            while arr.ndim > 2:
                arr = arr[0]
            if arr.ndim != 2:
                continue
            try:
                w = WCS(hdu.header).celestial
                x, y = w.world_to_pixel(coord)
            except Exception:
                continue
            ny, nx = arr.shape
            if np.isfinite(x) and np.isfinite(y) and 0 <= x < nx and 0 <= y < ny:
                return idx, arr.astype(float), w, float(x), float(y)
    return None


def _render_png(data, x, y, out_path: Path, title: str, p_lo: float, p_hi: float):
    import matplotlib.pyplot as plt
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return False
    lo, hi = np.nanpercentile(finite, [p_lo, p_hi])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = np.nanmin(finite), np.nanmax(finite)
    fig = plt.figure(figsize=(8, 8), dpi=160)
    ax = fig.add_subplot(111)
    ax.imshow(data, origin="lower", cmap="gray", vmin=lo, vmax=hi)
    ax.plot([x], [y], marker="+", markersize=18, markeredgewidth=2)
    ax.set_title(title)
    ax.set_xlabel("HST FITS quicklook — exact catalog coordinate marked +")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def run(prereg_path: Path, output_dir: Path) -> dict:
    from astroquery.mast import Observations

    cfg = json.loads(prereg_path.read_text(encoding="utf-8"))
    output_dir.mkdir(parents=True, exist_ok=True)
    qrad = float(cfg["mast"]["query_radius_arcmin"]) * u.arcmin
    preferred = list(cfg["mast"]["preferred_subgroups"])
    max_dl = int(cfg["mast"]["max_downloaded_products_per_target"])
    p_lo, p_hi = map(float, cfg["render"]["quicklook_percentiles"])

    results = []
    for target in cfg["targets"]:
        label = target["label"]
        coord = SkyCoord(float(target["ra_deg_icrs"]) * u.deg,
                         float(target["dec_deg_icrs"]) * u.deg,
                         frame="icrs")
        entry = {
            "label": label,
            "ra_deg_icrs": float(target["ra_deg_icrs"]),
            "dec_deg_icrs": float(target["dec_deg_icrs"]),
            "query_radius_arcmin": float(qrad.to_value(u.arcmin)),
            "observation_count": 0,
            "candidate_product_count": 0,
            "download_attempts": 0,
            "wcs_contained_products": [],
            "status": "NO_HST_COVERAGE"
        }
        try:
            obs = Observations.query_region(coord, radius=qrad)
            keep = []
            for row in obs:
                if str(_v(row, "obs_collection", "")).upper() != "HST":
                    continue
                if str(_v(row, "dataproduct_type", "")).lower() != "image":
                    continue
                if str(_v(row, "intentType", "science")).lower() != "science":
                    continue
                if str(_v(row, "dataRights", "PUBLIC")).upper() != "PUBLIC":
                    continue
                keep.append(row)
            entry["observation_count"] = len(keep)
            if not keep:
                results.append(entry)
                continue

            from astropy.table import Table
            products = Observations.get_product_list(Table(rows=keep, names=obs.colnames))
            cand = []
            seen = set()
            for row in products:
                fn = str(_v(row, "productFilename", ""))
                uri = str(_v(row, "dataURI", ""))
                if not fn.lower().endswith((".fits", ".fits.gz")):
                    continue
                if str(_v(row, "productType", "SCIENCE")).upper() != "SCIENCE":
                    continue
                key = (uri, fn)
                if key in seen:
                    continue
                seen.add(key)
                cand.append(row)
            cand.sort(key=lambda r: _score_product(r, preferred), reverse=True)
            entry["candidate_product_count"] = len(cand)

            tdir = output_dir / label.lower()
            tdir.mkdir(parents=True, exist_ok=True)
            for row in cand[:max_dl]:
                uri = str(_v(row, "dataURI", ""))
                fn = str(_v(row, "productFilename", "")) or uri.rsplit("/", 1)[-1]
                local = tdir / fn
                entry["download_attempts"] += 1
                try:
                    Observations.download_file(uri, local_path=str(local), cache=False)
                    hit = _find_containing_hdu(local, coord)
                    if hit is None:
                        local.unlink(missing_ok=True)
                        continue
                    hdu_idx, data, _w, x, y = hit
                    png = tdir / (fn.replace(".fits.gz", "").replace(".fits", "") + "__REAL_HST.png")
                    _render_png(data, x, y, png,
                                f"{label} — real HST archive image ({fn})",
                                p_lo, p_hi)
                    entry["wcs_contained_products"].append({
                        "filename": fn,
                        "dataURI": uri,
                        "hdu_index": hdu_idx,
                        "target_pixel": [x, y],
                        "quicklook_png": str(png.relative_to(output_dir)),
                        "fits_file": str(local.relative_to(output_dir)),
                        "instrument": str(_v(row, "instrument_name", "")),
                        "filter": str(_v(row, "filters", "")),
                        "subgroup": str(_v(row, "productSubGroupDescription", ""))
                    })
                except Exception as exc:
                    entry.setdefault("download_errors", []).append(f"{type(exc).__name__}: {exc}")
            if entry["wcs_contained_products"]:
                entry["status"] = "REAL_HST_FITS_AND_QUICKLOOK_READY"
            else:
                entry["status"] = "NO_EXACT_WCS_CONTAINED_HST_PRODUCT_IN_FROZEN_DOWNLOAD_SET"
        except Exception as exc:
            entry["status"] = "MAST_ERROR"
            entry["error"] = f"{type(exc).__name__}: {exc}"
        results.append(entry)

    receipt = {
        "schema": "janus.cosmos.love_edem.hst_direct.result.v1",
        "experiment_id": cfg["experiment_id"],
        "source": "MAST/STScI Hubble archive",
        "synthetic_generation_used": False,
        "targets": results,
        "claim_firewall": cfg["claim_firewall"]
    }
    (output_dir / "love-edem-hst-direct-result.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prereg", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()
    result = run(Path(args.prereg), Path(args.output_dir))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
