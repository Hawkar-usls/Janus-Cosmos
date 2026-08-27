#!/usr/bin/env python3
"""File-backed I/O shim for the Gen4 GMRT live pass.

netCDF4 does not accept BytesIO through xarray. Keep the scientific contract
and metrics unchanged; only materialize the HTTP payload to a temporary file
before parsing it.
"""
from __future__ import annotations

import tempfile

import numpy as np
import requests
import xarray as xr

import kusto_frozen_point_morphology_gen4 as core


def fetch_gmrt_mask(halfwidth_km: float = 46.0):
    west, east, south, north = core.bbox_for_halfwidth_km(core.TARGET_LAT, core.TARGET_LON, halfwidth_km)
    params = {
        "west": f"{west:.8f}",
        "east": f"{east:.8f}",
        "south": f"{south:.8f}",
        "north": f"{north:.8f}",
        "layer": "topo-mask",
        "format": "coards",
        "resolution": "max",
    }
    r = requests.get(core.GMRT_URL, params=params, timeout=180)
    r.raise_for_status()
    payload_sha = core.sha256_bytes(r.content)

    with tempfile.NamedTemporaryFile(suffix=".nc") as tmp:
        tmp.write(r.content)
        tmp.flush()
        with xr.open_dataset(tmp.name, engine="netcdf4") as ds:
            lat, lon, z, variable = core._pick_grid(ds)
            attrs = {k: str(v) for k, v in ds.attrs.items()}

    mlat, mlon = core.meters_per_degree(core.TARGET_LAT)
    dy = abs(float(np.median(np.diff(lat)))) * mlat if lat.size > 1 else float("inf")
    dx = abs(float(np.median(np.diff(lon)))) * mlon if lon.size > 1 else float("inf")
    resolution_m = max(dx, dy)
    return {
        "lat": lat,
        "lon": lon,
        "z": z,
        "resolution_m": resolution_m,
        "grid_dx_m": dx,
        "grid_dy_m": dy,
        "variable": variable,
        "request_url": r.url,
        "http_status": r.status_code,
        "payload_sha256": payload_sha,
        "attrs": attrs,
        "finite_fraction_entire_request": float(np.isfinite(z).mean()),
    }


core.fetch_gmrt_mask = fetch_gmrt_mask

if __name__ == "__main__":
    raise SystemExit(core.main())
