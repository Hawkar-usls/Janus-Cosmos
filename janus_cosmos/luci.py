from __future__ import annotations

from pathlib import Path
from typing import Mapping

import numpy as np


class LuciProvenanceError(RuntimeError):
    """Raised when a FITS file cannot prove LUCI/LUCIFER provenance."""


_IMAGING_GRATING_TOKENS = {"", "NONE", "N/A", "NA", "OPEN", "OUT", "CLEAR", "MIRROR", "IMAGING", "IMAGE"}


def _norm(value: object) -> str:
    return "".join(ch for ch in str(value or "").upper().strip() if ch.isalnum())


def canonical_luci_instrument(value: object) -> str:
    token = _norm(value)
    if token in {"LUCI1", "LUCIFER1"}:
        return "LUCI1"
    if token in {"LUCI2", "LUCIFER2"}:
        return "LUCI2"
    if token in {"LUCI", "LUCIFER"}:
        return "LUCI_LEGACY_UNNUMBERED"
    raise LuciProvenanceError(f"FITS INSTRUME is not LUCI/LUCIFER: {value!r}")


def _first(headers: list[Mapping[str, object]], *keys: str) -> str:
    for header in headers:
        for key in keys:
            value = header.get(key)
            if value not in (None, ""):
                return str(value).strip()
    return ""


def inspect_luci_headers(primary: Mapping[str, object], image_header: Mapping[str, object] | None = None, *, require_imaging: bool = True) -> dict:
    headers = [primary]
    if image_header is not None:
        headers.append(image_header)

    instrument_raw = _first(headers, "INSTRUME", "INSTRUMENT")
    instrument = canonical_luci_instrument(instrument_raw)

    telescope = _first(headers, "TELESCOP", "TELESCOPE")
    telescope_norm = _norm(telescope)
    if telescope and not ("LBT" in telescope_norm or "LARGEBINOCULAR" in telescope_norm):
        raise LuciProvenanceError(f"TELESCOP is incompatible with LBT: {telescope!r}")

    grating = _first(headers, "GRATING", "GRATNAME", "GRISM", "GRAT")
    grating_norm = str(grating or "").upper().strip()
    inferred_mode = "imaging" if grating_norm in _IMAGING_GRATING_TOKENS else "spectroscopy"
    if require_imaging and inferred_mode != "imaging":
        raise LuciProvenanceError(
            f"LUCI frame appears spectroscopic (grating={grating!r}); Janus-Cosmos geometry gate accepts imaging only"
        )

    filt = _first(headers, "FILTER", "FILTER1", "FILTER2", "FILTNAME")
    return {
        "instrument": instrument,
        "instrument_raw": instrument_raw,
        "telescope": telescope or "UNRECORDED",
        "mode": inferred_mode,
        "grating": grating,
        "filter_header": filt,
        "object_header": _first(headers, "OBJECT", "OBJNAME", "TARGNAME"),
        "date_obs": _first(headers, "DATE-OBS", "DATEOBS"),
        "exptime": _first(headers, "EXPTIME", "DIT"),
        "provenance_gate": "LUCI_ONLY_PASS",
    }


def _collapse_image(data: np.ndarray) -> tuple[np.ndarray, str]:
    arr = np.asarray(data, dtype=np.float32)
    if arr.ndim == 2:
        return arr, "native_2d"
    if arr.ndim == 3 and min(arr.shape) >= 1:
        # LUCI can store repeated integrations as a cube. Collapse only the
        # integration axis; no spectral cube interpretation is attempted.
        axis = int(np.argmin(arr.shape)) if min(arr.shape) < 16 else 0
        return np.nanmedian(arr, axis=axis).astype(np.float32), f"median_axis_{axis}"
    raise LuciProvenanceError(f"unsupported LUCI science plane dimensions: {arr.shape}")


def read_luci_fits_image(path: Path, *, require_imaging: bool = True, expected_instrument: str | None = None) -> tuple[np.ndarray, dict]:
    from astropy.io import fits

    with fits.open(path, memmap=False, lazy_load_hdus=True) as hdul:
        primary_header = hdul[0].header
        candidates: list[tuple[int, int, int, str, np.ndarray, Mapping[str, object], str]] = []
        for idx, hdu in enumerate(hdul):
            data = getattr(hdu, "data", None)
            if data is None or np.ndim(data) not in (2, 3):
                continue
            try:
                image, collapse = _collapse_image(np.asarray(data))
            except LuciProvenanceError:
                continue
            if image.size < 1024:
                continue
            extname = str(hdu.header.get("EXTNAME", ""))
            candidates.append((image.size, 1 if extname.upper() == "SCI" else 0, idx, extname, image, hdu.header, collapse))

        if not candidates:
            raise LuciProvenanceError(f"no usable 2-D/3-D LUCI image plane in {path}")

        _, _, idx, extname, image, image_header, collapse = max(candidates, key=lambda x: (x[1], x[0]))
        provenance = inspect_luci_headers(primary_header, image_header, require_imaging=require_imaging)

    if expected_instrument:
        expected = canonical_luci_instrument(expected_instrument)
        actual = provenance["instrument"]
        legacy_wildcard = expected == "LUCI_LEGACY_UNNUMBERED" or actual == "LUCI_LEGACY_UNNUMBERED"
        if expected != actual and not legacy_wildcard:
            raise LuciProvenanceError(f"manifest instrument {expected!r} does not match FITS header {actual!r}")

    meta = {
        **provenance,
        "selected_hdu": int(idx),
        "selected_extname": extname,
        "native_shape": [int(x) for x in np.shape(image)],
        "cube_collapse": collapse,
        "nan_fraction": float(np.mean(~np.isfinite(image))),
    }
    return image, meta
