from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from janus_cosmos.luci import LuciProvenanceError, canonical_luci_instrument, read_luci_fits_image


def _write(path: Path, *, instrument="LUCI1", telescope="LBT", grating="MIRROR", cube=False):
    data = np.arange(3 * 64 * 64 if cube else 64 * 64, dtype=np.float32)
    data = data.reshape((3, 64, 64) if cube else (64, 64))
    hdu = fits.PrimaryHDU(data=data)
    hdu.header["INSTRUME"] = instrument
    hdu.header["TELESCOP"] = telescope
    hdu.header["FILTER"] = "K"
    hdu.header["GRATING"] = grating
    hdu.header["OBJECT"] = "SYNTH-LUCI"
    hdu.writeto(path)


def test_legacy_name_maps_to_luci():
    assert canonical_luci_instrument("LUCIFER1") == "LUCI1"
    assert canonical_luci_instrument("luci2") == "LUCI2"


def test_luci1_imaging_fits_passes(tmp_path: Path):
    p = tmp_path / "luci1.fits"
    _write(p)
    image, meta = read_luci_fits_image(p, expected_instrument="LUCI1")
    assert image.shape == (64, 64)
    assert meta["instrument"] == "LUCI1"
    assert meta["mode"] == "imaging"
    assert meta["provenance_gate"] == "LUCI_ONLY_PASS"


def test_non_luci_instrument_fails_closed(tmp_path: Path):
    p = tmp_path / "hst.fits"
    _write(p, instrument="ACS")
    with pytest.raises(LuciProvenanceError):
        read_luci_fits_image(p)


def test_non_lbt_telescope_fails_closed(tmp_path: Path):
    p = tmp_path / "wrong_scope.fits"
    _write(p, telescope="HST")
    with pytest.raises(LuciProvenanceError):
        read_luci_fits_image(p)


def test_spectroscopic_frame_rejected_by_geometry_runner(tmp_path: Path):
    p = tmp_path / "spec.fits"
    _write(p, grating="G210")
    with pytest.raises(LuciProvenanceError):
        read_luci_fits_image(p, require_imaging=True)


def test_luci_cube_collapses_to_2d(tmp_path: Path):
    p = tmp_path / "cube.fits"
    _write(p, cube=True)
    image, meta = read_luci_fits_image(p)
    assert image.shape == (64, 64)
    assert meta["cube_collapse"].startswith("median_axis_")


def test_manifest_instrument_mismatch_rejected(tmp_path: Path):
    p = tmp_path / "luci1.fits"
    _write(p, instrument="LUCI1")
    with pytest.raises(LuciProvenanceError):
        read_luci_fits_image(p, expected_instrument="LUCI2")
