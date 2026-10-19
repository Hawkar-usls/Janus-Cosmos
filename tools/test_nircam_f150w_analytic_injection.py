#!/usr/bin/env python3
"""Standalone optical-profile invariants; the full replay tests MAST pixels."""
import math
import unittest

import numpy as np

from inject_nircam_f150w_empirical_psf_brackets import (
    AP_ARCSEC, EE50_ARCSEC, EE80_ARCSEC, core_halo_fit, moffat_fit,
    psf_ee, radial_density,
)


class TestF150WProfileModels(unittest.TestCase):
    def test_both_profiles_reproduce_official_i2d_ee_anchors(self):
        for model, parameters in [("moffat", moffat_fit()),
                                  ("two_gaussian", core_halo_fit())]:
            with self.subTest(profile=model):
                self.assertAlmostEqual(psf_ee(EE50_ARCSEC, model, parameters), .5, places=10)
                self.assertAlmostEqual(psf_ee(EE80_ARCSEC, model, parameters), .8, places=10)
                self.assertTrue(.7 < psf_ee(AP_ARCSEC, model, parameters) < .8)

    def test_density_integrates_to_expected_encircled_energy(self):
        # Riemann integral of I(r) 2*pi*r dr out to the official EE80 radius.
        radii=(np.arange(20000)+.5)*(EE80_ARCSEC/20000)
        for model,parameters in [("moffat",moffat_fit()),
                                 ("two_gaussian",core_halo_fit())]:
            with self.subTest(profile=model):
                integral=np.sum(radial_density(radii**2,model,parameters)
                                *2*math.pi*radii)*EE80_ARCSEC/len(radii)
                self.assertAlmostEqual(integral,.8,places=6)

    def test_extended_gaussian_source_has_lower_aperture_encircled_energy(self):
        core,halo,weight=core_halo_fit()
        energies=[]
        for intrinsic_sigma in [0,.03,.06,.10,.15,.20]:
            blurred=(math.hypot(core,intrinsic_sigma),
                     math.hypot(halo,intrinsic_sigma),weight)
            energies.append(psf_ee(AP_ARCSEC,"two_gaussian",blurred))
        self.assertTrue(all(a>b for a,b in zip(energies,energies[1:])))


if __name__ == "__main__":
    unittest.main()
