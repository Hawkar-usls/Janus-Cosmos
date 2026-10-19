"""Offline unit, synthetic guard, and real-output nonpromotion checks for V6-R1."""

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workspace"))

from janus_cosmos_fermi_spectral_v6r1 import (
    AuditError, BACKGROUND, band_shape, band_plus_pl_shape,
    exponential_time_kernel, integrate_photon_bins,
)
from janus_cosmos_fermi_lat_irf_v6r1 import (
    audit_physical_parameter_bounds, nearest_irf_parameter,
)
from janus_cosmos_fermi_source_lag_v6r1 import (
    exact_lag_profile_degeneracy, naive_eta_ci, simulate,
)
from janus_cosmos_simulation_audit import digest

BASE = ROOT / "data/cosmos"
RECEIPTS = (
    "COSMOS-ALI-V6R1-GBM-BAND-POWERLAW-DRIFT-2026-10-19.json",
    "COSMOS-ALI-V6R1-LAT-FULL-IRF-READINESS-2026-10-19.json",
    "COSMOS-ALI-V6R1-SOURCE-LAG-TOY-INJECTION-COVERAGE-2026-10-19.json",
)


class TestV6R1Band(unittest.TestCase):
    def test_band_continuity_and_positive_hard_component(self):
        alpha, beta, epeak = -1., -2.5, 400.
        break_keV = (alpha - beta) * epeak / (alpha + 2.)
        left = band_shape(break_keV * (1. - 1e-8), alpha, beta, epeak)
        right = band_shape(break_keV * (1. + 1e-8), alpha, beta, epeak)
        self.assertAlmostEqual(left / right, 1., places=6)
        self.assertAlmostEqual(band_plus_pl_shape(1000., alpha, beta, epeak, 2.) /
                               band_shape(1000., alpha, beta, epeak), 3.)
        self.assertEqual(band_plus_pl_shape(200., alpha, beta, epeak, 0.),
                         band_shape(200., alpha, beta, epeak))

    def test_reject_invalid_band_parameters(self):
        for a, b, e in ((-2.1, -2.5, 200.), (-1., -0.5, 400.), (-1., -3., 0.)):
            with self.assertRaises(AuditError):
                band_shape(100., a, b, e)
        with self.assertRaises(AuditError):
            band_plus_pl_shape(100., -1., -2.4, 400., -1.)

    def test_energy_bin_quadrature_constant_and_linear(self):
        self.assertAlmostEqual(integrate_photon_bins([(10., 20.)], lambda x: 1.)[0], 10.)
        self.assertAlmostEqual(integrate_photon_bins([(10., 20.)], lambda x: x)[0], 150.)

    def test_time_exposure_kernel_is_exact_for_zero_slope(self):
        self.assertAlmostEqual(exponential_time_kernel(5., 9., 10., 0.), 4.)
        slope = .013
        self.assertAlmostEqual(exponential_time_kernel(1., 9., 0., slope),
                               (math.exp(9 * slope) - math.exp(slope)) / slope, places=9)
        with self.assertRaises(AuditError):
            exponential_time_kernel(9., 5., 0., 0.)

    def test_background_train_holdout_time_windows_disjoint(self):
        for grb, windows in BACKGROUND.items():
            self.assertLessEqual(windows[0][1], windows[1][0], grb)
            self.assertLessEqual(windows[1][1], windows[2][0], grb)
            self.assertLess(windows[2][1], 0., grb)


class TestV6R1CALDB(unittest.TestCase):
    @staticmethod
    def example_grid():
        return {"ENERG_LO": [10. + 10 * i for i in range(23)],
                "ENERG_HI": [20. + 10 * i for i in range(23)],
                "CTHETA_LO": [i / 8 for i in range(8)],
                "CTHETA_HI": [(i + 1) / 8 for i in range(8)],
                "F": [0.65] * 184, "S1": [1.] * 184, "K1": [1.] * 184,
                "BIAS": [0.] * 184, "BIAS2": [0.] * 184,
                "S2": [2.] * 184, "K2": [1.] * 184,
                "PINDEX1": [2.] * 184, "PINDEX2": [2.] * 184}

    def test_edisp_mixture_probability_bounds_and_grid_index(self):
        grid = self.example_grid()
        audit_physical_parameter_bounds("edisp", grid)
        self.assertAlmostEqual(nearest_irf_parameter(grid, "F", 100., 40.), .65)
        grid["F"][0] = 1.2
        with self.assertRaises(AuditError):
            audit_physical_parameter_bounds("edisp", grid)

    def test_reject_invalid_psf_width_and_tail(self):
        grid = {key: [1.5] * 184 for key in ("NCORE", "NTAIL", "SCORE", "STAIL", "GCORE", "GTAIL")}
        audit_physical_parameter_bounds("psf", grid)
        grid["GTAIL"][3] = 1.0
        with self.assertRaises(AuditError):
            audit_physical_parameter_bounds("psf", grid)


class TestV6R1SourceLag(unittest.TestCase):
    def test_intrinsic_beta_exactly_absorbs_eta_for_any_eta(self):
        k, y, sigma = (1., 2., .5), (1., -.4, 2.), (.2, .1, .4)
        for eta in (-100., 0., 2., 100.):
            self.assertEqual(exact_lag_profile_degeneracy(k, y, sigma, eta)["profile_chisquare"], 0.)

    def test_naive_eta_interval_fails_adversarial_source_lag(self):
        k = (.8, 1.5, .3)
        beta = tuple(.8 * x for x in k)
        est, lo, hi = naive_eta_ci(beta, k, (.2, .2, .3))
        self.assertAlmostEqual(est, .8)
        self.assertFalse(lo <= 0. <= hi)

    def test_seed_reproduces_toy_coverage_and_zero_is_not_discovery(self):
        a = simulate(seed=3, n=120)
        b = simulate(seed=3, n=120)
        self.assertEqual(a, b)
        self.assertLess(a["scenarios"][2]["naive_95pct_interval_empirical_coverage"], .1)
        self.assertIsNone(a["observed_multiburst_LIV_p_value"])


class TestV6R1FrozenReceipts(unittest.TestCase):
    def test_all_three_reproducible_digests_and_nonpromotion(self):
        for filename in RECEIPTS:
            r = json.loads((BASE / filename).read_text())
            d = dict(r)
            d.pop("result_sha256")
            self.assertEqual(digest(d), r["result_sha256"])
            self.assertIsNone(r.get("observed_liv_p_value", r.get("observed_multiburst_LIV_p_value")))
            self.assertEqual(r["scientific_gates"].get("BHB_BOOM"), "NOT_ESTABLISHED")

    def test_actual_r1_gbm_spectral_segments_and_independent_holdout(self):
        r = json.loads((BASE / RECEIPTS[0]).read_text())
        self.assertEqual(len(r["gbm_bursts"]), 3)
        self.assertEqual(sum(len(b["spectral_segments"]) for b in r["gbm_bursts"]), 17)
        for b in r["gbm_bursts"]:
            self.assertEqual(len(b["detector_pretrigger_train_and_holdout_background"]), 3)
            for s in b["spectral_segments"]:
                self.assertEqual(s["three_detector_fit"]["num_predeclared_shape_templates"], 96)
                self.assertTrue(s["nai_train_bgo_holdout"]["bgo_predeclared_independent_detector_evaluation"])
                self.assertIsNone(s["three_detector_fit"]["model_selection_significance"])

    def test_real_caldb_audit_covers_nine_families_and_no_kernels_claimed(self):
        r = json.loads((BASE / RECEIPTS[1]).read_text())
        self.assertEqual(len(r["response_families"]), 9)
        self.assertEqual(sum(len(x["event_type_parameter_audits"]) for x in r["response_families"]), 30)
        self.assertEqual(r["scientific_gates"]["LAT_NORMALIZED_PSF_AND_EDISP"], "BLOCKED_NOT_EVALUATED")


if __name__ == "__main__":
    unittest.main()
