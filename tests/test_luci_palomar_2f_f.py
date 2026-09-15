import unittest
import numpy as np
from experiments.luci.run_palomar_2f_f import gaussian_retained_fraction, edge_aware_local_injection, build_untested_rows
from experiments.luci.run_palomar_2f_e import local_coordinate_injection
class Test2FF(unittest.TestCase):
    def test_retained_flux_center_and_edge(self):
        self.assertGreater(gaussian_retained_fraction((128,128),64,64,4),0.999999); self.assertLess(gaussian_retained_fraction((128,128),64,126,6),0.995)
    def test_edge_recovery_repairs_only_margin(self):
        rng=np.random.default_rng(7); a=rng.normal(0,1,(128,128)); strict=local_coordinate_injection(a,119.0,64.0,4.0); self.assertEqual(strict.get('reason'),'COORDINATE_TOO_CLOSE_TO_EDGE'); edge=edge_aware_local_injection(a,119.0,64.0,4.0); self.assertGreaterEqual(edge['retained_psf_flux_fraction'],0.995); self.assertTrue(edge['passed'])
    def test_build_untested_excludes_representatives(self):
        import experiments.luci.run_palomar_2f_f as m
        oldp,oldf=m.EXPECTED_UNTESTED_PAIRS,m.EXPECTED_UNTESTED_FILES; m.EXPECTED_UNTESTED_PAIRS=m.EXPECTED_UNTESTED_FILES=1
        try:
            e=[{'src_id':'s','file_name':'a'},{'src_id':'s','file_name':'b'}]; r=[{'src_id':'s','file_name':'a'}]; self.assertEqual(build_untested_rows(e,r,('s',))[0]['file_name'],'b')
        finally: m.EXPECTED_UNTESTED_PAIRS,m.EXPECTED_UNTESTED_FILES=oldp,oldf
if __name__=='__main__': unittest.main()
