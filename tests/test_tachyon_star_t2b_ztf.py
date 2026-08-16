import unittest
from experiments.luci.run_tachyon_star_t2b_ztf import _product_url,_quality_ok,_classify_trial

class T(unittest.TestCase):
 def row(self):
  return {
   'ra_deg':'10.5','dec_deg':'20.25','field':'123','ccdid':'4','qid':'2','filtercode':'zr','exptime_s':'30',
   'a_filefracday':'20180411467847','a_imgtypecode':'o','a_infobits':'0',
   'b_filefracday':'20180411467999','b_imgtypecode':'o','b_infobits':'0',
   'c_filefracday':'20180411468111','c_imgtypecode':'o','c_infobits':'0'
  }
 def test_product_url(self):
  u=_product_url(self.row(),'a','sciimg.fits')
  self.assertIn('/2018/0411/467847/',u)
  self.assertIn('ztf_20180411467847_000123_zr_c04_o_q2_sciimg.fits',u)
  self.assertIn('center=10.5%2C20.25',u)
 def test_quality_threshold(self):
  r=self.row();self.assertTrue(_quality_ok(r))
  r['b_infobits']='33554432';self.assertFalse(_quality_ok(r))
 def test_isolated_b(self):
  a={'status':'QUALIFIED_ABSENCE'};b={'status':'SOURCE_PRESENT'};c={'status':'QUALIFIED_ABSENCE'}
  self.assertEqual(_classify_trial(a,b,c),'ISOLATED_B_L0')
 def test_null(self):
  a={'status':'QUALIFIED_ABSENCE'};b={'status':'QUALIFIED_ABSENCE'};c={'status':'QUALIFIED_ABSENCE'}
  self.assertEqual(_classify_trial(a,b,c),'NO_ISOLATED_B_EVENT')
 def test_unresolved(self):
  a={'status':'QUALIFIED_ABSENCE'};b={'status':'BLOCKED_LOCAL_SENSITIVITY'};c={'status':'QUALIFIED_ABSENCE'}
  self.assertEqual(_classify_trial(a,b,c),'UNRESOLVED_TRIAL')

if __name__=='__main__': unittest.main()
