import unittest
from experiments.luci.run_tachyon_star_t2a import eligible_ztf_triples,primary_ztf
class T(unittest.TestCase):
 def row(self,t,ff='1',sid='s',ex='30',field='1',ccd='2',qid='3',f='zr'):
  return {'src_id':sid,'obsdate':t,'filefracday':ff,'field':field,'ccdid':ccd,'qid':qid,'filtercode':f,'exptime':ex,'pid':'1','ipac_pub_date':'2020-01-01'}
 def test_eligible(self):
  r=[self.row('2024-01-01T00:00:00','1'),self.row('2024-01-01T00:04:00','2'),self.row('2024-01-01T00:11:00','3')]
  x=eligible_ztf_triples(r);self.assertEqual(len(x),1);self.assertEqual(x[0]['delta_pre_s'],240);self.assertEqual(x[0]['delta_post_s'],420)
 def test_gap_reject(self):
  r=[self.row('2024-01-01T00:00:00','1'),self.row('2024-01-01T01:00:00','2'),self.row('2024-01-01T01:05:00','3')]
  self.assertEqual(eligible_ztf_triples(r),[])
 def test_primary_per_source(self):
  a={'src_id':'a','timing_distance_s':10,'b_obsdate':'2024-01-02','b_filefracday':'2','field':'1','ccdid':'1','qid':'1'}
  b={**a,'timing_distance_s':2,'b_obsdate':'2024-01-01','b_filefracday':'1'}
  c={**a,'src_id':'b','timing_distance_s':4}
  x=primary_ztf(sorted([a,b,c],key=lambda r:(r['src_id'],r['timing_distance_s'],r['b_obsdate'])));self.assertEqual(len(x),2);self.assertEqual(x[0]['timing_distance_s'],2)
if __name__=='__main__':unittest.main()
