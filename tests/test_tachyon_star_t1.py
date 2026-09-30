import unittest
from datetime import datetime
from experiments.luci.run_tachyon_star_t1 import classify_trial

def q(status,date='2024-01-01T00:00:00',exptime='60'):
    return {'status':status,'image_meta':{'date_obs':date,'exptime':exptime}}

class T(unittest.TestCase):
    def test_isolated(self):
        x=classify_trial(q('QUALIFIED_NO_COUNTERPART_INHERITED_R1'),q('COUNTERPART_CANDIDATE'),q('QUALIFIED_NO_COUNTERPART_INHERITED_R1'))
        self.assertEqual(x['class'],'ISOLATED_B_EVENT_CANDIDATE')
    def test_no_event(self):
        x=classify_trial(q('QUALIFIED_NO_COUNTERPART_INHERITED_R1'),q('QUALIFIED_NO_COUNTERPART_INHERITED_R1'),q('QUALIFIED_NO_COUNTERPART_INHERITED_R1'))
        self.assertEqual(x['class'],'NO_ISOLATED_B_EVENT')
    def test_blocked(self):
        x=classify_trial(q('BLOCKED_LOCAL_SENSITIVITY'),q('COUNTERPART_CANDIDATE'),q('QUALIFIED_NO_COUNTERPART_INHERITED_R1'))
        self.assertEqual(x['class'],'UNRESOLVED_PRIMARY_TRIAL')
    def test_dual_overlap(self):
        a=q('QUALIFIED_NO_COUNTERPART_INHERITED_R1')
        b=q('COUNTERPART_CANDIDATE','2024-01-01T00:00:00','60')
        c=q('QUALIFIED_NO_COUNTERPART_INHERITED_R1')
        p=q('COUNTERPART_CANDIDATE','2024-01-01T00:00:05','60')
        x=classify_trial(a,b,c,p)
        self.assertEqual(x['paired_class'],'DUAL_INSTRUMENT_COINCIDENCE_CANDIDATE')
        self.assertTrue(x['paired_overlap']['adjudicating'])
    def test_paired_sensitive_veto(self):
        a=q('QUALIFIED_NO_COUNTERPART_INHERITED_R1')
        b=q('COUNTERPART_CANDIDATE','2024-01-01T00:00:00','60')
        c=q('QUALIFIED_NO_COUNTERPART_INHERITED_R1')
        p=q('QUALIFIED_NO_COUNTERPART_INHERITED_R1','2024-01-01T00:00:01','60')
        x=classify_trial(a,b,c,p)
        self.assertEqual(x['paired_class'],'PRIMARY_ONLY_WITH_PAIRED_SENSITIVE_ABSENCE')

if __name__=='__main__':unittest.main()
