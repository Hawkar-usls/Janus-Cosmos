"""SIGPROC range-header parsing and no power access unit checks for SIT2E."""
import io
import struct
import unittest
import janus_sit2e_parkes_range_headers as p

def ss(s):
    b=s.encode("utf8")
    return struct.pack("<i",len(b))+b

class RangeHeaderContracts(unittest.TestCase):
    def test_july_real_source_format_is_filterbank(self):
        self.assertEqual(p.q.P["archive_sessions"][0]["date"],"2025-07-31")

    def test_sigproc_header_contract(self):
        x=ss("HEADER_START")
        for k,v in [("source_name","3IATLAS_S"),("nchans",8192),("nifs",2),
                    ("fch1",3000.),("foff",-0.002),("tstart",60887.6),("tsamp",1.5)]:
            x+=ss(k)
            if isinstance(v,str):x+=ss(v)
            elif isinstance(v,int):x+=struct.pack("<i",v)
            else:x+=struct.pack("<d",v)
        x+=ss("HEADER_END")
        r=p.parse_sigproc_header(x)
        self.assertEqual(r["nchans"],8192)
        self.assertAlmostEqual(r["channel_width_hz"],2000)
        self.assertTrue(r["no_radio_power_decoded"])

    def test_unknown_sigproc_parameter_fails_closed(self):
        with self.assertRaises(ValueError):
            p.parse_sigproc_header(ss("HEADER_START")+ss("nonexistent_param"))

    def test_truncated_fails_closed(self):
        with self.assertRaises(ValueError):
            p.parse_sigproc_header(ss("HEADER_START")+ss("nchans")+b"\x00")

    def test_role_requires_archive_suffix(self):
        self.assertEqual(p.q.role_hint("guppi_60930_19193_3IATLAS_R_0001.0002.h5"),"REFERENCE_NAME_ONLY")
        self.assertEqual(p.q.role_hint("guppi_60930_17044_3IATLAS_S_0001.0002.h5"),"SOURCE_NAME_ONLY")

if __name__=="__main__":unittest.main()
