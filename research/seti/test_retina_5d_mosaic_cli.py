import hashlib
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from retina_5d_mosaic_cli import from_topa_frame_manifest
try:
    from PIL import Image
except ImportError:
    Image = None


@unittest.skipUnless(Image is not None, "Pillow needed for optional TOPA QC tests")
class TestTOPAManifest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        frames=self.root/"frames"
        frames.mkdir()
        data=[]
        for i in range(3):
            pix=np.zeros((16,16),dtype=np.uint8)
            pix[:4,:4]=180+i
            path=frames/f"frame-{i:08d}.jpg"
            Image.fromarray(pix).save(path)
            data.append({"file":path.name,
                         "sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                         "pts_time_s":float(i)})
        self.payload={"schema":"janus.topa.retina.frame-manifest.v1",
                      "source_sha256":"a"*64,"frames":data}
        self.manifest=self.root/"frame_manifest.json"

    def run_manifest(self):
        self.manifest.write_text(json.dumps(self.payload))
        return from_topa_frame_manifest(self.manifest)

    def test_valid_hashed_frames(self):
        out=self.run_manifest()
        self.assertEqual(out["sampled_frames"],3)
        self.assertEqual(out["topa_frame_pts_time_s"],[0.,1.,2.])
        self.assertFalse(out["original_video_rescanned"])

    def test_tamper_is_rejected(self):
        self.run_manifest()
        (self.root/"frames"/self.payload["frames"][0]["file"]).write_bytes(b"corrupt")
        with self.assertRaises(ValueError):
            self.run_manifest()

    def test_path_escape_is_rejected(self):
        self.payload["frames"][0]["file"]="../steal.jpg"
        with self.assertRaises(ValueError):
            self.run_manifest()


if __name__ == "__main__":
    unittest.main()
