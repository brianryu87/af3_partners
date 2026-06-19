import os
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import make_inputs


@unittest.skipUnless(os.environ.get("AF3_LIVE") == "1",
                     "set AF3_LIVE=1 to run live network smoke test")
class TestLiveSmoke(unittest.TestCase):
    def test_rps24_real(self):
        with tempfile.TemporaryDirectory() as d:
            zip_path = make_inputs.build("RPS24", d)
            with zipfile.ZipFile(zip_path) as z:
                jsons = [n for n in z.namelist() if n.endswith(".json")]
                self.assertGreater(len(jsons), 0)
                # RPS24 should find ribosomal partners (e.g. RPS15) via STRING
                self.assertTrue(any("ribosomal_protein/" in n for n in jsons))
