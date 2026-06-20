import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from af3partners import af3
from af3partners.models import InputSeq, Partner


class TestAf3(unittest.TestCase):
    def test_to_rna_converts_t_to_u(self):
        self.assertEqual(af3.to_rna("ACGTacgt"), "ACGUacgu")

    def test_protein_pair_shape(self):
        inp = InputSeq("P62847", "P62847", "MND", True)
        partner = Partner(gene="RPS15", partner_id="P62841", sequence="AAA")
        job = af3.af3_protein_pair("P62847_P62841", inp, partner)
        self.assertEqual(job["dialect"], "alphafold3")
        self.assertEqual(job["version"], 4)
        self.assertEqual(job["name"], "P62847_P62841")
        self.assertEqual(len(job["sequences"]), 2)
        self.assertEqual(job["sequences"][0]["protein"]["id"], "A")
        self.assertEqual(job["sequences"][0]["protein"]["sequence"], "MND")
        self.assertEqual(job["sequences"][1]["protein"]["id"], "B")
        self.assertEqual(job["sequences"][1]["protein"]["sequence"], "AAA")
        self.assertEqual(len(job["modelSeeds"]), 1)

    def test_rna_pair_uses_rna_key_and_uracil(self):
        inp = InputSeq("P62847", "P62847", "MND", True)
        partner = Partner(gene="ZFAS1", partner_id="ZFAS1", sequence="ACGT")
        job = af3.af3_rna_pair("P62847_ZFAS1", inp, partner)
        self.assertIn("rna", job["sequences"][1])
        self.assertEqual(job["sequences"][1]["rna"]["sequence"], "ACGU")
        self.assertEqual(job["sequences"][1]["rna"]["description"], "ZFAS1")

    def test_af3_max_tokens_constant(self):
        self.assertEqual(af3.AF3_MAX_TOKENS, 5000)


if __name__ == "__main__":
    unittest.main()
