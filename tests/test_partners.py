import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from af3partners import partners

STRING_TSV = (
    "stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\t"
    "nscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore\n"
    "9606.ENSP1\t9606.ENSP2\tRPS24\tRPS15\t9606\t0.999\t0\t0\t0\t0.991\t0.994\t0.9\t0.727\n"
    "9606.ENSP1\t9606.ENSP3\tRPS24\tTP53\t9606\t0.45\t0\t0\t0\t0\t0\t0\t0.45\n"
)

INTACT_JSON = """
{"content": [
  {"moleculeA": "RPS24", "moleculeB": "RPS15", "uniqueIdA": "P62847",
   "uniqueIdB": "P62841", "intactMiscore": 0.7, "taxIdA": 9606, "taxIdB": 9606,
   "typeA": "protein", "typeB": "protein"},
  {"moleculeA": "RPS24", "moleculeB": "YEASTX", "uniqueIdA": "P62847",
   "uniqueIdB": "Q00001", "intactMiscore": 0.9, "taxIdA": 9606, "taxIdB": 559292,
   "typeA": "protein", "typeB": "protein"},
  {"moleculeA": "OTHER", "moleculeB": "NOISE", "uniqueIdA": "Z9", "uniqueIdB": "Z8",
   "intactMiscore": 0.9, "taxIdA": 9606, "taxIdB": 9606, "typeA": "protein", "typeB": "protein"}
]}
"""


class TestParseString(unittest.TestCase):
    def setUp(self):
        self.d = partners.parse_string(STRING_TSV)

    def test_keys_uppercase_gene(self):
        self.assertIn("RPS15", self.d)
        self.assertIn("TP53", self.d)

    def test_scores(self):
        p = self.d["RPS15"]
        self.assertEqual(p.string_score, 0.999)
        self.assertEqual(p.string_experiments, 0.994)
        self.assertEqual(p.string_textmining, 0.727)
        self.assertEqual(p.sources, ["string"])

    def test_empty(self):
        self.assertEqual(partners.parse_string(""), {})


class TestParseIntact(unittest.TestCase):
    def setUp(self):
        self.d = partners.parse_intact(INTACT_JSON, {"P62847"})

    def test_keeps_human_partner_of_input(self):
        self.assertIn("RPS15", self.d)
        self.assertEqual(self.d["RPS15"].intact_mi, 0.7)
        self.assertEqual(self.d["RPS15"].partner_id, "P62841")
        self.assertEqual(self.d["RPS15"].sources, ["intact"])

    def test_drops_non_human_partner(self):
        self.assertNotIn("YEASTX", self.d)

    def test_drops_interactions_not_involving_input(self):
        self.assertNotIn("NOISE", self.d)


class TestMerge(unittest.TestCase):
    def test_union_merges_scores_and_sources(self):
        s = partners.parse_string(STRING_TSV)
        i = partners.parse_intact(INTACT_JSON, {"P62847"})
        m = partners.merge_partners(s, i)
        rps15 = m["RPS15"]
        self.assertEqual(rps15.string_score, 0.999)
        self.assertEqual(rps15.intact_mi, 0.7)
        self.assertEqual(sorted(rps15.sources), ["intact", "string"])
        self.assertEqual(rps15.partner_id, "P62841")  # accession carried from intact
