import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from af3partners import classify
from af3partners.models import Partner


class TestClassifyKind(unittest.TestCase):
    def test_rps_is_ribosomal(self):
        self.assertEqual(classify.classify_kind("RPS15"), "ribosomal")

    def test_mrpl_is_ribosomal(self):
        self.assertEqual(classify.classify_kind("MRPL12"), "ribosomal")

    def test_explicit_set_is_ribosomal(self):
        self.assertEqual(classify.classify_kind("RACK1"), "ribosomal")
        self.assertEqual(classify.classify_kind("FAU"), "ribosomal")

    def test_other_is_nonribosomal(self):
        self.assertEqual(classify.classify_kind("TP53"), "nonribosomal")

    def test_rps24_itself_is_ribosomal_but_case_insensitive(self):
        self.assertEqual(classify.classify_kind("rps24"), "ribosomal")


class TestDeriveTier(unittest.TestCase):
    def test_string_experiments_makes_high(self):
        p = Partner(gene="X", string_score=0.2, string_experiments=0.5)
        self.assertEqual(classify.derive_tier(p), "high")

    def test_intact_above_threshold_makes_high(self):
        p = Partner(gene="X", intact_mi=0.5)
        self.assertEqual(classify.derive_tier(p), "high")

    def test_high_combined_score_makes_high(self):
        p = Partner(gene="X", string_score=0.8)
        self.assertEqual(classify.derive_tier(p), "high")

    def test_mid_score_makes_medium(self):
        p = Partner(gene="X", string_score=0.5, string_textmining=0.5)
        self.assertEqual(classify.derive_tier(p), "medium")

    def test_low_score_textmining_only_makes_low(self):
        p = Partner(gene="X", string_score=0.2, string_textmining=0.2)
        self.assertEqual(classify.derive_tier(p), "low")

    def test_encori_evidence_makes_high(self):
        p = Partner(gene="X", encori_evidence="12 CLIP experiments")
        self.assertEqual(classify.derive_tier(p), "high")
