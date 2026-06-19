import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import rna

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")

ENCORI_TEXT = (
    "#please cite ...\n"
    "RBP\tgeneID\tgeneName\tgeneType\tclusterNum\ttotalClipExpNum\ttotalClipSiteNum\n"
    "ELAVL1\tENSG1\tTSPAN6\tprotein_coding\t8\t12\t20\n"
    "ELAVL1\tENSG1\tTSPAN6\tprotein_coding\t8\t12\t21\n"
    "ELAVL1\tENSG2\tMALAT1\tlincRNA\t3\t30\t40\n"
)


class TestParseEncori(unittest.TestCase):
    def setUp(self):
        self.d = rna.parse_encori(ENCORI_TEXT)

    def test_dedups_by_gene(self):
        self.assertIn("TSPAN6", self.d)
        self.assertIn("MALAT1", self.d)
        self.assertEqual(len(self.d), 2)

    def test_kind_is_rna_and_evidence_recorded(self):
        p = self.d["MALAT1"]
        self.assertEqual(p.kind, "rna")
        self.assertEqual(p.partner_id, "MALAT1")
        self.assertIn("30", p.encori_evidence)
        self.assertEqual(p.sources, ["encori"])

    def test_empty_text(self):
        self.assertEqual(rna.parse_encori("#comment only\n"), {})


class TestLoadRnaTsv(unittest.TestCase):
    def setUp(self):
        self.parts = rna.load_rna_tsv(os.path.join(FIXTURES, "rna_example.tsv"))

    def test_loads_two(self):
        self.assertEqual(len(self.parts), 2)

    def test_fields(self):
        zfas1 = [p for p in self.parts if p.gene == "ZFAS1"][0]
        self.assertEqual(zfas1.kind, "rna")
        self.assertEqual(zfas1.partner_id, "ZFAS1")
        self.assertEqual(zfas1.sequence, "ACGUACGUACGUAAAUUUGGGCCC")
        self.assertEqual(zfas1.sources, ["curated"])

    def test_dna_normalized_to_rna(self):
        rrna = [p for p in self.parts if p.gene == "RNA18SN5"][0]
        self.assertNotIn("T", rrna.sequence)
        self.assertTrue(rrna.sequence.startswith("UACCUGG"))


class TestDiscoverRnaPartners(unittest.TestCase):
    def test_merges_encori_and_curated_curated_wins_sequence(self):
        parts = rna.discover_rna_partners(
            "ELAVL1",
            os.path.join(FIXTURES, "rna_example.tsv"),
            http=lambda url: ENCORI_TEXT,
        )
        genes = {p.gene for p in parts}
        self.assertIn("MALAT1", genes)   # from ENCORI (no sequence)
        self.assertIn("ZFAS1", genes)    # from curated (has sequence)
        malat1 = [p for p in parts if p.gene == "MALAT1"][0]
        self.assertIsNone(malat1.sequence)


class TestEncoriFailure(unittest.TestCase):
    def test_encori_failure_falls_back_to_curated(self):
        def boom(url):
            raise OSError("encori down")
        parts = rna.discover_rna_partners(
            "RPS24",
            os.path.join(FIXTURES, "rna_example.tsv"),
            http=boom,
        )
        genes = {p.gene for p in parts}
        self.assertIn("ZFAS1", genes)    # curated rows still loaded
        self.assertEqual(len(parts), 2)  # ENCORI error swallowed; only the 2 curated rows
