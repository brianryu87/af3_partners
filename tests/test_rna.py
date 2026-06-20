import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from af3partners import rna

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


class TestFetchRnacentralSequence(unittest.TestCase):
    def test_strips_taxid_and_normalizes_to_rna(self):
        captured = {}
        def fake(url):
            captured["url"] = url
            return '{"sequence": "CAAAGTGCTG", "length": 10}'
        seq = rna.fetch_rnacentral_sequence("URS0000149452_9606", http=fake)
        self.assertIn("/URS0000149452.json", captured["url"])  # taxid suffix stripped
        self.assertEqual(seq, "CAAAGUGCUG")                    # T -> U

    def test_failure_returns_none(self):
        def boom(url):
            raise OSError("rnacentral down")
        self.assertIsNone(rna.fetch_rnacentral_sequence("URS0000149452_9606", http=boom))

    def test_missing_sequence_field_returns_none(self):
        self.assertIsNone(rna.fetch_rnacentral_sequence("URS1", http=lambda u: '{"length": 0}'))


INTACT_RNA_JSON = """
{"content": [
  {"moleculeA": "ELAVL1", "moleculeB": "hsa-mir-93", "uniqueIdA": "Q15717",
   "uniqueIdB": "URS0000149452_9606", "intactMiscore": 0.56, "taxIdA": 9606,
   "taxIdB": 9606, "typeA": "protein", "typeB": "mirna"},
  {"moleculeA": "ELAVL1", "moleculeB": "lncrna_h19", "uniqueIdA": "Q15717",
   "uniqueIdB": "URS0000767B73_10090", "intactMiscore": 0.9, "taxIdA": 9606,
   "taxIdB": 10090, "typeA": "protein", "typeB": "lncrna"},
  {"moleculeA": "ELAVL1", "moleculeB": "mrna_x", "uniqueIdA": "Q15717",
   "uniqueIdB": "ENST00000407627", "intactMiscore": 0.5, "taxIdA": 9606,
   "taxIdB": 9606, "typeA": "protein", "typeB": "mrna"},
  {"moleculeA": "ELAVL1", "moleculeB": "PABPC1", "uniqueIdA": "Q15717",
   "uniqueIdB": "P11940", "intactMiscore": 0.8, "taxIdA": 9606, "taxIdB": 9606,
   "typeA": "protein", "typeB": "protein"}
]}
"""


class TestParseIntactRna(unittest.TestCase):
    def setUp(self):
        self.d = rna.parse_intact_rna(INTACT_RNA_JSON, {"Q15717"})

    def test_keeps_human_urs_rna(self):
        self.assertIn("URS0000149452", self.d)
        p = self.d["URS0000149452"]
        self.assertEqual(p.kind, "rna")
        self.assertEqual(p.partner_id, "URS0000149452")  # taxid stripped
        self.assertEqual(p.gene, "hsa-mir-93")
        self.assertEqual(p.intact_mi, 0.56)
        self.assertEqual(p.sources, ["intact"])
        self.assertIsNone(p.sequence)

    def test_drops_non_human_urs(self):
        self.assertNotIn("URS0000767B73", self.d)

    def test_drops_mrna_enst_and_protein(self):
        self.assertTrue(all(k.startswith("URS") for k in self.d))
        self.assertNotIn("P11940", self.d)

    def test_only_one_kept(self):
        self.assertEqual(len(self.d), 1)


class TestDiscoverRnaWithIntact(unittest.TestCase):
    def test_intact_rna_resolved_and_appended(self):
        def fake(url):
            if "rnacentral" in url:
                return '{"sequence": "CAAAGUGCUGUUCGUGCAGGUAG", "length": 23}'
            if "encori" in url:
                return "#cite\nRBP\tgeneName\tgeneType\ttotalClipExpNum\n"
            raise AssertionError(f"unexpected url: {url}")
        parts = rna.discover_rna_partners(
            "ELAVL1", None, http=fake,
            intact_text=INTACT_RNA_JSON, input_accessions={"Q15717"})
        mirs = [p for p in parts if p.partner_id == "URS0000149452"]
        self.assertEqual(len(mirs), 1)
        self.assertEqual(mirs[0].sequence, "CAAAGUGCUGUUCGUGCAGGUAG")
        self.assertEqual(mirs[0].kind, "rna")
        self.assertEqual(mirs[0].sources, ["intact"])

    def test_no_intact_text_means_no_intact_rna(self):
        parts = rna.discover_rna_partners(
            "ELAVL1", None, http=lambda u: "#c\nRBP\tgeneName\tgeneType\ttotalClipExpNum\n")
        self.assertEqual(parts, [])
