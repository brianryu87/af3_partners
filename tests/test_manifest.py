import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import manifest
from models import InputSeq, Partner


class TestManifestRow(unittest.TestCase):
    def test_row_has_all_columns(self):
        inp = InputSeq("P62847", "P62847-2", "MND", False)
        p = Partner(gene="RPS15", partner_id="P62841", kind="ribosomal",
                    sources=["string", "intact"], string_score=0.99,
                    string_experiments=0.9, intact_mi=0.7, tier="high")
        row = manifest.manifest_row("RPS24", inp, p, "AF3_inputs/ribosomal_protein/high/P62847-2_P62841.json")
        for col in manifest.MANIFEST_COLUMNS:
            self.assertIn(col, row)
        self.assertEqual(row["input_gene"], "RPS24")
        self.assertEqual(row["input_isoform"], "P62847-2")
        self.assertEqual(row["partner_gene"], "RPS15")
        self.assertEqual(row["partner_kind"], "ribosomal")
        self.assertEqual(row["sources"], "string;intact")
        self.assertEqual(row["tier"], "high")
        self.assertEqual(row["input_reviewed"], "False")


class TestWriteManifest(unittest.TestCase):
    def test_roundtrip(self):
        inp = InputSeq("P62847", "P62847", "MND", True)
        p = Partner(gene="RPS15", partner_id="P62841", kind="ribosomal", tier="high")
        rows = [manifest.manifest_row("RPS24", inp, p, "x.json")]
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "manifest.tsv")
            manifest.write_manifest(path, rows)
            with open(path, newline="") as f:
                read = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(read[0]["partner_gene"], "RPS15")
        self.assertEqual(list(read[0].keys()), manifest.MANIFEST_COLUMNS)


class TestWriteReadme(unittest.TestCase):
    def test_mentions_symbol_and_counts(self):
        partners = [Partner(gene="RPS15", kind="ribosomal", tier="high"),
                    Partner(gene="ZFAS1", kind="rna", tier="high")]
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "README.txt")
            manifest.write_readme(path, "RPS24", partners, "2026-06-19")
            text = open(path).read()
        self.assertIn("RPS24", text)
        self.assertIn("ribosomal", text)
        self.assertIn("2026-06-19", text)


class TestReadmeNoRnaCaveat(unittest.TestCase):
    def test_caveat_present_when_no_rna(self):
        partners = [Partner(gene="RPS15", kind="ribosomal", tier="high")]
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "README.txt")
            manifest.write_readme(path, "RPS24", partners, "2026-06-19")
            text = open(path, encoding="utf-8").read()
        self.assertIn("no RNA partners", text)
