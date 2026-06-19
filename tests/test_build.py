import json
import os
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import make_inputs

FASTA = (
    ">sp|P62847|RS24_HUMAN OS=Homo sapiens\nMNDTVTIR\n"
    ">sp|P62847-2|RS24_HUMAN Isoform 2\nMNDTVT\n"
)
STRING_TSV = (
    "stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\t"
    "nscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore\n"
    "9606.E1\t9606.E2\tRPS24\tRPS15\t9606\t0.999\t0\t0\t0\t0.9\t0.99\t0.9\t0.7\n"
)
INTACT_JSON = '{"content": []}'
ENCORI_TEXT = "#cite\nRBP\tgeneName\tgeneType\ttotalClipExpNum\n"  # no RNA hits
ACC_TSV = "Entry\tReviewed\nP62841\treviewed\n"
PARTNER_FASTA = ">sp|P62841|RS15\nMAAAA\n"
UNIPROT_INT_TSV = "Entry\tInteracts with\nP62847\t\n"


def fake_http(url):
    if "/stream?" in url:
        return FASTA
    if "string-db.org" in url:
        return STRING_TSV
    if "intact/ws" in url:
        return INTACT_JSON
    if "encori" in url:
        return ENCORI_TEXT
    if "cc_interaction" in url or "fields=accession,cc_interaction" in url:
        return UNIPROT_INT_TSV
    if "/search?" in url and "accession,reviewed" in url:
        return ACC_TSV
    if url.endswith("P62841.fasta"):
        return PARTNER_FASTA
    raise AssertionError(f"unexpected url: {url}")


class TestBuild(unittest.TestCase):
    def test_produces_zip_with_json_and_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            zip_path = make_inputs.build("RPS24", d, rna_tsv=None, http=fake_http)
            self.assertTrue(zip_path.endswith("RPS24.zip"))
            with zipfile.ZipFile(zip_path) as z:
                names = z.namelist()
                jsons = [n for n in names if n.endswith(".json")]
                # 2 isoforms x 1 ribosomal partner (RPS15) = 2 JSONs
                self.assertEqual(len(jsons), 2)
                self.assertTrue(any("ribosomal_protein/high/" in n for n in jsons))
                self.assertTrue(any(n.endswith("manifest.tsv") for n in names))
                self.assertTrue(any(n.endswith("README.txt") for n in names))
                sample = json.loads(z.read(jsons[0]))
                self.assertEqual(sample[0]["dialect"], "alphafold3")
                self.assertEqual(len(sample[0]["sequences"]), 2)


def fake_http_empty(url):
    if "/stream?" in url:
        return FASTA
    if "string-db.org" in url:
        return ("stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\t"
                "nscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore\n")
    if "intact/ws" in url:
        return '{"content": []}'
    if "encori" in url:
        return "#cite\nRBP\tgeneName\tgeneType\ttotalClipExpNum\n"
    if "cc_interaction" in url:
        return "Entry\tInteracts with\nP62847\t\n"
    raise AssertionError(f"unexpected url: {url}")


class TestBuildZeroPartners(unittest.TestCase):
    def test_zero_partners_still_valid_zip(self):
        with tempfile.TemporaryDirectory() as d:
            zip_path = make_inputs.build("RPS24", d, rna_tsv=None, http=fake_http_empty)
            with zipfile.ZipFile(zip_path) as z:
                names = z.namelist()
                self.assertEqual([n for n in names if n.endswith(".json")], [])
                self.assertTrue(any(n.endswith("manifest.tsv") for n in names))
                self.assertTrue(any(n.endswith("README.txt") for n in names))


INTACT_RNA_BUILD = ('{"content": [{"moleculeA": "RPS24", "moleculeB": "mir-x",'
                    ' "uniqueIdA": "P62847", "uniqueIdB": "URS0000ABCDEF_9606",'
                    ' "intactMiscore": 0.6, "taxIdA": 9606, "taxIdB": 9606,'
                    ' "typeA": "protein", "typeB": "mirna"}]}')


def fake_http_rna(url):
    if "/stream?" in url:
        return FASTA
    if "string-db.org" in url:
        return ("stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\t"
                "nscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore\n")
    if "intact/ws" in url:
        return INTACT_RNA_BUILD
    if "encori" in url:
        return "#cite\nRBP\tgeneName\tgeneType\ttotalClipExpNum\n"
    if "rnacentral" in url:
        return '{"sequence": "ACGUACGUACGU", "length": 12}'
    if "cc_interaction" in url:
        return "Entry\tInteracts with\nP62847\t\n"
    raise AssertionError(f"unexpected url: {url}")


class TestBuildWithIntactRna(unittest.TestCase):
    def test_intact_rna_produces_rna_json(self):
        with tempfile.TemporaryDirectory() as d:
            zip_path = make_inputs.build("RPS24", d, rna_tsv=None, http=fake_http_rna)
            with zipfile.ZipFile(zip_path) as z:
                jsons = [n for n in z.namelist() if n.endswith(".json")]
                rna_jsons = [n for n in jsons if "/rna/" in n]
                self.assertTrue(rna_jsons, f"expected an rna/ JSON, got {jsons}")
                self.assertTrue(any("URS0000ABCDEF" in n for n in rna_jsons))
                sample = json.loads(z.read(rna_jsons[0]))[0]
                self.assertIn("rna", sample["sequences"][1])
                self.assertEqual(sample["sequences"][1]["rna"]["sequence"], "ACGUACGUACGU")
