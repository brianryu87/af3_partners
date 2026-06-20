import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from af3partners import uniprot

FASTA = (
    ">sp|P62847|RS24_HUMAN 40S ribosomal protein S24 OS=Homo sapiens\n"
    "MNDTVTIRTRKFMTNRLLQRK\n"
    "QMVIDVLHPGKATVPK\n"
    ">sp|P62847-2|RS24_HUMAN Isoform 2 OS=Homo sapiens\n"
    "MNDTVTIRT\n"
    ">tr|A0A2R8Y849|A0A2R8Y849_HUMAN Ribosomal protein S24 OS=Homo sapiens\n"
    "MNDTVTIRT\n"
)

ACC_TSV = "Entry\tReviewed\nP62841\treviewed\nA0A123\tunreviewed\n"
INTERACTION_TSV = 'Entry\tInteracts with\nP62847\tP62841; Q9Y2B4\n'


class TestParseInputFasta(unittest.TestCase):
    def setUp(self):
        self.seqs = uniprot.parse_input_fasta(FASTA)

    def test_count(self):
        self.assertEqual(len(self.seqs), 3)

    def test_canonical(self):
        c = self.seqs[0]
        self.assertEqual(c.accession, "P62847")
        self.assertEqual(c.isoform_id, "P62847")
        self.assertTrue(c.reviewed)
        self.assertEqual(c.sequence, "MNDTVTIRTRKFMTNRLLQRKQMVIDVLHPGKATVPK")

    def test_isoform_id_and_base_accession(self):
        iso = self.seqs[1]
        self.assertEqual(iso.accession, "P62847")
        self.assertEqual(iso.isoform_id, "P62847-2")

    def test_unreviewed_flag(self):
        tr = self.seqs[2]
        self.assertEqual(tr.isoform_id, "A0A2R8Y849")
        self.assertFalse(tr.reviewed)


class TestResolveInputSequences(unittest.TestCase):
    def test_uses_http_and_parses(self):
        seqs = uniprot.resolve_input_sequences("RPS24", http=lambda url: FASTA)
        self.assertEqual(len(seqs), 3)

    def test_empty_raises(self):
        with self.assertRaises(SystemExit):
            uniprot.resolve_input_sequences("NOTAGENE", http=lambda url: "")


class TestResolveAccession(unittest.TestCase):
    def test_prefers_reviewed(self):
        acc = uniprot.resolve_accession("RPS15", http=lambda url: ACC_TSV)
        self.assertEqual(acc, "P62841")

    def test_none_when_empty(self):
        self.assertIsNone(uniprot.resolve_accession("X", http=lambda url: "Entry\tReviewed\n"))


class TestFetchSequence(unittest.TestCase):
    def test_concatenates_fasta_body(self):
        single = ">sp|P62841|RS15\nAAAA\nBBBB\n"
        self.assertEqual(uniprot.fetch_sequence("P62841", http=lambda url: single), "AAAABBBB")


class TestUniprotInteractions(unittest.TestCase):
    def test_parses_interactant_accessions(self):
        accs = uniprot.parse_uniprot_interactions(INTERACTION_TSV)
        self.assertEqual(accs, {"P62841", "Q9Y2B4"})

    def test_empty_interactions(self):
        self.assertEqual(uniprot.parse_uniprot_interactions("Entry\tInteracts with\nP62847\t\n"), set())
