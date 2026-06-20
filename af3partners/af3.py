import random

from .models import InputSeq, Partner

AF3_DIALECT = "alphafold3"
AF3_VERSION = 4  # latest AF3 input JSON schema version
AF3_MAX_TOKENS = 5000  # AlphaFold Server per-job total: 1 residue or 1 nucleotide = 1 token


def to_rna(seq):
    return seq.replace("T", "U").replace("t", "u")


def af3_protein_pair(name, inp: InputSeq, partner: Partner):
    return {
        "name": name,
        "modelSeeds": [random.randint(1, 999)],
        "sequences": [
            {"protein": {"id": "A", "sequence": inp.sequence, "description": inp.isoform_id}},
            {"protein": {"id": "B", "sequence": partner.sequence, "description": partner.partner_id}},
        ],
        "dialect": AF3_DIALECT,
        "version": AF3_VERSION,
    }


def af3_rna_pair(name, inp: InputSeq, partner: Partner):
    return {
        "name": name,
        "modelSeeds": [random.randint(1, 999)],
        "sequences": [
            {"protein": {"id": "A", "sequence": inp.sequence, "description": inp.isoform_id}},
            {"rna": {"id": "B", "sequence": to_rna(partner.sequence), "description": partner.gene}},
        ],
        "dialect": AF3_DIALECT,
        "version": AF3_VERSION,
    }
