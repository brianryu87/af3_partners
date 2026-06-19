import re

from models import Partner

RIBOSOMAL_RE = re.compile(r"^(RPS|RPL|MRPS|MRPL)\d")
RIBOSOMAL_EXPLICIT = {"RPSA", "FAU", "RACK1", "UBA52", "RPLP0", "RPLP1", "RPLP2"}

STRING_HIGH = 0.7
STRING_MED = 0.4
INTACT_EXPERIMENTAL = 0.45


def classify_kind(gene):
    g = gene.upper()
    if RIBOSOMAL_RE.match(g) or g in RIBOSOMAL_EXPLICIT:
        return "ribosomal"
    return "nonribosomal"


def derive_tier(p: Partner):
    experimental = (
        (p.string_experiments or 0) > 0
        or (p.intact_mi or 0) >= INTACT_EXPERIMENTAL
        or p.uniprot_curated
        or bool(p.encori_evidence)
    )
    score = p.string_score or 0
    if experimental or score >= STRING_HIGH:
        return "high"
    if score >= STRING_MED:
        return "medium"
    return "low"
