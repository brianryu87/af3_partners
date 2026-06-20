from dataclasses import dataclass, field
from typing import Optional


@dataclass
class InputSeq:
    accession: str        # base accession, e.g. "P62847"
    isoform_id: str       # canonical or isoform id, e.g. "P62847" or "P62847-2"
    sequence: str
    reviewed: bool


@dataclass
class Partner:
    gene: str
    partner_id: str = ""                       # UniProt accession (protein) or gene (rna)
    kind: str = ""                             # "ribosomal" | "nonribosomal" | "rna"
    sequence: Optional[str] = None
    sources: list = field(default_factory=list)
    string_score: Optional[float] = None
    string_experiments: Optional[float] = None
    string_textmining: Optional[float] = None
    intact_mi: Optional[float] = None
    uniprot_curated: bool = False
    encori_evidence: Optional[str] = None
    tier: str = ""
