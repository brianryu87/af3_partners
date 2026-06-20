import re
import urllib.parse

from . import httpget
from .models import InputSeq

ORGANISM = 9606
BASE = "https://rest.uniprot.org/uniprotkb"


def _fasta_records(text):
    header, seq = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def parse_input_fasta(text):
    out = []
    for header, seq in _fasta_records(text):
        parts = header.split("|")
        db, acc = parts[0], parts[1]            # 'sp'/'tr', 'P62847' or 'P62847-2'
        out.append(InputSeq(
            accession=acc.split("-")[0],
            isoform_id=acc,
            sequence=seq,
            reviewed=(db == "sp"),
        ))
    return out


def resolve_input_sequences(symbol, http=httpget.get):
    q = f"gene:{symbol} AND organism_id:{ORGANISM}"
    url = (f"{BASE}/stream?query={urllib.parse.quote(q)}"
           "&format=fasta&includeIsoform=true")
    seqs = parse_input_fasta(http(url))
    if not seqs:
        raise SystemExit(f"No UniProt entries found for human gene '{symbol}'.")
    return seqs


def fetch_sequence(accession, http=httpget.get):
    text = http(f"{BASE}/{accession}.fasta")
    recs = list(_fasta_records(text))
    return recs[0][1] if recs else ""


def resolve_accession(gene, http=httpget.get):
    q = f"gene_exact:{gene} AND organism_id:{ORGANISM}"
    url = (f"{BASE}/search?query={urllib.parse.quote(q)}"
           "&format=tsv&fields=accession,reviewed&size=10")
    rows = [l.split("\t") for l in http(url).splitlines()[1:] if l]
    if not rows:
        return None
    for acc, reviewed in rows:
        if reviewed == "reviewed":
            return acc
    return rows[0][0]


def fetch_protein_by_gene(gene, http=httpget.get):
    acc = resolve_accession(gene, http)
    if not acc:
        return None, None
    return acc, fetch_sequence(acc, http)


def fetch_uniprot_interactions(symbol, http=httpget.get):
    q = f"gene_exact:{symbol} AND organism_id:{ORGANISM} AND reviewed:true"
    url = (f"{BASE}/search?query={urllib.parse.quote(q)}"
           "&format=tsv&fields=accession,cc_interaction&size=5")
    return http(url)


def parse_uniprot_interactions(tsv_text):
    accs = set()
    for line in tsv_text.splitlines()[1:]:
        if "\t" not in line:
            continue
        field = line.split("\t", 1)[1].strip()
        for tok in field.replace(",", ";").split(";"):
            tok = tok.strip()
            base = tok.split("-")[0]
            if re.fullmatch(r"[A-Z0-9]{6,10}", base):
                accs.add(base)
    return accs
