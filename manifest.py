import csv
from collections import Counter

from models import InputSeq, Partner

MANIFEST_COLUMNS = [
    "input_gene", "input_accession", "input_isoform", "input_reviewed",
    "partner_id", "partner_gene", "partner_kind", "sources",
    "string_score", "string_experiments", "string_textmining",
    "intact_mi", "uniprot_curated", "encori_evidence", "tier", "json_path",
]


def manifest_row(symbol, inp: InputSeq, partner: Partner, json_path):
    return {
        "input_gene": symbol,
        "input_accession": inp.accession,
        "input_isoform": inp.isoform_id,
        "input_reviewed": str(inp.reviewed),
        "partner_id": partner.partner_id,
        "partner_gene": partner.gene,
        "partner_kind": partner.kind,
        "sources": ";".join(partner.sources),
        "string_score": "" if partner.string_score is None else partner.string_score,
        "string_experiments": "" if partner.string_experiments is None else partner.string_experiments,
        "string_textmining": "" if partner.string_textmining is None else partner.string_textmining,
        "intact_mi": "" if partner.intact_mi is None else partner.intact_mi,
        "uniprot_curated": str(partner.uniprot_curated),
        "encori_evidence": partner.encori_evidence or "",
        "tier": partner.tier,
        "json_path": json_path,
    }


def write_manifest(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_readme(path, symbol, partners, date):
    by_kind = Counter(p.kind for p in partners)
    by_tier = Counter(p.tier for p in partners)
    lines = [
        f"AF3 partner inputs for {symbol}",
        f"Generated: {date}",
        "",
        "Sources queried: UniProt, STRING, IntAct, ENCORI (+ curated --rna-tsv if supplied).",
        "All external references use their latest release.",
        "",
        f"Total partners: {len(partners)}",
        f"  By kind: " + ", ".join(f"{k}={v}" for k, v in sorted(by_kind.items())),
        f"  By tier: " + ", ".join(f"{k}={v}" for k, v in sorted(by_tier.items())),
        "",
        "Each JSON pairs one input-protein isoform with one partner (AF3-local format).",
        "Files are organized by partner kind and confidence tier; see manifest.tsv for scores.",
    ]
    if by_kind.get("rna", 0) == 0:
        lines.append("")
        lines.append("Note: no RNA partners — gene is not a catalogued RBP in ENCORI and "
                     "no --rna-tsv was supplied (or supplied rows were empty).")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
