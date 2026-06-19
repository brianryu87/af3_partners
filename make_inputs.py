import argparse
import json
import os
import shutil
import sys
import tempfile
import zipfile
from datetime import date

import httpget
import uniprot
import partners as partners_mod
import rna as rna_mod
import manifest as manifest_mod
from af3 import af3_protein_pair, af3_rna_pair
from classify import classify_kind, derive_tier

KIND_DIR = {"ribosomal": "ribosomal_protein", "nonribosomal": "nonribosomal_protein", "rna": "rna"}


def _collect_partners(symbol, input_seqs, rna_tsv, http):
    input_accs = {s.accession for s in input_seqs}
    result = []

    proteins = partners_mod.discover_protein_partners(symbol, input_accs, http)
    seen_acc = {}
    for p in proteins.values():
        p.kind = classify_kind(p.gene)
        acc, seq = uniprot.fetch_protein_by_gene(p.gene, http)
        if not seq:
            print(f"  WARN: no sequence for {p.gene}, skipping")
            continue
        if acc in seen_acc:
            existing = seen_acc[acc]
            for s in p.sources:
                if s not in existing.sources:
                    existing.sources.append(s)
            print(f"  WARN: {p.gene} resolves to {acc}, already added; merged sources")
            continue
        p.partner_id = acc
        p.sequence = seq
        seen_acc[acc] = p
        result.append(p)

    for p in rna_mod.discover_rna_partners(symbol, rna_tsv, http):
        if not p.sequence:
            print(f"  WARN: no sequence for RNA {p.gene} (ENCORI-only), skipping")
            continue
        if len(p.sequence) > rna_mod.RNA_MAX_LEN:
            print(f"  WARN: RNA {p.gene} sequence {len(p.sequence)} nt exceeds RNA_MAX_LEN={rna_mod.RNA_MAX_LEN}, skipping")
            continue
        result.append(p)

    try:
        curated = uniprot.parse_uniprot_interactions(
            uniprot.fetch_uniprot_interactions(symbol, http))
    except Exception as e:
        print(f"  WARN: UniProt interactions failed: {e}")
        curated = set()
    for p in result:
        if p.kind == "rna":
            continue
        if p.partner_id in curated:
            p.uniprot_curated = True
            if "uniprot" not in p.sources:
                p.sources.append("uniprot")

    for p in result:
        p.tier = derive_tier(p)
    return result


def build(symbol, out_dir, rna_tsv=None, http=httpget.get):
    symbol = symbol.upper()
    input_seqs = uniprot.resolve_input_sequences(symbol, http)
    partners = _collect_partners(symbol, input_seqs, rna_tsv, http)

    work = tempfile.mkdtemp()
    try:
        root = os.path.join(work, symbol)
        rows = []
        for inp in input_seqs:
            for p in partners:
                name = f"{inp.isoform_id}_{p.partner_id}"
                rel = os.path.join("AF3_inputs", KIND_DIR[p.kind], p.tier, f"{name}.json")
                abspath = os.path.join(root, rel)
                os.makedirs(os.path.dirname(abspath), exist_ok=True)
                job = (af3_rna_pair if p.kind == "rna" else af3_protein_pair)(name, inp, p)
                with open(abspath, "w") as f:
                    json.dump([job], f, indent=2)
                rows.append(manifest_mod.manifest_row(symbol, inp, p, rel.replace(os.sep, "/")))

        os.makedirs(root, exist_ok=True)
        manifest_mod.write_manifest(os.path.join(root, "manifest.tsv"), rows)
        manifest_mod.write_readme(os.path.join(root, "README.txt"),
                                  symbol, partners, date.today().isoformat())

        os.makedirs(out_dir, exist_ok=True)
        zip_path = os.path.join(out_dir, f"{symbol}.zip")
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
            for dirpath, _, files in os.walk(root):
                for fn in files:
                    full = os.path.join(dirpath, fn)
                    z.write(full, os.path.relpath(full, work))
        print(f"Wrote {zip_path}: {len(rows)} JSONs, {len(partners)} partners")
        return zip_path
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate AF3 input JSONs for a gene's interaction partners.")
    parser.add_argument("symbol", help="Human gene symbol, e.g. RPS24")
    parser.add_argument("--out", default=".", help="Output directory for SYMBOL.zip")
    parser.add_argument("--rna-tsv", default=None, help="Optional curated RNA partners TSV (gene, sequence)")
    args = parser.parse_args(argv)
    build(args.symbol, args.out, args.rna_tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
