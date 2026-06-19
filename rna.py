import csv
import urllib.parse

import httpget
from af3 import to_rna
from models import Partner

ENCORI_ASSEMBLY = "hg38"
ENCORI_MAX = 25
RNA_MAX_LEN = 1500


def fetch_encori(symbol, http=httpget.get):
    url = ("https://rnasysu.com/encori/api/RBPTarget/"
           f"?assembly={ENCORI_ASSEMBLY}&geneType=all&RBP={urllib.parse.quote(symbol)}"
           "&clipExpNum=1&pancancerNum=0&target=all&cellType=all")
    return http(url)


def parse_encori(text):
    lines = [l for l in text.splitlines() if l and not l.startswith("#")]
    if not lines:
        return {}
    idx = {name: i for i, name in enumerate(lines[0].split("\t"))}
    if "geneName" not in idx:
        return {}
    best = {}
    for line in lines[1:]:
        c = line.split("\t")
        gene = c[idx["geneName"]]
        try:
            evidence = int(c[idx["totalClipExpNum"]])
        except (ValueError, KeyError, IndexError):
            evidence = 0
        key = gene.upper()
        if key not in best or evidence > best[key][0]:
            best[key] = (evidence, gene)
    top = sorted(best.items(), key=lambda kv: kv[1][0], reverse=True)[:ENCORI_MAX]
    out = {}
    for key, (evidence, gene) in top:
        out[key] = Partner(gene=gene, partner_id=gene, kind="rna",
                           sources=["encori"],
                           encori_evidence=f"{evidence} CLIP experiments")
    return out


def load_rna_tsv(path):
    out = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = (row.get("gene") or "").strip()
            seq = (row.get("sequence") or "").strip()
            if not gene or not seq:
                continue
            out.append(Partner(gene=gene, partner_id=gene, kind="rna",
                               sequence=to_rna(seq), sources=["curated"]))
    return out


def discover_rna_partners(symbol, rna_tsv, http=httpget.get):
    by_gene = {}
    try:
        for key, p in parse_encori(fetch_encori(symbol, http)).items():
            by_gene[key] = p
    except Exception as e:
        print(f"  WARN: ENCORI failed: {e}")
    if rna_tsv:
        for p in load_rna_tsv(rna_tsv):
            key = p.gene.upper()
            if key in by_gene:
                # curated sequence enriches an ENCORI hit
                by_gene[key].sequence = p.sequence
                if "curated" not in by_gene[key].sources:
                    by_gene[key].sources.append("curated")
            else:
                by_gene[key] = p
    return list(by_gene.values())
