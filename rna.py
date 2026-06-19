import csv
import json
import urllib.parse

import httpget
from af3 import to_rna
from models import Partner

ENCORI_ASSEMBLY = "hg38"
ENCORI_MAX = 25
RNACENTRAL = "https://rnacentral.org/api/v1/rna"


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


def fetch_rnacentral_sequence(urs, http=httpget.get):
    base = urs.split("_")[0]  # strip any _<taxid> suffix
    try:
        data = json.loads(http(f"{RNACENTRAL}/{base}.json"))
    except Exception:
        return None
    seq = data.get("sequence")
    if not seq:
        return None
    return to_rna(seq)


def parse_intact_rna(json_text, input_accessions):
    data = json.loads(json_text)
    out = {}
    for item in data.get("content", []):
        if item.get("taxIdA") != 9606 or item.get("taxIdB") != 9606:
            continue
        ida, idb = item.get("uniqueIdA"), item.get("uniqueIdB")
        if ida in input_accessions:
            other_id, other_name = idb, item.get("moleculeB")
        elif idb in input_accessions:
            other_id, other_name = ida, item.get("moleculeA")
        else:
            continue
        if not other_id or not other_id.startswith("URS"):
            continue
        urs = other_id.split("_")[0]
        mi = item.get("intactMiscore")
        existing = out.get(urs)
        if existing is None or (mi if mi is not None else 0) > (existing.intact_mi if existing.intact_mi is not None else 0):
            out[urs] = Partner(gene=other_name or urs, partner_id=urs, kind="rna",
                               sources=["intact"], intact_mi=mi)
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


def discover_rna_partners(symbol, rna_tsv, http=httpget.get, intact_text=None, input_accessions=None):
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
                by_gene[key].sequence = p.sequence
                if "curated" not in by_gene[key].sources:
                    by_gene[key].sources.append("curated")
            else:
                by_gene[key] = p
    result = list(by_gene.values())
    if intact_text and input_accessions:
        for urs, p in parse_intact_rna(intact_text, input_accessions).items():
            p.sequence = fetch_rnacentral_sequence(p.partner_id, http)
            result.append(p)
    return result
