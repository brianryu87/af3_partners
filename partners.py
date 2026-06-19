import json
import urllib.parse

import httpget
from models import Partner

ORGANISM = 9606


def fetch_string(symbol, http=httpget.get):
    url = ("https://string-db.org/api/tsv/interaction_partners"
           f"?identifiers={urllib.parse.quote(symbol)}&species={ORGANISM}")
    return http(url)


def parse_string(tsv_text):
    lines = [l for l in tsv_text.splitlines() if l]
    if not lines:
        return {}
    idx = {name: i for i, name in enumerate(lines[0].split("\t"))}
    out = {}
    for line in lines[1:]:
        c = line.split("\t")
        gene = c[idx["preferredName_B"]]
        out[gene.upper()] = Partner(
            gene=gene,
            sources=["string"],
            string_score=float(c[idx["score"]]),
            string_experiments=float(c[idx["escore"]]),
            string_textmining=float(c[idx["tscore"]]),
        )
    return out


def fetch_intact(symbol, http=httpget.get):
    url = ("https://www.ebi.ac.uk/intact/ws/interaction/findInteractions/"
           f"{urllib.parse.quote(symbol)}?page=0&pageSize=500")
    return http(url)


def parse_intact(json_text, input_accessions):
    data = json.loads(json_text)
    out = {}
    for item in data.get("content", []):
        if item.get("taxIdA") != ORGANISM or item.get("taxIdB") != ORGANISM:
            continue
        ida, idb = item.get("uniqueIdA"), item.get("uniqueIdB")
        if ida in input_accessions:
            other_id, other_gene, other_type = idb, item.get("moleculeB"), item.get("typeB")
        elif idb in input_accessions:
            other_id, other_gene, other_type = ida, item.get("moleculeA"), item.get("typeA")
        else:
            continue
        if other_type != "protein" or not other_gene or other_id in input_accessions:
            continue
        key = other_gene.upper()
        mi = item.get("intactMiscore")
        existing = out.get(key)
        if existing is None or (mi if mi is not None else 0) > (existing.intact_mi if existing.intact_mi is not None else 0):
            out[key] = Partner(gene=other_gene, partner_id=other_id,
                               sources=["intact"], intact_mi=mi)
    return out


def merge_partners(*dicts):
    merged = {}
    for d in dicts:
        for key, p in d.items():
            m = merged.get(key)
            if m is None:
                merged[key] = Partner(gene=p.gene, partner_id=p.partner_id,
                                      sources=list(p.sources),
                                      string_score=p.string_score,
                                      string_experiments=p.string_experiments,
                                      string_textmining=p.string_textmining,
                                      intact_mi=p.intact_mi)
                continue
            if not m.partner_id and p.partner_id:
                m.partner_id = p.partner_id
            for src in p.sources:
                if src not in m.sources:
                    m.sources.append(src)
            for attr in ("string_score", "string_experiments", "string_textmining", "intact_mi"):
                if getattr(m, attr) is None and getattr(p, attr) is not None:
                    setattr(m, attr, getattr(p, attr))
    return merged


def discover_protein_partners(symbol, input_accessions, http=httpget.get):
    parts = []
    try:
        parts.append(parse_string(fetch_string(symbol, http)))
    except Exception as e:
        print(f"  WARN: STRING failed: {e}")
    try:
        parts.append(parse_intact(fetch_intact(symbol, http), input_accessions))
    except Exception as e:
        print(f"  WARN: IntAct failed: {e}")
    return merge_partners(*parts)
