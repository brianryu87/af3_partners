# CLAUDE.md

Generalizes the RPS24-specific AF3 input pipeline to any human gene.

## Run

    python make_inputs.py SYMBOL [--out DIR] [--rna-tsv FILE]   # writes SYMBOL.zip

Pure stdlib. Network is isolated in `httpget.get` and injected into every source
function as `http=` for testability.

## Modules

- `httpget.py` — single HTTP helper.
- `models.py` — `InputSeq`, `Partner` dataclasses.
- `uniprot.py` — gene → isoform sequences; accession/sequence lookup; curated interactions.
- `partners.py` — STRING + IntAct protein-partner discovery + union.
- `rna.py` — ENCORI + IntAct-RNA (RNAcentral) + curated-TSV RNA partners.
- `classify.py` — ribosomal/non-ribosomal classification; confidence tiering.
- `af3.py` — AF3-local JSON builders.
- `manifest.py` — manifest.tsv + README.txt.
- `make_inputs.py` — orchestration, zip, CLI.

## Notes

- Human only (`ORGANISM = 9606`).
- All external references pinned to latest via named constants (`AF3_VERSION`,
  `ENCORI_ASSEMBLY`); never hard-pin an old release.
- IntAct `findInteractions` substring-matches names — always filter by input
  UniProt accessions and `taxId == 9606`.
- ENCORI returns binding sites (many rows/target) and no sequences — deduped and
  bounded (`ENCORI_MAX`); usable RNA sequences come from `--rna-tsv`.
- Each emitted JSON is one pairwise job; a pair is skipped when
  `len(input isoform) + len(partner) > AF3_MAX_TOKENS` (5000, the AlphaFold
  Server per-job token limit). This replaced the old per-RNA `RNA_MAX_LEN` cap.
- IntAct also reports RNA interactors (miRNA/lncRNA/ncRNA) carrying RNAcentral
  URS ids; these are kept when `uniqueId` starts with `URS` and both sides are
  human, with sequences resolved from RNAcentral (`RNACENTRAL`). IntAct is
  fetched once per run and shared between the protein and RNA paths.
- Tests: stdlib `unittest`, fixtures + injected fake HTTP; live smoke gated by `AF3_LIVE=1`.
