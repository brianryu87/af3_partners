# CLAUDE.md

Generalizes the RPS24-specific AF3 input pipeline to any human gene.

## Run

    af3partners SYMBOL [--out DIR] [--rna-tsv FILE]    # after pip install
    python -m af3partners SYMBOL [--out DIR] [--rna-tsv FILE]   # from repo root

Pure stdlib. Installs via `pip install .` (hatchling, pure stdlib, no dependencies).
Network is isolated in `httpget.get` and injected into every source
function as `http=` for testability.

## Modules

- `af3partners/httpget.py` — single HTTP helper.
- `af3partners/models.py` — `InputSeq`, `Partner` dataclasses.
- `af3partners/uniprot.py` — gene → isoform sequences; accession/sequence lookup; curated interactions.
- `af3partners/partners.py` — STRING + IntAct protein-partner discovery + union.
- `af3partners/rna.py` — ENCORI + IntAct-RNA (RNAcentral) + curated-TSV RNA partners.
- `af3partners/classify.py` — ribosomal/non-ribosomal classification; confidence tiering.
- `af3partners/af3.py` — AF3-local JSON builders.
- `af3partners/manifest.py` — manifest.tsv + README.txt.
- `af3partners/make_inputs.py` — orchestration, zip, CLI.

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
