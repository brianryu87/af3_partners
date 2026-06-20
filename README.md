# af3_partners

Turn a human gene symbol into a zip of AlphaFold3 (AF3) local-format input JSONs,
one per (input-protein isoform × interacting partner).

## Install

    git clone https://github.com/brianryu87/af3_partners.git
    cd af3_partners
    pip install .          # or: pip install -e .   (for development)

## Usage

    af3partners RPS24 --out .
    af3partners RPS24 --rna-tsv my_rna_partners.tsv

Without installing, run it from the repo root with:

    python -m af3partners RPS24

Output: `RPS24.zip` containing:

    RPS24/
      AF3_inputs/{ribosomal_protein,nonribosomal_protein,rna}/{high,medium,low}/*.json
      manifest.tsv   # per-pair sources, scores, tier
      README.txt     # run summary and caveats

## Partner sources (latest release of each)

- Sequences & isoforms: UniProt (reviewed + computational/TrEMBL).
- Protein partners: STRING (combined + per-channel scores) and IntAct (MI score).
- Curated UniProt binary interactions mark `uniprot_curated`.
- RNA partners: ENCORI (best-effort; bounded), curated RNA from IntAct (RNAcentral
  URS interactors, sequences resolved via RNAcentral), and optional `--rna-tsv`
  (gene, sequence).

## Job size limit

Each JSON is one pairwise job. Pairs whose total length (input isoform + partner,
counting 1 residue or 1 nucleotide as 1 token) exceeds `AF3_MAX_TOKENS` (5000,
the AlphaFold Server limit) are skipped with a warning.

## Confidence tiers

- high: experimental evidence (STRING experiments channel, IntAct MI ≥ 0.45,
  UniProt-curated, or ENCORI CLIP) or STRING combined ≥ 0.7
- medium: STRING combined in [0.4, 0.7)
- low: below 0.4 / text-mining-only

## Tests

    python -m unittest discover -s tests -t .
    AF3_LIVE=1 python -m unittest tests.test_live_smoke   # opt-in, hits live APIs

Pure Python standard library; no third-party dependencies.
