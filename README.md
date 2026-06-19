# af3_partners

Turn a human gene symbol into a zip of AlphaFold3 (AF3) local-format input JSONs,
one per (input-protein isoform × interacting partner).

## Usage

    python make_inputs.py RPS24 --out .
    python make_inputs.py RPS24 --rna-tsv my_rna_partners.tsv

Output: `RPS24.zip` containing:

    RPS24/
      AF3_inputs/{ribosomal_protein,nonribosomal_protein,rna}/{high,medium,low}/*.json
      manifest.tsv   # per-pair sources, scores, tier
      README.txt     # run summary and caveats

## Partner sources (latest release of each)

- Sequences & isoforms: UniProt (reviewed + computational/TrEMBL).
- Protein partners: STRING (combined + per-channel scores) and IntAct (MI score).
- Curated UniProt binary interactions mark `uniprot_curated`.
- RNA partners: ENCORI (best-effort; bounded) + optional `--rna-tsv` (gene, sequence).

## Confidence tiers

- high: experimental evidence (STRING experiments channel, IntAct MI ≥ 0.45,
  UniProt-curated, or ENCORI CLIP) or STRING combined ≥ 0.7
- medium: STRING combined in [0.4, 0.7)
- low: below 0.4 / text-mining-only

## Tests

    python -m unittest discover -s tests -t .
    AF3_LIVE=1 python -m unittest tests.test_live_smoke   # opt-in, hits live APIs

Pure Python standard library; no third-party dependencies.
