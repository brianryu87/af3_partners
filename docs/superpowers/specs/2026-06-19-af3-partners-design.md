# af3_partners — Design Spec

**Date:** 2026-06-19
**Status:** Approved (pending spec review)

## Goal

A standalone CLI tool that takes a **gene symbol** and produces a **zip** of
AlphaFold3 (AF3) local-format input JSON files. Each JSON pairs the input
protein (one isoform) with one interacting partner — a ribosomal protein, a
non-ribosomal protein, or an RNA. Partners are discovered automatically from
public databases; every partner is included and **marked by confidence** so
low-quality pairs are visible but never silently dropped.

This generalizes the RPS24-specific pipeline in `projects/rps24/` to any human
gene.

## Scope

**In scope (v1):**
- Input: a single human gene symbol.
- Output: `SYMBOL.zip` of AF3-local-format JSONs + a manifest + a README.
- Partner discovery from STRING, IntAct, UniProt, and (best-effort) ENCORI.
- Optional curated RNA partners via `--rna-tsv`.
- Confidence tiering (high/medium/low) via manifest + tiered folders.

**Out of scope (v1):**
- Non-human organisms (fixed to taxon 9606).
- AF-Server JSON format (AF3-local only).
- RNA transcript trimming heuristics beyond a simple length cap.

## Interface

```
python make_inputs.py SYMBOL [--out DIR] [--rna-tsv FILE]
```

- `SYMBOL` — human gene symbol (e.g. `RPS24`).
- `--out DIR` — directory to write `SYMBOL.zip` into (default: current dir).
- `--rna-tsv FILE` — optional curated RNA partners: TSV with columns
  `gene`, `sequence` (DNA or RNA alphabet; T→U normalized). Always included.

Pure Python standard library only (`urllib`, `json`, `zipfile`, `csv`,
`argparse`), matching the workspace convention. No third-party dependencies.

## External references — always latest release

Every external data source is queried at its **current/latest release**; no
source is hard-pinned to an old version. Where a version or assembly must be
named, it lives in a single named constant so it is trivial to bump:

| Source | Endpoint | Latest pin (constant) |
|---|---|---|
| UniProt | `rest.uniprot.org/uniprotkb` | current release (unversioned REST) |
| STRING | `string-db.org/api` | latest stable API (`STRING_API_VERSION`) |
| IntAct | `ebi.ac.uk/intact/ws` | current ws |
| ENCORI/starBase | `rnasysu.com/encori/api` | assembly `ENCORI_ASSEMBLY = "hg38"` (GRCh38) |
| AF3 JSON schema | n/a | `AF3_VERSION` = latest input schema version (currently 4), dialect `alphafold3` |

## Pipeline

Each stage is an isolated, independently testable unit.

### 1. Resolve input gene → sequences
UniProt REST query for **all human entries** of the symbol — reviewed
(Swiss-Prot) **and** unreviewed/computational (TrEMBL, e.g. `A0A2R8…`) — plus
the isoforms of reviewed entries. Returns a list of input "chain A" sequences:
`{accession, isoform_id, sequence, reviewed: bool}`. If no entry is found, exit
with a clear error.

### 2. Discover protein partners
Union of three sources, keyed by **UniProt accession**:
- **STRING** `interaction_partners` (tsv): combined `score` plus per-channel
  scores (`escore` experiments, `dscore` database, `tscore` textmining, etc.).
  STRING returns Ensembl protein ids (`ENSP…`) and `preferredName`; map to
  UniProt accession via UniProt ID-mapping (gene name → accession is sufficient
  for human, with ENSP as fallback).
- **IntAct** `ws/interaction/findInteractions/<symbol>` (JSON): MI confidence
  score. Results mix organisms — **filter to human UniProt accessions** using
  the `idA`/`idB` / `uniqueId` fields.
- **UniProt** `cc_interaction` field — included when present (frequently empty;
  not relied upon as the primary source).

Each partner accumulates a record:
`{accession, gene, kind, sources: {string, intact, uniprot}, scores: {...}}`.

### 3. Discover RNA partners (hybrid)
- **ENCORI** `RBPTarget` (best-effort): returns RNA targets only when the input
  gene is a catalogued RNA-binding protein from CLIP-seq; otherwise empty. The
  README records why RNA is empty when it is.
- **`--rna-tsv`** (always honored): curated `gene`, `sequence` rows. Reliable
  source of usable RNA sequences.

### 4. Classify protein partners
- **ribosomal** if the HGNC symbol matches `^(RPS|RPL|MRPS|MRPL)\d` OR is in the
  explicit set `{RPSA, FAU, RACK1, UBA52, RPLP0, RPLP1, RPLP2}`.
- **non-ribosomal** otherwise.
- RNA partners are their own kind.

### 5. Fetch partner sequences
- Protein partners: UniProt sequence by accession; fetches deduplicated/cached.
- RNA partners: sequence from `--rna-tsv` (preferred) or ENCORI/Ensembl
  best-effort; normalized T→U; capped at a max length constant
  (`RNA_MAX_LEN`). Over-length RNA is recorded in the manifest as skipped.

### 6. Derive confidence tier
Named-constant thresholds:
- **high** — any experimental evidence (STRING `escore` > 0, IntAct MI ≥ 0.45,
  UniProt-curated, or ENCORI CLIP) **or** STRING combined score ≥ 0.7.
- **medium** — STRING combined score in [0.4, 0.7) without strong experimental
  evidence.
- **low** — STRING combined < 0.4, or text-mining-only / prediction-only.

### 7. Emit AF3-local JSON
One JSON per (input isoform × partner), reusing the structure from
`projects/rps24/gen/make_inputs.py` (`af3_protein_pair` / `af3_rna_pair`):
dialect `alphafold3`, version `AF3_VERSION`, chain A = input isoform, chain B =
partner (`protein` or `rna` key). Filename `{inputAcc}_{partnerId}.json`.

## Output: zip layout

```
SYMBOL/
  AF3_inputs/
    ribosomal_protein/{high,medium,low}/{inputAcc}_{partnerAcc}.json
    nonribosomal_protein/{high,medium,low}/{inputAcc}_{partnerAcc}.json
    rna/{high,medium,low}/{inputAcc}_{rnaGene}.json
  manifest.tsv
  README.txt
```

**`manifest.tsv` columns:**
`input_gene, input_accession, input_isoform, input_reviewed, partner_id,
partner_gene, partner_kind, sources, string_score, string_experiments,
string_textmining, intact_mi, uniprot_curated, encori_evidence, tier,
json_path`

**`README.txt`:** symbol, run date, sources queried + their resolved
release/assembly, partner counts per kind/tier, and caveats (e.g. "RNA empty:
gene is not a catalogued RBP and no --rna-tsv supplied").

## Error handling

- Each external source is wrapped: on network error / non-200 / parse failure,
  emit a warning and continue. A STRING outage or an empty RNA result never
  aborts the run.
- Gene-not-found in UniProt → clear error, non-zero exit.
- Zero partners discovered → still produce a valid zip with an empty tree and a
  manifest/README noting zero partners.

## Module decomposition

- `uniprot.py` — gene → input sequences; accession → sequence; ID mapping.
- `partners.py` — STRING + IntAct + UniProt protein-partner discovery; union.
- `rna.py` — ENCORI query + `--rna-tsv` loader.
- `classify.py` — ribosomal/non-ribosomal classification; tier derivation.
- `af3.py` — AF3-local JSON builders (ported from rps24).
- `manifest.py` — manifest row assembly + TSV/README writers.
- `make_inputs.py` — CLI, orchestration, zip assembly.

## Testing

- **Unit tests against recorded (fixture) API responses** — no live network:
  STRING parse, IntAct parse + human filter, UniProt isoform resolution,
  classification rule, tier derivation, AF3 JSON shape, manifest rows, zip
  assembly.
- **One opt-in live smoke test** on `RPS24` (guarded by an env flag), since live
  network in CI is fragile.

## Success criteria

- `python make_inputs.py RPS24` produces `RPS24.zip` containing AF3-local JSONs
  partitioned by partner kind and confidence tier, a populated `manifest.tsv`,
  and a `README.txt`.
- Each JSON validates as AF3-local (dialect `alphafold3`, latest version) with
  exactly two chains (input isoform + partner).
- Ribosomal partners (e.g. RPS15, RPS2, RPS8 for RPS24) land in
  `ribosomal_protein/`; high-confidence experimental partners land in `high/`.
- All unit tests pass with no live network.
- A source outage degrades gracefully (warning, partial output) rather than
  crashing.
