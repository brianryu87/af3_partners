# IntAct-RNA + RNAcentral RNA Partners — Design Spec

**Date:** 2026-06-19
**Status:** Approved (pending spec review)
**Extends:** [2026-06-19-af3-partners-design.md](2026-06-19-af3-partners-design.md)

## Goal

Add curated, experimentally-scored **RNA partners from IntAct** — which the
tool already queries for protein partners — resolving their sequences via
**RNAcentral**. This augments the existing RNA sources (ENCORI best-effort +
optional curated `--rna-tsv`); it does not replace them.

## Motivation

IntAct's `findInteractions` response includes RNA interactors (miRNA, lncRNA,
generic `rna`) that the current parser discards via its `typeB == "protein"`
filter. These RNA interactors carry **RNAcentral URS identifiers** (e.g.
`URS0000149452_9606`) and come with curated **MI confidence scores** — i.e.
high-quality RNA partners are already flowing into the tool and being dropped.
RNAcentral resolves the bare URS id (taxid suffix stripped) to a sequence in RNA
alphabet via `https://rnacentral.org/api/v1/rna/{URS}.json`.

## Scope

**In scope:**
- RNA interactors from IntAct whose `uniqueId` is an RNAcentral URS id.
- Human only (`taxId 9606` on both sides).
- Sequence resolution via RNAcentral only.

**Out of scope (v1):**
- mRNA interactors (Ensembl ENST ids) — usually exceed `RNA_MAX_LEN` anyway.
- DNA (`ss dna`/`ds dna`), small molecules, peptides.
- Cross-source dedup between IntAct-RNA and ENCORI/curated (kept separate).

## Single IntAct Fetch, Shared

Today `partners.discover_protein_partners` fetches IntAct internally. To avoid a
duplicate call, `_collect_partners` fetches IntAct **once** and passes the JSON
text to both paths:

- `partners.discover_protein_partners(symbol, input_accessions, http=httpget.get, intact_text=None)`
  — if `intact_text` is `None`, fetch via `fetch_intact` (backward compatible);
  otherwise use the supplied text.
- `rna.discover_rna_partners(symbol, rna_tsv, http=httpget.get, intact_text=None)`
  — uses the supplied text for IntAct-RNA discovery; if `None`, skips IntAct-RNA.

## New Units (`rna.py`)

### `parse_intact_rna(json_text, input_accessions) -> dict[str, Partner]`
From the IntAct JSON `content[]`:
- Keep an interaction only if one side's `uniqueId` is in `input_accessions`
  (the input protein) **and** the other side's `uniqueId` starts with `"URS"`
  (an RNAcentral identifier). This single rule selects all RNAcentral-resolvable
  RNA (miRNA / lncRNA / generic `rna`) and excludes mRNA (ENST), DNA, protein,
  and small molecules.
- Require `taxIdA == 9606 and taxIdB == 9606` (drops non-human RNA such as the
  observed mouse `URS..._10090`).
- Build a `Partner`: `gene` = the RNA side's molecule name, `partner_id` = the
  URS with any `_<taxid>` suffix stripped, `kind = "rna"`, `sources = ["intact"]`,
  `intact_mi = intactMiscore`, `sequence = None` (resolved later).
- Dedup by URS, keeping the row with the highest `intact_mi`.

### `fetch_rnacentral_sequence(urs, http=httpget.get) -> str | None`
- Strip any `_<taxid>` suffix from `urs`.
- GET `https://rnacentral.org/api/v1/rna/{URS}.json`, parse the `sequence` field,
  return `to_rna(sequence)` (idempotent — RNAcentral already returns U alphabet).
- On any failure (network, non-200, missing field) return `None`.

Constant: `RNACENTRAL = "https://rnacentral.org/api/v1/rna"`.

## Data Flow

`discover_rna_partners` merges three sources:
1. ENCORI (existing, best-effort).
2. Curated `--rna-tsv` (existing).
3. **IntAct-RNA (new):** if `intact_text` is supplied, `parse_intact_rna` →
   for each partner, `fetch_rnacentral_sequence(partner_id)`; set `sequence` if
   resolved, else leave `None`.

IntAct-RNA partners are deduped among themselves by URS and appended to the
result list (no cross-source reconciliation with ENCORI/curated in v1). Partners
with `sequence is None` are skipped downstream by the existing `_collect_partners`
loop, which also applies the existing `RNA_MAX_LEN` cap — so short ncRNAs
(miRNAs ~22 nt) pass and oversized lncRNAs are dropped with a warning.

## Classification & Tiering

No change to `classify.py`. IntAct-RNA partners are `kind = "rna"`; `derive_tier`
already promotes `intact_mi >= 0.45` to `high`, so well-scored IntAct-RNA
partners land in `rna/high`. Lower-scored ones fall to `rna/low` (no STRING
score). The manifest already carries `intact_mi` and `sources`.

## Naming

- Output filename: `{isoform}_{URS}.json`.
- Manifest `partner_gene` = IntAct molecule name; `partner_id` = URS.

## Error Handling

- IntAct fetch: wrapped (warn + continue), as today.
- Each RNAcentral fetch: wrapped per-URS → `None` → that partner is skipped.
- A source failure never aborts the run; zero RNA partners remains valid.

## Testing

- `parse_intact_rna`: fixture IntAct JSON with the input protein interacting with
  (a) a human URS RNA — kept; (b) a mouse URS RNA (`taxId 10090`) — dropped;
  (c) an mRNA with an ENST id — dropped; (d) a protein — dropped; (e) a self
  edge — dropped. Verify URS taxid-suffix stripping and max-MI dedup.
- `fetch_rnacentral_sequence`: fake http returning RNAcentral JSON → strips
  taxid, returns U-alphabet sequence; failure path → `None`.
- `discover_rna_partners` with `intact_text`: ENCORI empty + curated row +
  IntAct-RNA row, with a fake http serving both IntAct and RNAcentral → result
  contains the curated partner and the IntAct-RNA partner with a resolved
  sequence.
- `test_build`: extend the fake http so an IntAct RNA interactor + RNAcentral
  resolution produces an `rna/<tier>/` JSON in the zip.

## Success Criteria

- For an RBP input with human IntAct RNA interactors, the zip contains
  `AF3_inputs/rna/<tier>/<isoform>_<URS>.json` files with the RNA sequence as
  chain B.
- Non-human RNA, mRNA (ENST), DNA, protein, and self interactions are excluded.
- All existing tests still pass; new unit tests pass; no live network in tests.
- IntAct is fetched once per run (not twice).
