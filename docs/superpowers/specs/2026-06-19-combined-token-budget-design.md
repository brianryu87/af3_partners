# Combined-Token Budget (5000) — Design Spec

**Date:** 2026-06-19
**Status:** Approved (pending spec review)
**Extends:** [2026-06-19-af3-partners-design.md](2026-06-19-af3-partners-design.md)

## Goal

Replace the per-RNA length cap (`RNA_MAX_LEN = 1500`) with a per-job
**combined-token budget of 5000**, applied to every emitted pair (protein+protein
and protein+RNA). This matches the AlphaFold Server limit, where a job's total
size is the sum of all chain lengths and 1 amino acid or 1 nucleotide = 1 token.

## Motivation

An AF3 job here is one (input-protein isoform × partner) pair. The real
constraint is the *total* tokens in that job, not the length of one chain. The
current `RNA_MAX_LEN = 1500` per-RNA check is both too low and the wrong shape:
for RPS24 (133 aa) it wrongly excludes 18S rRNA (1,869 nt) even though the job
totals only ~2,002 tokens — well under 5,000. It also never checks
protein+protein pairs, which can exceed the budget for a long input isoform
paired with a very large protein.

## Design

### Constant
- Add `AF3_MAX_TOKENS = 5000` to `af3.py` (an AF3-level limit; 1 residue or
  1 nucleotide counts as 1 token).
- Remove `RNA_MAX_LEN` from `rna.py` (no longer used).

### Move the cap from partner-level to pair-level
- In `make_inputs._collect_partners`: **remove** the `RNA_MAX_LEN` length check
  on RNA partners. **Keep** the existing "no sequence → skip" filter (sequence-less
  ENCORI-only or unresolved-RNAcentral partners cannot be folded). After this,
  every partner remaining in the list has a non-empty `sequence`.
- In `make_inputs.build()`'s emit loop over `(inp × partner)`: before writing a
  pair's JSON, skip it when
  `len(inp.sequence) + len(partner.sequence) > AF3_MAX_TOKENS`. A skipped pair
  produces **no JSON file and no manifest row**, and emits a warning naming the
  isoform, partner, and total token count.

This applies uniformly to all pair kinds. Because the budget depends on the
input isoform's length, the same partner may be emitted against a short isoform
and skipped against a longer one — correct, since each pair is its own job.

## Behavior change

- RPS24 (133) + 18S rRNA (1,869) = 2,002 tokens < 5,000 → 18S rRNA is now
  emitted under `rna/<tier>/`.
- A pair whose input isoform + partner exceeds 5,000 tokens (e.g. a long input
  isoform paired with a very large protein) is now skipped with a warning,
  whereas previously protein pairs were emitted unchecked.

## Error handling

Unchanged. The over-budget skip is a warning, never an abort. Zero emitted
pairs still produces a valid zip with manifest + README (existing behavior).

## Testing

- `tests/test_af3.py`: assert `af3.AF3_MAX_TOKENS == 5000`.
- `tests/test_build.py`: a fake-http build where the input isoform paired with a
  deliberately long partner sequence exceeds 5,000 tokens is skipped — its JSON
  is absent from the zip and it has no manifest row — while an under-budget pair
  (including a long-but-within-budget RNA) is emitted. Existing tests use short
  sequences and are unaffected.
- No live network in tests (injected fake http).

## Documentation

- README.md and CLAUDE.md: describe the combined 5,000-token budget (per job =
  input isoform + partner) in place of the `RNA_MAX_LEN` description.

## Out of scope

- The README "no RNA partners" note wording (does not distinguish "none found"
  from "found but skipped") — tracked separately.
- Counting tokens as anything other than raw sequence length (no modified
  residues / ligand accounting; all partners here are protein or RNA chains).

## Success criteria

- `af3.AF3_MAX_TOKENS == 5000`; `RNA_MAX_LEN` no longer present in the codebase.
- A pair exceeding 5,000 combined tokens is skipped (no JSON, no manifest row,
  with a warning); under-budget pairs are emitted.
- A live RPS24 run produces `rna/<tier>/<isoform>_<URS>.json` for 18S rRNA.
- All existing tests still pass; new tests pass; no live network in tests.
