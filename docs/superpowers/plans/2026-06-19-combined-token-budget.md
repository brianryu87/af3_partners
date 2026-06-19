# Combined-Token Budget (5000) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the per-RNA `RNA_MAX_LEN=1500` cap with a per-job combined-token budget of 5000 applied to every emitted (input isoform × partner) pair.

**Architecture:** Add `AF3_MAX_TOKENS = 5000` to `af3.py`. Move the length check out of `_collect_partners` (which checked one RNA sequence) into `build()`'s pair-emit loop, where it can sum the input isoform length and the partner length. Remove the now-unused `RNA_MAX_LEN` from `rna.py`.

**Tech Stack:** Python 3 standard library only; stdlib `unittest`; injected fake HTTP for tests.

## Global Constraints

- Python standard library ONLY — no third-party imports in any module or test.
- `AF3_MAX_TOKENS = 5000`, defined in `af3.py` (1 residue or 1 nucleotide = 1 token).
- The budget is per job = per (input isoform × partner) pair: skip when `len(inp.sequence) + len(partner.sequence) > AF3_MAX_TOKENS`. A skipped pair produces NO JSON file and NO manifest row, and emits a warning.
- Applies uniformly to all pair kinds (protein+protein and protein+RNA).
- `RNA_MAX_LEN` must no longer exist anywhere in the codebase.
- The existing "no sequence → skip" filter for RNA partners in `_collect_partners` stays.
- Tests use stdlib `unittest`, run from project root with `python3 -m unittest ...`; no live network (injected fake http).
- Commit after every task with `git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit`.

---

### Task 1: Combined-token budget (constant + pair-level check + remove RNA_MAX_LEN)

**Files:**
- Modify: `af3.py` (add `AF3_MAX_TOKENS`)
- Modify: `make_inputs.py` (import `AF3_MAX_TOKENS`; remove RNA_MAX_LEN check in `_collect_partners`; add pair-budget check in `build()` emit loop)
- Modify: `rna.py` (remove `RNA_MAX_LEN` constant)
- Test: `tests/test_af3.py` (constant), `tests/test_build.py` (budget behavior)

**Interfaces:**
- Produces: `af3.AF3_MAX_TOKENS = 5000`
- The pair-emit loop in `make_inputs.build()` skips any pair where `len(inp.sequence) + len(p.sequence) > AF3_MAX_TOKENS`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_af3.py` inside the existing `TestAf3` class (a new method):
```python
    def test_af3_max_tokens_constant(self):
        self.assertEqual(af3.AF3_MAX_TOKENS, 5000)
```

Append to `tests/test_build.py` (reuses the module-level `FASTA`; note `FASTA` defines input RPS24 isoforms with short sequences — confirm by reading the top of the file):
```python
def fake_http_budget(url):
    # RPS24 input + two STRING partners: RPS15 (short -> emitted), BIGPROT (long -> over budget)
    if "/stream?" in url:
        return FASTA
    if "string-db.org" in url:
        return ("stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\t"
                "nscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore\n"
                "9606.E1\t9606.E2\tRPS24\tRPS15\t9606\t0.99\t0\t0\t0\t0.9\t0.99\t0.9\t0.7\n"
                "9606.E1\t9606.E3\tRPS24\tBIGPROT\t9606\t0.99\t0\t0\t0\t0.9\t0.99\t0.9\t0.7\n")
    if "intact/ws" in url:
        return '{"content": []}'
    if "encori" in url:
        return "#cite\nRBP\tgeneName\tgeneType\ttotalClipExpNum\n"
    if "cc_interaction" in url:
        return "Entry\tInteracts with\nP62847\t\n"
    if "gene_exact" in url and "RPS15" in url:
        return "Entry\tReviewed\nP62841\treviewed\n"
    if "gene_exact" in url and "BIGPROT" in url:
        return "Entry\tReviewed\nQ99999\treviewed\n"
    if url.endswith("P62841.fasta"):
        return ">sp|P62841|RS15\nMAAAA\n"            # 5 aa -> under budget
    if url.endswith("Q99999.fasta"):
        return ">sp|Q99999|BIG\n" + ("A" * 6000) + "\n"  # 6000 aa -> over budget
    raise AssertionError(f"unexpected url: {url}")


class TestBuildTokenBudget(unittest.TestCase):
    def test_over_budget_pair_skipped_under_budget_emitted(self):
        with tempfile.TemporaryDirectory() as d:
            zip_path = make_inputs.build("RPS24", d, rna_tsv=None, http=fake_http_budget)
            with zipfile.ZipFile(zip_path) as z:
                names = z.namelist()
                jsons = [n for n in names if n.endswith(".json")]
                # RPS15 emitted for both input isoforms; BIGPROT skipped for both
                self.assertTrue(jsons, "expected RPS15 JSONs")
                self.assertTrue(all("P62841" in n for n in jsons),
                                f"only RPS15 (P62841) should be emitted, got {jsons}")
                self.assertFalse(any("Q99999" in n for n in jsons),
                                 "BIGPROT (Q99999) should be skipped (over budget)")
                # manifest has no row pointing at the skipped partner
                manifest = z.read("RPS24/manifest.tsv").decode()
                self.assertNotIn("Q99999", manifest)
                self.assertIn("P62841", manifest)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m unittest tests.test_af3.TestAf3.test_af3_max_tokens_constant tests.test_build.TestBuildTokenBudget -v`
Expected: FAIL — `AF3_MAX_TOKENS` missing (AttributeError), and the build test would emit BIGPROT JSONs (no budget check yet).

- [ ] **Step 3: Add the constant to `af3.py`**

After the existing constants (currently `af3.py:5-6`):
```python
AF3_DIALECT = "alphafold3"
AF3_VERSION = 4  # latest AF3 input JSON schema version
```
add:
```python
AF3_MAX_TOKENS = 5000  # AlphaFold Server per-job total: 1 residue or 1 nucleotide = 1 token
```

- [ ] **Step 4: Remove the RNA_MAX_LEN check in `make_inputs._collect_partners`**

The current RNA loop reads:
```python
    for p in rna_mod.discover_rna_partners(symbol, rna_tsv, http,
                                           intact_text=intact_text, input_accessions=input_accs):
        if not p.sequence:
            print(f"  WARN: no sequence for RNA {p.gene}, skipping")
            continue
        if len(p.sequence) > rna_mod.RNA_MAX_LEN:
            print(f"  WARN: RNA {p.gene} sequence {len(p.sequence)} nt exceeds RNA_MAX_LEN={rna_mod.RNA_MAX_LEN}, skipping")
            continue
        result.append(p)
```
Remove the two `if len(p.sequence) > rna_mod.RNA_MAX_LEN:` lines and their `continue`, leaving:
```python
    for p in rna_mod.discover_rna_partners(symbol, rna_tsv, http,
                                           intact_text=intact_text, input_accessions=input_accs):
        if not p.sequence:
            print(f"  WARN: no sequence for RNA {p.gene}, skipping")
            continue
        result.append(p)
```

- [ ] **Step 5: Add the pair-budget check in `make_inputs.build()` and import the constant**

Change the af3 import at the top of `make_inputs.py` from:
```python
from af3 import af3_protein_pair, af3_rna_pair
```
to:
```python
from af3 import af3_protein_pair, af3_rna_pair, AF3_MAX_TOKENS
```

In `build()`, the emit loop currently begins:
```python
        for inp in input_seqs:
            for p in partners:
                name = f"{inp.isoform_id}_{p.partner_id}"
                rel = os.path.join("AF3_inputs", KIND_DIR[p.kind], p.tier, f"{name}.json")
```
Insert the budget check as the first statement inside `for p in partners:`, before `name = ...`:
```python
        for inp in input_seqs:
            for p in partners:
                if len(inp.sequence) + len(p.sequence) > AF3_MAX_TOKENS:
                    print(f"  WARN: skipping {inp.isoform_id} x {p.partner_id}: "
                          f"{len(inp.sequence) + len(p.sequence)} tokens > {AF3_MAX_TOKENS}")
                    continue
                name = f"{inp.isoform_id}_{p.partner_id}"
                rel = os.path.join("AF3_inputs", KIND_DIR[p.kind], p.tier, f"{name}.json")
```
(Leave the rest of the loop body — abspath, makedirs, job build, write, `rows.append(...)` — unchanged.)

- [ ] **Step 6: Remove the `RNA_MAX_LEN` constant from `rna.py`**

Delete the line (currently `rna.py:11`):
```python
RNA_MAX_LEN = 1500
```
Leave `ENCORI_ASSEMBLY`, `ENCORI_MAX`, and `RNACENTRAL` in place. (Confirm no other reference to `RNA_MAX_LEN` remains anywhere: `grep -rn RNA_MAX_LEN .` over the project should return nothing after this task.)

- [ ] **Step 7: Run the focused tests, then the full suite**

Run: `python3 -m unittest tests.test_af3.TestAf3.test_af3_max_tokens_constant tests.test_build.TestBuildTokenBudget -v`
Expected: PASS.

Run: `python3 -m unittest discover -s tests -t . -v`
Expected: ALL pass, pristine output. Report the total count.

Also run: `grep -rn "RNA_MAX_LEN" . --include='*.py'`
Expected: no matches (constant fully removed).

- [ ] **Step 8: Commit**

```bash
git add af3.py make_inputs.py rna.py tests/test_af3.py tests/test_build.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: cap pairs by combined 5000-token AF3 budget; drop RNA_MAX_LEN"
```

---

### Task 2: Documentation

**Files:**
- Modify: `README.md`, `CLAUDE.md`

**Interfaces:** none.

- [ ] **Step 1: Update `README.md`**

The RNA partners bullet currently reads (it mentions ENCORI/IntAct/`--rna-tsv`). Find the line describing RNA partners and confirm there is no stale `RNA_MAX_LEN` wording. Then add a sentence to the README's confidence/size description (immediately after the "Partner sources" bulleted list, before the "Confidence tiers" heading) so the size limit is documented:
```markdown
## Job size limit

Each JSON is one pairwise job. Pairs whose total length (input isoform + partner,
counting 1 residue or 1 nucleotide as 1 token) exceeds `AF3_MAX_TOKENS` (5000,
the AlphaFold Server limit) are skipped with a warning.
```
If the README already contains any `RNA_MAX_LEN` reference, replace it with the above.

- [ ] **Step 2: Update `CLAUDE.md`**

Find the note describing ENCORI/RNA sizing. The current ENCORI note reads:
```markdown
- ENCORI returns binding sites (many rows/target) and no sequences — deduped and
  bounded (`ENCORI_MAX`); usable RNA sequences come from `--rna-tsv`.
```
Append a new note after it:
```markdown
- Each emitted JSON is one pairwise job; a pair is skipped when
  `len(input isoform) + len(partner) > AF3_MAX_TOKENS` (5000, the AlphaFold
  Server per-job token limit). This replaced the old per-RNA `RNA_MAX_LEN` cap.
```
If `CLAUDE.md` contains any `RNA_MAX_LEN` reference, update it to reflect the combined-token budget.

- [ ] **Step 3: Verify suite still green**

Run: `python3 -m unittest discover -s tests -t .`
Expected: OK, same count as end of Task 1.

- [ ] **Step 4: Commit**

```bash
git add README.md CLAUDE.md
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "docs: document combined 5000-token AF3 job-size budget"
```

---

## Self-Review

**Spec coverage:**
- `AF3_MAX_TOKENS = 5000` in `af3.py` → Task 1 Step 3. ✓
- Remove `RNA_MAX_LEN` everywhere → Task 1 Steps 4 & 6 + grep check Step 7. ✓
- Pair-level combined-token skip (no JSON, no manifest row, warning) → Task 1 Step 5. ✓
- Applies to all pair kinds → the check is in the generic `for p in partners` loop, before kind-specific emit. ✓
- Keep the no-sequence RNA filter → Task 1 Step 4 leaves it. ✓
- Behavior change verified (over-budget skipped, under-budget emitted, manifest row absent for skipped) → `TestBuildTokenBudget`. ✓
- Docs updated → Task 2. ✓
- Tests stdlib unittest, no live network → fake http throughout. ✓

**Placeholder scan:** No TBD/TODO; complete code in every code step; exact run commands + expected outputs. ✓

**Type consistency:** `AF3_MAX_TOKENS` defined in `af3.py`, imported into `make_inputs.py`, referenced in the build loop and the test. `inp.sequence` / `p.sequence` are existing `InputSeq` / `Partner` fields. No new signatures. ✓

**Regression note for implementer:** existing `tests/test_build.py` tests use short sequences (well under 5000 combined), so the new budget check never trips for them — they must still pass unchanged. If any existing test fails, the budget check is mis-placed (it must be inside the `for p in partners` loop, after which the emit proceeds for under-budget pairs).
