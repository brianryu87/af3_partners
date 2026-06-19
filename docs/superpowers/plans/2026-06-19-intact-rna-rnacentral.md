# IntAct-RNA + RNAcentral Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add curated, experimentally-scored RNA partners from IntAct (already queried for proteins), resolving their sequences via RNAcentral, into the existing `af3_partners` RNA pipeline.

**Architecture:** Two new pure-ish functions in `rna.py` (an IntAct-RNA parser and an RNAcentral sequence fetch), wired into `discover_rna_partners`. IntAct is fetched once in the orchestrator and shared with both the protein and RNA paths via a new optional `intact_text` parameter. Network stays isolated in `httpget.get` and injected for tests.

**Tech Stack:** Python 3 standard library only (`json`, `urllib`, `csv`); stdlib `unittest`; injected fake HTTP for tests.

## Global Constraints

- Python standard library ONLY — no third-party imports in any module or test.
- Human only: keep only interactions with `taxIdA == 9606 and taxIdB == 9606`.
- RNA scope: resolve ONLY interactors whose IntAct `uniqueId` starts with `"URS"` (RNAcentral). mRNA (ENST), DNA, protein, small molecules are excluded by that rule.
- Sequence source for IntAct-RNA is RNAcentral only: `RNACENTRAL = "https://rnacentral.org/api/v1/rna"`; fetch `/{URS}.json` with the `_<taxid>` suffix stripped; sequence normalized with `to_rna` (idempotent).
- IntAct must be fetched ONCE per run and shared (not fetched twice).
- All existing tests must still pass; no live network in tests (run with `python3 -m unittest ...` from the project root).
- Graceful degradation: an IntAct or RNAcentral failure warns/returns None and never aborts the run.
- Partner identity for IntAct-RNA: `partner_id` = URS (taxid stripped), `gene` = IntAct molecule name, `kind = "rna"`, `sources = ["intact"]`.
- Commit after every task with `git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit`.

---

### Task 1: RNAcentral sequence fetch

**Files:**
- Modify: `rna.py` (add `import json`, `RNACENTRAL` constant, `fetch_rnacentral_sequence`)
- Test: `tests/test_rna.py` (add `TestFetchRnacentralSequence`)

**Interfaces:**
- Consumes: `httpget.get`, `af3.to_rna` (already imported in `rna.py`)
- Produces: `rna.fetch_rnacentral_sequence(urs: str, http=httpget.get) -> str | None`; `rna.RNACENTRAL`

- [ ] **Step 1: Write the failing test** — append to `tests/test_rna.py`:

```python
class TestFetchRnacentralSequence(unittest.TestCase):
    def test_strips_taxid_and_normalizes_to_rna(self):
        captured = {}
        def fake(url):
            captured["url"] = url
            return '{"sequence": "CAAAGTGCTG", "length": 10}'
        seq = rna.fetch_rnacentral_sequence("URS0000149452_9606", http=fake)
        self.assertIn("/URS0000149452.json", captured["url"])  # taxid suffix stripped
        self.assertEqual(seq, "CAAAGUGCUG")                    # T -> U

    def test_failure_returns_none(self):
        def boom(url):
            raise OSError("rnacentral down")
        self.assertIsNone(rna.fetch_rnacentral_sequence("URS0000149452_9606", http=boom))

    def test_missing_sequence_field_returns_none(self):
        self.assertIsNone(rna.fetch_rnacentral_sequence("URS1", http=lambda u: '{"length": 0}'))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m unittest tests.test_rna.TestFetchRnacentralSequence -v`
Expected: FAIL/ERROR — `module 'rna' has no attribute 'fetch_rnacentral_sequence'`.

- [ ] **Step 3: Implement** — in `rna.py`, change the import line `import csv` block to also import `json`, add the constant near the others, and add the function. Specifically:

Change the top imports from:
```python
import csv
import urllib.parse
```
to:
```python
import csv
import json
import urllib.parse
```

Add after the existing constants (`RNA_MAX_LEN = 1500`):
```python
RNACENTRAL = "https://rnacentral.org/api/v1/rna"
```

Add this function (place it after `parse_encori`):
```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m unittest tests.test_rna.TestFetchRnacentralSequence -v`
Expected: PASS (3 tests OK).

- [ ] **Step 5: Commit**

```bash
git add rna.py tests/test_rna.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: resolve RNA sequences from RNAcentral by URS id"
```

---

### Task 2: IntAct-RNA parser

**Files:**
- Modify: `rna.py` (add `parse_intact_rna`)
- Test: `tests/test_rna.py` (add `INTACT_RNA_JSON` fixture + `TestParseIntactRna`)

**Interfaces:**
- Consumes: `models.Partner` (already imported in `rna.py`)
- Produces: `rna.parse_intact_rna(json_text: str, input_accessions: set[str]) -> dict[str, Partner]` (keyed by URS)

**IntAct JSON shape (verified):** `content[]` items carry `moleculeA`/`moleculeB` (names), `uniqueIdA`/`uniqueIdB` (accession for proteins, `URS<...>_<taxid>` for RNAcentral RNA, `ENST...` for mRNA), `intactMiscore` (float), `taxIdA`/`taxIdB` (int), `typeA`/`typeB`.

- [ ] **Step 1: Write the failing test** — append to `tests/test_rna.py` (fixture at module level near the other fixtures, class with the others):

```python
INTACT_RNA_JSON = """
{"content": [
  {"moleculeA": "ELAVL1", "moleculeB": "hsa-mir-93", "uniqueIdA": "Q15717",
   "uniqueIdB": "URS0000149452_9606", "intactMiscore": 0.56, "taxIdA": 9606,
   "taxIdB": 9606, "typeA": "protein", "typeB": "mirna"},
  {"moleculeA": "ELAVL1", "moleculeB": "lncrna_h19", "uniqueIdA": "Q15717",
   "uniqueIdB": "URS0000767B73_10090", "intactMiscore": 0.9, "taxIdA": 9606,
   "taxIdB": 10090, "typeA": "protein", "typeB": "lncrna"},
  {"moleculeA": "ELAVL1", "moleculeB": "mrna_x", "uniqueIdA": "Q15717",
   "uniqueIdB": "ENST00000407627", "intactMiscore": 0.5, "taxIdA": 9606,
   "taxIdB": 9606, "typeA": "protein", "typeB": "mrna"},
  {"moleculeA": "ELAVL1", "moleculeB": "PABPC1", "uniqueIdA": "Q15717",
   "uniqueIdB": "P11940", "intactMiscore": 0.8, "taxIdA": 9606, "taxIdB": 9606,
   "typeA": "protein", "typeB": "protein"}
]}
"""


class TestParseIntactRna(unittest.TestCase):
    def setUp(self):
        self.d = rna.parse_intact_rna(INTACT_RNA_JSON, {"Q15717"})

    def test_keeps_human_urs_rna(self):
        self.assertIn("URS0000149452", self.d)
        p = self.d["URS0000149452"]
        self.assertEqual(p.kind, "rna")
        self.assertEqual(p.partner_id, "URS0000149452")  # taxid stripped
        self.assertEqual(p.gene, "hsa-mir-93")
        self.assertEqual(p.intact_mi, 0.56)
        self.assertEqual(p.sources, ["intact"])
        self.assertIsNone(p.sequence)

    def test_drops_non_human_urs(self):
        self.assertNotIn("URS0000767B73", self.d)

    def test_drops_mrna_enst_and_protein(self):
        self.assertTrue(all(k.startswith("URS") for k in self.d))
        self.assertNotIn("P11940", self.d)

    def test_only_one_kept(self):
        self.assertEqual(len(self.d), 1)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python3 -m unittest tests.test_rna.TestParseIntactRna -v`
Expected: FAIL/ERROR — `module 'rna' has no attribute 'parse_intact_rna'`.

- [ ] **Step 3: Implement** — add to `rna.py` (after `fetch_rnacentral_sequence`):

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python3 -m unittest tests.test_rna.TestParseIntactRna -v`
Expected: PASS (4 tests OK).

- [ ] **Step 5: Commit**

```bash
git add rna.py tests/test_rna.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: parse RNAcentral RNA interactors from IntAct payload"
```

---

### Task 3: Wire IntAct-RNA into discovery + single shared IntAct fetch

**Files:**
- Modify: `rna.py` (`discover_rna_partners` — new params + IntAct-RNA branch)
- Modify: `partners.py` (`discover_protein_partners` — optional `intact_text`)
- Modify: `make_inputs.py` (`_collect_partners` — fetch IntAct once, pass to both)
- Test: `tests/test_rna.py` (add `TestDiscoverRnaWithIntact`), `tests/test_build.py` (add `TestBuildWithIntactRna`)

**Interfaces:**
- Consumes: `rna.parse_intact_rna`, `rna.fetch_rnacentral_sequence` (Tasks 1-2); `partners.fetch_intact`, `partners.parse_intact`, `partners.parse_string` (existing).
- Produces:
  - `rna.discover_rna_partners(symbol, rna_tsv, http=httpget.get, intact_text=None, input_accessions=None) -> list[Partner]`
  - `partners.discover_protein_partners(symbol, input_accessions, http=httpget.get, intact_text=None) -> dict[str, Partner]`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_rna.py` (uses `INTACT_RNA_JSON` from Task 2):
```python
class TestDiscoverRnaWithIntact(unittest.TestCase):
    def test_intact_rna_resolved_and_appended(self):
        def fake(url):
            if "rnacentral" in url:
                return '{"sequence": "CAAAGUGCUGUUCGUGCAGGUAG", "length": 23}'
            if "encori" in url:
                return "#cite\nRBP\tgeneName\tgeneType\ttotalClipExpNum\n"
            raise AssertionError(f"unexpected url: {url}")
        parts = rna.discover_rna_partners(
            "ELAVL1", None, http=fake,
            intact_text=INTACT_RNA_JSON, input_accessions={"Q15717"})
        mirs = [p for p in parts if p.partner_id == "URS0000149452"]
        self.assertEqual(len(mirs), 1)
        self.assertEqual(mirs[0].sequence, "CAAAGUGCUGUUCGUGCAGGUAG")
        self.assertEqual(mirs[0].kind, "rna")
        self.assertEqual(mirs[0].sources, ["intact"])

    def test_no_intact_text_means_no_intact_rna(self):
        parts = rna.discover_rna_partners(
            "ELAVL1", None, http=lambda u: "#c\nRBP\tgeneName\tgeneType\ttotalClipExpNum\n")
        self.assertEqual(parts, [])
```

Append to `tests/test_build.py` (reuses the module-level `FASTA`):
```python
INTACT_RNA_BUILD = ('{"content": [{"moleculeA": "RPS24", "moleculeB": "mir-x",'
                    ' "uniqueIdA": "P62847", "uniqueIdB": "URS0000ABCDEF_9606",'
                    ' "intactMiscore": 0.6, "taxIdA": 9606, "taxIdB": 9606,'
                    ' "typeA": "protein", "typeB": "mirna"}]}')


def fake_http_rna(url):
    if "/stream?" in url:
        return FASTA
    if "string-db.org" in url:
        return ("stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\t"
                "nscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore\n")
    if "intact/ws" in url:
        return INTACT_RNA_BUILD
    if "encori" in url:
        return "#cite\nRBP\tgeneName\tgeneType\ttotalClipExpNum\n"
    if "rnacentral" in url:
        return '{"sequence": "ACGUACGUACGU", "length": 12}'
    if "cc_interaction" in url:
        return "Entry\tInteracts with\nP62847\t\n"
    raise AssertionError(f"unexpected url: {url}")


class TestBuildWithIntactRna(unittest.TestCase):
    def test_intact_rna_produces_rna_json(self):
        with tempfile.TemporaryDirectory() as d:
            zip_path = make_inputs.build("RPS24", d, rna_tsv=None, http=fake_http_rna)
            with zipfile.ZipFile(zip_path) as z:
                jsons = [n for n in z.namelist() if n.endswith(".json")]
                rna_jsons = [n for n in jsons if "/rna/" in n]
                self.assertTrue(rna_jsons, f"expected an rna/ JSON, got {jsons}")
                self.assertTrue(any("URS0000ABCDEF" in n for n in rna_jsons))
                sample = json.loads(z.read(rna_jsons[0]))[0]
                self.assertIn("rna", sample["sequences"][1])
                self.assertEqual(sample["sequences"][1]["rna"]["sequence"], "ACGUACGUACGU")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python3 -m unittest tests.test_rna.TestDiscoverRnaWithIntact tests.test_build.TestBuildWithIntactRna -v`
Expected: FAIL — `discover_rna_partners` got unexpected keyword argument `intact_text` (and the build test produces no `rna/` JSON).

- [ ] **Step 3a: Implement `rna.discover_rna_partners`** — replace the existing function body. Change the signature and append the IntAct-RNA branch:

```python
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
```

- [ ] **Step 3b: Implement `partners.discover_protein_partners`** — replace the existing function (`partners.py:89-99`) so IntAct text can be supplied (fetched only if not):

```python
def discover_protein_partners(symbol, input_accessions, http=httpget.get, intact_text=None):
    parts = []
    try:
        parts.append(parse_string(fetch_string(symbol, http)))
    except Exception as e:
        print(f"  WARN: STRING failed: {e}")
    try:
        if intact_text is None:
            intact_text = fetch_intact(symbol, http)
        parts.append(parse_intact(intact_text, input_accessions))
    except Exception as e:
        print(f"  WARN: IntAct failed: {e}")
    return merge_partners(*parts)
```

- [ ] **Step 3c: Implement the `make_inputs._collect_partners` wiring** — in `make_inputs.py`, fetch IntAct once and pass it to both paths.

Replace this (currently `make_inputs.py:25`):
```python
    proteins = partners_mod.discover_protein_partners(symbol, input_accs, http)
```
with:
```python
    try:
        intact_text = partners_mod.fetch_intact(symbol, http)
    except Exception as e:
        print(f"  WARN: IntAct fetch failed: {e}")
        intact_text = None

    proteins = partners_mod.discover_protein_partners(symbol, input_accs, http, intact_text=intact_text)
```

And replace the RNA discovery call (currently `make_inputs.py:45`):
```python
    for p in rna_mod.discover_rna_partners(symbol, rna_tsv, http):
```
with:
```python
    for p in rna_mod.discover_rna_partners(symbol, rna_tsv, http,
                                           intact_text=intact_text, input_accessions=input_accs):
```
(Leave the rest of that loop — the no-sequence skip and the `RNA_MAX_LEN` cap — unchanged.)

- [ ] **Step 4: Run the focused tests, then the full suite**

Run: `python3 -m unittest tests.test_rna.TestDiscoverRnaWithIntact tests.test_build.TestBuildWithIntactRna -v`
Expected: PASS.

Run: `python3 -m unittest discover -s tests -t . -v`
Expected: ALL tests pass, output pristine. Report the total count (should be the prior count + the new tests).

- [ ] **Step 5: Commit**

```bash
git add rna.py partners.py make_inputs.py tests/test_rna.py tests/test_build.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: discover IntAct RNA partners via RNAcentral; fetch IntAct once and share"
```

---

### Task 4: Documentation

**Files:**
- Modify: `README.md` (RNA partners bullet)
- Modify: `CLAUDE.md` (rna.py description + notes)

**Interfaces:** none.

- [ ] **Step 1: Update `README.md`** — replace the line (`README.md:23`):
```markdown
- RNA partners: ENCORI (best-effort; bounded) + optional `--rna-tsv` (gene, sequence).
```
with:
```markdown
- RNA partners: ENCORI (best-effort; bounded), curated RNA from IntAct (RNAcentral
  URS interactors, sequences resolved via RNAcentral), and optional `--rna-tsv`
  (gene, sequence).
```

- [ ] **Step 2: Update `CLAUDE.md`** — replace the line (`CLAUDE.md:18`):
```markdown
- `rna.py` — ENCORI + curated-TSV RNA partners.
```
with:
```markdown
- `rna.py` — ENCORI + IntAct-RNA (RNAcentral) + curated-TSV RNA partners.
```

Then append after the ENCORI note (`CLAUDE.md:31-32`):
```markdown
- IntAct also reports RNA interactors (miRNA/lncRNA/ncRNA) carrying RNAcentral
  URS ids; these are kept when `uniqueId` starts with `URS` and both sides are
  human, with sequences resolved from RNAcentral (`RNACENTRAL`). IntAct is
  fetched once per run and shared between the protein and RNA paths.
```

- [ ] **Step 3: Verify suite still green** (docs only, but confirm nothing broke)

Run: `python3 -m unittest discover -s tests -t .`
Expected: OK (same count as end of Task 3).

- [ ] **Step 4: Commit**

```bash
git add README.md CLAUDE.md
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "docs: document IntAct-RNA + RNAcentral RNA source"
```

---

## Self-Review

**Spec coverage:**
- IntAct RNA interactors kept by `uniqueId` startswith `URS` → Task 2 `parse_intact_rna`. ✓
- Human-only `taxId 9606` filter → Task 2. ✓
- URS taxid-suffix stripping → Tasks 1 & 2. ✓
- RNAcentral sequence resolution, `to_rna`, failure→None → Task 1. ✓
- Single IntAct fetch shared between protein + RNA paths → Task 3 (orchestrator fetch + optional `intact_text` on both discovery functions). ✓
- Merge as third RNA source; sequence-less skipped + `RNA_MAX_LEN` cap reused → Task 3 (append to result; existing `_collect_partners` loop unchanged). ✓
- Tiering unchanged (`intact_mi ≥ 0.45 → high` already in `derive_tier`) → no classify task needed. ✓
- Naming: filename `{isoform}_{URS}`, manifest `partner_gene`=molecule, `partner_id`=URS → Task 2 sets fields; existing build/manifest use them. ✓
- Docs updated → Task 4. ✓
- Tests w/o live network → all tasks use injected fake http. ✓

**Placeholder scan:** No TBD/TODO; every code step shows complete code; exact run commands + expected outcomes given. ✓

**Type consistency:** `parse_intact_rna(json_text, input_accessions) -> dict[URS, Partner]`; `fetch_rnacentral_sequence(urs, http) -> str|None`; `discover_rna_partners(..., intact_text=None, input_accessions=None)`; `discover_protein_partners(..., intact_text=None)`. The orchestrator passes `intact_text=intact_text` and `input_accessions=input_accs` (the existing local `input_accs` set at `make_inputs.py:22`). `Partner` fields used (`gene`, `partner_id`, `kind`, `sources`, `intact_mi`, `sequence`) all exist. Consistent across tasks. ✓

**Regression note for implementer:** the existing `test_build.TestBuild` and `TestBuildZeroPartners` use a `fake_http` whose `intact/ws` branch returns `{"content": []}` — after Task 3 the orchestrator fetches IntAct once and passes it to `discover_rna_partners`; with empty content `parse_intact_rna` returns `{}` and `fetch_rnacentral_sequence` is never called, so those tests need no change and must still pass.
