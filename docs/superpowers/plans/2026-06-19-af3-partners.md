# af3_partners Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A zero-dependency CLI that turns a human gene symbol into a zip of AlphaFold3 local-format input JSONs — one per (input-protein isoform × interacting partner) — with a confidence manifest and tiered folders.

**Architecture:** A linear pipeline of small, single-responsibility modules. Network access is isolated in one `httpget` helper and injected into every source function as a default argument, so all parsing logic is unit-tested against recorded fixtures with a fake HTTP. The orchestrator wires the stages and assembles the zip.

**Tech Stack:** Python 3 standard library only (`urllib`, `json`, `csv`, `zipfile`, `argparse`, `dataclasses`, `re`). Tests use stdlib `unittest`. No third-party packages at runtime or test time.

## Global Constraints

- Python standard library only — no third-party imports in any module or test.
- Human only: `ORGANISM = 9606` everywhere.
- Every external reference uses its latest release; versions/assemblies live in named constants, never hard-pinned to an old release: `AF3_VERSION` (latest AF3 input schema, currently 4), `ENCORI_ASSEMBLY = "hg38"`, UniProt/STRING/IntAct unversioned/current endpoints.
- AF3-local format only: dialect `"alphafold3"`, key `"protein"`/`"rna"`, top-level JSON is a single-element list `[job]`.
- All discovered partners are included; none are dropped for low confidence — they are marked via manifest + `high`/`medium`/`low` folders.
- Run from the project root `/home/brianryu87/projects/af3_partners`. Tests run with `python -m unittest tests.test_<name>`.
- Commit after every task with `git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit`.

---

### Task 1: Scaffold, data models, HTTP helper, AF3 builders

**Files:**
- Create: `httpget.py`, `models.py`, `af3.py`, `tests/__init__.py`, `tests/test_af3.py`

**Interfaces:**
- Produces:
  - `httpget.get(url: str, accept: str | None = None) -> str`
  - `models.InputSeq(accession: str, isoform_id: str, sequence: str, reviewed: bool)`
  - `models.Partner(gene: str, partner_id: str = "", kind: str = "", sequence: str | None = None, sources: list = [], string_score: float | None = None, string_experiments: float | None = None, string_textmining: float | None = None, intact_mi: float | None = None, uniprot_curated: bool = False, encori_evidence: str | None = None, tier: str = "")`
  - `af3.AF3_DIALECT = "alphafold3"`, `af3.AF3_VERSION = 4`
  - `af3.to_rna(seq: str) -> str`
  - `af3.af3_protein_pair(name: str, inp: InputSeq, partner: Partner) -> dict`
  - `af3.af3_rna_pair(name: str, inp: InputSeq, partner: Partner) -> dict`

- [ ] **Step 1: Create `httpget.py`**

```python
import urllib.request

TIMEOUT = 30


def get(url, accept=None):
    """GET a URL, return decoded body text. Raises urllib.error.URLError on failure."""
    req = urllib.request.Request(url)
    if accept:
        req.add_header("Accept", accept)
    with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
        return resp.read().decode("utf-8")
```

- [ ] **Step 2: Create `models.py`**

```python
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class InputSeq:
    accession: str        # base accession, e.g. "P62847"
    isoform_id: str       # canonical or isoform id, e.g. "P62847" or "P62847-2"
    sequence: str
    reviewed: bool


@dataclass
class Partner:
    gene: str
    partner_id: str = ""                       # UniProt accession (protein) or gene (rna)
    kind: str = ""                             # "ribosomal" | "nonribosomal" | "rna"
    sequence: Optional[str] = None
    sources: list = field(default_factory=list)
    string_score: Optional[float] = None
    string_experiments: Optional[float] = None
    string_textmining: Optional[float] = None
    intact_mi: Optional[float] = None
    uniprot_curated: bool = False
    encori_evidence: Optional[str] = None
    tier: str = ""
```

- [ ] **Step 3: Write the failing test `tests/test_af3.py`**

```python
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import af3
from models import InputSeq, Partner


class TestAf3(unittest.TestCase):
    def test_to_rna_converts_t_to_u(self):
        self.assertEqual(af3.to_rna("ACGTacgt"), "ACGUacgu")

    def test_protein_pair_shape(self):
        inp = InputSeq("P62847", "P62847", "MND", True)
        partner = Partner(gene="RPS15", partner_id="P62841", sequence="AAA")
        job = af3.af3_protein_pair("P62847_P62841", inp, partner)
        self.assertEqual(job["dialect"], "alphafold3")
        self.assertEqual(job["version"], 4)
        self.assertEqual(job["name"], "P62847_P62841")
        self.assertEqual(len(job["sequences"]), 2)
        self.assertEqual(job["sequences"][0]["protein"]["id"], "A")
        self.assertEqual(job["sequences"][0]["protein"]["sequence"], "MND")
        self.assertEqual(job["sequences"][1]["protein"]["id"], "B")
        self.assertEqual(job["sequences"][1]["protein"]["sequence"], "AAA")
        self.assertEqual(len(job["modelSeeds"]), 1)

    def test_rna_pair_uses_rna_key_and_uracil(self):
        inp = InputSeq("P62847", "P62847", "MND", True)
        partner = Partner(gene="ZFAS1", partner_id="ZFAS1", sequence="ACGT")
        job = af3.af3_rna_pair("P62847_ZFAS1", inp, partner)
        self.assertIn("rna", job["sequences"][1])
        self.assertEqual(job["sequences"][1]["rna"]["sequence"], "ACGU")
        self.assertEqual(job["sequences"][1]["rna"]["description"], "ZFAS1")


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 4: Run test to verify it fails**

Run: `python -m unittest tests.test_af3 -v`
Expected: FAIL / ERROR with "No module named 'af3'".

- [ ] **Step 5: Create `af3.py`**

```python
import random

from models import InputSeq, Partner

AF3_DIALECT = "alphafold3"
AF3_VERSION = 4  # latest AF3 input JSON schema version


def to_rna(seq):
    return seq.replace("T", "U").replace("t", "u")


def af3_protein_pair(name, inp: InputSeq, partner: Partner):
    return {
        "name": name,
        "modelSeeds": [random.randint(1, 999)],
        "sequences": [
            {"protein": {"id": "A", "sequence": inp.sequence, "description": inp.isoform_id}},
            {"protein": {"id": "B", "sequence": partner.sequence, "description": partner.partner_id}},
        ],
        "dialect": AF3_DIALECT,
        "version": AF3_VERSION,
    }


def af3_rna_pair(name, inp: InputSeq, partner: Partner):
    return {
        "name": name,
        "modelSeeds": [random.randint(1, 999)],
        "sequences": [
            {"protein": {"id": "A", "sequence": inp.sequence, "description": inp.isoform_id}},
            {"rna": {"id": "B", "sequence": to_rna(partner.sequence), "description": partner.gene}},
        ],
        "dialect": AF3_DIALECT,
        "version": AF3_VERSION,
    }
```

- [ ] **Step 6: Create empty `tests/__init__.py`** (so `tests` is a package)

- [ ] **Step 7: Run test to verify it passes**

Run: `python -m unittest tests.test_af3 -v`
Expected: PASS (3 tests OK).

- [ ] **Step 8: Commit**

```bash
git add httpget.py models.py af3.py tests/__init__.py tests/test_af3.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: scaffold models, httpget, AF3 builders"
```

---

### Task 2: Classification and confidence tiering

**Files:**
- Create: `classify.py`, `tests/test_classify.py`

**Interfaces:**
- Consumes: `models.Partner`
- Produces:
  - `classify.classify_kind(gene: str) -> str` → `"ribosomal"` | `"nonribosomal"`
  - `classify.derive_tier(p: Partner) -> str` → `"high"` | `"medium"` | `"low"`
  - Constants `STRING_HIGH = 0.7`, `STRING_MED = 0.4`, `INTACT_EXPERIMENTAL = 0.45`

- [ ] **Step 1: Write the failing test `tests/test_classify.py`**

```python
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import classify
from models import Partner


class TestClassifyKind(unittest.TestCase):
    def test_rps_is_ribosomal(self):
        self.assertEqual(classify.classify_kind("RPS15"), "ribosomal")

    def test_mrpl_is_ribosomal(self):
        self.assertEqual(classify.classify_kind("MRPL12"), "ribosomal")

    def test_explicit_set_is_ribosomal(self):
        self.assertEqual(classify.classify_kind("RACK1"), "ribosomal")
        self.assertEqual(classify.classify_kind("FAU"), "ribosomal")

    def test_other_is_nonribosomal(self):
        self.assertEqual(classify.classify_kind("TP53"), "nonribosomal")

    def test_rps24_itself_is_ribosomal_but_case_insensitive(self):
        self.assertEqual(classify.classify_kind("rps24"), "ribosomal")


class TestDeriveTier(unittest.TestCase):
    def test_string_experiments_makes_high(self):
        p = Partner(gene="X", string_score=0.2, string_experiments=0.5)
        self.assertEqual(classify.derive_tier(p), "high")

    def test_intact_above_threshold_makes_high(self):
        p = Partner(gene="X", intact_mi=0.5)
        self.assertEqual(classify.derive_tier(p), "high")

    def test_high_combined_score_makes_high(self):
        p = Partner(gene="X", string_score=0.8)
        self.assertEqual(classify.derive_tier(p), "high")

    def test_mid_score_makes_medium(self):
        p = Partner(gene="X", string_score=0.5, string_textmining=0.5)
        self.assertEqual(classify.derive_tier(p), "medium")

    def test_low_score_textmining_only_makes_low(self):
        p = Partner(gene="X", string_score=0.2, string_textmining=0.2)
        self.assertEqual(classify.derive_tier(p), "low")

    def test_encori_evidence_makes_high(self):
        p = Partner(gene="X", encori_evidence="12 CLIP experiments")
        self.assertEqual(classify.derive_tier(p), "high")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m unittest tests.test_classify -v`
Expected: FAIL / ERROR with "No module named 'classify'".

- [ ] **Step 3: Create `classify.py`**

```python
import re

from models import Partner

RIBOSOMAL_RE = re.compile(r"^(RPS|RPL|MRPS|MRPL)\d")
RIBOSOMAL_EXPLICIT = {"RPSA", "FAU", "RACK1", "UBA52", "RPLP0", "RPLP1", "RPLP2"}

STRING_HIGH = 0.7
STRING_MED = 0.4
INTACT_EXPERIMENTAL = 0.45


def classify_kind(gene):
    g = gene.upper()
    if RIBOSOMAL_RE.match(g) or g in RIBOSOMAL_EXPLICIT:
        return "ribosomal"
    return "nonribosomal"


def derive_tier(p: Partner):
    experimental = (
        (p.string_experiments or 0) > 0
        or (p.intact_mi or 0) >= INTACT_EXPERIMENTAL
        or p.uniprot_curated
        or bool(p.encori_evidence)
    )
    score = p.string_score or 0
    if experimental or score >= STRING_HIGH:
        return "high"
    if score >= STRING_MED:
        return "medium"
    return "low"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m unittest tests.test_classify -v`
Expected: PASS (all tests OK).

- [ ] **Step 5: Commit**

```bash
git add classify.py tests/test_classify.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: partner classification and confidence tiering"
```

---

### Task 3: UniProt input resolution and sequence fetch

**Files:**
- Create: `uniprot.py`, `tests/test_uniprot.py`

**Interfaces:**
- Consumes: `httpget.get`, `models.InputSeq`
- Produces:
  - `uniprot.ORGANISM = 9606`
  - `uniprot.parse_input_fasta(text: str) -> list[InputSeq]`
  - `uniprot.resolve_input_sequences(symbol: str, http=httpget.get) -> list[InputSeq]` (raises `SystemExit` if none)
  - `uniprot.fetch_sequence(accession: str, http=httpget.get) -> str`
  - `uniprot.resolve_accession(gene: str, http=httpget.get) -> str | None`
  - `uniprot.fetch_protein_by_gene(gene: str, http=httpget.get) -> tuple[str | None, str | None]`  (accession, sequence)
  - `uniprot.parse_uniprot_interactions(tsv_text: str) -> set[str]`  (interactant accessions)
  - `uniprot.fetch_uniprot_interactions(symbol: str, http=httpget.get) -> str`

- [ ] **Step 1: Write the failing test `tests/test_uniprot.py`**

```python
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import uniprot

FASTA = (
    ">sp|P62847|RS24_HUMAN 40S ribosomal protein S24 OS=Homo sapiens\n"
    "MNDTVTIRTRKFMTNRLLQRK\n"
    "QMVIDVLHPGKATVPK\n"
    ">sp|P62847-2|RS24_HUMAN Isoform 2 OS=Homo sapiens\n"
    "MNDTVTIRTRKFMT\n"
    ">tr|A0A2R8Y849|A0A2R8Y849_HUMAN Ribosomal protein S24 OS=Homo sapiens\n"
    "MNDTVTIRT\n"
)

ACC_TSV = "Entry\tReviewed\nP62841\treviewed\nA0A123\tunreviewed\n"
INTERACTION_TSV = 'Entry\tInteracts with\nP62847\tP62841; Q9Y2B4\n'


class TestParseInputFasta(unittest.TestCase):
    def setUp(self):
        self.seqs = uniprot.parse_input_fasta(FASTA)

    def test_count(self):
        self.assertEqual(len(self.seqs), 3)

    def test_canonical(self):
        c = self.seqs[0]
        self.assertEqual(c.accession, "P62847")
        self.assertEqual(c.isoform_id, "P62847")
        self.assertTrue(c.reviewed)
        self.assertEqual(c.sequence, "MNDTVTIRTRKFMTNRLLQRKQMVIDVLHPGKATVPK")

    def test_isoform_id_and_base_accession(self):
        iso = self.seqs[1]
        self.assertEqual(iso.accession, "P62847")
        self.assertEqual(iso.isoform_id, "P62847-2")

    def test_unreviewed_flag(self):
        tr = self.seqs[2]
        self.assertEqual(tr.isoform_id, "A0A2R8Y849")
        self.assertFalse(tr.reviewed)


class TestResolveInputSequences(unittest.TestCase):
    def test_uses_http_and_parses(self):
        seqs = uniprot.resolve_input_sequences("RPS24", http=lambda url: FASTA)
        self.assertEqual(len(seqs), 3)

    def test_empty_raises(self):
        with self.assertRaises(SystemExit):
            uniprot.resolve_input_sequences("NOTAGENE", http=lambda url: "")


class TestResolveAccession(unittest.TestCase):
    def test_prefers_reviewed(self):
        acc = uniprot.resolve_accession("RPS15", http=lambda url: ACC_TSV)
        self.assertEqual(acc, "P62841")

    def test_none_when_empty(self):
        self.assertIsNone(uniprot.resolve_accession("X", http=lambda url: "Entry\tReviewed\n"))


class TestFetchSequence(unittest.TestCase):
    def test_concatenates_fasta_body(self):
        single = ">sp|P62841|RS15\nAAAA\nBBBB\n"
        self.assertEqual(uniprot.fetch_sequence("P62841", http=lambda url: single), "AAAABBBB")


class TestUniprotInteractions(unittest.TestCase):
    def test_parses_interactant_accessions(self):
        accs = uniprot.parse_uniprot_interactions(INTERACTION_TSV)
        self.assertEqual(accs, {"P62841", "Q9Y2B4"})

    def test_empty_interactions(self):
        self.assertEqual(uniprot.parse_uniprot_interactions("Entry\tInteracts with\nP62847\t\n"), set())
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m unittest tests.test_uniprot -v`
Expected: FAIL / ERROR with "No module named 'uniprot'".

- [ ] **Step 3: Create `uniprot.py`**

```python
import urllib.parse

import httpget
from models import InputSeq

ORGANISM = 9606
BASE = "https://rest.uniprot.org/uniprotkb"


def _fasta_records(text):
    header, seq = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header, seq = line[1:], []
        elif line.strip():
            seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def parse_input_fasta(text):
    out = []
    for header, seq in _fasta_records(text):
        parts = header.split("|")
        db, acc = parts[0], parts[1]            # 'sp'/'tr', 'P62847' or 'P62847-2'
        out.append(InputSeq(
            accession=acc.split("-")[0],
            isoform_id=acc,
            sequence=seq,
            reviewed=(db == "sp"),
        ))
    return out


def resolve_input_sequences(symbol, http=httpget.get):
    q = f"gene:{symbol} AND organism_id:{ORGANISM}"
    url = (f"{BASE}/stream?query={urllib.parse.quote(q)}"
           "&format=fasta&includeIsoform=true")
    seqs = parse_input_fasta(http(url))
    if not seqs:
        raise SystemExit(f"No UniProt entries found for human gene '{symbol}'.")
    return seqs


def fetch_sequence(accession, http=httpget.get):
    text = http(f"{BASE}/{accession}.fasta")
    recs = list(_fasta_records(text))
    return recs[0][1] if recs else ""


def resolve_accession(gene, http=httpget.get):
    q = f"gene_exact:{gene} AND organism_id:{ORGANISM}"
    url = (f"{BASE}/search?query={urllib.parse.quote(q)}"
           "&format=tsv&fields=accession,reviewed&size=10")
    rows = [l.split("\t") for l in http(url).splitlines()[1:] if l]
    if not rows:
        return None
    for acc, reviewed in rows:
        if reviewed == "reviewed":
            return acc
    return rows[0][0]


def fetch_protein_by_gene(gene, http=httpget.get):
    acc = resolve_accession(gene, http)
    if not acc:
        return None, None
    return acc, fetch_sequence(acc, http)


def fetch_uniprot_interactions(symbol, http=httpget.get):
    q = f"gene_exact:{symbol} AND organism_id:{ORGANISM} AND reviewed:true"
    url = (f"{BASE}/search?query={urllib.parse.quote(q)}"
           "&format=tsv&fields=accession,cc_interaction&size=5")
    return http(url)


def parse_uniprot_interactions(tsv_text):
    accs = set()
    for line in tsv_text.splitlines()[1:]:
        if "\t" not in line:
            continue
        field = line.split("\t", 1)[1].strip()
        for tok in field.replace(",", ";").split(";"):
            tok = tok.strip()
            # interactant tokens look like UniProt accessions; ignore isoform suffix
            if tok and tok[0].isalpha() and tok[1:2].isdigit() is False and "-" not in tok[:1]:
                base = tok.split("-")[0]
                if base.isalnum():
                    accs.add(base)
    return accs
```

Note: `cc_interaction` values are space/semicolon-separated accession tokens (e.g. `P62841; Q9Y2B4`). The parser keeps alphanumeric accession-shaped tokens and strips isoform suffixes. The test fixture `P62841; Q9Y2B4` exercises this.

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m unittest tests.test_uniprot -v`
Expected: PASS. If `parse_uniprot_interactions` mis-tokenizes, simplify the keep-rule to: keep tokens matching `re.fullmatch(r"[A-Z0-9]{6,10}", base)`. Adjust until the two interaction tests pass.

- [ ] **Step 5: Commit**

```bash
git add uniprot.py tests/test_uniprot.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: UniProt input resolution, sequence fetch, interactions"
```

---

### Task 4: Protein partner discovery (STRING + IntAct + UniProt union)

**Files:**
- Create: `partners.py`, `tests/test_partners.py`

**Interfaces:**
- Consumes: `httpget.get`, `models.Partner`
- Produces:
  - `partners.fetch_string(symbol: str, http=httpget.get) -> str`
  - `partners.parse_string(tsv_text: str) -> dict[str, Partner]`  (keyed by UPPERCASE gene)
  - `partners.fetch_intact(symbol: str, http=httpget.get) -> str`
  - `partners.parse_intact(json_text: str, input_accessions: set[str]) -> dict[str, Partner]`
  - `partners.merge_partners(*dicts) -> dict[str, Partner]`
  - `partners.discover_protein_partners(symbol: str, input_accessions: set[str], http=httpget.get) -> dict[str, Partner]`

**Real API shapes (verified 2026-06-19):**
- STRING `interaction_partners` TSV columns include: `preferredName_A`, `preferredName_B`, `score`, `escore` (experiments), `tscore` (textmining).
- IntAct `findInteractions` JSON: `content[]` items with `moleculeA`/`moleculeB` (gene names), `uniqueIdA`/`uniqueIdB` (UniProt accessions), `intactMiscore` (float), `taxIdA`/`taxIdB` (int), `typeA`/`typeB` (`"protein"`). **`findInteractions` does substring name matching**, so filter by input accessions, not by name.

- [ ] **Step 1: Write the failing test `tests/test_partners.py`**

```python
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import partners

STRING_TSV = (
    "stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\t"
    "nscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore\n"
    "9606.ENSP1\t9606.ENSP2\tRPS24\tRPS15\t9606\t0.999\t0\t0\t0\t0.991\t0.994\t0.9\t0.727\n"
    "9606.ENSP1\t9606.ENSP3\tRPS24\tTP53\t9606\t0.45\t0\t0\t0\t0\t0\t0\t0.45\n"
)

INTACT_JSON = """
{"content": [
  {"moleculeA": "RPS24", "moleculeB": "RPS15", "uniqueIdA": "P62847",
   "uniqueIdB": "P62841", "intactMiscore": 0.7, "taxIdA": 9606, "taxIdB": 9606,
   "typeA": "protein", "typeB": "protein"},
  {"moleculeA": "RPS24", "moleculeB": "YEASTX", "uniqueIdA": "P62847",
   "uniqueIdB": "Q00001", "intactMiscore": 0.9, "taxIdA": 9606, "taxIdB": 559292,
   "typeA": "protein", "typeB": "protein"},
  {"moleculeA": "OTHER", "moleculeB": "NOISE", "uniqueIdA": "Z9", "uniqueIdB": "Z8",
   "intactMiscore": 0.9, "taxIdA": 9606, "taxIdB": 9606, "typeA": "protein", "typeB": "protein"}
]}
"""


class TestParseString(unittest.TestCase):
    def setUp(self):
        self.d = partners.parse_string(STRING_TSV)

    def test_keys_uppercase_gene(self):
        self.assertIn("RPS15", self.d)
        self.assertIn("TP53", self.d)

    def test_scores(self):
        p = self.d["RPS15"]
        self.assertEqual(p.string_score, 0.999)
        self.assertEqual(p.string_experiments, 0.994)
        self.assertEqual(p.string_textmining, 0.727)
        self.assertEqual(p.sources, ["string"])

    def test_empty(self):
        self.assertEqual(partners.parse_string(""), {})


class TestParseIntact(unittest.TestCase):
    def setUp(self):
        self.d = partners.parse_intact(INTACT_JSON, {"P62847"})

    def test_keeps_human_partner_of_input(self):
        self.assertIn("RPS15", self.d)
        self.assertEqual(self.d["RPS15"].intact_mi, 0.7)
        self.assertEqual(self.d["RPS15"].partner_id, "P62841")
        self.assertEqual(self.d["RPS15"].sources, ["intact"])

    def test_drops_non_human_partner(self):
        self.assertNotIn("YEASTX", self.d)

    def test_drops_interactions_not_involving_input(self):
        self.assertNotIn("NOISE", self.d)


class TestMerge(unittest.TestCase):
    def test_union_merges_scores_and_sources(self):
        s = partners.parse_string(STRING_TSV)
        i = partners.parse_intact(INTACT_JSON, {"P62847"})
        m = partners.merge_partners(s, i)
        rps15 = m["RPS15"]
        self.assertEqual(rps15.string_score, 0.999)
        self.assertEqual(rps15.intact_mi, 0.7)
        self.assertEqual(sorted(rps15.sources), ["intact", "string"])
        self.assertEqual(rps15.partner_id, "P62841")  # accession carried from intact
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m unittest tests.test_partners -v`
Expected: FAIL / ERROR with "No module named 'partners'".

- [ ] **Step 3: Create `partners.py`**

```python
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
        if existing is None or (mi or 0) > (existing.intact_mi or 0):
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m unittest tests.test_partners -v`
Expected: PASS (all tests OK).

- [ ] **Step 5: Commit**

```bash
git add partners.py tests/test_partners.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: protein partner discovery from STRING and IntAct"
```

---

### Task 5: RNA partner discovery (ENCORI + curated TSV)

**Files:**
- Create: `rna.py`, `tests/test_rna.py`, `tests/fixtures/rna_example.tsv`

**Interfaces:**
- Consumes: `httpget.get`, `models.Partner`
- Produces:
  - `rna.ENCORI_ASSEMBLY = "hg38"`, `rna.ENCORI_MAX = 25`, `rna.RNA_MAX_LEN = 1500`
  - `rna.fetch_encori(symbol: str, http=httpget.get) -> str`
  - `rna.parse_encori(text: str) -> dict[str, Partner]`  (keyed by UPPERCASE gene; top ENCORI_MAX by evidence)
  - `rna.load_rna_tsv(path: str) -> list[Partner]`
  - `rna.discover_rna_partners(symbol: str, rna_tsv: str | None, http=httpget.get) -> list[Partner]`

**ENCORI shape (verified 2026-06-19):** TSV after `#`-comment lines, header columns include `RBP`, `geneName` (RNA target), `geneType`, `totalClipExpNum`. Many rows per target (binding sites) → dedup by `geneName`, keep max `totalClipExpNum`. Sequences are NOT in this response; ENCORI RNA partners have no usable sequence in v1 and are recorded but produce a JSON only if a sequence is supplied via `--rna-tsv` for the same gene. Curated TSV is the reliable sequence source.

- [ ] **Step 1: Create fixture `tests/fixtures/rna_example.tsv`**

```
gene	sequence
ZFAS1	ACGUACGUACGUAAAUUUGGGCCC
RNA18SN5	TACCTGGTTGATCCTGCCAGTAGCATATGCTTG
```

- [ ] **Step 2: Write the failing test `tests/test_rna.py`**

```python
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import rna

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")

ENCORI_TEXT = (
    "#please cite ...\n"
    "RBP\tgeneID\tgeneName\tgeneType\tclusterNum\ttotalClipExpNum\ttotalClipSiteNum\n"
    "ELAVL1\tENSG1\tTSPAN6\tprotein_coding\t8\t12\t20\n"
    "ELAVL1\tENSG1\tTSPAN6\tprotein_coding\t8\t12\t21\n"
    "ELAVL1\tENSG2\tMALAT1\tlincRNA\t3\t30\t40\n"
)


class TestParseEncori(unittest.TestCase):
    def setUp(self):
        self.d = rna.parse_encori(ENCORI_TEXT)

    def test_dedups_by_gene(self):
        self.assertIn("TSPAN6", self.d)
        self.assertIn("MALAT1", self.d)
        self.assertEqual(len(self.d), 2)

    def test_kind_is_rna_and_evidence_recorded(self):
        p = self.d["MALAT1"]
        self.assertEqual(p.kind, "rna")
        self.assertEqual(p.partner_id, "MALAT1")
        self.assertIn("30", p.encori_evidence)
        self.assertEqual(p.sources, ["encori"])

    def test_empty_text(self):
        self.assertEqual(rna.parse_encori("#comment only\n"), {})


class TestLoadRnaTsv(unittest.TestCase):
    def setUp(self):
        self.parts = rna.load_rna_tsv(os.path.join(FIXTURES, "rna_example.tsv"))

    def test_loads_two(self):
        self.assertEqual(len(self.parts), 2)

    def test_fields(self):
        zfas1 = [p for p in self.parts if p.gene == "ZFAS1"][0]
        self.assertEqual(zfas1.kind, "rna")
        self.assertEqual(zfas1.partner_id, "ZFAS1")
        self.assertEqual(zfas1.sequence, "ACGUACGUACGUAAAUUUGGGCCC")
        self.assertEqual(zfas1.sources, ["curated"])

    def test_dna_normalized_to_rna(self):
        rrna = [p for p in self.parts if p.gene == "RNA18SN5"][0]
        self.assertNotIn("T", rrna.sequence)
        self.assertTrue(rrna.sequence.startswith("UACCUGG"))


class TestDiscoverRnaPartners(unittest.TestCase):
    def test_merges_encori_and_curated_curated_wins_sequence(self):
        parts = rna.discover_rna_partners(
            "ELAVL1",
            os.path.join(FIXTURES, "rna_example.tsv"),
            http=lambda url: ENCORI_TEXT,
        )
        genes = {p.gene for p in parts}
        self.assertIn("MALAT1", genes)   # from ENCORI (no sequence)
        self.assertIn("ZFAS1", genes)    # from curated (has sequence)
        malat1 = [p for p in parts if p.gene == "MALAT1"][0]
        self.assertIsNone(malat1.sequence)
```

- [ ] **Step 3: Run test to verify it fails**

Run: `python -m unittest tests.test_rna -v`
Expected: FAIL / ERROR with "No module named 'rna'".

- [ ] **Step 4: Create `rna.py`**

```python
import csv
import urllib.parse

import httpget
from af3 import to_rna
from models import Partner

ENCORI_ASSEMBLY = "hg38"
ENCORI_MAX = 25
RNA_MAX_LEN = 1500
ORGANISM = 9606


def fetch_encori(symbol, http=httpget.get):
    url = ("https://rnasysu.com/encori/api/RBPTarget/"
           f"?assembly={ENCORI_ASSEMBLY}&geneType=all&RBP={urllib.parse.quote(symbol)}"
           "&clipExpNum=1&pancancerNum=0&target=all&cellType=all")
    return http(url)


def parse_encori(text):
    lines = [l for l in text.splitlines() if l and not l.startswith("#")]
    if not lines:
        return {}
    idx = {name: i for i, name in enumerate(lines[0].split("\t"))}
    if "geneName" not in idx:
        return {}
    best = {}
    for line in lines[1:]:
        c = line.split("\t")
        gene = c[idx["geneName"]]
        try:
            evidence = int(c[idx["totalClipExpNum"]])
        except (ValueError, KeyError, IndexError):
            evidence = 0
        key = gene.upper()
        if key not in best or evidence > best[key][0]:
            best[key] = (evidence, gene)
    top = sorted(best.items(), key=lambda kv: kv[1][0], reverse=True)[:ENCORI_MAX]
    out = {}
    for key, (evidence, gene) in top:
        out[key] = Partner(gene=gene, partner_id=gene, kind="rna",
                           sources=["encori"],
                           encori_evidence=f"{evidence} CLIP experiments")
    return out


def load_rna_tsv(path):
    out = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = (row.get("gene") or "").strip()
            seq = (row.get("sequence") or "").strip()
            if not gene or not seq:
                continue
            out.append(Partner(gene=gene, partner_id=gene, kind="rna",
                               sequence=to_rna(seq), sources=["curated"]))
    return out


def discover_rna_partners(symbol, rna_tsv, http=httpget.get):
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
                # curated sequence enriches an ENCORI hit
                by_gene[key].sequence = p.sequence
                if "curated" not in by_gene[key].sources:
                    by_gene[key].sources.append("curated")
            else:
                by_gene[key] = p
    return list(by_gene.values())
```

- [ ] **Step 5: Run test to verify it passes**

Run: `python -m unittest tests.test_rna -v`
Expected: PASS (all tests OK).

- [ ] **Step 6: Commit**

```bash
git add rna.py tests/test_rna.py tests/fixtures/rna_example.tsv
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: RNA partner discovery from ENCORI and curated TSV"
```

---

### Task 6: Manifest and README writers

**Files:**
- Create: `manifest.py`, `tests/test_manifest.py`

**Interfaces:**
- Consumes: `models.InputSeq`, `models.Partner`
- Produces:
  - `manifest.MANIFEST_COLUMNS: list[str]`
  - `manifest.manifest_row(inp: InputSeq, partner: Partner, json_path: str) -> dict`
  - `manifest.write_manifest(path: str, rows: list[dict]) -> None`
  - `manifest.write_readme(path: str, symbol: str, partners: list[Partner], date: str) -> None`

- [ ] **Step 1: Write the failing test `tests/test_manifest.py`**

```python
import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import manifest
from models import InputSeq, Partner


class TestManifestRow(unittest.TestCase):
    def test_row_has_all_columns(self):
        inp = InputSeq("P62847", "P62847-2", "MND", False)
        p = Partner(gene="RPS15", partner_id="P62841", kind="ribosomal",
                    sources=["string", "intact"], string_score=0.99,
                    string_experiments=0.9, intact_mi=0.7, tier="high")
        row = manifest.manifest_row(inp, p, "AF3_inputs/ribosomal_protein/high/P62847-2_P62841.json")
        for col in manifest.MANIFEST_COLUMNS:
            self.assertIn(col, row)
        self.assertEqual(row["input_isoform"], "P62847-2")
        self.assertEqual(row["partner_gene"], "RPS15")
        self.assertEqual(row["partner_kind"], "ribosomal")
        self.assertEqual(row["sources"], "string;intact")
        self.assertEqual(row["tier"], "high")
        self.assertEqual(row["input_reviewed"], "False")


class TestWriteManifest(unittest.TestCase):
    def test_roundtrip(self):
        inp = InputSeq("P62847", "P62847", "MND", True)
        p = Partner(gene="RPS15", partner_id="P62841", kind="ribosomal", tier="high")
        rows = [manifest.manifest_row(inp, p, "x.json")]
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "manifest.tsv")
            manifest.write_manifest(path, rows)
            with open(path, newline="") as f:
                read = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(read[0]["partner_gene"], "RPS15")
        self.assertEqual(list(read[0].keys()), manifest.MANIFEST_COLUMNS)


class TestWriteReadme(unittest.TestCase):
    def test_mentions_symbol_and_counts(self):
        partners = [Partner(gene="RPS15", kind="ribosomal", tier="high"),
                    Partner(gene="ZFAS1", kind="rna", tier="high")]
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "README.txt")
            manifest.write_readme(path, "RPS24", partners, "2026-06-19")
            text = open(path).read()
        self.assertIn("RPS24", text)
        self.assertIn("ribosomal", text)
        self.assertIn("2026-06-19", text)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m unittest tests.test_manifest -v`
Expected: FAIL / ERROR with "No module named 'manifest'".

- [ ] **Step 3: Create `manifest.py`**

```python
import csv
from collections import Counter

from models import InputSeq, Partner

MANIFEST_COLUMNS = [
    "input_gene", "input_accession", "input_isoform", "input_reviewed",
    "partner_id", "partner_gene", "partner_kind", "sources",
    "string_score", "string_experiments", "string_textmining",
    "intact_mi", "uniprot_curated", "encori_evidence", "tier", "json_path",
]


def manifest_row(inp: InputSeq, partner: Partner, json_path):
    return {
        "input_gene": inp.accession,
        "input_accession": inp.accession,
        "input_isoform": inp.isoform_id,
        "input_reviewed": str(inp.reviewed),
        "partner_id": partner.partner_id,
        "partner_gene": partner.gene,
        "partner_kind": partner.kind,
        "sources": ";".join(partner.sources),
        "string_score": "" if partner.string_score is None else partner.string_score,
        "string_experiments": "" if partner.string_experiments is None else partner.string_experiments,
        "string_textmining": "" if partner.string_textmining is None else partner.string_textmining,
        "intact_mi": "" if partner.intact_mi is None else partner.intact_mi,
        "uniprot_curated": str(partner.uniprot_curated),
        "encori_evidence": partner.encori_evidence or "",
        "tier": partner.tier,
        "json_path": json_path,
    }


def write_manifest(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_readme(path, symbol, partners, date):
    by_kind = Counter(p.kind for p in partners)
    by_tier = Counter(p.tier for p in partners)
    lines = [
        f"AF3 partner inputs for {symbol}",
        f"Generated: {date}",
        "",
        "Sources queried: UniProt, STRING, IntAct, ENCORI (+ curated --rna-tsv if supplied).",
        "All external references use their latest release.",
        "",
        f"Total partners: {len(partners)}",
        f"  By kind: " + ", ".join(f"{k}={v}" for k, v in sorted(by_kind.items())),
        f"  By tier: " + ", ".join(f"{k}={v}" for k, v in sorted(by_tier.items())),
        "",
        "Each JSON pairs one input-protein isoform with one partner (AF3-local format).",
        "Files are organized by partner kind and confidence tier; see manifest.tsv for scores.",
    ]
    if by_kind.get("rna", 0) == 0:
        lines.append("")
        lines.append("Note: no RNA partners — gene is not a catalogued RBP in ENCORI and "
                     "no --rna-tsv was supplied (or supplied rows were empty).")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m unittest tests.test_manifest -v`
Expected: PASS (all tests OK).

- [ ] **Step 5: Commit**

```bash
git add manifest.py tests/test_manifest.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: manifest and README writers"
```

---

### Task 7: Orchestration, zip assembly, and CLI

**Files:**
- Create: `make_inputs.py`, `tests/test_build.py`

**Interfaces:**
- Consumes: all prior modules.
- Produces:
  - `make_inputs.build(symbol, out_dir, rna_tsv=None, http=httpget.get) -> str`  (returns zip path)
  - `make_inputs.main(argv=None) -> int`

**Orchestration order (in `build`):**
1. `resolve_input_sequences(symbol)` → input seqs.
2. `discover_protein_partners(symbol, {s.accession for inputs})` → protein partners (dict).
3. For each protein partner: `classify_kind(gene)`; resolve sequence+accession via `fetch_protein_by_gene` (skip partner if no sequence).
4. `discover_rna_partners(symbol, rna_tsv)` → rna partners (skip those with no sequence — ENCORI-only hits).
5. Mark `uniprot_curated` on protein partners whose accession ∈ `parse_uniprot_interactions(fetch_uniprot_interactions(symbol))`.
6. `derive_tier(p)` for every partner.
7. For each (input seq × partner with a sequence): build AF3 JSON, choose path `AF3_inputs/<kind>_protein|rna/<tier>/<isoform>_<partner_id>.json`, write into a temp tree; collect manifest rows.
8. Write `manifest.tsv`, `README.txt`; zip the `SYMBOL/` tree to `<out_dir>/SYMBOL.zip`.

- [ ] **Step 1: Write the failing integration test `tests/test_build.py`**

```python
import json
import os
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import make_inputs

FASTA = (
    ">sp|P62847|RS24_HUMAN OS=Homo sapiens\nMNDTVTIR\n"
    ">sp|P62847-2|RS24_HUMAN Isoform 2\nMNDTVT\n"
)
STRING_TSV = (
    "stringId_A\tstringId_B\tpreferredName_A\tpreferredName_B\tncbiTaxonId\tscore\t"
    "nscore\tfscore\tpscore\tascore\tescore\tdscore\ttscore\n"
    "9606.E1\t9606.E2\tRPS24\tRPS15\t9606\t0.999\t0\t0\t0\t0.9\t0.99\t0.9\t0.7\n"
)
INTACT_JSON = '{"content": []}'
ENCORI_TEXT = "#cite\nRBP\tgeneName\tgeneType\ttotalClipExpNum\n"  # no RNA hits
ACC_TSV = "Entry\tReviewed\nP62841\treviewed\n"
PARTNER_FASTA = ">sp|P62841|RS15\nMAAAA\n"
UNIPROT_INT_TSV = "Entry\tInteracts with\nP62847\t\n"


def fake_http(url):
    if "/stream?" in url:
        return FASTA
    if "string-db.org" in url:
        return STRING_TSV
    if "intact/ws" in url:
        return INTACT_JSON
    if "encori" in url:
        return ENCORI_TEXT
    if "cc_interaction" in url or "fields=accession,cc_interaction" in url:
        return UNIPROT_INT_TSV
    if "/search?" in url and "accession,reviewed" in url:
        return ACC_TSV
    if url.endswith("P62841.fasta"):
        return PARTNER_FASTA
    raise AssertionError(f"unexpected url: {url}")


class TestBuild(unittest.TestCase):
    def test_produces_zip_with_json_and_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            zip_path = make_inputs.build("RPS24", d, rna_tsv=None, http=fake_http)
            self.assertTrue(zip_path.endswith("RPS24.zip"))
            with zipfile.ZipFile(zip_path) as z:
                names = z.namelist()
                jsons = [n for n in names if n.endswith(".json")]
                # 2 isoforms x 1 ribosomal partner (RPS15) = 2 JSONs
                self.assertEqual(len(jsons), 2)
                self.assertTrue(any("ribosomal_protein/high/" in n for n in jsons))
                self.assertTrue(any(n.endswith("manifest.tsv") for n in names))
                self.assertTrue(any(n.endswith("README.txt") for n in names))
                sample = json.loads(z.read(jsons[0]))
                self.assertEqual(sample[0]["dialect"], "alphafold3")
                self.assertEqual(len(sample[0]["sequences"]), 2)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m unittest tests.test_build -v`
Expected: FAIL / ERROR with "No module named 'make_inputs'".

- [ ] **Step 3: Create `make_inputs.py`**

```python
import argparse
import json
import os
import shutil
import sys
import tempfile
import zipfile
from datetime import date

import httpget
import uniprot
import partners as partners_mod
import rna as rna_mod
import manifest as manifest_mod
from af3 import af3_protein_pair, af3_rna_pair
from classify import classify_kind, derive_tier

KIND_DIR = {"ribosomal": "ribosomal_protein", "nonribosomal": "nonribosomal_protein", "rna": "rna"}


def _collect_partners(symbol, input_seqs, rna_tsv, http):
    input_accs = {s.accession for s in input_seqs}
    result = []

    proteins = partners_mod.discover_protein_partners(symbol, input_accs, http)
    for p in proteins.values():
        p.kind = classify_kind(p.gene)
        acc, seq = uniprot.fetch_protein_by_gene(p.gene, http)
        if not seq:
            print(f"  WARN: no sequence for {p.gene}, skipping")
            continue
        p.partner_id = acc
        p.sequence = seq
        result.append(p)

    for p in rna_mod.discover_rna_partners(symbol, rna_tsv, http):
        if not p.sequence:
            print(f"  WARN: no sequence for RNA {p.gene} (ENCORI-only), skipping")
            continue
        result.append(p)

    try:
        curated = uniprot.parse_uniprot_interactions(
            uniprot.fetch_uniprot_interactions(symbol, http))
    except Exception as e:
        print(f"  WARN: UniProt interactions failed: {e}")
        curated = set()
    for p in result:
        if p.partner_id in curated:
            p.uniprot_curated = True
            if "uniprot" not in p.sources:
                p.sources.append("uniprot")

    for p in result:
        p.tier = derive_tier(p)
    return result


def build(symbol, out_dir, rna_tsv=None, http=httpget.get):
    symbol = symbol.upper()
    input_seqs = uniprot.resolve_input_sequences(symbol, http)
    partners = _collect_partners(symbol, input_seqs, rna_tsv, http)

    work = tempfile.mkdtemp()
    try:
        root = os.path.join(work, symbol)
        rows = []
        for inp in input_seqs:
            for p in partners:
                name = f"{inp.isoform_id}_{p.partner_id}"
                rel = os.path.join("AF3_inputs", KIND_DIR[p.kind], p.tier, f"{name}.json")
                abspath = os.path.join(root, rel)
                os.makedirs(os.path.dirname(abspath), exist_ok=True)
                job = (af3_rna_pair if p.kind == "rna" else af3_protein_pair)(name, inp, p)
                with open(abspath, "w") as f:
                    json.dump([job], f, indent=2)
                rows.append(manifest_mod.manifest_row(inp, p, rel.replace(os.sep, "/")))

        os.makedirs(root, exist_ok=True)
        manifest_mod.write_manifest(os.path.join(root, "manifest.tsv"), rows)
        manifest_mod.write_readme(os.path.join(root, "README.txt"),
                                  symbol, partners, date.today().isoformat())

        os.makedirs(out_dir, exist_ok=True)
        zip_path = os.path.join(out_dir, f"{symbol}.zip")
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
            for dirpath, _, files in os.walk(root):
                for fn in files:
                    full = os.path.join(dirpath, fn)
                    z.write(full, os.path.relpath(full, work))
        print(f"Wrote {zip_path}: {len(rows)} JSONs, {len(partners)} partners")
        return zip_path
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate AF3 input JSONs for a gene's interaction partners.")
    parser.add_argument("symbol", help="Human gene symbol, e.g. RPS24")
    parser.add_argument("--out", default=".", help="Output directory for SYMBOL.zip")
    parser.add_argument("--rna-tsv", default=None, help="Optional curated RNA partners TSV (gene, sequence)")
    args = parser.parse_args(argv)
    build(args.symbol, args.out, args.rna_tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m unittest tests.test_build -v`
Expected: PASS.

- [ ] **Step 5: Run the full suite**

Run: `python -m unittest discover -s tests -t . -v`
Expected: all tests across the 6 test files PASS.

- [ ] **Step 6: Commit**

```bash
git add make_inputs.py tests/test_build.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "feat: orchestration, zip assembly, and CLI"
```

---

### Task 8: README, project CLAUDE.md, and opt-in live smoke test

**Files:**
- Create: `README.md`, `CLAUDE.md`, `tests/test_live_smoke.py`

**Interfaces:** none new.

- [ ] **Step 1: Create opt-in live smoke test `tests/test_live_smoke.py`**

```python
import os
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import make_inputs


@unittest.skipUnless(os.environ.get("AF3_LIVE") == "1",
                     "set AF3_LIVE=1 to run live network smoke test")
class TestLiveSmoke(unittest.TestCase):
    def test_rps24_real(self):
        with tempfile.TemporaryDirectory() as d:
            zip_path = make_inputs.build("RPS24", d)
            with zipfile.ZipFile(zip_path) as z:
                jsons = [n for n in z.namelist() if n.endswith(".json")]
                self.assertGreater(len(jsons), 0)
                # RPS24 should find ribosomal partners (e.g. RPS15) via STRING
                self.assertTrue(any("ribosomal_protein/" in n for n in jsons))
```

- [ ] **Step 2: Run the live smoke test**

Run: `AF3_LIVE=1 python -m unittest tests.test_live_smoke -v`
Expected: PASS (real network). If a transient API error occurs, re-run; graceful degradation means partial sources still produce a zip, but RPS24 via STRING should yield ribosomal partners.

- [ ] **Step 3: Create `README.md`**

```markdown
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
```

- [ ] **Step 4: Create `CLAUDE.md`** (project guidance, mirroring the workspace style)

```markdown
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
- `rna.py` — ENCORI + curated-TSV RNA partners.
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
- Tests: stdlib `unittest`, fixtures + injected fake HTTP; live smoke gated by `AF3_LIVE=1`.
```

- [ ] **Step 5: Commit**

```bash
git add README.md CLAUDE.md tests/test_live_smoke.py
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit -m "docs: README, CLAUDE.md, and opt-in live smoke test"
```

---

## Self-Review

**Spec coverage:**
- Gene → zip of AF3-local JSONs → Tasks 1, 7. ✓
- All isoforms incl. computational → Task 3 (`includeIsoform=true`, reviewed flag). ✓
- STRING + IntAct + UniProt protein discovery → Task 4 + Task 7 step 5 (curated marking). ✓
- ENCORI + curated RNA hybrid → Task 5. ✓
- Ribosomal/non-ribosomal/RNA classification → Task 2 + Task 7. ✓
- All partners included, marked by confidence → tiers in Task 2, folders + manifest in Tasks 6/7. ✓
- Manifest + tiered folders + README → Tasks 6, 7. ✓
- Latest-release constants → Global Constraints, used in Tasks 1/3/4/5. ✓
- Human only → Global Constraints. ✓
- Graceful degradation on source failure → try/except in Tasks 4, 5, 7. ✓
- Zero partners still produces a valid zip → Task 7 (`os.makedirs(root)` before writing manifest even if no rows). ✓
- Unit tests w/o live network + opt-in live smoke → all tasks + Task 8. ✓

**Placeholder scan:** No TBD/TODO; every code step shows complete code; expected outputs given for each run. ✓

**Type consistency:** `Partner`/`InputSeq` field names match across Tasks 1–7; `partner_id` carries accession (protein) or gene (rna); `discover_protein_partners` returns a dict, `discover_rna_partners` returns a list, both consumed correctly in Task 7. `manifest_row(inp, partner, json_path)` signature matches its call. ✓

**Known approximation flagged for the implementer:** `parse_uniprot_interactions` token rule (Task 3) is the one spot likely to need a small tweak; Task 3 Step 4 gives the fallback regex `re.fullmatch(r"[A-Z0-9]{6,10}", base)` if the first rule mis-tokenizes.
