# Package af3_partners as an Installable Package — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the flat-layout `af3_partners` tool into an installable Python package (`af3partners`) that colleagues clone and `pip install`, with an `af3partners` CLI command and importable API.

**Architecture:** Move all nine root modules into a single `af3partners/` package directory, convert cross-module imports to package-relative, add `__init__.py`/`__main__.py`, and add a hatchling `pyproject.toml` with a console-script entry point. No behavior changes. Distribution is clone-and-install; no PyPI.

**Tech Stack:** Python ≥3.9, hatchling build backend, stdlib only (no runtime dependencies), stdlib `unittest`.

## Global Constraints

- Distribution name / import name / CLI command: all `af3partners` (one word).
- Build backend: hatchling.
- Layout: flat package dir `af3partners/` (NOT `src/`).
- `requires-python = ">=3.9"`.
- Runtime dependencies: none (pure stdlib).
- License: MIT, `Copyright (c) 2026 Brian Ryu`.
- Version: `0.1.0`, stated statically in both `pyproject.toml` and `af3partners/__init__.py`.
- No changes to tool behavior, data sources, or output format — imports and packaging only.
- Commit with inline git identity (this environment has none configured):
  `git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' commit ...`
- Run all commands from the repo root `/home/brianryu87/projects/af3_partners`.

---

### Task 1: Create the `af3partners` package

Move the nine modules into a package directory, convert their internal imports to
package-relative, add `__init__.py` and `__main__.py`, and update the test
imports. The existing `unittest` suite is the regression test for this refactor:
it is green before and must be green after.

**Files:**
- Move (via `git mv`): `af3.py`, `classify.py`, `httpget.py`, `make_inputs.py`, `manifest.py`, `models.py`, `partners.py`, `rna.py`, `uniprot.py` → `af3partners/`
- Create: `af3partners/__init__.py`, `af3partners/__main__.py`
- Modify (relative imports): `af3partners/{af3,classify,make_inputs,manifest,partners,rna,uniprot}.py`
- Modify (test imports): `tests/{test_af3,test_build,test_classify,test_live_smoke,test_manifest,test_partners,test_rna,test_uniprot}.py`

**Interfaces:**
- Produces: package `af3partners` importable from repo root; `af3partners.make_inputs:main` (returns `int`) and `af3partners.make_inputs:build(symbol, out_dir, rna_tsv=None, http=httpget.get)`; `af3partners.__version__` (str); `af3partners.build` re-export.

- [ ] **Step 1: Confirm the baseline test suite is green**

Run: `python -m unittest discover -s tests -t . -v`
Expected: all tests PASS (the `AF3_LIVE` smoke test is skipped). Record this as the green baseline.

- [ ] **Step 2: Move the nine modules into the package directory**

```bash
cd /home/brianryu87/projects/af3_partners
mkdir af3partners
git mv af3.py classify.py httpget.py make_inputs.py manifest.py models.py partners.py rna.py uniprot.py af3partners/
```

- [ ] **Step 3: Rewrite internal imports to package-relative**

In `af3partners/af3.py`:
`from models import InputSeq, Partner` → `from .models import InputSeq, Partner`

In `af3partners/classify.py`:
`from models import Partner` → `from .models import Partner`

In `af3partners/manifest.py`:
`from models import InputSeq, Partner` → `from .models import InputSeq, Partner`

In `af3partners/partners.py`:
`import httpget` → `from . import httpget`
`from models import Partner` → `from .models import Partner`

In `af3partners/uniprot.py`:
`import httpget` → `from . import httpget`
`from models import InputSeq` → `from .models import InputSeq`

In `af3partners/rna.py`:
`import httpget` → `from . import httpget`
`from af3 import to_rna` → `from .af3 import to_rna`
`from models import Partner` → `from .models import Partner`

In `af3partners/make_inputs.py`:
`import httpget` → `from . import httpget`
`import uniprot` → `from . import uniprot`
`import partners as partners_mod` → `from . import partners as partners_mod`
`import rna as rna_mod` → `from . import rna as rna_mod`
`import manifest as manifest_mod` → `from . import manifest as manifest_mod`
`from af3 import af3_protein_pair, af3_rna_pair, AF3_MAX_TOKENS` → `from .af3 import af3_protein_pair, af3_rna_pair, AF3_MAX_TOKENS`
`from classify import classify_kind, derive_tier` → `from .classify import classify_kind, derive_tier`

(`af3partners/models.py` and `af3partners/httpget.py` have no internal imports — leave them.)

- [ ] **Step 4: Create `af3partners/__init__.py`**

```python
from .make_inputs import build

__version__ = "0.1.0"
__all__ = ["build"]
```

- [ ] **Step 5: Create `af3partners/__main__.py`**

```python
import sys

from .make_inputs import main

sys.exit(main())
```

- [ ] **Step 6: Update test imports**

In each test file, leave the `sys.path.insert(0, ...)` line as-is (it points at
the repo root, where `af3partners/` now lives) and change only the module imports:

`tests/test_af3.py`:
`import af3` → `from af3partners import af3`
`from models import InputSeq, Partner` → `from af3partners.models import InputSeq, Partner`

`tests/test_classify.py`:
`import classify` → `from af3partners import classify`
`from models import Partner` → `from af3partners.models import Partner`

`tests/test_manifest.py`:
`import manifest` → `from af3partners import manifest`
`from models import InputSeq, Partner` → `from af3partners.models import InputSeq, Partner`

`tests/test_partners.py`:
`import partners` → `from af3partners import partners`

`tests/test_rna.py`:
`import rna` → `from af3partners import rna`

`tests/test_uniprot.py`:
`import uniprot` → `from af3partners import uniprot`

`tests/test_build.py`:
`import make_inputs` → `from af3partners import make_inputs`

`tests/test_live_smoke.py`:
`import make_inputs` → `from af3partners import make_inputs`

- [ ] **Step 7: Run the test suite — verify still green**

Run: `python -m unittest discover -s tests -t . -v`
Expected: same passing set as the Step 1 baseline (live smoke skipped). Any `ImportError` here means a rewrite was missed.

- [ ] **Step 8: Verify the module entry point works**

Run: `python -m af3partners --help`
Expected: argparse usage text for `af3partners` (positional `symbol`, `--out`, `--rna-tsv`), exit code 0.

- [ ] **Step 9: Commit**

```bash
git add -A
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' \
  commit -m "refactor: move modules into af3partners package with relative imports"
```

---

### Task 2: Add build configuration and install

Add the package metadata so it installs, plus the LICENSE referenced by it, and
ignore build artifacts. Verify a real editable install in a throwaway venv (so the
user's conda env is untouched).

**Files:**
- Create: `pyproject.toml`, `LICENSE`
- Modify: `.gitignore`

**Interfaces:**
- Consumes: `af3partners` package and `af3partners.make_inputs:main` from Task 1.
- Produces: console script `af3partners`; installable distribution `af3partners` 0.1.0.

- [ ] **Step 1: Create `LICENSE` (MIT)**

```
MIT License

Copyright (c) 2026 Brian Ryu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

- [ ] **Step 2: Create `pyproject.toml`**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "af3partners"
version = "0.1.0"
description = "Generate AlphaFold3 input JSONs for a human gene's interaction partners"
readme = "README.md"
requires-python = ">=3.9"
license = "MIT"
license-files = ["LICENSE"]
authors = [{ name = "Brian Ryu", email = "brianryu87@gmail.com" }]
dependencies = []

[project.scripts]
af3partners = "af3partners.make_inputs:main"

[tool.hatch.build.targets.wheel]
packages = ["af3partners"]
```

(If the installed hatchling is too old to accept the SPDX string `license = "MIT"`,
fall back to `license = { text = "MIT" }` and drop the `license-files` line. The
Step 4 install will reveal this.)

- [ ] **Step 3: Update `.gitignore`**

Append these lines to the existing `.gitignore`:

```
build/
dist/
*.egg-info/
```

- [ ] **Step 4: Verify an editable install in a throwaway venv**

```bash
rm -rf /tmp/af3pkgtest
python -m venv /tmp/af3pkgtest
/tmp/af3pkgtest/bin/pip install -e . -q
/tmp/af3pkgtest/bin/af3partners --help
/tmp/af3pkgtest/bin/python -c "from af3partners import build, __version__; print(__version__)"
```
Expected: install succeeds with no error; `af3partners --help` prints usage (exit 0); the last command prints `0.1.0`.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml LICENSE .gitignore
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' \
  commit -m "build: add hatchling pyproject, MIT license, ignore build artifacts"
```

---

### Task 3: Update documentation

Point README and CLAUDE.md at the new install flow and invocation, replacing the
old `python make_inputs.py` script-style usage.

**Files:**
- Modify: `README.md`, `CLAUDE.md`

**Interfaces:**
- Consumes: install/usage behavior established in Tasks 1–2.

- [ ] **Step 1: Update `README.md` usage**

Replace the `## Usage` body (currently `python make_inputs.py RPS24 --out .` etc.)
with an install step followed by the command invocation:

```markdown
## Install

    git clone <repo-url>
    cd af3_partners
    pip install .          # or: pip install -e .   (for development)

## Usage

    af3partners RPS24 --out .
    af3partners RPS24 --rna-tsv my_rna_partners.tsv

Without installing, run it from the repo root with:

    python -m af3partners RPS24
```

Leave the rest of the README (partner sources, job-size limit, tiers, tests) unchanged.

- [ ] **Step 2: Update the Tests section of `README.md`**

The test command is already run from the repo root and is unchanged:
`python -m unittest discover -s tests -t .`. Confirm it still reads correctly; no edit needed unless wording references the old flat layout.

- [ ] **Step 3: Update `CLAUDE.md`**

- In the `## Run` section, replace `python make_inputs.py SYMBOL [...]` with:

```
    af3partners SYMBOL [--out DIR] [--rna-tsv FILE]    # after pip install
    python -m af3partners SYMBOL [--out DIR] [--rna-tsv FILE]   # from repo root
```

- In the `## Modules` section, prefix each module with `af3partners/`
  (e.g. `af3partners/httpget.py`, `af3partners/make_inputs.py`, …), since the
  modules now live under the package directory.
- Add a one-line note that the package is installed with `pip install .` (hatchling, pure stdlib, no dependencies).

- [ ] **Step 4: Verify no stale script-style references remain**

Run: `grep -rn "python make_inputs.py" README.md CLAUDE.md`
Expected: no matches.

- [ ] **Step 5: Commit**

```bash
git add README.md CLAUDE.md
git -c user.name='Brian Ryu' -c user.email='brianryu87@gmail.com' \
  commit -m "docs: document pip install and af3partners command"
```

---

## Self-Review

**Spec coverage:**
- Target structure (package dir, `__init__`, `__main__`, moved modules) → Task 1.
- Import rewrites (7 source modules) → Task 1, Step 3; behavioral consequence (no direct script run) documented in README/CLAUDE → Task 3.
- `__init__.py` exposes `build` + `__version__` → Task 1, Step 4.
- `__main__.py` → Task 1, Step 5.
- pyproject.toml (hatchling, entry point, requires-python, no deps) → Task 2, Step 2.
- LICENSE (MIT) → Task 2, Step 1.
- Tests updated, both run modes work → Task 1, Steps 6–8.
- `.gitignore` build artifacts → Task 2, Step 3.
- README + CLAUDE.md → Task 3.
- All five verification criteria covered: install (T2 S4), `af3partners --help` + `python -m af3partners --help` (T2 S4, T1 S8), unittest suite (T1 S7), `from af3partners import build` (T2 S4).

**Placeholder scan:** none — every step has exact paths, exact import edits, and concrete commands with expected output.

**Type consistency:** `main()`→`int` and `build(symbol, out_dir, rna_tsv=None, http=...)` signatures are reused consistently between the Task 1 interface block, `__init__.py`, `__main__.py`, and the pyproject entry point `af3partners.make_inputs:main`.

## Post-implementation follow-up (not a code task)

After the branch is merged, update the `af3-partners-tool` memory: invocation is
now `af3partners SYMBOL` / `python -m af3partners SYMBOL`; the old
`python3 .../make_inputs.py` absolute-path script run no longer works because the
modules use package-relative imports.
