# Design: Package af3_partners as an installable Python package

Date: 2026-06-20

## Goal

Let colleagues download the tool and use it. They will `git clone` the repo and
`pip install .` (or `-e .`) into their own environment, getting an importable
`af3partners` package and an `af3partners` command-line entry point. No PyPI
publishing.

## Decisions (from brainstorming)

- **Distribution:** clone-and-install. No public index, no PyPI release process.
- **Packaging depth:** a real installable package (pyproject.toml + importable
  package + console entry point), not just scripts.
- **Layout:** Approach A — flat package directory `af3partners/` (not `src/`).
- **Build backend:** hatchling.
- **Names:** distribution `af3partners`, import `af3partners`, CLI `af3partners`
  (one word everywhere).
- **License:** MIT, `Copyright (c) 2026 Brian Ryu`.
- **Python floor:** `>=3.9`.

## Why this is needed

The current layout is flat: all nine modules sit at the repo root and import each
other by bare top-level name (`import httpget`, `from models import Partner`).
Installed as-is, those generic names (`models`, `httpget`, `rna`, `classify`,
`manifest`, …) would occupy the global import namespace and collide with other
packages. Packaging therefore requires moving the modules under a single
`af3partners/` namespace and converting the cross-imports to package-relative.

## Target structure

```
af3_partners/                  # repo root (directory name on disk unchanged)
├── pyproject.toml             # NEW
├── LICENSE                    # NEW (MIT)
├── README.md                  # updated
├── CLAUDE.md                  # updated
├── .gitignore                 # updated (build/, dist/, *.egg-info/)
├── af3partners/               # NEW package directory
│   ├── __init__.py            # NEW
│   ├── __main__.py            # NEW
│   ├── af3.py
│   ├── classify.py
│   ├── httpget.py
│   ├── make_inputs.py
│   ├── manifest.py
│   ├── models.py
│   ├── partners.py
│   ├── rna.py
│   └── uniprot.py
├── tests/                     # imports updated to af3partners.*
└── docs/                      # unchanged
```

## Component details

### 1. Module move + import rewrites

Move all nine `.py` modules into `af3partners/`. Convert every internal import to
package-relative:

- `import httpget` → `from . import httpget`
- `import uniprot` → `from . import uniprot`
- `import partners as partners_mod` → `from . import partners as partners_mod`
- `import rna as rna_mod` → `from . import rna as rna_mod`
- `import manifest as manifest_mod` → `from . import manifest as manifest_mod`
- `from models import Partner` → `from .models import Partner`
- `from af3 import ...` → `from .af3 import ...`
- `from classify import ...` → `from .classify import ...`

Modules with internal imports to fix: `make_inputs.py`, `partners.py`,
`uniprot.py`, `rna.py`, `classify.py`, `af3.py`, `manifest.py`.
`models.py` and `httpget.py` have no internal imports.

No logic changes — imports only.

### 2. `__init__.py`

Exposes the library API and version:

```python
from .make_inputs import build

__version__ = "0.1.0"
__all__ = ["build"]
```

### 3. `__main__.py`

Enables `python -m af3partners`:

```python
import sys

from .make_inputs import main

sys.exit(main())
```

`make_inputs.py`'s existing `main()` (returns an int) and its
`if __name__ == "__main__"` block are kept unchanged.

### 4. pyproject.toml

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
authors = [{ name = "Brian Ryu", email = "brianryu87@gmail.com" }]
dependencies = []

[project.scripts]
af3partners = "af3partners.make_inputs:main"

[tool.hatch.build.targets.wheel]
packages = ["af3partners"]
```

### 5. LICENSE

Standard MIT license text, `Copyright (c) 2026 Brian Ryu`.

### 6. Tests

Update imports in each test module from bare names to `af3partners.*`
(e.g. `from rna import parse_encori` → `from af3partners.rna import parse_encori`).
No test logic changes. Tests resolve both ways:

- from the repo root without installing (the `af3partners/` directory is
  importable from cwd), and
- after `pip install -e .`.

Live smoke test stays gated behind `AF3_LIVE=1`.

### 7. Docs

- **README.md:** add an Install section (`git clone … && cd af3_partners &&
  pip install .` — or `-e .` for development) and switch usage examples from
  `python make_inputs.py RPS24` to `af3partners RPS24`.
- **CLAUDE.md:** update the module list to the `af3partners/` paths and the run
  command to `af3partners SYMBOL` / `python -m af3partners SYMBOL`.

### 8. .gitignore

Add `build/`, `dist/`, `*.egg-info/`.

## Behavioral consequence

With package-relative imports, running the module file directly
(`python af3partners/make_inputs.py`) no longer works — Python rejects relative
imports in a top-level script. Supported invocations become:

- `af3partners SYMBOL` (after install), or
- `python -m af3partners SYMBOL` (from repo root, no install).

This supersedes the current absolute-path script-invocation habit recorded in the
`af3-partners-tool` memory; that memory will be updated after implementation.

## Out of scope

- PyPI publishing / release automation.
- `src/` layout.
- Single-sourcing the version (static `0.1.0` in both pyproject and `__init__`
  is acceptable for now).
- Any change to tool behavior, data sources, or output format.

## Verification criteria

1. `pip install -e .` succeeds in a clean environment.
2. `af3partners --help` works (console entry point resolves).
3. `python -m af3partners --help` works.
4. `python -m unittest discover -s tests -t .` — all existing tests pass.
5. `python -c "from af3partners import build, __version__"` succeeds.
