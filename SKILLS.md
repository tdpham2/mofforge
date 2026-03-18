# mofforge - AI Agent Skills

## Project Overview

mofforge is a Python find-and-replace tool for atomistic crystal structure models, especially Metal-Organic Frameworks (MOFs). It features VF2 subgraph isomorphism for pattern matching, SVD-based Procrustes alignment for fragment replacement, SMARTS-like pattern matching, multi-step pipelines, batch processing, structure validation, and provenance tracking.

**Domain**: Computational chemistry / materials science / crystallography.

## Architecture

```
src/mofforge/
├── core/              # Core data structures
│   ├── crystal.py     # Crystal class (pymatgen Structure + NetworkX Graph)
│   ├── bonding.py     # Bond inference via covalent radii
│   └── moiety.py      # Fragment loading, anchor atom handling
├── search/            # Pattern matching
│   ├── search.py      # MatchResult API, result grouping
│   └── isomorphism.py # VF2 subgraph isomorphism (NetworkX wrapper)
├── replace/           # Fragment replacement
│   ├── replace.py     # Main replacement pipeline (5 modes)
│   ├── alignment.py   # SVD Procrustes alignment
│   └── conglomerate.py # Periodic boundary reassembly (reassemble())
├── io/                # File I/O
│   ├── cif.py         # CIF read/write (pymatgen)
│   └── xyz.py         # XYZ read/write with anchor tags
├── utils/             # Utilities
│   ├── config.py      # Global config, covalent/vdW radii tables
│   └── periodic.py    # PBC utilities (wrap, min image)
├── smarts.py          # SMARTS-like string pattern matching
├── pipeline.py        # Multi-step Pipeline API (fluent interface)
├── validation.py      # Structure validation (clashes, bonds, coordination)
├── batch.py           # YAML-driven batch processing
├── provenance.py      # Provenance/audit trail tracking
└── cli.py             # CLI via Click (search, replace, remove, validate, batch)
```

## Key Concepts

- **Crystal**: Central data structure wrapping a pymatgen `Structure` (lattice + atoms + fractional coords) with a NetworkX `Graph` (bonding network). Supports indexing, slicing, addition, CIF/XYZ I/O.
- **Fragment**: A molecular fragment loaded from XYZ files via `fragment()`. Atoms tagged with `!` (e.g., `H!`, `C!`) are anchor sites used for alignment and removal during replacement.
- **MatchResult**: VF2 subgraph isomorphism via `find_pattern()` finds query fragments in a parent crystal's bond graph. Results are grouped by location (unique atom sets) and orientation (different mappings at same location).
- **Replace**: SVD Procrustes alignment via `replace_pattern()` positions replacement fragments onto matched locations. Handles bond reconstruction, periodic boundary crossings, and coordinate wrapping.
- **Anchor atom convention**: `!`-tagged atoms in XYZ files mark anchor sites. During search, tags are stripped so `H!` matches `H`. During replace, tagged atoms define where the fragment connects to the parent.

## Primary API Names

| Name | Purpose |
|---|---|
| `fragment()` | Load fragment from XYZ |
| `find_pattern()` | Pattern matching |
| `replace_pattern()` | Replace matched patterns |
| `swap()` | One-step search+replace |
| `MatchResult` | Search result class |
| `reassemble()` | PB reassembly |
| `anchor_indices()` | Find anchor atoms |
| `untag_anchor()` | Strip anchor tags |
| `subtract_anchor()` | Remove anchor atoms |

## Tech Stack

- **Python >= 3.10** (uses `from __future__ import annotations`, `X | Y` unions, `list[str]` generics)
- **pymatgen** - crystallographic data structures, CIF I/O, neighbor search
- **NetworkX** - bond graph representation, VF2 subgraph isomorphism
- **NumPy/SciPy** - numerical operations, SVD alignment
- **Click** - CLI framework
- **PyYAML** - batch config parsing
- **tqdm** - progress bars
- **Hatchling** - build system (PEP 517/621)
- **pytest** - testing (73 tests across 11 files)
- **ruff** - linting and formatting

## Development Commands

```bash
# Install for development
pip install -e ".[dev]"

# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_crystal.py -v

# Run with coverage
pytest tests/ --cov=mofforge --cov-report=term-missing

# Lint
ruff check src/ tests/

# Format
ruff format src/ tests/

# Lint + fix auto-fixable issues
ruff check src/ tests/ --fix

# Run CLI
mofforge --help
mofforge search -p crystal.cif -q query.xyz
mofforge replace -p crystal.cif -q query.xyz -r replacement.xyz -o output.cif
mofforge validate structure.cif
mofforge batch -c config.yaml
```

## Code Conventions

- All modules use `from __future__ import annotations` for forward references.
- Type hints are used throughout (Python 3.10+ syntax).
- `typing.TYPE_CHECKING` guards prevent circular imports.
- Dataclasses are used for data containers (`Crystal`, `MatchResult`, `Alignment`, `BondingRule`, etc.).
- Logging via `logging.getLogger("mofforge")` in every module.
- Line length: 100 characters.
- Import sorting: isort via ruff.
- Docstrings follow Google style with `Args:` and `Returns:` sections.
- Version is defined in `pyproject.toml` and read at runtime via `importlib.metadata`.

## File Formats

- **CIF files** (`.cif`): Crystallographic Information Files for periodic crystal structures. Read/written via pymatgen.
- **XYZ files** (`.xyz`): Simple atomic coordinate files for molecular fragments. Extended with `!` suffix on element symbols to mark anchor atoms.

## Test Data

- `tests/data/crystals/`: CIF files (IRMOF-1, UiO-66, Ti-MIL-125, MOF-74, etc.)
- `tests/data/moieties/`: XYZ fragment files
- `examples/data/crystals/`: CIF files for example scripts
- `examples/data/moieties/`: XYZ fragment files for example scripts

## Common Workflows

### Adding a new replacement mode
1. Add mode logic in `src/mofforge/replace/replace.py` in the `replace_pattern()` function's mode selection block.
2. Add tests in `tests/test_replace.py`.
3. Update CLI options in `src/mofforge/cli.py` if applicable.
4. Update `docs/python_api.md` and `docs/cli.md`.

### Adding a new validation check
1. Add the check function `_check_*()` in `src/mofforge/validation.py`.
2. Add a parameter to `validate_structure()` to enable/disable it.
3. Add results field to `ValidationReport` dataclass.
4. Add tests in `tests/test_validation.py`.

### Adding a new CLI subcommand
1. Add the command function with `@main.command()` in `src/mofforge/cli.py`.
2. Use lazy imports inside the function for faster startup.
3. Add tests in `tests/test_cli.py`.
4. Update `docs/cli.md`.

### Adding support for a new file format
1. Create a new module in `src/mofforge/io/` (e.g., `poscar.py`).
2. Export read/write functions from `src/mofforge/io/__init__.py`.
3. Add to `src/mofforge/__init__.py` exports and `__all__`.
4. Add tests in `tests/test_io.py`.
5. Update `docs/python_api.md`.

### Working with Crystal objects
```python
from mofforge import Crystal, infer_bonds, fragment, find_pattern, replace_pattern

# Load and bond
parent = Crystal.from_cif("structure.cif")
parent = infer_bonds(parent, periodic=True)

# Fragment operations
query = fragment("fragment.xyz", fragment_path="./moieties")

# Search
match = find_pattern(query, parent)
print(f"{match.nb_locations()} locations, {match.nb_isomorphisms()} total matches")

# Replace
child = replace_pattern(match, replacement, nb_loc=6)
child.write_cif("output.cif")
```

## Known Limitations

- Version `0.1.0` (Alpha) -- API may change.
- Charge balance validation is a placeholder (requires oxidation state data).
- The `MofforgeConfig` singleton is not thread-safe.
- Large structures (>10,000 atoms) may be slow during bond inference due to pairwise distance calculations.
- SMARTS parser supports a subset of the full SMARTS specification (elements, bonds, rings, wildcards -- no recursive SMARTS or logical operators).

## Performance Notes

- Bond inference is the most expensive operation. For large crystals, the neighbor search via pymatgen (periodic) or pairwise distance (non-periodic) dominates.
- Bonding rules are cached via `functools.lru_cache` to avoid regenerating ~12,000 rules per call.
- The bonding rule lookup uses a dictionary keyed by `(species_i, species_j)` for O(1) lookups instead of linear scan.
- The VF2 algorithm scales well for typical MOF queries (small fragments in moderate-sized unit cells).
- Batch processing supports parallel execution via `ProcessPoolExecutor`.
