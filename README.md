# mofforge

A Python find-and-replace tool for atomistic crystal structure models, especially metal-organic frameworks (MOFs). Built on [pymatgen](https://pymatgen.org/), [NetworkX](https://networkx.org/), and [SciPy](https://scipy.org/).

## Features

**Core**
- Pattern matching via VF2 graph isomorphism
- Fragment replacement with SVD-based Procrustes alignment
- Periodic boundary condition handling (cross-boundary bonds, reassembly)
- Anchor atom masking (`!`-tagged atoms) for targeted replacement
- Multiple replacement modes (all locations, N random, specific, random orientation)
- CIF and XYZ file I/O

**Advanced**
- SMARTS-like pattern matching — search with strings like `"[Zn]-[O]-C"` instead of XYZ files
- Multi-step Pipeline API — chain operations fluently
- Structure validation — detect steric clashes, unusual bonds, coordination issues
- Batch processing — YAML-driven bulk operations with parallel execution
- Provenance tracking — record what modifications were made and when
- CLI — command-line interface for all core operations

## Installation

```bash
pip install mofforge
```

For development:

```bash
git clone <repo-url>
cd mofforge
pip install -e ".[dev]"
```

### Dependencies

- Python >= 3.10
- pymatgen >= 2024.1.1
- networkx >= 3.0
- numpy >= 1.24
- scipy >= 1.10
- click >= 8.0
- pyyaml >= 6.0
- tqdm >= 4.60

## Quick Start

```python
from mofforge import Crystal, infer_bonds, fragment, find_pattern, replace_pattern

# Load parent crystal and infer bonds
parent = Crystal.from_cif("IRMOF-1.cif")
parent = infer_bonds(parent, periodic=True)

# Load query and replacement fragments from XYZ files
query = fragment("2-!-p-phenylene.xyz", fragment_path="./moieties")
replacement = fragment("2-acetylamido-p-phenylene.xyz", fragment_path="./moieties")

# Search and replace
match = find_pattern(query, parent)
child = replace_pattern(match, replacement, nb_loc=6)
child.write_cif("output.cif")
```

See the full [Python API Manual](docs/python_api.md) and [CLI Reference](docs/cli.md) for details.

## Documentation

- **[Python API Manual](docs/python_api.md)** — Complete guide to using mofforge as a Python library
- **[CLI Reference](docs/cli.md)** — Command-line interface reference

## Examples

All examples are in the `examples/` directory with sample data in `examples/data/`.

### Core examples

| Script | Description |
|--------|-------------|
| [`pattern_matching.py`](examples/pattern_matching.py) | Find p-phenylene linkers in IRMOF-1 |
| [`linker_functionalization.py`](examples/linker_functionalization.py) | Functionalize BDC linkers with acetylamido groups |
| [`selective_modification.py`](examples/selective_modification.py) | All 5 replacement modes demonstrated |
| [`structure_repair.py`](examples/structure_repair.py) | Repair missing hydrogen atoms from X-ray data |
| [`cleanup_and_repair.py`](examples/cleanup_and_repair.py) | Fix disordered rings and remove guest molecules |
| [`defect_engineering.py`](examples/defect_engineering.py) | Engineer missing-linker defects in UiO-66 |
| [`symmetry_analysis.py`](examples/symmetry_analysis.py) | Replacement with symmetry analysis and supercell construction |

### Advanced examples

| Script | Description |
|--------|-------------|
| [`string_pattern_search.py`](examples/string_pattern_search.py) | String-based pattern matching: `"[Zn]-[O]"`, rings, wildcards |
| [`multi_step_pipeline.py`](examples/multi_step_pipeline.py) | Chain multiple operations with provenance tracking |
| [`structure_validation.py`](examples/structure_validation.py) | Post-modification structure validation |

Run any example:

```bash
cd examples
python pattern_matching.py
python linker_functionalization.py --nb-loc 6
python structure_validation.py path/to/structure.cif
```

## How It Works

mofforge implements a two-stage pipeline for crystal structure modification:

### 1. Find (Pattern Matching)

The parent crystal's bonding network is represented as a labeled graph (atoms = nodes, bonds = edges). A query fragment is matched against this graph using the **VF2 subgraph isomorphism algorithm** from NetworkX.

Matches are grouped by **location** (which set of parent atoms) and **orientation** (which mapping within that set).

### 2. Replace (Fragment Installation)

For each matched location, the replacement fragment is aligned to the parent geometry using **SVD-based Procrustes alignment** (orthogonal point cloud registration). The algorithm:

1. Identifies corresponding atoms between replacement and parent
2. Centers both point clouds and solves for the optimal rotation via SVD
3. Applies the rigid-body transformation to position the replacement
4. Reconstructs bonds between replacement atoms and parent neighbors
5. Removes obsolete (matched query) atoms
6. Handles periodic boundary crossings via reassembly

### Anchor Atom Masking Convention

Atoms tagged with `!` in XYZ files (e.g., `H!`, `C!`) mark anchor sites:

- During **search**: tags are stripped so `H!` matches `H` in the parent
- During **replace**: tagged atoms and their parent counterparts are removed; untagged atoms are used for alignment

```
# 2-!-p-phenylene.xyz
10
C   -1.71069   0.96969  -0.46280
C   -0.48337   1.30874   0.11690
...
H!   1.06706   0.70670   1.48683   <-- this H will be replaced
H    0.00122   2.23972  -0.14750
```

## Project Structure

```
mofforge/
├── src/mofforge/
│   ├── core/          # Crystal, bonding, fragment loading
│   ├── search/        # VF2 isomorphism, MatchResult API
│   ├── replace/       # Alignment, reassembly, replacement pipeline
│   ├── io/            # CIF and XYZ I/O
│   ├── utils/         # Config, periodic boundary utilities
│   ├── smarts.py      # SMARTS-like pattern matching
│   ├── pipeline.py    # Multi-step Pipeline API
│   ├── validation.py  # Structure validation
│   ├── batch.py       # Batch processing
│   ├── provenance.py  # Provenance tracking
│   └── cli.py         # CLI (Click)
├── tests/             # pytest tests
├── examples/          # 10 example scripts + data
└── docs/              # Manuals
```

## Testing

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## License

MIT
