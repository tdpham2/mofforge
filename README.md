# mofforge

A Python toolkit for building, modifying, and analyzing atomistic crystal structure models, especially metal-organic frameworks (MOFs).

mofforge provides three core capabilities:

1. **Build** MOFs from scratch using topology templates and building blocks, powered by [TOBACCO 3.0](https://github.com/tobacco-mof/tobacco_3.0) and [Pormake](https://github.com/Sangwon91/PORMAKE)
2. **Modify** existing structures via substructure search and fragment replacement, adapted partially from [PoreMatMod.jl](https://github.com/SimonEnsemble/PoreMatMod.jl)
3. **Analyze** structures with adsorbate placement, validation, and visualization

Built on [pymatgen](https://pymatgen.org/), [NetworkX](https://networkx.org/), and [SciPy](https://scipy.org/).

> **Alpha Software** -- mofforge is under active development. Parts of this
> codebase were developed with AI assistance and have not been exhaustively
> tested. Expect rough edges in the build and adsorbate subsystems.
> Always validate generated structures before using them in production workflows.
> If you encounter bugs or unexpected behavior, please open an issue on GitHub.

## Features

### Substructure Search & Replacement

- Pattern matching via VF2 graph isomorphism
- Fragment replacement with SVD-based Procrustes alignment
- Periodic boundary condition handling (cross-boundary bonds, reassembly)
- Anchor atom masking (`!`-tagged atoms) for targeted replacement
- Multiple replacement modes (all locations, N random, specific, random orientation)
- SMARTS-like pattern matching -- search with strings like `"[Zn]-[O]-C"` instead of XYZ files
- CIF and XYZ file I/O

### MOF Construction

- Build MOFs from topology + building blocks using two backends:
  - **TOBACCO 3.0** -- template-based assembly (requires separate installation)
  - **Pormake** -- topology-driven construction (`pip install pormake`)
- SMILES-to-building-block conversion with automatic connection point detection
- Unified `MOFBuilder` API across both backends
- Support for carboxylate, direct, and carboxylic connection modes

### Adsorbate Placement

- Automatic adsorption site detection (void sites via 3D grid sampling, open-metal sites via coordination analysis)
- 20 built-in adsorbate molecules (H2, N2, CO2, H2O, CH4, NH3, and more)
- Random orientation and intermolecular distance filtering
- Integration with the validation system for clash detection

### Visualization & AI Integration

- Structure rendering to PNG via 3Dmol.js + Playwright (ball-stick, stick, sphere representations)
- MCP (Model Context Protocol) server exposing 8 tools for AI agent integration
- Atom labels, unit cell edges, chemical formula overlay

### Pipeline & Automation

- Multi-step Pipeline API -- chain operations fluently
- Structure validation -- detect steric clashes, unusual bonds, coordination issues
- Batch processing -- YAML-driven bulk operations with parallel execution
- Provenance tracking -- record what modifications were made and when
- CLI -- command-line interface for all operations

## Installation

```bash
pip install mofforge
```

Optional dependencies for specific features:

```bash
pip install mofforge[build]   # Pormake backend
pip install mofforge[vis]     # PNG rendering (Playwright)
pip install mofforge[chem]    # SMILES conversion (RDKit)
pip install mofforge[mcp]     # MCP server for AI agents
pip install mofforge[all]     # everything above
```

For development:

```bash
git clone <repo-url>
cd mofforge
pip install -e ".[dev]"
```

### Dependencies

**Required:**

- Python >= 3.10
- pymatgen
- networkx >= 3.0
- numpy >= 1.24
- scipy >= 1.10
- click >= 8.0
- pyyaml >= 6.0

**Optional:**

- `pormake` -- Pormake construction backend
- `rdkit` -- SMILES parsing and 3D coordinate generation
- `playwright` -- headless browser for PNG rendering
- `mcp[cli]` -- Model Context Protocol server

**External (not pip-installable):**

- [TOBACCO 3.0](https://github.com/tobacco-mof/tobacco_3.0) -- must be cloned separately and configured via `mofforge.toml` or the `MOFFORGE_TOBACCO_PATH` environment variable
- [Architector](https://github.com/lanl/Architector) -- used for molecule visualization (future: general inorganic complex construction)

## Quick Start

### Modify an existing MOF

```python
from mofforge import Crystal, infer_bonds, fragment, find_pattern, replace_pattern

# Load parent crystal and infer bonds
parent = Crystal.from_cif("IRMOF-1.cif")
parent = infer_bonds(parent, periodic=True)

# Load query and replacement fragments from XYZ files
query = fragment("2-!-p-phenylene.xyz", fragment_path="./moieties")
replacement = fragment("2-acetylamido-p-phenylene.xyz", fragment_path="./moieties")

# Search and replace at 6 random locations
match = find_pattern(query, parent)
child = replace_pattern(match, replacement, nb_loc=6)
child.write_cif("output.cif")
```

### Build a MOF from scratch

```python
from mofforge.build import MOFBuilder

# Build using the Pormake backend
builder = MOFBuilder(backend="pormake")

# Add building blocks
builder.add_node("Zn_node.cif")
builder.add_edge("BDC_linker.cif")

# Build with a specific topology
result = builder.build(topology="pcu", output_dir="./output")
print(result.output_paths)
```

Or from SMILES:

```python
from mofforge.build.smiles_to_bb import smiles_to_tobacco_edge_cif

# Convert a SMILES linker to a TOBACCO-compatible building block
smiles_to_tobacco_edge_cif(
    smiles="O=C(O)c1ccc(C(=O)O)cc1",  # BDC linker
    output_path="BDC_edge.cif",
)
```

### Place adsorbates

```python
from mofforge import Crystal, infer_bonds
from mofforge.adsorbate import place_adsorbate

crystal = Crystal.from_cif("MOF-5.cif")
crystal = infer_bonds(crystal, periodic=True)

# Place 10 CO2 molecules at void sites
result = place_adsorbate(crystal, "CO2", n_adsorbates=10, strategy="void")
result.crystal.write_cif("MOF-5_loaded.cif")
print(f"Placed {result.n_adsorbates} adsorbates, {result.clashes} clashes")
```

## Documentation

- **[Python API Manual](docs/python_api.md)** -- complete guide to the search-and-replace system
- **[MOF Construction Guide](docs/build.md)** -- building MOFs with TOBACCO and Pormake
- **[CLI Reference](docs/cli.md)** -- command-line interface reference

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

### Build (MOF Construction)

mofforge wraps two established MOF construction codes behind a unified `MOFBuilder` API:

- **TOBACCO 3.0**: reads topology templates and building block CIF files, assembles them into periodic crystal structures. mofforge manages TOBACCO's filesystem layout and module isolation automatically.
- **Pormake**: uses the RCSR topology database and an in-memory building block registry. Building blocks can come from files, the Pormake database, or SMILES strings (converted automatically via RDKit).

Both backends produce CIF files. The output can be fed directly into the find-and-replace pipeline for post-synthetic modification.

See the [MOF Construction Guide](docs/build.md) for full details.

### Find (Pattern Matching)

The parent crystal's bonding network is represented as a labeled graph (atoms = nodes, bonds = edges). A query fragment is matched against this graph using the **VF2 subgraph isomorphism algorithm** from NetworkX.

Matches are grouped by **location** (which set of parent atoms) and **orientation** (which mapping within that set).

### Replace (Fragment Installation)

For each matched location, the replacement fragment is aligned to the parent geometry using **SVD-based Procrustes alignment** (orthogonal point cloud registration). The algorithm:

1. Identifies corresponding atoms between replacement and parent
2. Centers both point clouds and solves for the optimal rotation via SVD
3. Applies the rigid-body transformation to position the replacement
4. Reconstructs bonds between replacement atoms and parent neighbors
5. Removes obsolete (matched query) atoms
6. Handles periodic boundary crossings via reassembly

The find-and-replace paradigm -- subgraph isomorphism search, anchor atom masking, SVD-based geometric alignment, and multi-mode replacement -- is adapted from the approach described in [PoreMatMod.jl](https://github.com/SimonEnsemble/PoreMatMod.jl). mofforge reimplements these ideas in Python on top of pymatgen and NetworkX.

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
│   ├── build/         # MOF construction (MOFBuilder, TOBACCO, Pormake, SMILES-to-BB)
│   ├── adsorbate/     # Adsorption site detection and adsorbate placement
│   ├── io/            # CIF and XYZ I/O
│   ├── vis/           # PNG rendering via 3Dmol.js + Playwright
│   ├── mcp/           # MCP server for AI agent integration
│   ├── utils/         # Config, periodic boundary utilities
│   ├── smarts.py      # SMARTS-like pattern matching
│   ├── pipeline.py    # Multi-step Pipeline API
│   ├── validation.py  # Structure validation
│   ├── batch.py       # Batch processing
│   ├── provenance.py  # Provenance tracking
│   └── cli.py         # CLI (Click)
├── tests/             # pytest tests
├── examples/          # Example scripts + data
└── docs/              # Manuals (Python API, MOF construction, CLI)
```

## Acknowledgments & Inspirations

mofforge builds on and is inspired by several projects:

- **[PoreMatMod.jl](https://github.com/SimonEnsemble/PoreMatMod.jl)** -- The substructure search-and-replace paradigm in mofforge (VF2 isomorphism matching, anchor atom masking, SVD Procrustes alignment, multi-mode replacement) is heavily inspired by and adapted from PoreMatMod.jl by the Ensemble lab at Oregon State University.
- **[TOBACCO 3.0](https://github.com/tobacco-mof/tobacco_3.0)** -- Topologically-based crystal constructor, used as a MOF assembly backend.
- **[Pormake](https://github.com/Sangwon91/PORMAKE)** -- Porous materials maker, used as a MOF assembly backend.
- **[Architector](https://github.com/lanl/Architector)** -- Automated molecular/inorganic complex construction. Currently used for molecule visualization; planned for deeper integration in future releases.
- **[pymatgen](https://pymatgen.org/)** -- Crystal structure representation, CIF I/O, and periodic neighbor search.
- **[NetworkX](https://networkx.org/)** -- Bond graph representation and VF2 subgraph isomorphism.

## Testing

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## License

MIT
