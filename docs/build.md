# MOF Construction Guide

Complete guide to building MOFs from topology templates and building blocks using mofforge.

> **Early Development** -- The MOF construction subsystem is in alpha. The
> SMILES-to-building-block conversion, backend integrations, and adsorbate
> placement were partially developed with AI assistance and may produce
> unexpected results for edge cases. Please validate generated structures
> carefully and report any issues you encounter.

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
  - [TOBACCO 3.0 Setup](#tobacco-30-setup)
  - [Pormake Setup](#pormake-setup)
- [Configuration](#configuration)
  - [Configuration File (mofforge.toml)](#configuration-file-mofforgetoml)
  - [Environment Variables](#environment-variables)
  - [Priority Order](#priority-order)
- [Quick Start](#quick-start)
  - [Building with TOBACCO](#building-with-tobacco)
  - [Building with Pormake](#building-with-pormake)
- [MOFBuilder API](#mofbuilder-api)
  - [Creating a Builder](#creating-a-builder)
  - [Adding Building Blocks](#adding-building-blocks)
  - [Listing Topologies](#listing-topologies)
  - [Building a MOF](#building-a-mof)
  - [Managing Building Blocks](#managing-building-blocks)
  - [Status and Configuration](#status-and-configuration)
- [SMILES-to-Building-Block Conversion](#smiles-to-building-block-conversion)
  - [Connection Point Detection](#connection-point-detection)
  - [TOBACCO Format (CIF)](#tobacco-format-cif)
  - [Pormake Format (XYZ)](#pormake-format-xyz)
  - [Connection Modes](#connection-modes)
  - [Common MOF Linker Examples](#common-mof-linker-examples)
- [TOBACCO Backend Details](#tobacco-backend-details)
  - [Directory Layout](#directory-layout)
  - [How TOBACCO Isolation Works](#how-tobacco-isolation-works)
  - [TOBACCO Configuration](#tobacco-configuration)
  - [Parallel Execution](#parallel-execution)
  - [Output Structure](#output-structure)
- [Pormake Backend Details](#pormake-backend-details)
  - [In-Memory Registry vs Database](#in-memory-registry-vs-database)
  - [Topology Node/Edge Type Mapping](#topology-nodeedge-type-mapping)
  - [Output Naming Convention](#output-naming-convention)
- [Build-then-Modify Workflow](#build-then-modify-workflow)
- [CLI Reference](#cli-reference)
  - [mofforge build](#mofforge-build)
  - [mofforge build-status](#mofforge-build-status)
  - [mofforge build-list](#mofforge-build-list)
- [API Reference](#api-reference)
  - [Classes](#classes)
  - [Functions](#functions)
- [Troubleshooting](#troubleshooting)

---

## Overview

mofforge can construct MOFs from scratch by combining **topology templates** (network blueprints like pcu, dia, etc.) with **building blocks** (metal clusters as nodes, organic linkers as edges). This is powered by two established construction codes:

- **[TOBACCO 3.0](https://github.com/tobacco-mof/tobacco_3.0)** -- a template-based MOF assembler. Building blocks are placed at template sites according to the topology. TOBACCO is not pip-installable; it must be cloned separately.
- **[Pormake](https://github.com/Sangwon91/PORMAKE)** -- a topology-driven porous materials maker that uses the RCSR database. Pormake is pip-installable.

mofforge wraps both behind a unified `MOFBuilder` API, so you can switch between backends with a single parameter change.

Additionally, mofforge can convert **SMILES strings** directly into building block files (CIF for TOBACCO, XYZ for Pormake) with automatic detection of connection points -- no manual building block preparation needed for common linker chemistries.

---

## Prerequisites

### TOBACCO 3.0 Setup

TOBACCO is not distributed via pip. You must clone it and tell mofforge where to find it.

**1. Clone TOBACCO:**

```bash
git clone https://github.com/tobacco-mof/tobacco_3.0.git
cd tobacco_3.0
```

No further installation is needed -- mofforge imports TOBACCO's modules directly from the cloned directory.

**2. Tell mofforge where TOBACCO lives** (choose one method):

Option A -- configuration file:

```toml
# mofforge.toml (place in your project directory or ~/.mofforge.toml)
[backends.tobacco]
path = "/path/to/tobacco_3.0"
```

Option B -- environment variable:

```bash
export MOFFORGE_TOBACCO_PATH="/path/to/tobacco_3.0"
```

Option C -- pass directly in code:

```python
builder = MOFBuilder(backend="tobacco", tobacco_path="/path/to/tobacco_3.0")
```

**3. Validate the installation:**

```python
from mofforge.build.config import validate_tobacco_path
from pathlib import Path

errors = validate_tobacco_path(Path("/path/to/tobacco_3.0"))
if errors:
    for e in errors:
        print(f"  ERROR: {e}")
else:
    print("TOBACCO installation is valid")
```

mofforge checks that the directory contains the required files (`tobacco.py`, `configuration.py`) and directories (`templates/`, `nodes/`, `edges/`).

### Pormake Setup

Pormake is pip-installable:

```bash
pip install pormake
# or install via mofforge's optional dependency group:
pip install mofforge[build]
```

No additional configuration is needed. Pormake ships with its own topology and building block databases.

---

## Configuration

### Configuration File (mofforge.toml)

Create a `mofforge.toml` file in your project directory or at `~/.mofforge.toml`:

```toml
[backends.tobacco]
path = "/absolute/path/to/tobacco_3.0"

[backends.pormake]
output_dir = "./pormake_output"
```

mofforge searches for this file in two locations (in order):

1. `./mofforge.toml` (current working directory)
2. `~/.mofforge.toml` (home directory)

The first file found is used.

### Environment Variables

| Variable | Description |
|----------|-------------|
| `MOFFORGE_TOBACCO_PATH` | Path to the TOBACCO 3.0 installation directory |
| `MOFFORGE_PORMAKE_OUTPUT_DIR` | Default output directory for Pormake builds |

### Priority Order

When the same setting is specified in multiple places, the highest-priority source wins:

1. **Explicit kwargs** (highest) -- passed to `MOFBuilder(...)` or `BuildConfig.load(...)`
2. **Environment variables** -- `MOFFORGE_TOBACCO_PATH`, etc.
3. **TOML file** (lowest) -- `mofforge.toml`

```python
# This tobacco_path will override both TOML and environment variable:
builder = MOFBuilder(backend="tobacco", tobacco_path="/my/override/path")
```

---

## Quick Start

### Building with TOBACCO

```python
from mofforge.build import MOFBuilder

# Create a builder with the TOBACCO backend
builder = MOFBuilder(backend="tobacco")

# Add building blocks (CIF files only for TOBACCO)
builder.add_node("my_node.cif")
builder.add_edge("my_edge.cif")

# List available topologies
topologies = builder.list_topologies()
print(topologies)  # ['pcu', 'dia', 'sra', ...]

# Build with a specific topology
result = builder.build(topology="pcu", output_dir="./output")

if result.success:
    print(f"Built in {result.elapsed_seconds:.1f}s")
    for p in result.output_paths:
        print(f"  -> {p}")
else:
    print(f"Build failed: {result.errors}")
```

From the command line:

```bash
mofforge build \
    -b tobacco \
    -t pcu \
    -n my_node.cif \
    -e my_edge.cif \
    -o ./output
```

### Building with Pormake

```python
from mofforge.build import MOFBuilder

# Create a builder with the Pormake backend
builder = MOFBuilder(backend="pormake")

# Add building blocks (files or pormake database names)
builder.add_node("Zn_paddle_wheel.cif")
builder.add_edge("BDC_linker.xyz")

# Build
result = builder.build(topology="pcu", output_dir="./output")

if result.success:
    print(f"Output: {result.output_paths}")
    # The Crystal object is also available directly:
    crystal = result.crystal
    crystal.write_cif("my_mof.cif")
```

From the command line:

```bash
mofforge build \
    -b pormake \
    -t pcu \
    -n Zn_paddle_wheel.cif \
    -e BDC_linker.xyz \
    -o ./output
```

---

## MOFBuilder API

`MOFBuilder` is the primary user-facing class. It provides a unified interface over both backends.

### Creating a Builder

```python
from mofforge.build import MOFBuilder

# TOBACCO backend (default)
builder = MOFBuilder(backend="tobacco")

# Pormake backend
builder = MOFBuilder(backend="pormake")

# With explicit configuration
builder = MOFBuilder(
    backend="tobacco",
    tobacco_path="/path/to/tobacco_3.0",
)

builder = MOFBuilder(
    backend="pormake",
    output_dir="./my_output",
    bb_dir="./my_building_blocks",
)
```

### Adding Building Blocks

Building blocks are categorized as **nodes** (metal clusters / SBUs) or **edges** (organic linkers):

```python
# Add a node from a CIF file
builder.add_node("Zn_paddle_wheel.cif")

# Add an edge from a file, with a custom name
builder.add_edge("BDC.cif", name="BDC_linker")

# Add with explicit connection points
builder.add_node("my_node.cif", connection_points=[0, 5, 10, 15])
```

**Arguments:**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `source` | `str \| Path` | required | Path to a CIF/XYZ/MOL2 file, or a SMILES string |
| `name` | `str \| None` | `None` | Name for the building block (auto-generated from filename if `None`) |
| `connection_points` | `list[int] \| None` | `None` | Atom indices marking connection sites |

If `name` is not provided:
- For files (`.cif`, `.xyz`, `.mol2`): the filename stem is used (e.g., `"BDC"` from `"BDC.cif"`)
- For SMILES strings: the first 20 characters are used

### Listing Topologies

```python
# List all available topologies
topos = builder.list_topologies()
print(topos)  # ['pcu', 'dia', 'sra', 'acs', ...]

# Get a description of a specific topology
desc = builder.describe_topology("pcu")
print(desc)
```

### Building a MOF

```python
result = builder.build(
    topology="pcu",
    output_dir="./output",
)
```

The `build()` method returns a `BuildResult` object:

```python
result.success           # bool: did the build succeed?
result.output_paths      # list[Path]: paths to generated CIF files
result.crystal           # Crystal | None: the first output as a Crystal object
result.errors            # list[str]: error messages (if any)
result.elapsed_seconds   # float: wall-clock time
result.backend           # str: which backend was used
result.metadata          # dict: backend-specific metadata
```

**Arguments:**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `topology` | `str` | required | Topology name (e.g., `"pcu"`, `"dia"`) |
| `output_dir` | `str \| Path` | `"."` | Directory for output files |
| `**options` | | | Backend-specific options (e.g., `parallel=True` for TOBACCO) |

### Managing Building Blocks

```python
# List what's registered
nodes = builder.list_nodes()
edges = builder.list_edges()

# Preview what would be removed (dry run, the default)
builder.remove_edges(["BDC"], dry_run=True)

# Actually remove
builder.remove_edges(["BDC"], dry_run=False)

# Clear all building blocks of a type
builder.clear_nodes(dry_run=False)
builder.clear_edges(dry_run=False)
```

All remove/clear operations default to `dry_run=True` for safety. You must explicitly pass `dry_run=False` to perform the actual deletion.

**Copying from a backend's database:**

```python
# List what's available in the database
available = builder.copy_from_database(role="node", names=None)
print(available)  # shows available database entries

# Copy specific entries (dry run first)
builder.copy_from_database(role="node", names=["Zn_pw"], dry_run=True)

# Actually copy
builder.copy_from_database(role="node", names=["Zn_pw"], dry_run=False)
```

### Status and Configuration

```python
# Get backend status (counts, paths, availability)
status = builder.status()
print(status)

# Get current configuration
config = builder.get_configuration()

# Set a configuration value
builder.set_configuration("output_dir", "./new_output")
```

---

## SMILES-to-Building-Block Conversion

mofforge can convert SMILES strings into building block files suitable for TOBACCO or Pormake. This requires RDKit:

```bash
pip install rdkit
# or:
pip install mofforge[chem]
```

### Connection Point Detection

Before generating a building block file, mofforge analyzes the SMILES to find connection points automatically:

```python
from mofforge.build.smiles_to_bb import detect_connection_points

# Detects carboxylate groups in BDC linker
info = detect_connection_points("O=C(O)c1ccc(C(=O)O)cc1")
print(info.mode)                     # "carboxylate"
print(info.connection_atom_indices)  # [carbon indices of carboxylate groups]
print(info.carboxylate_groups)       # [CarboxylateGroup(...), ...]
```

Detection strategy:
1. **Carboxylate mode** (tried first): searches for `C(=O)[O,OH]` SMARTS patterns. If exactly `n_points` matches are found, this mode is used.
2. **Direct mode** (fallback): if no carboxylate groups are found, the two atoms at the ends of the molecule's graph diameter are used as connection points.

For molecules with -COOH groups that should be fully stripped during assembly:

```python
from mofforge.build.smiles_to_bb import detect_carboxylic_groups

info = detect_carboxylic_groups("OC(=O)c1ccc(C(=O)O)cc1")
print(info.mode)  # "carboxylic"
```

**Arguments:**

| Function | Arguments | Returns |
|----------|-----------|---------|
| `detect_connection_points(smiles, n_points=2)` | SMILES string, expected number of connection points | `ConnectionInfo` |
| `detect_carboxylic_groups(smiles)` | SMILES string (must have exactly 2 -COOH groups) | `ConnectionInfo` |

### TOBACCO Format (CIF)

```python
from mofforge.build.smiles_to_bb import smiles_to_tobacco_edge_cif

# BDC linker (1,4-benzenedicarboxylate)
path = smiles_to_tobacco_edge_cif(
    smiles="O=C(O)c1ccc(C(=O)O)cc1",
    output_path="BDC_edge.cif",
    name="BDC",
)
```

**Arguments:**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `smiles` | `str` | required | SMILES string for the linker molecule |
| `output_path` | `str \| Path` | required | Where to write the CIF file |
| `name` | `str` | `"edge"` | Data block name in the CIF |
| `cell_length` | `float` | `40.0` | Cubic unit cell side length in Angstrom |
| `uff_max_iters` | `int` | `2000` | Max iterations for UFF force-field optimization |
| `mode` | `str` | `"auto"` | Connection mode: `"auto"`, `"carboxylate"`, `"direct"`, or `"carboxylic"` |

The pipeline:
1. Parse SMILES with RDKit
2. Detect connection points (carboxylate or direct)
3. Add explicit hydrogens
4. Generate 3D coordinates (ETKDGv3)
5. Optimize geometry (UFF force field)
6. Build the CIF with connection atoms labeled as `X` (francium)
7. Write the TOBACCO-format CIF

### Pormake Format (XYZ)

```python
from mofforge.build.smiles_to_bb import smiles_to_pormake_edge_xyz

path = smiles_to_pormake_edge_xyz(
    smiles="O=C(O)c1ccc(C(=O)O)cc1",
    output_path="BDC_edge.xyz",
)
```

**Arguments:**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `smiles` | `str` | required | SMILES string for the linker molecule |
| `output_path` | `str \| Path` | required | Where to write the XYZ file |
| `uff_max_iters` | `int` | `2000` | Max iterations for UFF force-field optimization |

The output is a Pormake extended XYZ file. Connection points are marked with `X` dummy atoms placed at a fixed distance (0.75 Angstrom) outward from the connection atom along the bond axis.

### Connection Modes

Three connection modes control how SMILES molecules are converted to building blocks:

| Mode | Description | When used |
|------|-------------|-----------|
| **carboxylate** | Keeps carboxylate atoms in the building block. Adds `X` dummy atoms at the centroid of the two carboxylate oxygens. | Auto-detected when the molecule has exactly N carboxylate groups matching `C(=O)[O,OH]`. |
| **direct** | Labels the connection-point atoms themselves as `X`. Removes their hydrogen atoms. | Auto-detected fallback when no carboxylate groups are found. Uses graph-diameter endpoints. |
| **carboxylic** | Strips the entire -COOH group. Labels the anchor atom (the atom bonded to the carboxylate carbon) as `X`. | Must be explicitly requested via `mode="carboxylic"`. |

For TOBACCO CIF output, `mode="auto"` (the default) tries carboxylate first, then falls back to direct. You can force a specific mode:

```python
# Force direct mode even if carboxylates are present
smiles_to_tobacco_edge_cif(
    smiles="O=C(O)c1ccc(C(=O)O)cc1",
    output_path="BDC_direct.cif",
    mode="direct",
)

# Strip -COOH groups entirely
smiles_to_tobacco_edge_cif(
    smiles="OC(=O)c1ccc(C(=O)O)cc1",
    output_path="BDC_carboxylic.cif",
    mode="carboxylic",
)
```

### Common MOF Linker Examples

```python
from mofforge.build.smiles_to_bb import smiles_to_tobacco_edge_cif

# BDC (1,4-benzenedicarboxylate) -- the linker in MOF-5, UiO-66, MIL-53
smiles_to_tobacco_edge_cif("O=C(O)c1ccc(C(=O)O)cc1", "BDC.cif", name="BDC")

# BPDC (biphenyl-4,4'-dicarboxylate) -- the linker in UiO-67
smiles_to_tobacco_edge_cif(
    "O=C(O)c1ccc(-c2ccc(C(=O)O)cc2)cc1",
    "BPDC.cif",
    name="BPDC",
)

# NDC (2,6-naphthalenedicarboxylate) -- the linker in IRMOF-8
smiles_to_tobacco_edge_cif(
    "O=C(O)c1ccc2cc(C(=O)O)ccc2c1",
    "NDC.cif",
    name="NDC",
)

# NH2-BDC (2-amino-1,4-benzenedicarboxylate) -- the linker in UiO-66-NH2
smiles_to_tobacco_edge_cif(
    "O=C(O)c1ccc(C(=O)O)c(N)c1",
    "NH2_BDC.cif",
    name="NH2_BDC",
)
```

---

## TOBACCO Backend Details

### Directory Layout

TOBACCO expects a specific filesystem layout:

```
tobacco_3.0/
├── tobacco.py           # main TOBACCO module
├── configuration.py     # TOBACCO settings (bond length, charges, etc.)
├── templates/           # topology template CIF files
│   ├── pcu.cif
│   ├── dia.cif
│   └── ...
├── nodes/               # node building block CIF files (active)
│   ├── Zn_paddle_wheel.cif
│   └── ...
├── edges/               # edge building block CIF files (active)
│   ├── BDC.cif
│   └── ...
├── nodes_database/      # additional nodes available for copy_from_database
├── edges_database/      # additional edges available for copy_from_database
├── template_database/   # additional templates
└── output_cifs/         # TOBACCO's default output directory
```

When you call `builder.add_node()` or `builder.add_edge()`, mofforge copies the CIF file into the appropriate TOBACCO directory (`nodes/` or `edges/`).

**Important**: TOBACCO only accepts CIF files as building blocks. XYZ or MOL2 files are not supported.

### How TOBACCO Isolation Works

TOBACCO was designed to run as a standalone script from its own directory. mofforge handles this transparently using a context manager that:

1. Saves the current working directory, `sys.path`, and any cached TOBACCO modules
2. Changes to the TOBACCO directory and inserts it into `sys.path`
3. Runs the TOBACCO build
4. Restores everything to its original state

This means TOBACCO never permanently pollutes your Python environment. You can use mofforge and TOBACCO in the same process without conflicts.

### TOBACCO Configuration

TOBACCO has its own `configuration.py` file that controls build parameters (bond lengths, charges, cell optimization, etc.). You can read and modify it through mofforge:

```python
# Read current configuration
config = builder.get_configuration()
print(config)

# Set a configuration value
builder.set_configuration("bond_length", 1.5)
```

The `get_configuration()` method parses the Python configuration file using `ast.literal_eval` and returns a dictionary. `set_configuration()` performs a regex replacement in the file.

### Parallel Execution

TOBACCO supports building from multiple templates simultaneously:

```python
# Build all templates in parallel
result = builder.build(topology="pcu", output_dir="./output", parallel=True)
```

### Output Structure

After a build, the `BuildResult` contains:

- `output_paths`: paths to all generated CIF files, copied from TOBACCO's `output_cifs/` to your specified `output_dir`
- `crystal`: the first output file loaded as a `Crystal` object (ready for further modification)
- `metadata`: backend-specific information about the build

You can also inspect TOBACCO's output directory directly:

```python
backend = builder.backend
outputs = backend.list_outputs()  # organized dict of CIF files
```

---

## Pormake Backend Details

### In-Memory Registry vs Database

Pormake maintains two sources of building blocks:

1. **In-memory registry**: building blocks added via `add_node()` / `add_edge()`. These are stored as `BuildingBlock` objects and persist for the lifetime of the `MOFBuilder` instance.
2. **Pormake database**: the built-in database shipped with pormake. You can list and use database entries directly.

When listing building blocks, both sources are merged:

```python
# Lists both registered and database building blocks
all_edges = builder.list_edges()
```

You can copy database entries into the registry:

```python
# See what's in the database
available = builder.copy_from_database(role="edge", names=None)

# Register database entries for use
builder.copy_from_database(role="edge", names=["E41"], dry_run=False)
```

If a `bb_dir` was specified at construction, file-based building blocks are also copied to that directory for persistence:

```python
builder = MOFBuilder(backend="pormake", bb_dir="./my_bbs")
builder.add_edge("BDC.xyz")  # file is also copied to ./my_bbs/
```

### Topology Node/Edge Type Mapping

Pormake topologies define specific node and edge **types** (e.g., a topology might have 2 node types and 1 edge type). mofforge handles the mapping automatically:

- If you provide **1 node**, it is broadcast to all node types in the topology
- If you provide **1 edge**, it is broadcast to all edge types
- If you provide **multiple** nodes/edges, they are mapped to types in order

```python
# Single node/edge -- automatically broadcast to all types
builder.add_node("my_node.cif")
builder.add_edge("my_edge.xyz")
result = builder.build(topology="pcu")

# Multiple nodes for topologies with multiple node types
builder.add_node("node_type_A.cif")
builder.add_node("node_type_B.cif")
builder.add_edge("edge.xyz")
result = builder.build(topology="some_multi_node_topology")
```

### Output Naming Convention

Pormake output CIF files are named using the pattern:

```
{topology}_{node_names}_{edge_names}.cif
```

For example: `pcu_Zn_pw_BDC.cif`

---

## Build-then-Modify Workflow

A common workflow is to build a MOF and then modify it using mofforge's search-and-replace system. The `Pipeline` class supports this directly:

```python
from mofforge import Pipeline

# Build a MOF, then functionalize its linkers
child = (
    Pipeline.build_mof(
        backend="pormake",
        topology="pcu",
        nodes=["Zn_node.cif"],
        edges=["BDC_edge.cif"],
        fragment_path="./moieties",
    )
    .replace(
        query="2-!-p-phenylene.xyz",
        replacement="2-amino-p-phenylene.xyz",
        nb_loc=6,
    )
    .validate()
    .build(name="functionalized_MOF")
)

child.write_cif("functionalized_MOF.cif")
```

Or manually:

```python
from mofforge import Crystal, infer_bonds, fragment, find_pattern, replace_pattern
from mofforge.build import MOFBuilder

# Step 1: Build
builder = MOFBuilder(backend="pormake")
builder.add_node("Zn_node.cif")
builder.add_edge("BDC_edge.cif")
result = builder.build(topology="pcu")

# Step 2: Modify
parent = result.crystal
parent = infer_bonds(parent, periodic=True)

query = fragment("2-!-p-phenylene.xyz", fragment_path="./moieties")
replacement = fragment("2-amino-p-phenylene.xyz", fragment_path="./moieties")

match = find_pattern(query, parent)
child = replace_pattern(match, replacement, nb_loc=6)
child.write_cif("modified_MOF.cif")
```

---

## CLI Reference

### mofforge build

Build a MOF from topology and building blocks.

```
mofforge build [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-b, --backend` | `tobacco \| pormake` | `tobacco` | Which backend to use |
| `-t, --topology` | `TEXT` | required | Topology name (e.g., `pcu`, `dia`) |
| `-n, --node` | `PATH` | (multiple) | Node building block file(s) |
| `-e, --edge` | `PATH` | (multiple) | Edge building block file(s) |
| `-o, --output` | `DIR` | `.` | Output directory |
| `--tobacco-path` | `PATH` | from config | Override TOBACCO directory |
| `--parallel` | flag | off | Run TOBACCO in parallel mode |
| `-v, --verbose` | flag | off | Verbose output |

**Examples:**

```bash
# Build with TOBACCO
mofforge build -b tobacco -t pcu -n node.cif -e edge.cif -o ./output

# Build with Pormake
mofforge build -b pormake -t pcu -n node.cif -e edge.xyz -o ./output

# Multiple building blocks
mofforge build -b tobacco -t pcu \
    -n node_A.cif -n node_B.cif \
    -e edge.cif \
    -o ./output

# With parallel execution (TOBACCO only)
mofforge build -b tobacco -t pcu -n node.cif -e edge.cif --parallel
```

### mofforge build-status

Show the status of a build backend.

```
mofforge build-status [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-b, --backend` | `tobacco \| pormake` | `tobacco` | Backend to check |
| `--tobacco-path` | `PATH` | from config | Override TOBACCO directory |
| `-v, --verbose` | flag | off | Verbose output |

**Example:**

```bash
mofforge build-status -b tobacco
```

Output:

```json
{
  "templates": 25,
  "nodes": 3,
  "edges": 2,
  "outputs": 10,
  "configuration": { ... }
}
```

### mofforge build-list

List available topologies, nodes, or edges.

```
mofforge build-list [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-b, --backend` | `tobacco \| pormake` | `tobacco` | Backend to query |
| `--type` | `topologies \| nodes \| edges` | required | What to list |
| `--tobacco-path` | `PATH` | from config | Override TOBACCO directory |
| `-v, --verbose` | flag | off | Verbose output |

**Examples:**

```bash
# List topologies
mofforge build-list -b pormake --type topologies

# List registered nodes
mofforge build-list -b tobacco --type nodes

# List registered edges
mofforge build-list -b pormake --type edges
```

---

## API Reference

### Classes

| Class | Module | Description |
|-------|--------|-------------|
| `MOFBuilder` | `mofforge.build.builder` | Unified facade for MOF construction |
| `BuildingBlock` | `mofforge.build.base` | A node or edge building block |
| `BuildResult` | `mofforge.build.base` | Outcome of a build operation |
| `Topology` | `mofforge.build.base` | A network topology descriptor |
| `BuildConfig` | `mofforge.build.config` | Configuration loaded from TOML/env/kwargs |
| `ConfigError` | `mofforge.build.config` | Exception for invalid build configuration |
| `ConnectionInfo` | `mofforge.build.smiles_to_bb` | Result of connection point detection |
| `CarboxylateGroup` | `mofforge.build.smiles_to_bb` | Atom indices for one carboxylate group |
| `TobaccoBackend` | `mofforge.build.tobacco_backend` | TOBACCO backend (implements `BuilderBackend`) |
| `PormakeBackend` | `mofforge.build.pormake_backend` | Pormake backend (implements `BuilderBackend`) |

### Functions

| Function | Module | Description |
|----------|--------|-------------|
| `detect_connection_points(smiles, n_points=2)` | `mofforge.build.smiles_to_bb` | Auto-detect connection points (carboxylate or direct) |
| `detect_carboxylic_groups(smiles)` | `mofforge.build.smiles_to_bb` | Detect -COOH groups specifically |
| `smiles_to_tobacco_edge_cif(smiles, output_path, ...)` | `mofforge.build.smiles_to_bb` | Convert SMILES to TOBACCO-format CIF |
| `smiles_to_pormake_edge_xyz(smiles, output_path, ...)` | `mofforge.build.smiles_to_bb` | Convert SMILES to Pormake-format XYZ |
| `validate_tobacco_path(path)` | `mofforge.build.config` | Validate a TOBACCO installation directory |
| `BuildConfig.load(**overrides)` | `mofforge.build.config` | Load merged configuration from all sources |

### BuildingBlock Properties

| Property | Type | Description |
|----------|------|-------------|
| `name` | `str` | Identifier for the building block |
| `role` | `"node" \| "edge"` | Whether this is a node (SBU) or edge (linker) |
| `source` | `Path \| str` | File path or SMILES string |
| `connection_points` | `list[int] \| None` | Atom indices marking connection sites |
| `is_smiles` | `bool` (property) | `True` if `source` looks like a SMILES string |

### BuildResult Properties

| Property | Type | Description |
|----------|------|-------------|
| `success` | `bool` | Whether the build succeeded |
| `output_paths` | `list[Path]` | Paths to generated CIF files |
| `crystal` | `Crystal \| None` | First output loaded as a Crystal object |
| `errors` | `list[str]` | Error messages (empty on success) |
| `elapsed_seconds` | `float` | Wall-clock build time |
| `backend` | `str` | Which backend was used |
| `metadata` | `dict` | Backend-specific metadata |

---

## Troubleshooting

### TOBACCO path not configured

```
ConfigError: TOBACCO path is not configured
```

You need to tell mofforge where TOBACCO is installed. See [TOBACCO 3.0 Setup](#tobacco-30-setup) for options (TOML file, environment variable, or direct kwarg).

### Missing required TOBACCO files

```
ConfigError: Invalid TOBACCO installation at /path/to/tobacco: Missing required file: tobacco.py
```

Your TOBACCO directory is missing required files. Make sure you cloned the full repository and haven't renamed any files. The directory must contain `tobacco.py`, `configuration.py`, and the `templates/`, `nodes/`, and `edges/` subdirectories.

### RDKit not installed

```
ImportError: rdkit is required for SMILES-to-building-block conversion.
```

Install RDKit:

```bash
pip install rdkit
# or:
pip install mofforge[chem]
```

### Pormake not installed

```
ModuleNotFoundError: No module named 'pormake'
```

Install Pormake:

```bash
pip install pormake
# or:
pip install mofforge[build]
```

### SMILES connection point detection fails

```
ValueError: Expected 2 carboxylate groups but found 0
```

The molecule doesn't have carboxylate groups detectable by the `C(=O)[O,OH]` SMARTS pattern. Options:

- Use `mode="direct"` to fall back to graph-diameter endpoint detection
- Use `mode="carboxylic"` if the molecule has -COOH groups in a different arrangement
- Provide connection points manually via `connection_points=[...]` when calling `add_edge()`

### RDKit cannot generate 3D coordinates

```
ValueError: RDKit could not generate 3-D coordinates for: '...'
```

Some SMILES strings are difficult for RDKit's ETKDGv3 embedding. Try simplifying the SMILES, or generate the building block CIF/XYZ manually using other tools (e.g., Avogadro, OpenBabel) and provide the file directly to `add_edge()`.

### TOBACCO build produces no output

Check:
1. That your building blocks are compatible with the chosen topology (correct number of connection points)
2. That `nodes/` and `edges/` in the TOBACCO directory contain your files (use `builder.status()` to verify)
3. TOBACCO's `configuration.py` settings (use `builder.get_configuration()`)
4. Run with verbose logging to see TOBACCO's output:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```
