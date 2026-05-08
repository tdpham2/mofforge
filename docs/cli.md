# CLI Reference

mofforge provides a command-line interface for all core operations. After installation, the `mofforge` command is available.

```bash
mofforge --help
mofforge --version
```

## Table of Contents

- [search](#search) — Find substructures in a crystal
- [replace](#replace) — Find and replace substructures
- [remove](#remove) — Remove guest molecules
- [validate](#validate) — Validate a crystal structure
- [batch](#batch) — Batch process multiple structures
- [build](#build) — Build a MOF from topology + building blocks
- [build-status](#build-status) — Show backend status
- [build-list](#build-list) — List topologies, nodes, or edges
- [csd](#csd) — Look up MOFs in the CSD database
- [render](#render) — Render a structure to PNG

---

## search

Search for a substructure (fragment) in a crystal structure.

```
mofforge search -p PARENT -q QUERY [OPTIONS]
```

**Required:**

| Option | Description |
|--------|-------------|
| `-p, --parent TEXT` | Path to the parent CIF file |
| `-q, --query TEXT` | Path to the query fragment XYZ file |

**Optional:**

| Option | Description |
|--------|-------------|
| `--disconnected` | Search for isolated components only (e.g., guest molecules) |
| `--fragment-path TEXT` | Directory containing fragment XYZ files |
| `-v, --verbose` | Enable verbose/debug output |

**Examples:**

```bash
# Search for p-phenylene in IRMOF-1
mofforge search -p IRMOF-1.cif -q moieties/p-phenylene.xyz

# Search for guest molecules (isolated components)
mofforge search -p loaded_MOF.cif -q moieties/acetylene.xyz --disconnected

# Verbose output
mofforge search -p IRMOF-1.cif -q moieties/p-phenylene.xyz -v
```

**Output:**

```
Loading parent: IRMOF-1.cif
Loading query: p-phenylene.xyz
Found 96 matches at 24 locations
  Isomorphisms: 96
  Locations: 24
  Orientations per location: [4, 4, 4, ...]
```

---

## replace

Find a substructure and replace it with a different fragment.

```
mofforge replace -p PARENT -q QUERY -r REPLACEMENT [OPTIONS]
```

**Required:**

| Option | Description |
|--------|-------------|
| `-p, --parent TEXT` | Path to the parent CIF file |
| `-q, --query TEXT` | Path to the query fragment XYZ file |
| `-r, --replacement TEXT` | Path to the replacement fragment XYZ file |

**Optional:**

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output TEXT` | `new_xtal.cif` | Output CIF file path |
| `--nb-loc INTEGER` | `0` (all) | Number of random locations to replace |
| `--random` | off | Use random orientations instead of optimal |
| `--validate` | off | Run structure validation on the output |
| `--fragment-path TEXT` | query's directory | Directory containing fragment XYZ files |
| `-v, --verbose` | off | Enable verbose output |

**Examples:**

```bash
# Replace all locations with optimal orientation (default)
mofforge replace \
    -p IRMOF-1.cif \
    -q moieties/2-!-p-phenylene.xyz \
    -r moieties/2-acetylamido-p-phenylene.xyz \
    -o functionalized.cif

# Replace at 6 random locations
mofforge replace \
    -p IRMOF-1.cif \
    -q moieties/2-!-p-phenylene.xyz \
    -r moieties/2-acetylamido-p-phenylene.xyz \
    --nb-loc 6 \
    -o partial.cif

# Replace with random orientations and validate
mofforge replace \
    -p IRMOF-1.cif \
    -q moieties/2-!-p-phenylene.xyz \
    -r moieties/2-nitro-p-phenylene.xyz \
    --random --validate \
    -o random_nitro.cif

# Repair missing hydrogens (replacement is superset of query)
mofforge replace \
    -p IRMOF-1_noH.cif \
    -q moieties/1,4-C-phenylene_noH.xyz \
    -r moieties/1,4-C-phenylene.xyz \
    -o repaired.cif
```

**Output:**

```
Loading parent: IRMOF-1.cif
Found 96 matches at 24 locations
Output written to: functionalized.cif
  Atoms: 592, Bonds: 680
```

---

## remove

Remove guest molecules from a crystal structure. This is equivalent to `replace` with `None` as the replacement, using `--disconnected` search mode.

```
mofforge remove -p PARENT -g GUEST [OPTIONS]
```

**Required:**

| Option | Description |
|--------|-------------|
| `-p, --parent TEXT` | Path to the parent CIF file |
| `-g, --guest TEXT` | Path to the guest fragment XYZ file |

**Optional:**

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output TEXT` | `clean.cif` | Output CIF file path |
| `--fragment-path TEXT` | guest's directory | Directory containing fragment XYZ files |
| `-v, --verbose` | off | Enable verbose output |

**Examples:**

```bash
# Remove acetylene guests
mofforge remove \
    -p SIFSIX-2-Cu-i.cif \
    -g moieties/acetylene.xyz \
    -o clean.cif

# Remove with verbose output
mofforge remove -p loaded_MOF.cif -g moieties/guest.xyz -v
```

**Output:**

```
Loading parent: SIFSIX-2-Cu-i.cif
Found 8 guest molecule(s)
Output written to: clean.cif
  Atoms: 168 (removed 32)
```

---

## validate

Validate a crystal structure for steric clashes, unusual bond lengths, and coordination geometry issues.

```
mofforge validate STRUCTURE [OPTIONS]
```

**Required:**

| Argument | Description |
|----------|-------------|
| `STRUCTURE` | Path to the CIF file to validate |

**Optional:**

| Option | Description |
|--------|-------------|
| `-v, --verbose` | Enable verbose/debug output |

**Examples:**

```bash
# Validate a structure
mofforge validate IRMOF-1.cif

# Validate a modified structure
mofforge validate functionalized_MOF.cif -v
```

**Output:**

```
Loading structure: IRMOF-1.cif
Validation Report (valid=True):
  All checks passed.
```

Or, if issues are found:

```
Validation Report (valid=False):
  Steric clashes: 3
    atoms 42-198: 1.523 A
    atoms 43-199: 1.523 A
    atoms 44-200: 1.523 A
  Coordination issues: 1
    atom 0 (Zn): CN=3 (expected 4-6)
```

---

## batch

Run batch processing on multiple crystal structures from a YAML configuration file.

```
mofforge batch -c CONFIG [OPTIONS]
```

**Required:**

| Option | Description |
|--------|-------------|
| `-c, --config TEXT` | Path to the YAML configuration file |

**Optional:**

| Option | Description |
|--------|-------------|
| `-v, --verbose` | Enable verbose output |

**YAML Configuration Format:**

```yaml
# batch_config.yaml

# Parent structures (supports glob patterns)
parents:
  - path: "structures/*.cif"

# Operations to apply (in order)
operations:
  - type: replace
    query: BDC.xyz
    replacement: NH2-BDC.xyz
    mode: all_optimal          # or: random, nb_loc_6
  - type: validate

# Output settings
output:
  directory: results/
  format: cif
  naming: "{parent_name}_functionalized"

# Parallelism (0 = sequential)
parallel: 4

# Path to fragment XYZ files
fragment_path: ./data/fragments
```

**Supported operation types:**

| Type | Description | Required fields |
|------|-------------|-----------------|
| `replace` | Find and replace | `query`, `replacement` |
| `remove` | Delete guest molecules | `guest` or `query` |
| `validate` | Validate structure | (none) |

**Replacement modes** (set via the `mode` field):

| Mode | Description |
|------|-------------|
| `all_optimal` (default) | All locations, optimal orientation |
| `random` | All locations, random orientation |
| `nb_loc_N` | N random locations (e.g., `nb_loc_6`) |

**Examples:**

```bash
# Run batch processing
mofforge batch -c batch_config.yaml

# Verbose output
mofforge batch -c batch_config.yaml -v
```

**Output:**

```
Batch processing 15 structures

Batch Results (15 structures):
  IRMOF-1: OK
    -> results/IRMOF-1_functionalized.cif
  UiO-66: OK
    -> results/UiO-66_functionalized.cif
  MOF-74: FAILED: No matches found
```

---

## build

Build a MOF from a topology template and building blocks. See the [MOF Construction Guide](build.md) for full details.

```
mofforge build [OPTIONS]
```

**Required:**

| Option | Description |
|--------|-------------|
| `-t, --topology TEXT` | Topology name (e.g., `pcu`, `dia`) |

**Optional:**

| Option | Default | Description |
|--------|---------|-------------|
| `-b, --backend [tobacco\|pormake]` | `tobacco` | Which construction backend to use |
| `-n, --node PATH` | | Node building block file (can be repeated) |
| `-e, --edge PATH` | | Edge building block file (can be repeated) |
| `-o, --output DIR` | `.` | Output directory |
| `--tobacco-path PATH` | from config | Override TOBACCO installation directory |
| `--parallel` | off | Run TOBACCO in parallel mode |
| `-v, --verbose` | off | Enable verbose output |

**Examples:**

```bash
# Build with TOBACCO
mofforge build -b tobacco -t pcu -n node.cif -e edge.cif -o ./output

# Build with Pormake
mofforge build -b pormake -t pcu -n node.cif -e edge.xyz -o ./output

# Multiple building blocks
mofforge build -b tobacco -t pcu -n node_A.cif -n node_B.cif -e edge.cif

# Parallel execution (TOBACCO only)
mofforge build -b tobacco -t pcu -n node.cif -e edge.cif --parallel
```

---

## build-status

Show the status of a build backend (number of templates, nodes, edges, outputs, and configuration).

```
mofforge build-status [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-b, --backend [tobacco\|pormake]` | `tobacco` | Backend to check |
| `--tobacco-path PATH` | from config | Override TOBACCO directory |
| `-v, --verbose` | off | Verbose output |

**Example:**

```bash
mofforge build-status -b tobacco
```

---

## build-list

List available topologies, nodes, or edges for a backend.

```
mofforge build-list [OPTIONS]
```

**Required:**

| Option | Description |
|--------|-------------|
| `--type [topologies\|nodes\|edges]` | What to list |

**Optional:**

| Option | Default | Description |
|--------|---------|-------------|
| `-b, --backend [tobacco\|pormake]` | `tobacco` | Backend to query |
| `--tobacco-path PATH` | from config | Override TOBACCO directory |
| `-v, --verbose` | off | Verbose output |

**Examples:**

```bash
# List topologies
mofforge build-list -b pormake --type topologies

# List registered nodes
mofforge build-list -b tobacco --type nodes
```

---

## csd

Look up MOF entries in the CSD (Cambridge Structural Database) by REFcode, DOI, CCDC deposition number, chemical name, or formula.

This command requires a CSD MOF subset data file exported as a **Tab Separated List** from the [ConQuest](https://www.ccdc.cam.ac.uk/solutions/software/conquest/) application. A CSD license from the CCDC is required. The data file is not distributed with mofforge.

```
mofforge csd QUERY [OPTIONS]
```

**Required:**

| Argument | Description |
|----------|-------------|
| `QUERY` | Search term (REFcode, DOI, CCDC number, or name) |

**Optional:**

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --field [auto\|refcode\|name\|doi\|formula\|ccdc]` | `auto` | Field to search (auto-detects by default) |
| `-n, --limit INT` | `10` | Maximum number of results to display |
| `--data-path TEXT` | from config | Path to the CSD Tab Separated List file |
| `-v, --verbose` | off | Show additional details (formula, journal, space group) |

**Auto-detection rules** (when `--field auto`):

| Query pattern | Detected as |
|---------------|-------------|
| Starts with `10.` | DOI |
| 5-8 digits | CCDC deposition number |
| 3-8 uppercase letters + optional digits | REFcode (falls back to name search if no match) |
| Anything else | Chemical name (substring search) |

**Configuration:**

The data file path can be set via (in order of priority):

1. `--data-path` CLI option
2. `set_paths(csd_data=...)` in Python
3. `MOFFORGE_CSD_DATA_PATH` environment variable
4. `[csd] data_path` in `mofforge.toml`

On first use, mofforge parses the Tab Separated List and builds a local SQLite cache (`.db` file alongside the `.tab` file) for fast subsequent lookups.

**Examples:**

```bash
# Look up a REFcode
mofforge csd ABACUF --data-path /path/to/MOF_subset.tab

# Search by common name
mofforge csd "HKUST-1" -n 5

# Search by DOI
mofforge csd "10.1038/46248" -v

# Search by CCDC deposition number
mofforge csd 1100034 --field ccdc

# Search by formula
mofforge csd "C18 H6 Cu3" --field formula
```

**Output:**

```
CSD lookup: 1 match(es) for 'ABACUF' (field: refcode)
  ABACUF: catena(Tetra-aqua-tetrakis(formato)-di-barium-copper)
    CCDC: 1100034
```

With `-v`:

```
CSD lookup: 1 match(es) for 'ABACUF' (field: refcode)
  ABACUF: catena(Tetra-aqua-tetrakis(formato)-di-barium-copper)
    DOI: 10.1016/j.molstruc.2004.03.051
    CCDC: 1100034
    Formula: (C6 H14 Ba2 Cu1 O16)n
    Journal: J.Mol.Struct. (2004)
    Space group: P-1
```

---

## render

Render a crystal structure (CIF or XYZ) to a PNG image using 3Dmol.js and Playwright.

```
mofforge render [OPTIONS]
```

**Required:**

| Option | Description |
|--------|-------------|
| `-i, --input TEXT` | Path to the CIF or XYZ file |

**Optional:**

| Option | Default | Description |
|--------|---------|-------------|
| `-o, --output TEXT` | `structure.png` | Output PNG file path |
| `--label-mode [sequential\|per_element\|none]` | `sequential` | Atom label style |
| `--representation [ball_stick\|stick\|sphere]` | `ball_stick` | Visual representation |
| `--width INT` | `800` | Image width in pixels |
| `--height INT` | `600` | Image height in pixels |
| `--show-unit-cell` | auto | Show unit cell edges (auto-enabled for CIF) |
| `--show-formula` | on | Show chemical formula overlay |
| `--bg-color TEXT` | `white` | Background color |
| `--label-size INT` | `14` | Atom label font size |

**Examples:**

```bash
# Render a CIF file
mofforge render -i IRMOF-1.cif -o irmof1.png

# Render with specific options
mofforge render -i structure.cif -o output.png \
    --representation sphere \
    --width 1200 --height 900 \
    --label-mode per_element

# Render an XYZ fragment
mofforge render -i fragment.xyz -o fragment.png --label-mode none
```

---

## Common Workflows

### Substructure search only

```bash
mofforge search -p parent.cif -q query.xyz
```

### Functionalize a MOF at N random sites

```bash
mofforge replace -p parent.cif -q query_with_mask.xyz -r replacement.xyz --nb-loc 6 -o output.cif
```

### Repair missing hydrogens

```bash
mofforge replace -p noH_structure.cif -q bare_ring.xyz -r ring_with_H.xyz -o repaired.cif
```

### Remove guests then validate

```bash
mofforge remove -p loaded.cif -g guest.xyz -o clean.cif
mofforge validate clean.cif
```

### Batch process a directory of MOFs

```bash
mofforge batch -c my_config.yaml -v
```

### Build a MOF and render it

```bash
mofforge build -b pormake -t pcu -n node.cif -e edge.cif -o ./output
mofforge render -i ./output/pcu_node_edge.cif -o mof_preview.png
```

### Look up a MOF and find related structures

```bash
# Find the REFcode for HKUST-1
mofforge csd "HKUST-1" --data-path /path/to/MOF_subset.tab

# Find all structures from a specific paper
mofforge csd "10.1126/science.283.5405.1148" --field doi -v
```

### Build, then functionalize

```bash
# Build the base MOF
mofforge build -b tobacco -t pcu -n Zn_pw.cif -e BDC.cif -o ./output

# Functionalize linkers in the built MOF
mofforge replace \
    -p ./output/pcu_MOF.cif \
    -q moieties/2-!-p-phenylene.xyz \
    -r moieties/2-amino-p-phenylene.xyz \
    --nb-loc 6 \
    -o functionalized.cif
```

---

## Notes

- All CIF files are read via pymatgen and support standard crystallographic formats.
- XYZ files follow the standard format with optional `!`-tagged species for R-group masking.
- The `--fragment-path` option is only needed if the query/replacement/guest file is given as a bare filename rather than a full path.
- When a full path is provided for `-q`/`-r`/`-g`, the parent directory is used as the fragment path automatically.
- Bond inference is performed automatically by all CLI commands.
- Output files are written in CIF format.
