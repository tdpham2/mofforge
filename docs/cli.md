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

---

## Notes

- All CIF files are read via pymatgen and support standard crystallographic formats.
- XYZ files follow the standard format with optional `!`-tagged species for R-group masking.
- The `--fragment-path` option is only needed if the query/replacement/guest file is given as a bare filename rather than a full path.
- When a full path is provided for `-q`/`-r`/`-g`, the parent directory is used as the fragment path automatically.
- Bond inference is performed automatically by all CLI commands.
- Output files are written in CIF format.
