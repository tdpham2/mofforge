# Python API Manual

Complete guide to using mofforge as a Python library.

## Table of Contents

- [Loading Data](#loading-data)
  - [Crystal Structures (CIF)](#crystal-structures-cif)
  - [Fragments (XYZ)](#fragments-xyz)
  - [Configuring Paths](#configuring-paths)
- [Bond Inference](#bond-inference)
- [Pattern Matching](#pattern-matching)
  - [Basic Search](#basic-search)
  - [Interpreting Results](#interpreting-results)
  - [Disconnected Component Mode](#disconnected-component-mode)
- [Pattern Replacement](#pattern-replacement)
  - [Basic Replacement](#basic-replacement)
  - [Replacement Modes](#replacement-modes)
  - [Replacing with Nothing (Deletion)](#replacing-with-nothing-deletion)
  - [Convenience Function](#convenience-function)
  - [Keyword Arguments](#keyword-arguments)
- [Anchor Atom Masking](#anchor-atom-masking)
- [SMARTS-like Pattern Matching](#smarts-like-pattern-matching)
- [Pipeline API](#pipeline-api)
- [Structure Validation](#structure-validation)
- [Batch Processing](#batch-processing)
- [Provenance Tracking](#provenance-tracking)
- [Database Lookups](#database-lookups)
  - [CSD Database](#csd-database)
  - [CoRE MOF Database](#core-mof-database)
  - [CSD-to-CoreMOF Bridge](#csd-to-coremof-bridge)
- [API Reference](#api-reference)

---

## Loading Data

### Crystal Structures (CIF)

Load a periodic crystal structure from a CIF file:

```python
from mofforge import Crystal

parent = Crystal.from_cif("IRMOF-1.cif")

print(parent.n_atoms)       # number of atoms
print(parent.species)       # list of species labels
print(parent.frac_coords)   # fractional coordinates (N, 3)
print(parent.cart_coords)   # Cartesian coordinates in Angstroms (N, 3)
print(parent.lattice)       # pymatgen Lattice object
```

You can also create a Crystal from a pymatgen Structure:

```python
from pymatgen.core import Structure

structure = Structure.from_file("POSCAR")
xtal = Crystal.from_structure(structure, name="my_crystal")
```

### Fragments (XYZ)

Load a molecular fragment from an XYZ file:

```python
from mofforge import fragment

query = fragment("p-phenylene.xyz", fragment_path="./moieties")
```

The `fragment()` function:
- Reads Cartesian coordinates from the XYZ file
- Places the fragment in an arbitrary cubic box
- Infers bonds (non-periodic)
- Sorts atoms by bond degree (highest connectivity first) for search efficiency
- Moves anchor atoms (`!`-tagged) to the end of the atom list

**Arguments:**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `name` | `str` or `None` | required | XYZ filename, or `None` for an empty crystal |
| `fragment_path` | `str` or `Path` | `config.moiety_path` | Directory containing XYZ files |
| `bonding_rules` | `list[BondingRule]` or `None` | `None` | Custom bonding rules |
| `presort` | `bool` | `True` | Sort atoms by bond degree |

### Configuring Paths

Set default directories for crystal and moiety files:

```python
from mofforge import set_paths

set_paths(
    crystals="./data/crystals",
    moieties="./data/moieties",
)
```

Or access the global config directly:

```python
from mofforge import config

config.crystal_path = Path("./my_crystals")
config.moiety_path = Path("./my_moieties")
config.bond_pad = 0.25  # padding added to covalent radii sum (Angstroms)
```

---

## Bond Inference

Before searching, bonds must be inferred on the parent crystal:

```python
from mofforge import infer_bonds, remove_bonds

# Infer bonds with periodic boundary conditions
parent = infer_bonds(parent, periodic=True)

print(parent.n_bonds)  # number of bonds

# Examine bond graph (NetworkX)
for u, v, data in parent.bonds.edges(data=True):
    print(f"  {parent.species[u]}-{parent.species[v]}: "
          f"{data['distance']:.3f} A, cross_PB={data['cross_boundary']}")

# Remove all bonds
parent = remove_bonds(parent)
```

Bond inference uses covalent radii-based distance rules. Two atoms are bonded if their distance is less than the sum of their covalent radii plus a padding value (default 0.25 A).

**Arguments:**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `crystal` | `Crystal` | required | Crystal to infer bonds for |
| `periodic` | `bool` | `True` | Consider periodic boundary conditions |
| `bonding_rules` | `list[BondingRule]` or `None` | `None` | Custom bonding rules |

---

## Pattern Matching

### Basic Search

```python
from mofforge import find_pattern

match = find_pattern(query, parent)
```

This returns a `MatchResult` object containing all subgraph isomorphisms of the query fragment within the parent's bonding network.

### Interpreting Results

```python
match.nb_isomorphisms()   # total number of isomorphisms
match.nb_locations()      # number of unique locations (sets of parent atoms)
match.nb_ori_at_loc()     # list: number of orientations at each location

# Access individual isomorphisms
# match.isomorphisms[location_idx][orientation_idx] -> {query_atom: parent_atom}
first_mapping = match.isomorphisms[0][0]

# Extract parent atoms involved in all matches
hits = match.matched_substructures()
```

**Isomorphism structure:**

```
match.isomorphisms = [
    [                                    # Location 0
        {0: 233, 1: 306, 2: 318, ...},  #   Orientation 0
        {0: 233, 1: 318, 2: 306, ...},  #   Orientation 1
        ...
    ],
    [                                    # Location 1
        ...
    ],
    ...
]
```

### Disconnected Component Mode

To match isolated molecules (e.g., guest molecules in pores) rather than substructures of the framework:

```python
match = find_pattern(guest, parent, disconnected_component=True)
```

This only matches fragments that are complete connected components of the parent graph (not bonded to anything else).

---

## Pattern Replacement

### Basic Replacement

```python
from mofforge import replace_pattern

child = replace_pattern(match, replacement)
child.write_cif("output.cif")
```

### Replacement Modes

Five modes control where and how replacements are applied:

```python
# Mode 1: All locations, optimal orientation (default)
child = replace_pattern(match, replacement)

# Mode 2: N random locations, optimal orientation
child = replace_pattern(match, replacement, nb_loc=8)

# Mode 3: Specific locations, optimal orientation
child = replace_pattern(match, replacement, loc=[0, 5, 10, 15])

# Mode 4: Specific locations AND orientations
child = replace_pattern(match, replacement, loc=[0, 1, 2], ori=[0, 1, 2])

# Mode 5: Random orientations (skip alignment optimization)
child = replace_pattern(match, replacement, random=True)

# Combined: N random locations with random orientations
child = replace_pattern(match, replacement, nb_loc=4, random=True)
```

### Replacing with Nothing (Deletion)

To delete matched substructures, pass `None` as the replacement:

```python
child = replace_pattern(match, None)
```

### Convenience Function

One-step search and replace:

```python
from mofforge import swap

child = swap(parent, query, replacement, nb_loc=6)
```

### Keyword Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `match` | `MatchResult` | required | Search results |
| `replacement` | `Crystal` or `None` | required | Replacement fragment, or `None` to delete |
| `random` | `bool` | `False` | Use random orientations |
| `nb_loc` | `int` | `0` | Number of random locations (0 = all) |
| `loc` | `list[int]` | `None` | Specific location indices |
| `ori` | `list[int]` | `None` | Specific orientation indices (one per loc) |
| `name` | `str` | `"new_xtal"` | Name for the output crystal |
| `verbose` | `bool` | `False` | Log replacement details |
| `remove_duplicates` | `bool` | `False` | Remove overlapping atoms of same species |
| `periodic_boundaries` | `bool` | `True` | Use PBC for duplicate checking |
| `reinfer_bonds` | `bool` | `False` | Re-infer bonds after replacement |
| `wrap` | `bool` | `True` | Wrap coordinates to [0, 1) |
| `auto_supercell` | `bool` | `False` | Auto-expand unit cell if replacement is too large |

---

## Anchor Atom Masking

Atoms in XYZ fragment files can be tagged with `!` to mark anchor atoms (replacement sites). These were previously referred to as "R-group" atoms; the concept is the same, but the preferred terminology is now "anchor atoms."

```
10
  comment
C   -1.710   0.970  -0.463
C   -0.483   1.309   0.117
...
H!   1.067   0.707   1.487    <-- anchor atom: will be replaced
H    0.001   2.240  -0.148
```

**How masking works:**

1. **During search**: `!` tags are stripped (`H!` matches `H` in the parent). The full fragment (including anchor atoms) participates in matching.

2. **During replacement**: The unmasked portion of the query is found inside the replacement fragment to establish correspondence. Anchor atoms and their parent counterparts are deleted; unmasked atoms are used for SVD alignment.

**Utility functions:**

```python
from mofforge import anchor_indices, untag_anchor, subtract_anchor

species = ["C", "H!", "O", "C!"]

anchor_indices(species)    # [1, 3]
untag_anchor(species)      # ["C", "H", "O", "C"]

# Remove anchor atoms from a crystal
clean = subtract_anchor(crystal)
```

---

## SMARTS-like Pattern Matching

Search using string patterns instead of XYZ files:

```python
from mofforge import smarts_search

# Simple patterns
result = smarts_search("[Zn]-[O]", parent)         # Zn-O bonds
result = smarts_search("O-C-O", parent)             # carboxylate groups
result = smarts_search("[Zn]-[*]", parent)           # anything bonded to Zn

# Ring patterns
result = smarts_search("C1-C-C-C-C-C-1", parent)    # 6-membered C ring

print(result.nb_isomorphisms())
print(result.nb_locations())
```

**Supported syntax:**

| Syntax | Meaning | Example |
|--------|---------|---------|
| `C`, `Zn`, `O` | Element symbol | `C-O` |
| `[Zn]`, `[Cu]` | Bracketed atom | `[Zn]-[O]` |
| `[*]` or `*` | Any atom (wildcard) | `C-[*]-C` |
| `-` | Bond connector | `Zn-O-C` |
| `1`..`9` | Ring closure digit | `C1-C-C-C-C-C-1` |

You can also parse patterns into NetworkX graphs directly:

```python
from mofforge import parse_smarts

graph = parse_smarts("C1-C-C-C-C-C-1")
print(graph.number_of_nodes())  # 6
print(graph.number_of_edges())  # 6
```

---

## Pipeline API

Chain multiple operations fluently:

```python
from mofforge import Pipeline

child = (
    Pipeline("IRMOF-1.cif", fragment_path="./moieties")
    .replace(query="2-!-p-phenylene.xyz",
             replacement="2-acetylamido-p-phenylene.xyz",
             nb_loc=6)
    .replace(query="2-!-p-phenylene.xyz",
             replacement="2-nitro-p-phenylene.xyz",
             nb_loc=3)
    .validate()
    .build(name="multi_functionalized_IRMOF-1")
)
```

**Pipeline methods:**

| Method | Description |
|--------|-------------|
| `.replace(query, replacement, **kwargs)` | Queue a find-and-replace step |
| `.remove(guest, **kwargs)` | Queue a guest removal step |
| `.validate(**kwargs)` | Queue a validation check |
| `.build(name)` | Execute all steps, return final Crystal |
| `.build_all(name)` | Execute all steps, return list of intermediates |

**Access validation reports:**

```python
pipeline = Pipeline(parent, fragment_path="./moieties")
pipeline.replace(query="...", replacement="...")
pipeline.validate()
child = pipeline.build()

for report in pipeline.validation_reports:
    print(report.summary())
```

**Provenance tracking** is automatic. The final crystal's `.provenance` attribute records the full chain of operations:

```python
print(child.provenance.summary())
# Step 1: pipeline_start (2026-03-17T...)
# Step 2: replace (2026-03-17T...)
# Current: replace (2026-03-17T...)
```

---

## Structure Validation

Check a crystal structure for problems after modification:

```python
from mofforge import validate_structure

report = validate_structure(crystal)

print(report.is_valid)          # True/False
print(report.steric_clashes)    # [(atom_i, atom_j, distance), ...]
print(report.unusual_bonds)     # [(atom_i, atom_j, actual, expected), ...]
print(report.coordination_issues)  # [(atom_i, species, CN, expected_range), ...]
print(report.summary())         # human-readable report
```

**Arguments:**

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `crystal` | `Crystal` | required | Crystal to validate |
| `check_clashes` | `bool` | `True` | Check for steric clashes |
| `check_bonds` | `bool` | `True` | Check for unusual bond lengths |
| `check_coordination` | `bool` | `True` | Check metal coordination numbers |
| `check_charges` | `bool` | `False` | Check charge balance |
| `clash_tolerance` | `float` | `0.5` | Tolerance below vdW sum (A) |
| `bond_tolerance` | `float` | `0.3` | Fractional deviation threshold |

---

## Batch Processing

Process multiple structures from a YAML configuration:

```python
from mofforge import run_batch

results = run_batch("batch_config.yaml")

for r in results:
    print(f"{r.parent_name}: {'OK' if r.success else r.error}")
```

**YAML config format:**

```yaml
parents:
  - path: "structures/*.cif"

operations:
  - type: replace
    query: BDC.xyz
    replacement: NH2-BDC.xyz
  - type: validate

output:
  directory: results/
  format: cif
  naming: "{parent_name}_functionalized"

parallel: 4
moiety_path: ./data/moieties
```

**Supported operation types:** `replace`, `remove`, `validate`.

---

## Provenance Tracking

Every Crystal can carry provenance metadata:

```python
from mofforge import Provenance

# Manually create provenance
prov = Provenance(
    parent="IRMOF-1.cif",
    query="BDC.xyz",
    replacement="NH2-BDC.xyz",
    operation="replace",
    parameters={"nb_loc": 6},
)

# Save/load
prov.to_json("provenance.json")
loaded = Provenance.from_json("provenance.json")

# Chain provenance records
new_prov = prov.chain(Provenance(operation="validate"))
print(new_prov.summary())
```

The Pipeline API attaches provenance automatically.

---

## Database Lookups

mofforge includes two database modules for searching MOF structures by metadata and properties.

### CSD Database

Search the Cambridge Structural Database MOF subset by REFcode, DOI, name, or formula. See the [CSD Lookup Guide](csd.md) for full details.

```python
from mofforge.csd import get_database

db = get_database(data_path="/path/to/MOF_subset.tab")
result = db.search("HKUST-1")
for rec in result.records:
    print(f"{rec.refcode}: {rec.chemical_name_common}")
```

### CoRE MOF Database

Query ~10,000 simulation-ready MOF structures from the CoRE MOF database. See the [CoRE MOF Guide](coremof.md) for full details.

```python
from mofforge.coremof import get_database

db = get_database(data_path="/path/to/CR_data_CSD_modified_20250227.csv")

# Search by metal type
for rec in db.search_metal("Cu", limit=5):
    print(f"{rec.coreid}: {rec.topology_single}")

# Screen by properties
candidates = db.screen(
    metal="Cu",
    lcd_min=8.0,
    water_stability_min=0.7,
    has_oms=True,
)
for rec in candidates:
    print(f"{rec.coreid}: LCD={rec.lcd:.1f}A")
```

### CSD-to-CoreMOF Bridge

Map CSD refcodes to simulation-ready CoreMOF identifiers, or search a MOF name across both databases:

```python
from mofforge.coremof import csd_to_coremof, search_csd_name

# Single refcode bridge
records = csd_to_coremof("ABAVIJ")
for rec in records:
    print(f"{rec.coreid} ({rec.extension})")

# Name search across CSD + CoreMOF
for br in search_csd_name("HKUST"):
    print(f"CSD: {br.csd_record.refcode}")
    for rec in br.coremof_records:
        print(f"  -> {rec.coreid}")
```

---

## API Reference

### Core Classes

| Class | Module | Description |
|-------|--------|-------------|
| `Crystal` | `mofforge.core.crystal` | Crystal structure with bond graph |
| `MatchResult` | `mofforge.search.search` | Pattern matching results |
| `Alignment` | `mofforge.replace.alignment` | SVD alignment parameters |
| `BondingRule` | `mofforge.core.bonding` | Species-pair bonding rule |
| `Provenance` | `mofforge.provenance` | Modification metadata |
| `ValidationReport` | `mofforge.validation` | Validation results |
| `Pipeline` | `mofforge.pipeline` | Multi-step operation chain |
| `CoreMOFDatabase` | `mofforge.coremof` | CoRE MOF database with search/screen methods |
| `CoreMOFRecord` | `mofforge.coremof` | Single CoRE MOF entry with properties |
| `CoreMOFSearchResult` | `mofforge.coremof` | Search results container |
| `BridgeResult` | `mofforge.coremof` | CSD record paired with CoreMOF matches |
| `CSDDatabase` | `mofforge.csd` | CSD database with search methods |
| `CSDRecord` | `mofforge.csd` | Single CSD entry |

### Core Functions

| Function | Module | Description |
|----------|--------|-------------|
| `fragment(name)` | `mofforge.core.moiety` | Load fragment from XYZ |
| `infer_bonds(crystal, periodic)` | `mofforge.core.bonding` | Infer bonding network |
| `remove_bonds(crystal)` | `mofforge.core.bonding` | Clear all bonds |
| `find_pattern(query, parent)` | `mofforge.search.search` | Find matching substructures |
| `replace_pattern(match, replacement)` | `mofforge.replace.replace` | Replace matched substructures |
| `swap(parent, query, replacement)` | `mofforge.replace.replace` | One-step search+replace |
| `validate_structure(crystal)` | `mofforge.validation` | Validate a structure |
| `smarts_search(pattern, parent)` | `mofforge.smarts` | SMARTS pattern search |
| `parse_smarts(pattern)` | `mofforge.smarts` | Parse pattern to graph |
| `run_batch(config_path)` | `mofforge.batch` | Run batch processing |
| `set_paths(crystals, moieties, csd_data, coremof_data)` | `mofforge.utils.config` | Set data directories |
| `get_database(data_path)` | `mofforge.csd` | Get CSD database singleton |
| `get_coremof_database(data_path)` | `mofforge.coremof` | Get CoRE MOF database singleton |
| `csd_to_coremof(refcode)` | `mofforge.coremof` | Map CSD refcode to CoreMOF entries |
| `search_csd_name(name)` | `mofforge.coremof` | Search CSD name, return CoreMOF bridges |
| `reassemble()` | `mofforge.core.moiety` | Reassemble fragments |
| `anchor_indices(species)` | `mofforge.core.moiety` | Get indices of anchor atoms |
| `untag_anchor(species)` | `mofforge.core.moiety` | Strip `!` tags from species |
| `subtract_anchor(crystal)` | `mofforge.core.moiety` | Remove anchor atoms from crystal |

### Crystal Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `Crystal.from_cif(path)` | `Crystal` | Load from CIF file |
| `Crystal.from_xyz(species, coords)` | `Crystal` | Create from arrays |
| `Crystal.from_structure(struct)` | `Crystal` | Wrap pymatgen Structure |
| `Crystal.empty()` | `Crystal` | Empty crystal |
| `crystal[indices]` | `Crystal` | Extract sub-crystal |
| `crystal + other` | `Crystal` | Combine crystals |
| `crystal.copy()` | `Crystal` | Deep copy |
| `crystal.wrap()` | `Crystal` | Wrap coords to [0,1) |
| `crystal.write_cif(path)` | `None` | Write CIF file |
| `crystal.write_xyz(path)` | `None` | Write XYZ file |
