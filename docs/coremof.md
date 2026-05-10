# CoRE MOF Database Guide

mofforge includes a module for searching the [CoRE MOF database](https://github.com/mtap-research/CoRE-MOF-Tools) — a curated dataset of ~10,000 computation-ready MOF structures with pre-computed properties. Given a CSD refcode, MOF name, metal type, topology, or property ranges, it returns matching CoRE MOF entries with their `coreid` identifiers for downstream simulation workflows.

## Data Source

Download the CoRE MOF dataset from Zenodo:

**https://zenodo.org/records/14510695**

The key file is `CR_data_CSD_modified_YYYYMMDD.csv` — a CSV containing ~9,800 MOF entries with 56 columns of identifiers, structural properties, stability data, and publication metadata.

## Configuration

Configure the path to your CSV file using any of these methods (in order of priority):

### 1. mofforge.toml (recommended for persistent use)

Create or edit `mofforge.toml` in your project directory or home directory (`~/.mofforge.toml`):

```toml
[coremof]
data_path = "/path/to/CR_data_CSD_modified_20250227.csv"
```

### 2. Environment variable

```bash
export MOFFORGE_COREMOF_DATA_PATH="/path/to/CR_data_CSD_modified_20250227.csv"
```

### 3. Python API

```python
from mofforge import set_paths

set_paths(coremof_data="/path/to/CR_data_CSD_modified_20250227.csv")
```

### 4. CLI flag (one-off use)

```bash
mofforge coremof QUERY --data-path /path/to/CR_data.csv
```

## SQLite Caching

On first access, the CSV is parsed and indexed into a SQLite database (placed alongside the CSV file, e.g. `CR_data_CSD_modified_20250227.db`). Subsequent lookups use this cache for fast queries. The cache is automatically rebuilt if the source CSV is newer than the cache.

## Python API

### Basic Usage

```python
from mofforge.coremof import get_database

db = get_database()

# Look up by CoreMOF ID
rec = db.lookup_coreid("2004[Co][rtl]3[ASR]1")
print(rec.coreid, rec.metal_types, rec.topology_single)

# Search by name
results = db.search_name("HKUST")
for r in results:
    print(r.summary())
```

### CSD Bridge — Refcode to CoreMOF ID

The primary use case: given a CSD refcode, find all corresponding CoRE MOF entries with their simulation-ready `coreid` identifiers.

```python
from mofforge.coremof import csd_to_coremof

# CSD refcode -> CoreMOF entries
records = csd_to_coremof("ABAVIJ")
for rec in records:
    print(f"{rec.coreid} ({rec.extension})")
# Output:
#   2004[Co][rtl]3[ASR]1 (All Solvent Removed)
#   2004[Co][rtl]3[FSR]1 (Free Solvent Removed)
```

Multiple entries may exist per CSD refcode, corresponding to different processing variants (All Solvent Removed, Free Solvent Removed, with ion).

### Unified Search

The `search()` method auto-detects the query type:

```python
db = get_database()

# Auto-detects as CSD refcode
result = db.search("HKUST1")

# Auto-detects as DOI
result = db.search("10.1038/46248")

# Auto-detects as metal element
result = db.search("Cu")

# Auto-detects as topology
result = db.search("pcu")

# Explicit field selection
result = db.search("HKUST-1", field="name")

print(result.summary())
```

### Individual Search Methods

```python
db = get_database()

# By full CoreMOF refcode
records = db.lookup_refcode("ABAVIJ_ASR_pacman")

# By CSD base refcode (bridge)
records = db.lookup_base_refcode("ABAVIJ")

# By DOI
records = db.search_doi("10.1039/b404485a")

# By metal type (uses normalized table — no false positives)
records = db.search_metal("Cu")

# By topology
records = db.search_topology("pcu")
records = db.search_topology("dia", nodes="all")  # AllNodes topology

# By KH gas storage class
records = db.search_kh_class("superstrong_high_loading")

# By open metal sites
records = db.search_oms(has_oms=True)
```

### Property-Based Screening

Screen MOFs by combining property ranges and categorical filters:

```python
db = get_database()

# Find Cu MOFs with large pores and high water stability
candidates = db.screen(
    metal="Cu",
    lcd_min=8.0,
    water_stability_min=0.7,
    has_oms=True,
    limit=20,
)

for rec in candidates:
    print(f"{rec.coreid}: LCD={rec.lcd:.1f}A, H2O_stab={rec.water_stability:.2f}")
```

Available screening parameters:

| Parameter | Type | Description |
|-----------|------|-------------|
| `lcd_min` / `lcd_max` | float | Largest Cavity Diameter (A) |
| `pld_min` / `pld_max` | float | Pore Limiting Diameter (A) |
| `density_min` / `density_max` | float | Density (g/cm3) |
| `asa_min` / `asa_max` | float | Accessible Surface Area (m2/g) |
| `pore_volume_min` / `pore_volume_max` | float | Pore Volume (cm3/g) |
| `void_fraction_min` / `void_fraction_max` | float | Void Fraction |
| `water_stability_min` | float | Water stability (0-1) |
| `solvent_stability_min` | float | Solvent stability (0-1) |
| `thermal_stability_min` | float | Thermal stability (C) |
| `metal` | str | Metal element symbol |
| `topology` | str | Topology name |
| `has_oms` | bool | Open metal sites |
| `kh_class` | str | KH gas storage class |
| `extension` | str | Processing variant |

## CoreMOFRecord Fields

| Field | Type | Description |
|-------|------|-------------|
| `coreid` | str | CoRE MOF unique identifier |
| `refcode` | str | Full refcode (e.g. `ABAVIJ_ASR_pacman`) |
| `base_refcode` | str | Base CSD refcode (e.g. `ABAVIJ`) |
| `name` | str | MOF common name (or `"-"`) |
| `lcd` | float? | Largest Cavity Diameter (A) |
| `pld` | float? | Pore Limiting Diameter (A) |
| `density` | float? | Density (g/cm3) |
| `asa` | float? | Accessible Surface Area (m2/g) |
| `pore_volume` | float? | Pore Volume (cm3/g) |
| `void_fraction` | float? | Void Fraction |
| `topology_single` | str | Topology (SingleNodes) |
| `topology_all` | str | Topology (AllNodes) |
| `metal_types` | str | Metal elements (comma-separated) |
| `has_oms` | bool | Has open metal sites |
| `thermal_stability` | float? | Thermal stability (C) |
| `water_stability` | float? | Water stability (0-1) |
| `solvent_stability` | float? | Solvent stability (0-1) |
| `kh_class` | str | KH classification |
| `doi` | str | Publication DOI |
| `year` | str | Publication year |
| `extension` | str | Processing type |
| `raw` | dict | All 56 original CSV columns |

Fields marked `float?` are `None` when the value is unknown or missing.

## CLI

### Search

```bash
# Auto-detect query type
mofforge coremof ABAVIJ

# Explicit field
mofforge coremof Cu --field metal

# Verbose output with properties
mofforge coremof "HKUST-1" --field name -v

# Limit results
mofforge coremof pcu --field topology --limit 20
```

### CSD Bridge

```bash
# Get CoreMOF IDs for a CSD refcode
mofforge coremof ABAVIJ --bridge
```

Output:
```
CoreMOF bridge: 2 match(es) for CSD refcode 'ABAVIJ'
  2004[Co][rtl]3[ASR]1
    Refcode:   ABAVIJ_ASR_pacman
    Extension: All Solvent Removed
    Metals:    Co
  2004[Co][rtl]3[FSR]1
    Refcode:   ABAVIJ_FSR_pacman
    Extension: Free Solvent Removed
    Metals:    Co
```
