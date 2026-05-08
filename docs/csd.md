# CSD Lookup Guide

mofforge includes a module for searching the Cambridge Structural Database (CSD) MOF subset. Given a MOF name, DOI, or CCDC deposition number, it returns the corresponding CSD REFcode(s) and associated metadata.

## Prerequisites

### CSD License

The CSD is a proprietary database maintained by the [Cambridge Crystallographic Data Centre (CCDC)](https://www.ccdc.cam.ac.uk/). A valid CSD license is required to access the data. mofforge does **not** include or distribute any CSD data.

### Exporting the MOF Subset from ConQuest

The CSD data must be exported as a **Tab Separated List** (`.tab` file) using the [ConQuest](https://www.ccdc.cam.ac.uk/solutions/software/conquest/) application:

1. Open ConQuest and load the MOF subset (or your desired subset of structures)
2. Select the entries you want to export
3. Go to **File > Export** (or use the export function)
4. Choose **Tab Separated List** as the output format
5. Save the file (e.g., `MOF_subset.tab`)

The resulting `.tab` file contains one row per crystal structure with columns for REFcode, chemical names, formula, publication metadata, crystallographic parameters, and CCDC remarks (which include DOI and CCDC deposition numbers).

## Configuration

Configure the path to your `.tab` file using any of these methods (in order of priority):

### 1. mofforge.toml (recommended for persistent use)

Create or edit `mofforge.toml` in your project directory or home directory (`~/.mofforge.toml`):

```toml
[csd]
data_path = "/path/to/MOF_subset.tab"
```

### 2. Environment variable

```bash
export MOFFORGE_CSD_DATA_PATH="/path/to/MOF_subset.tab"
```

### 3. Python API

```python
from mofforge import set_paths

set_paths(csd_data="/path/to/MOF_subset.tab")
```

### 4. CLI flag (one-off use)

```bash
mofforge csd ABACUF --data-path /path/to/MOF_subset.tab
```

## SQLite Cache

On first use, mofforge parses the Tab Separated List and builds a SQLite cache (`.db` file) alongside the source file. This makes subsequent lookups near-instant. The cache is automatically rebuilt if the source `.tab` file is modified.

You can safely delete the `.db` file at any time -- it will be regenerated on the next query.

## Python API

### Basic Usage

```python
from mofforge.csd import get_database

db = get_database()  # uses configured data path

# Look up by REFcode
rec = db.lookup_refcode("ABACUF")
print(rec.refcode)                    # "ABACUF"
print(rec.chemical_name_common)       # common name
print(rec.chemical_formula_moiety)    # chemical formula
print(rec.doi)                        # DOI (if available)
print(rec.ccdc_number)                # CCDC deposition number
print(rec.space_group)                # space group
print(rec.year)                       # publication year
```

### Unified Search

The `search()` method auto-detects the query type:

```python
# Auto-detects REFcode
result = db.search("ABACUF")

# Auto-detects DOI (starts with "10.")
result = db.search("10.1038/46248")

# Auto-detects CCDC number (5-8 digits)
result = db.search("1100034")

# Falls back to name search
result = db.search("HKUST-1")
result = db.search("UiO-66")
```

You can also specify the field explicitly:

```python
result = db.search("C18 H6 Cu3", field="formula")
result = db.search("copper trimesate", field="name")
```

### Search Result

```python
result = db.search("HKUST-1")

print(result.n_matches)   # number of matches
print(result.field)        # which field was searched
print(result.query)        # the original query

for rec in result.records:
    print(rec.summary())   # one-line summary
```

### Individual Search Methods

```python
# Exact REFcode lookup (returns CSDRecord or None)
rec = db.lookup_refcode("ABACUF")

# Substring search on chemical names
records = db.search_name("UiO", limit=20)

# DOI search (exact or partial)
records = db.search_doi("10.1038/46248")

# Formula search (substring match)
records = db.search_formula("C18 H6 Cu3")

# CCDC deposition number
records = db.search_ccdc("1100034")
```

### CSDRecord Fields

Each `CSDRecord` contains:

| Field | Description |
|-------|-------------|
| `refcode` | CSD reference code (e.g., `ABACUF`) |
| `chemical_name_systematic` | IUPAC systematic name |
| `chemical_name_common` | Common/trivial name |
| `chemical_formula_moiety` | Chemical formula |
| `doi` | DOI (extracted from CCDC remarks, may be `None`) |
| `ccdc_number` | CCDC deposition number (may be `None`) |
| `authors` | Publication authors |
| `journal` | Journal name |
| `volume`, `pages`, `year` | Publication details |
| `space_group` | Crystallographic space group |
| `cell_a`, `cell_b`, `cell_c` | Unit cell lengths |
| `cell_alpha`, `cell_beta`, `cell_gamma` | Unit cell angles |
| `cell_volume` | Unit cell volume |
| `temperature` | Measurement temperature |
| `r_factor` | Crystallographic R-factor |
| `raw` | Dict of all 61 original columns for advanced use |

## CLI

```bash
# REFcode lookup
mofforge csd ABACUF

# Name search with limit
mofforge csd "HKUST-1" -n 5

# DOI search with verbose output
mofforge csd "10.1038/46248" -v

# CCDC number lookup
mofforge csd 1100034 --field ccdc

# Formula search
mofforge csd "C18 H6 Cu3" --field formula

# One-off data path
mofforge csd ABACUF --data-path /path/to/MOF_subset.tab
```

See the [CLI Reference](cli.md#csd) for the full option list.

## Data Privacy

The CSD data is proprietary and must not be committed to version control or distributed publicly. mofforge's `.gitignore` excludes the `external_data/` directory by default. If you place your `.tab` file elsewhere, ensure it is also excluded from version control.
