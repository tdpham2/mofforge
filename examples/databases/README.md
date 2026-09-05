# Database discovery and screening

These lessons use core dependencies. Their default mode stages **fictional**
CSD-format and CoRE-format metadata plus explicitly mapped copies of bundled
structures. Read the [fixture provenance](../data/databases/README.md):
the demonstration property values are not measurements or predictions for the
mapped IRMOF-1/UiO-66 geometries.

The database classes create SQLite caches beside metadata. Lessons copy metadata
into their output directories first, keeping source inputs separate from caches.

## 1. Search CSD metadata

```bash
python examples/databases/csd_lookup.py
python examples/databases/csd_lookup.py --query "Fictional zinc" --field name
```

[csd_lookup.py](csd_lookup.py) defaults to fictional refcode DEMOZN and also
demonstrates name, DOI, CCDC-number, and formula searches. Expect DEMOZN for the
specific lookup examples; the broad Zn formula search also finds DEMOMIS.
The summary records the auto-detected field and returned records.

For a real dataset:

```bash
python examples/databases/csd_lookup.py --data-path /data/csd/MOF_subset.tab --query HKUST --field name
```

Supply your own licensed ConQuest tab-separated export. mofforge searches the
exported metadata; this API does not connect to CCDC to download structures.
An explicit field avoids ambiguity between short chemical names and refcodes.

**Try:** query the intentionally unmatched DEMONONE record, then use the bridge
lesson to see that a CSD match need not have a CoRE counterpart. Empty search
results are valid observations. Missing or malformed metadata remains an error.

CLI equivalent after staging demo data:

```bash
mofforge csd DEMOZN --data-path examples/_outputs/csd_lookup/demo_data/csd_demo.tab
```

## 2. Screen CoRE MOF records

```bash
python examples/databases/coremof_screening.py --metal Zn --pld-min 3.8 --water-stability-min 0.7
```

[coremof_screening.py](coremof_screening.py) first searches a metal field, then
combines pore limiting diameter, a water-stability score, and the
`All Solvent Removed` extension. Default screening selects two fictional Zn
records, one with a deliberately missing CIF. A 10,000 Å PLD threshold produces
an intentional empty set.

PLD/LCD are in Å, density is in g/cm³, surface area is in m²/g, and pore volume
is in cm³/g. Stability fields must be interpreted according to the source
dataset: water/solvent values are scores, while thermal stability is in °C.
Missing numerical values become `None`; they do not satisfy a numeric bound.

Use `--data-path /data/coremof/metadata.csv` for real metadata and adjust filters.
CIFs are distributed separately; a metadata hit alone is not a usable structure.

**Try:** remove a bound and inspect how missing-property or alternative-extension
records affect the candidate list. Avoid treating ASR and FSR variants as
independent discoveries of the same framework.

CLI lookup equivalent:

```bash
mofforge coremof Zn --field metal --data-path examples/_outputs/coremof_screening/demo_data/coremof_demo.csv
```

The CLI exposes lookup; combined property screening is demonstrated through the
Python API and MCP tools.

## 3. Bridge names and resolve local structures

```bash
python examples/databases/bridge_and_resolve.py
```

[bridge_and_resolve.py](bridge_and_resolve.py) searches CSD names, associates
CoRE records with each base refcode, then calls `resolve_structure_path` for
each processing variant. The default report contains:

- DEMOZN with two variants, mapped to copies of the same teaching geometry.
- DEMOZR with a CIF but incomplete numerical properties.
- DEMOMIS with metadata but no CIF.
- DEMONONE with no CoRE counterpart.

Real mode requires all three paths together:

```bash
python examples/databases/bridge_and_resolve.py --csd-data /data/csd/MOF_subset.tab --coremof-data /data/coremof/metadata.csv --structures-dir /data/coremof/structures --name HKUST
```

The resolver searches local filenames and processing suffixes; it does not fetch
remote CIFs. Prefer the full coreid and inspect the resolved path/variant,
especially when a broad refcode can match several files.

**Try:** continue with [screen and place](../workflows/README.md), retaining
missing-structure outcomes instead of assuming every screened record is ready.
For licensed/live integration tests, set `MOFFORGE_RUN_DATABASE_TESTS=1`,
`MOFFORGE_CSD_DATA_PATH`, `MOFFORGE_COREMOF_DATA_PATH`, and
`MOFFORGE_COREMOF_STRUCTURES_PATH`; optionally set
`MOFFORGE_EXAMPLE_MOF_NAME` to a known name in your export.
